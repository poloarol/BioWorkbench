
import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd

from scipy import sparse

_EXCLUDED_COLOR_COLUMNS = {"_cell_id", "cellid"}
_EXCLUDED_COLOR_PREFIXES = ("center_", "polygon_center_")


def color_options(adata: sc.AnnData) -> list[str]:
    """obs columns suitable for coloring, without coordinates/polygon/cell ids."""
    return [
        col
        for col in adata.obs.columns
        if all(excluded not in str(col).lower() for excluded in _EXCLUDED_COLOR_COLUMNS)
        and not str(col).lower().startswith(_EXCLUDED_COLOR_PREFIXES)
    ]


MODALITIES = ("single_cell", "visium_hd", "merfish")

# Per-modality QC contract. "required" metrics must exist before filtering,
# "optional" ones are used when present, and "unavailable" ones are never
# produced for the modality, so nothing may depend on them.
QC_METRIC_CONTRACT = {
    "single_cell": {
        "required_obs": ("total_counts", "n_genes_by_counts"),
        "optional_obs": ("pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"),
        "unavailable_obs": ("blank_counts", "pct_counts_blank", "is_blank"),
        "required_obsm": (),
    },
    "visium_hd": {
        "required_obs": ("total_counts", "n_genes_by_counts"),
        "optional_obs": ("pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"),
        "unavailable_obs": ("blank_counts", "pct_counts_blank", "is_blank", "volume"),
        # center_x/center_y are derived from this by the loader, not required.
        "required_obsm": ("spatial",),
    },
    "merfish": {
        "required_obs": ("total_counts", "blank_counts", "pct_counts_blank"),
        "optional_obs": ("n_genes_by_counts", "volume"),
        "unavailable_obs": ("pct_counts_mt", "predicted_doublet"),
        "required_obsm": ("spatial",),
    },
}


def qc_metric_availability(adata: sc.AnnData, modality: str) -> dict[str, list[str]]:
    """Report which QC metrics of the modality contract are present or missing."""
    if modality not in QC_METRIC_CONTRACT:
        raise ValueError(
            f"Unknown modality {modality!r}; expected one of {MODALITIES}."
        )
    contract = QC_METRIC_CONTRACT[modality]
    return {
        "missing_required": [
            f"obs['{c}']" for c in contract["required_obs"] if c not in adata.obs.columns
        ] + [
            f"obsm['{k}']" for k in contract["required_obsm"] if k not in adata.obsm
        ],
        "missing_optional": [
            c for c in contract["optional_obs"] if c not in adata.obs.columns
        ],
        "unexpected_unavailable": [
            c for c in contract["unavailable_obs"] if c in adata.obs.columns
        ],
    }


def load_data(filepath: str, 
            cell_columns: list[str] = None, 
            genes: list[str] = None,
            is_spatial: bool = False,
            platform: str = None) -> dict[str, sc.AnnData]:
    """
    Load a dataset from a file and return a dictionary of AnnData objects.

    Parameters:
    - filepath: str
        The path to the file to load.
    - cell_columns: list[str], optional
        A list of column names to use as cell metadata. If None, all columns will be used.
    - genes: list[str], optional
        A list of column names to use as gene metadata. If None, all columns will be used.
    - is_spatial: bool, optional
        Whether the dataset is spatial. If True, the function will attempt to load spatial data.

    Returns:
    - data: dict[str, sc.AnnData]
        A dictionary containing the loaded data.
    """
    
    sdata = None

    if genes is None and cell_columns is None:
        adata = sc.read(filepath)
        if is_spatial:
            if platform is not None and platform.lower() == "merfish":
                validate_merfish_adata(adata)
                _calculate_qc_metrics(adata, is_spatial=is_spatial)
                _add_merfish_qc(adata)
            elif platform is not None and "visium" in platform.lower():
                validate_visium_hd_adata(adata)
                _calculate_qc_metrics(adata, is_spatial=False)
                _add_visium_center_columns(adata)
    # else:
    #     adata = sc.read(filepath, backed='r')
    #     obs_index = slice(None) if cell_columns is None else cell_columns
    #     var_index = slice(None) if genes is None else genes
    #     sdata = adata[obs_index, var_index].to_memory()
    #     if is_spatial:
    #         validate_spatial_adata(sdata, platform=platform)
    #     _calculate_qc_metrics(sdata, is_spatial=is_spatial)
    #     sdata.layers['raw'] = sdata.X.copy()
    #     if platform is not None and platform.lower() == "merfish":
    #         _add_merfish_qc(sdata)
    
    return {'raw': adata, 'subset': sdata}



def _validate_identifiers(names: pd.Index, label: str) -> None:
    """Reject identifiers that cannot safely key observations or genes."""
    if not names.is_unique:
        duplicated = names[names.duplicated()].unique()[:3].tolist()
        raise ValueError(
            f"{label} identifiers must be unique; duplicates include "
            f"{duplicated}. Make them unique explicitly rather than relying "
            "on automatic renaming."
        )

    as_text = names.astype(str).str.strip()
    blank = (as_text == "") | as_text.str.lower().isin({"nan", "none", "<na>"})
    if blank.any():
        raise ValueError(
            f"{label} identifiers must be non-empty; found {int(blank.sum())} "
            "empty or missing identifier(s)."
        )


def _add_visium_center_columns(adata: sc.AnnData) -> None:
    """Expose obsm['spatial'] as center_x/center_y for the plotting helpers.

    Existing columns are kept when they agree with obsm['spatial'] and are
    never overwritten; a disagreement is ambiguous, so it is reported.
    """
    coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    for column, values in (("center_x", coords[:, 0]), ("center_y", coords[:, 1])):
        if column not in adata.obs.columns:
            adata.obs[column] = values
            continue
        existing = pd.to_numeric(adata.obs[column], errors="coerce").to_numpy(dtype=float)
        if not np.allclose(existing, values, equal_nan=False):
            raise ValueError(
                f"adata.obs['{column}'] already exists and differs from "
                "adata.obsm['spatial']. Remove or rename the column, or make "
                "the two coordinate sources consistent."
            )


def validate_visium_hd_adata(
    adata: sc.AnnData,
    *,
    require_spatial_metadata: bool = False,
    require_raw_counts: bool = False,
) -> None:
    """
    Validate an AnnData object for Visium HD spatial analysis.

    Raises ValueError with an actionable message when validation fails.

    Parameters
    ----------
    adata
        AnnData object containing Visium HD data.
    require_spatial_metadata
        If True, require non-empty uns['spatial'] metadata.
    require_raw_counts
        If True, require a non-negative integer-like expression matrix.
        Use this only when adata.X is expected to contain raw counts.
    """

    # 1. Basic AnnData integrity
    if not isinstance(adata, sc.AnnData):
        raise ValueError(
            f"Expected an AnnData object; received {type(adata).__name__}."
        )

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(
            "The dataset has no observations or no genes "
            f"(shape {adata.n_obs} x {adata.n_vars})."
        )

    _validate_identifiers(adata.obs_names, "Observation")
    _validate_identifiers(adata.var_names, "Gene")

    # 2. Spatial coordinates
    if "spatial" not in adata.obsm:
        raise ValueError(
            "Missing adata.obsm['spatial']. Load the Visium HD data "
            "with a spatial-aware reader or provide valid coordinates."
        )

    raw_coords = adata.obsm["spatial"]
    if isinstance(raw_coords, pd.DataFrame):
        # AnnData already guarantees the index matches adata.obs_names.
        raw_coords = raw_coords.to_numpy()
    if sparse.issparse(raw_coords):
        raise ValueError("adata.obsm['spatial'] must be a dense numeric array.")

    coords = np.asarray(raw_coords)

    if coords.shape != (adata.n_obs, 2):
        raise ValueError(
            "Expected spatial coordinates with shape "
            f"({adata.n_obs}, 2); received {coords.shape}."
        )

    if coords.dtype == bool:
        raise ValueError("Spatial coordinates must be numeric, not boolean.")

    try:
        coords = coords.astype(np.float64, copy=False)
    except (TypeError, ValueError) as exc:
        raise ValueError("Spatial coordinates must be numeric.") from exc

    if not np.isfinite(coords).all():
        n_bad = int((~np.isfinite(coords).all(axis=1)).sum())
        raise ValueError(
            f"Spatial coordinates contain non-finite values in {n_bad} "
            "observation(s)."
        )

    if adata.n_obs > 1 and (coords == coords[0]).all():
        raise ValueError(
            "All spatial coordinates are identical, so the dataset has no "
            "spatial extent. Check that the coordinates were exported correctly."
        )

    # 3. Expression matrix
    X = adata.X

    if X is None:
        raise ValueError("The AnnData object has no expression matrix.")

    if not (sparse.issparse(X) or isinstance(X, np.ndarray)):
        raise ValueError(
            "The expression matrix must be a NumPy array or SciPy sparse "
            f"matrix; received {type(X).__name__}. Backed or lazy matrices "
            "are not supported here."
        )

    if X.ndim != 2 or tuple(X.shape) != (adata.n_obs, adata.n_vars):
        raise ValueError(
            "The expression matrix shape "
            f"{tuple(X.shape)} does not match the dataset dimensions "
            f"({adata.n_obs}, {adata.n_vars})."
        )

    values = X.data if sparse.issparse(X) else np.asarray(X)

    if values.dtype == bool or not np.issubdtype(values.dtype, np.number):
        raise ValueError("The expression matrix must contain numeric values.")

    if not np.isfinite(values).all():
        raise ValueError("The expression matrix contains non-finite values.")

    if np.any(values < 0):
        raise ValueError(
            "The expression matrix contains negative values, which are "
            "invalid for QC and filtering. Check whether it holds scaled or "
            "transformed expression instead of counts."
        )

    # Stored zeros in a sparse matrix are legitimate; only a matrix with no
    # signal at all is unusable.
    if not np.any(values):
        raise ValueError("The expression matrix is all zeros.")

    if require_raw_counts:
        if np.any(values < 0):
            raise ValueError(
                "Raw counts cannot be negative. Check whether the matrix "
                "contains normalized or transformed expression values."
            )

        if not np.allclose(values, np.round(values)):
            raise ValueError(
                "Raw counts should be integer-like, but non-integer values "
                "were found. Check which matrix or layer contains raw counts."
            )

    # 4. Visium spatial metadata
    spatial_metadata = adata.uns.get("spatial")

    if require_spatial_metadata:
        if not isinstance(spatial_metadata, dict) or not spatial_metadata:
            raise ValueError(
                "Missing or empty adata.uns['spatial']. "
                "Image overlays may not be available."
            )

    # 5. Optional image metadata sanity check
    if isinstance(spatial_metadata, dict):
        for library_id, library in spatial_metadata.items():
            if not isinstance(library, dict):
                raise ValueError(
                    f"Spatial metadata for '{library_id}' is malformed."
                )

            # Images and scale factors can be absent in some valid objects.
            images = library.get("images", {})
            if images is not None and not isinstance(images, dict):
                raise ValueError(
                    f"Image metadata for '{library_id}' must be a dictionary."
                )

            scalefactors = library.get("scalefactors", {})
            if scalefactors is not None and not isinstance(
                scalefactors, dict
            ):
                raise ValueError(
                    f"Scale-factor metadata for '{library_id}' must "
                    "be a dictionary."
                )

            # Optional numeric scale factors must be usable when supplied.
            for factor in (
                "spot_diameter_fullres",
                "tissue_hires_scalef",
                "tissue_lowres_scalef",
                "bin_size_um",
            ):
                value = (scalefactors or {}).get(factor)
                if value is None:
                    continue
                if (
                    isinstance(value, bool)
                    or not isinstance(value, (int, float, np.number))
                    or not np.isfinite(value)
                    or value <= 0
                ):
                    raise ValueError(
                        f"Scale factor '{factor}' for '{library_id}' must be "
                        f"a positive finite number; received {value!r}."
                    )

def validate_merfish_adata(adata: sc.AnnData) -> None:
    """Raise a ValueError with an actionable message if a MERFISH dataset is unusable."""

    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("The dataset is empty (no cells or no genes).")

    missing = [c for c in ("center_x", "center_y") if c not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"Spatial datasets require {', '.join(missing)} in adata.obs "
            "(cell centre coordinates). Add them or select the single-cell technology."
        )

    for column in ("center_x", "center_y"):
        if not pd.api.types.is_numeric_dtype(adata.obs[column]):
            raise ValueError(f"adata.obs['{column}'] must be numeric.")
        if not np.isfinite(adata.obs[column].to_numpy(dtype=float)).all():
            raise ValueError(
                f"adata.obs['{column}'] contains missing or non-finite values."
            )

    if not adata.var_names.str.startswith("Blank-").any():
        raise ValueError(
            "No blank genes found. MERFISH blank probes must be named with the "
            "'Blank-' prefix."
        )

    if np.asarray(adata.X.sum(axis=1)).ravel().max() <= 0:
        raise ValueError("All cells have zero total counts.")


def _add_merfish_qc(adata: sc.AnnData) -> None:
    """Add blank-count QC columns and the spatial embedding to a validated dataset."""

    blank_mask = adata.var_names.str.startswith("Blank-")
    blank_counts = np.asarray(adata[:, blank_mask].X.sum(axis=1)).ravel()
    total = adata.obs["total_counts"].to_numpy(dtype=float)
    adata.obs["blank_counts"] = blank_counts
    # Cells without counts have no defined blank fraction; report 0 rather than NaN.
    adata.obs["pct_counts_blank"] = np.divide(
        blank_counts * 100.0,
        total,
        out=np.zeros_like(total),
        where=total > 0,
    )
    adata.obsm["spatial"] = adata.obs[["center_x", "center_y"]].to_numpy(dtype=float)

# def load_data(
#     filepath: str,
#     cell_columns: list[str] | None = None,
#     genes: list[str] | None = None,
#     is_spatial: bool = False,
# ) -> dict[str, sc.AnnData | None]:
#     """
#     Load an AnnData dataset and optionally create an in-memory subset.

#     Returns
#     -------
#     dict
#         {
#             "raw": original AnnData,
#             "subset": selected AnnData or None
#         }
#     """

#     # Use backed mode for the original dataset
#     adata = sc.read(
#         filepath,
#         backed="r",
#     )

#     subset = None

#     # Create subset only when requested
#     if cell_columns is not None or genes is not None:

#         obs_indices = cell_columns if cell_columns is not None else slice(None)
#         var_indices = genes if genes is not None else slice(None)

#         subset = adata[obs_indices, var_indices].copy()

#         _prepare_adata(
#             subset,
#             is_spatial=is_spatial,
#         )

#     else:
#         # If no subset was requested, prepare the loaded dataset
#         # only if you actually want to modify the raw object.
#         _prepare_adata(
#             adata,
#             is_spatial=is_spatial,
#         )

#     return {
#         "raw": adata,
#         "subset": subset,
#     }


def _calculate_qc_metrics(adata, is_spatial=False):
    """
    Calculate QC metrics for an AnnData object.

    For spatial/MERFISH data, calculate the standard cell/spot
    metrics without top-N gene percentages.

    For single-cell data, calculate additional QC metrics based
    on available gene annotations.
    """

    if is_spatial:
        sc.pp.calculate_qc_metrics(
            adata,
            percent_top=None,
            inplace=True,
            log1p=True,
        )
        return

    # --------------------------------------------------------
    # Single-cell QC
    # --------------------------------------------------------

    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "Mt-"))
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=None,
        inplace=True,
        log1p=True,
    )


def write_to_disk(adata: sc.AnnData, filepath: str) -> None:
    """
    Write an AnnData object to disk.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to write to disk.
    - filepath: str
        The path to the file to write.

    Returns:
    - None
        The function writes the AnnData object to disk.
    """
    adata.write(filepath)


def calculate_module_score(
    adata: sc.AnnData,
    genes: list[str],
    score_name: str = "module_score",
):
    """
    Calculate a module score for a set of genes in an AnnData object.

    Parameters:
    - adata: sc.AnnData
        The AnnData object containing the gene expression data.
    - genes: list[str]
        A list of gene names to include in the module score.
    - score_name: str, optional (default="module_score")
        The name of the column in adata.obs to store the module score.

    Returns:
    - adata: sc.AnnData
        The AnnData object with the calculated module score added to adata.obs.
    """
    genes = [gene for gene in genes if gene in adata.var_names]

    if not genes:
        raise ValueError("None of the selected genes were found in the dataset.")

    sc.tl.score_genes(
        adata,
        gene_list=genes,
        score_name=score_name,
    )

    return adata


def get_spatially_variable_genes(
    adata: sc.AnnData,
    n_top_genes: int = 2000,
) -> list[str]:
    """
    Identify spatially variable genes using Moran's I in Squidpy.

    Parameters
    ----------
    adata
        AnnData object with spatial coordinates in `obsm["spatial"]`.
    n_top_genes
        Maximum number of top spatially autocorrelated genes to return.

    Returns
    -------
    list[str]
        Gene names ranked by Moran's I, excluding blank probes.

    Notes
    -----
    Uses the existing `spatial_connectivities` graph if available.
    Otherwise, constructs a Delaunay spatial graph.

    This function returns genes ranked by Moran's I, not necessarily
    statistically significant genes. For inferential use, examine the
    adjusted p-values in the Squidpy results.
    """
    if n_top_genes < 1:
        raise ValueError("n_top_genes must be at least 1.")

    if adata.n_obs < 3:
        raise ValueError("At least 3 observations are required.")

    if "spatial" not in adata.obsm:
        raise ValueError(
            "Spatial coordinates not found in adata.obsm['spatial']."
        )

    if adata.n_vars == 0:
        return []

    # Build the spatial graph if it doesn't already exist.
    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors_delaunay(
            adata,
            spatial_key="spatial",
        )

    # Exclude blank probes before testing.
    genes = [
        gene
        for gene in adata.var_names
        if "Blank-" not in str(gene)
    ]

    if not genes:
        return []

    # Calculate Moran's I and return the results without storing
    # the result table in adata.uns.
    results = sq.gr.spatial_autocorr(
        adata,
        mode="moran",
        genes=genes,
        attr="X",
        corr_method="fdr_bh",
        copy=True,
        n_jobs=1,
    )

    if results is None or results.empty:
        return []

    # Rank positive spatial autocorrelation first.
    results = results.sort_values("I", ascending=False)
    results = results.loc[results["I"] > 0]

    return results.index[:n_top_genes].astype(str).tolist()


# def _prepare_adata(
#     adata: sc.AnnData,
#     is_spatial: bool = False,
# ) -> sc.AnnData:

#     _calculate_qc_metrics(
#         adata,
#         is_spatial=is_spatial,
#     )

#     adata.layers["raw"] = adata.X.copy()

#     if is_spatial:

#         # Blank genes
#         blank_genes = adata.var_names[
#             adata.var_names.str.startswith("Blank-")
#         ].tolist()

#         if blank_genes:
#             blank_counts = np.asarray(
#                 adata[:, blank_genes].X.sum(axis=1)
#             ).ravel()

#             adata.obs["blank_counts"] = blank_counts

#             adata.obs["pct_counts_blank"] = (
#                 adata.obs["blank_counts"]
#                 / adata.obs["total_counts"]
#             ) * 100
#         else:
#             adata.obs["blank_counts"] = 0
#             adata.obs["pct_counts_blank"] = 0.0

#         # Spatial coordinates
#         if {"center_x", "center_y"}.issubset(adata.obs.columns):
#             adata.obsm["spatial"] = adata.obs[
#                 ["center_x", "center_y"]
#             ].to_numpy()

#     return adata
