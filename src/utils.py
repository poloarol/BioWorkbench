
import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd

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


def load_data(filepath: str, 
            cell_columns: list[str] = None, 
            genes: list[str] = None,
            is_spatial: bool = False) -> dict[str, sc.AnnData]:
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
            validate_spatial_adata(adata)
        _calculate_qc_metrics(adata, is_spatial=is_spatial)
        adata.layers['raw'] = adata.X.copy()
        
        if is_spatial:
            _add_spatial_qc(adata)
    else:
        adata = sc.read(filepath, backed='r')
        obs_index = slice(None) if cell_columns is None else cell_columns
        var_index = slice(None) if genes is None else genes
        sdata = adata[obs_index, var_index].to_memory()
        if is_spatial:
            validate_spatial_adata(sdata)
        _calculate_qc_metrics(sdata, is_spatial=is_spatial)
        sdata.layers['raw'] = sdata.X.copy()
        if is_spatial:
            _add_spatial_qc(sdata)
    
    return {'raw': adata, 'subset': sdata}


def validate_spatial_adata(adata: sc.AnnData) -> None:
    """Raise a ValueError with an actionable message if a spatial dataset is unusable."""

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


def _add_spatial_qc(adata: sc.AnnData) -> None:
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
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
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


def get_spatially_variable_genes(adata: sc.AnnData, n_top_genes: int = 2000) -> list[str]:
    """
    Identify spatially variable genes in an AnnData object.

    Parameters:
    - adata: sc.AnnData
        The AnnData object containing the gene expression data.
    - n_top_genes: int, optional (default=2000)
        The number of top spatially variable genes to return.

    Returns:
    - list[str]
        A list of the top spatially variable genes.
    """
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=n_top_genes,
        subset=False,
        inplace=True,
    )
    
    items = adata.var_names[adata.var["highly_variable"]].tolist()
    items = [x for x in items if "Blank-" not in x]

    
    return items


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
