import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from src.utils import MODALITIES, qc_metric_availability

REPORT_KEY = "filtering_report"
_REPORT_COLUMNS = [
    "step", "action", "n_obs_before", "n_obs_after",
    "n_vars_before", "n_vars_after", "n_flagged", "note",
]


def _record_step(report, step, action, before, after, n_flagged=0, note=""):
    """Append counts for one step; `action` is removed, flagged or skipped."""
    report.append({
        "step": step,
        "action": action,
        "n_obs_before": int(before[0]),
        "n_obs_after": int(after.n_obs),
        "n_vars_before": int(before[1]),
        "n_vars_after": int(after.n_vars),
        "n_flagged": int(n_flagged),
        "note": note,
    })


def _count_true(adata, column):
    if column not in adata.obs.columns:
        return 0
    return int(adata.obs[column].fillna(False).astype(bool).sum())


def _filter_cells(adata: sc.AnnData, min_genes: int = 200) -> sc.AnnData:
    """
    Filter cells in an AnnData object based on quality control metrics.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to filter.
    - min_genes: int, optional
        Minimum number of genes expressed in a cell to keep it. Default is 200.

    Returns:
    - filtered_adata: sc.AnnData
        The filtered AnnData object.
    """
    
    # Filter cells based on the specified criteria
    cell_subset, number_per_cell = sc.pp.filter_cells(
        adata,
        min_genes=min_genes,
        inplace=False
    )

    return adata[cell_subset, :].copy()

def _filter_genes(adata: sc.AnnData, min_cells: int = 3) -> sc.AnnData:
    """
    Filter genes in an AnnData object based on quality control metrics.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to filter.
    - min_cells: int, optional
        Minimum number of cells expressing a gene to keep it. Default is 3.

    Returns:
    - filtered_adata: sc.AnnData
        The filtered AnnData object.
    """


    gene_subset, number_per_gene = sc.pp.filter_genes(
        adata,
        min_cells=min_cells,
        inplace=False
    )

    return adata[:, gene_subset].copy()


def _filter_mitochondrial_genes(adata: sc.AnnData, max_mt_percentage: float = 20) -> sc.AnnData:
    """
    Filter cells in an AnnData object based on the fraction of mitochondrial genes.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to filter.
    - max_mt_percentage: float, optional
        Maximum percentage of mitochondrial genes allowed in a cell. Default is 20%.

    Returns:
    - filtered_adata: sc.AnnData
        The filtered AnnData object.
    """
    
    # Filter cells based on the specified criteria
    filtered_adata = adata[adata.obs["pct_counts_mt"] < max_mt_percentage].copy()
    
    return filtered_adata


def _identify_doublets(adata: sc.AnnData, expected_doublet_rate: float = 0.05) -> sc.AnnData:
    """
    Identify potential doublets in an AnnData object using Scrublet.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to analyze.
    - expected_doublet_rate: float, optional
        Expected doublet rate in the dataset. Default is 0.1 (10%).

    Returns:
    - adata_with_doublets: sc.AnnData
        The AnnData object with an additional column indicating potential doublets.
    """
    
    sc.pp.scrublet(adata, expected_doublet_rate=expected_doublet_rate)
    
    return adata

def _identify_blanks(adata: sc.AnnData, max_blank_percentage: float = 5) -> sc.AnnData:
    """
    Identify potential blank droplets in an AnnData object based on the fraction of blank counts.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to analyze.
    - max_blank_percentage: float, optional
        Maximum percentage of blank counts allowed in a cell. Default is 5%.

    Returns:
    - filtered_adata: sc.AnnData
        The filtered AnnData object.
    """

    # Filter cells based on the specified criteria
    adata.obs['is_blank'] = adata.obs["pct_counts_blank"] > max_blank_percentage
    return adata


THRESHOLD_RANGES = {
    "min_genes": (0, 1_000_000),
    "min_cells": (0, 1_000_000),
    "max_mt_percentage": (0.0, 100.0),
    "expected_doublet_rate": (0.0, 1.0),
    "blank_max_percentage": (0.0, 100.0),
    "negative_probe_max_percentage": (0.0, 100.0),
    "negative_codeword_max_percentage": (0.0, 100.0),
}


def validate_filter_thresholds(**thresholds):
    """Return a list of error messages for thresholds that are missing, non-finite or out of range."""
    errors = []
    for name, value in thresholds.items():
        if name not in THRESHOLD_RANGES:
            errors.append(f"Unknown threshold {name!r}.")
            continue
        low, high = THRESHOLD_RANGES[name]
        try:
            number = float(value)
        except (TypeError, ValueError):
            errors.append(f"{name} must be a number, got {value!r}.")
            continue
        if number != number or number in (float("inf"), float("-inf")):
            errors.append(f"{name} must be finite, got {value!r}.")
        elif not low <= number <= high:
            errors.append(f"{name} must be between {low} and {high}, got {value!r}.")
    return errors


def preview_filter_impact(
    adata,
    min_genes,
    min_cells,
    max_mt_percentage=None,
    negative_probe_max_percentage=None,
    negative_codeword_max_percentage=None,
):
    """
    Count what `data_filtering` would remove, without modifying `adata`.

    Mirrors the order and boundaries of the real steps: observations are kept
    when n_genes >= min_genes, genes are kept when detected in >= min_cells of the
    remaining observations, and observations are kept when pct_counts_mt is
    strictly below max_mt_percentage. `max_mt_percentage=None`, or a missing
    obs['pct_counts_mt'], skips the mitochondrial step (reported as None).
    """
    from scipy import sparse

    X = adata.X
    detected = (X > 0).tocsr() if sparse.issparse(X) else (np.asarray(X) > 0)
    n_genes = np.asarray(detected.sum(axis=1)).ravel()
    keep_obs = n_genes >= min_genes
    detected_kept = detected[keep_obs]
    per_gene = np.asarray(detected_kept.sum(axis=0)).ravel()
    keep_vars = per_gene >= min_cells

    mt_removed = None
    if max_mt_percentage is not None and "pct_counts_mt" in adata.obs.columns:
        mt = adata.obs["pct_counts_mt"].to_numpy(dtype=float)[keep_obs]
        mt_removed = int((~(mt < max_mt_percentage)).sum())

    controls_flagged = None
    if (
        negative_probe_max_percentage is not None
        and negative_codeword_max_percentage is not None
        and {"negative_probe_pct", "negative_codeword_pct"}.issubset(
            adata.obs.columns
        )
    ):
        probe_pct = adata.obs["negative_probe_pct"].to_numpy(dtype=float)[keep_obs]
        codeword_pct = adata.obs["negative_codeword_pct"].to_numpy(dtype=float)[keep_obs]
        controls_flagged = int(
            (
                (probe_pct > negative_probe_max_percentage)
                | (codeword_pct > negative_codeword_max_percentage)
            ).sum()
        )

    n_after_cells = int(keep_obs.sum())
    return {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "obs_removed_min_genes": int(adata.n_obs - n_after_cells),
        "vars_removed_min_cells": int(adata.n_vars - keep_vars.sum()),
        "obs_removed_mito": mt_removed,
        "obs_flagged_controls": controls_flagged,
        "obs_remaining": n_after_cells - (mt_removed or 0),
        "vars_remaining": int(keep_vars.sum()),
    }


def data_filtering(adata, 
           min_genes: int = 200, 
           min_cells: int = 3, 
           max_mt_percentage: float = 20, 
           expected_doublet_rate: float = 0.05,
           blank_max_percentage: float = 5,
           is_spatial: bool = False,
           modality: str = None,
           negative_probe_max_percentage: float = 5,
           negative_codeword_max_percentage: float = 5) -> sc.AnnData:
    """
    Filter an AnnData object based on quality control metrics.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to filter.
    - min_genes: int, optional
        Minimum number of genes expressed in a cell to keep it. Default is 200.
    - min_cells: int, optional
        Minimum number of cells expressing a gene to keep it. Default is 3.
    - max_mt_percentage: float, optional
        Maximum percentage of mitochondrial genes allowed in a cell. Default is 20%.
    - expected_doublet_rate: float, optional
        Expected doublet rate in the dataset. Default is 0.05 (5%).
    - is_spatial: bool, optional
        Legacy flag: True selects the MERFISH (blank-probe) path. Ignored when
        `modality` is given.
    - modality: str, optional
        One of "single_cell", "visium_hd", "merfish" or "xenium". Visium HD
        uses the single-cell QC path. Xenium flags observations exceeding either
        negative-control percentage threshold without removing them.
    - negative_probe_max_percentage: float, optional
        Maximum negative-probe percentage before a Xenium observation is flagged.
    - negative_codeword_max_percentage: float, optional
        Maximum negative-codeword percentage before a Xenium observation is flagged.
        Per-step counts are stored in `adata.uns["filtering_report"]`.
    - blank_max_percentage: float, optional
        Maximum percentage of blank counts allowed in a cell for spatial datasets. Default is 5%.

    Returns:
    - filtered_adata: sc.AnnData
        The filtered AnnData object.
    """
    
    # Filter cells based on the specified criteria
    if modality is None:
        modality = "merfish" if is_spatial else "single_cell"
    if modality not in MODALITIES:
        raise ValueError(f"Unknown modality {modality!r}; expected one of {MODALITIES}.")

    if modality == "xenium":
        threshold_errors = validate_filter_thresholds(
            negative_probe_max_percentage=negative_probe_max_percentage,
            negative_codeword_max_percentage=negative_codeword_max_percentage,
        )
        if threshold_errors:
            raise ValueError("Invalid Xenium filtering thresholds: " + " ".join(threshold_errors))

    if modality == "visium_hd":
        availability = qc_metric_availability(adata, modality)
        if availability["missing_required"]:
            raise ValueError(
                "Visium HD filtering requires "
                + ", ".join(availability["missing_required"])
                + ". Load the dataset with the Visium HD technology selected."
            )
    elif modality == "xenium":
        availability = qc_metric_availability(adata, modality)
        if availability["missing_required"]:
            raise ValueError(
                "Xenium filtering requires "
                + ", ".join(availability["missing_required"])
                + ". Load the dataset with the Xenium technology selected."
            )
        
        print(availability)

    report = []

    before = adata.shape
    adata = _filter_cells(adata, min_genes)
    _record_step(report, "min_genes", "removed", before, adata)

    before = adata.shape
    adata = _filter_genes(adata, min_cells)
    _record_step(report, "min_cells", "removed", before, adata)

    if modality == "merfish":
        before = adata.shape
        adata = _identify_blanks(adata, max_blank_percentage=blank_max_percentage)
        _record_step(
            report, "blank_probes", "flagged", before, adata,
            n_flagged=_count_true(adata, "is_blank"),
            note="Flagged in obs['is_blank']; observations are not removed.",
        )
    elif modality == "xenium":
        before = adata.shape
        high_negative_probes = (
            adata.obs["negative_probe_pct"] > negative_probe_max_percentage
        ).fillna(False)
        high_negative_codewords = (
            adata.obs["negative_codeword_pct"] > negative_codeword_max_percentage
        ).fillna(False)
        adata.obs["is_xenium_control"] = (
            high_negative_probes | high_negative_codewords
        )
        _record_step(
            report,
            "negative_controls",
            "flagged",
            before,
            adata,
            n_flagged=_count_true(adata, "is_xenium_control"),
            note=(
                "Flagged in obs['is_xenium_control'] when negative-probe "
                f"percentage exceeds {negative_probe_max_percentage}% or "
                "negative-codeword percentage exceeds "
                f"{negative_codeword_max_percentage}%; observations are not removed."
            ),
        )
    else:
        before = adata.shape
        if "pct_counts_mt" not in adata.obs.columns:
            _record_step(
                report, "max_mt_percentage", "skipped", before, adata,
                note="obs['pct_counts_mt'] is unavailable; no mitochondrial filtering was applied.",
            )
        elif "mt" in adata.var.columns and not adata.var["mt"].any():
            _record_step(
                report, "max_mt_percentage", "skipped", before, adata,
                note="No mitochondrial genes (prefix 'MT-') were annotated; no mitochondrial filtering was applied.",
            )
        else:
            adata = _filter_mitochondrial_genes(adata, max_mt_percentage)
            _record_step(report, "max_mt_percentage", "removed", before, adata)

        before = adata.shape
        # Scrublet requires enough cells and genes to build a model.
        if adata.n_obs > 1 and adata.n_vars > 0:
            adata = _identify_doublets(adata, expected_doublet_rate)
            _record_step(
                report, "doublets", "flagged", before, adata,
                n_flagged=_count_true(adata, "predicted_doublet"),
                note="Flagged in obs['predicted_doublet']; observations are not removed.",
            )
        else:
            _record_step(
                report, "doublets", "skipped", before, adata,
                note="Too few observations or genes for doublet detection.",
            )

    adata.uns[REPORT_KEY] = pd.DataFrame(report, columns=_REPORT_COLUMNS)

    return adata