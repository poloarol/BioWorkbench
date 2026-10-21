import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from src.filtering import (
    THRESHOLD_RANGES,
    data_filtering,
    preview_filter_impact,
    validate_filter_thresholds,
)


def make_adata(sparse_x=True):
    # Genes detected per observation: 0, 2, 3, 3; mito percent: 0, 10, 20, 50.
    X = np.array(
        [
            [0, 0, 0, 0],
            [1, 1, 0, 0],
            [1, 1, 1, 0],
            [1, 1, 1, 0],
        ],
        dtype=np.float32,
    )
    adata = ad.AnnData(
        sparse.csr_matrix(X) if sparse_x else X,
        obs=pd.DataFrame(
            {"pct_counts_mt": [0.0, 10.0, 20.0, 50.0]},
            index=[f"bin{i}" for i in range(4)],
        ),
        var=pd.DataFrame(index=[f"g{i}" for i in range(4)]),
    )
    adata.obsm["spatial"] = np.arange(8, dtype=float).reshape(4, 2)
    adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes_by_counts"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    return adata


@pytest.mark.parametrize("sparse_x", [True, False])
@pytest.mark.parametrize(
    "min_genes, expected_removed",
    [(0, 0), (2, 1), (3, 2), (4, 4)],
)
def test_min_genes_boundary_is_inclusive(min_genes, expected_removed, sparse_x):
    impact = preview_filter_impact(make_adata(sparse_x), min_genes, 0)
    assert impact["obs_removed_min_genes"] == expected_removed


def test_min_genes_preview_matches_real_filter(monkeypatch):
    monkeypatch.setattr(
        "src.filtering._identify_doublets", lambda adata, rate: adata
    )
    adata = make_adata()
    for min_genes in (0, 2, 3):
        impact = preview_filter_impact(adata, min_genes, 0, 20.0)
        result = data_filtering(
            adata, min_genes=min_genes, min_cells=0, max_mt_percentage=20.0,
            modality="single_cell",
        )
        assert result.n_obs == impact["obs_remaining"]


def test_min_cells_counts_only_observations_kept_by_min_genes():
    adata = make_adata()
    # min_genes=3 keeps two bins, where g0-g2 are each detected twice and g3 never.
    assert preview_filter_impact(adata, 3, 2)["vars_removed_min_cells"] == 1
    assert preview_filter_impact(adata, 3, 3)["vars_removed_min_cells"] == 4


def test_mito_boundary_is_strict():
    adata = make_adata()
    impact = preview_filter_impact(adata, 0, 0, max_mt_percentage=20.0)
    assert impact["obs_removed_mito"] == 2  # exactly 20 and 50 are removed
    impact = preview_filter_impact(adata, 0, 0, max_mt_percentage=20.1)
    assert impact["obs_removed_mito"] == 1


def test_mito_not_previewed_when_unavailable_or_disabled():
    adata = make_adata()
    assert preview_filter_impact(adata, 0, 0)["obs_removed_mito"] is None
    del adata.obs["pct_counts_mt"]
    assert preview_filter_impact(adata, 0, 0, 20.0)["obs_removed_mito"] is None


def test_preview_does_not_modify_input():
    adata = make_adata()
    before = adata.X.copy()
    preview_filter_impact(adata, 3, 3, 10.0)
    assert adata.shape == (4, 4)
    assert (adata.X != before).nnz == 0


def test_empty_result_is_reported_not_raised():
    impact = preview_filter_impact(make_adata(), 10, 1)
    assert impact["obs_remaining"] == 0
    assert impact["vars_remaining"] == 0


def test_data_filtering_empty_after_min_genes():
    result = data_filtering(
        make_adata(), min_genes=10, min_cells=0, modality="visium_hd",
    )
    assert result.n_obs == 0
    report = result.uns["filtering_report"]
    assert report.loc[report["step"] == "min_genes", "n_obs_after"].iloc[0] == 0
    assert report.loc[report["step"] == "doublets", "action"].iloc[0] == "skipped"


def test_validate_accepts_range_boundaries():
    for name, (low, high) in THRESHOLD_RANGES.items():
        assert validate_filter_thresholds(**{name: low}) == []
        assert validate_filter_thresholds(**{name: high}) == []


@pytest.mark.parametrize(
    "name, value",
    [
        ("min_genes", -1),
        ("min_cells", 1_000_001),
        ("max_mt_percentage", 100.01),
        ("max_mt_percentage", -0.1),
        ("expected_doublet_rate", 1.01),
        ("blank_max_percentage", float("nan")),
        ("max_mt_percentage", float("inf")),
        ("min_genes", "abc"),
        ("min_genes", None),
    ],
)
def test_validate_rejects_invalid_values(name, value):
    assert len(validate_filter_thresholds(**{name: value})) == 1


def test_validate_rejects_unknown_threshold():
    assert validate_filter_thresholds(bogus=1)
