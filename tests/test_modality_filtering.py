import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from unittest.mock import patch

from src.filtering import REPORT_KEY, data_filtering
from src.utils import QC_METRIC_CONTRACT, _calculate_qc_metrics, qc_metric_availability


def _counts(n_obs=40, n_vars=10, seed=0):
    rng = np.random.default_rng(seed)
    return rng.poisson(3.0, size=(n_obs, n_vars)).astype(np.int32) + 1


def make_visium(counts=None, genes=None, with_center=False):
    counts = _counts() if counts is None else counts
    genes = genes or ["MT-CO1", "MT-ND1"] + [f"G{i}" for i in range(counts.shape[1] - 2)]
    adata = ad.AnnData(
        sparse.csr_matrix(counts),
        obs=pd.DataFrame(index=[f"s_{i}-1" for i in range(counts.shape[0])]),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm["spatial"] = np.column_stack(
        [np.arange(counts.shape[0]) * 8.0, np.zeros(counts.shape[0])]
    )
    _calculate_qc_metrics(adata, is_spatial=False)
    if with_center:
        adata.obs["center_x"] = adata.obsm["spatial"][:, 0]
        adata.obs["center_y"] = adata.obsm["spatial"][:, 1]
    return adata


def make_merfish(counts=None):
    counts = _counts(n_vars=8) if counts is None else counts
    genes = [f"G{i}" for i in range(counts.shape[1] - 2)] + ["Blank-1", "Blank-2"]
    adata = ad.AnnData(
        counts.astype(np.float32),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(counts.shape[0])]),
        var=pd.DataFrame(index=genes),
    )
    _calculate_qc_metrics(adata, is_spatial=True)
    blank = counts[:, -2:].sum(axis=1)
    adata.obs["blank_counts"] = blank
    adata.obs["pct_counts_blank"] = blank * 100.0 / adata.obs["total_counts"].to_numpy()
    adata.obsm["spatial"] = np.random.default_rng(0).uniform(0, 100, (counts.shape[0], 2))
    return adata


def _fake_doublets(adata, rate):
    adata.obs["predicted_doublet"] = np.arange(adata.n_obs) % 5 == 0
    return adata


def run(adata, **kwargs):
    # Scrublet is not under test and needs realistically sized data.
    with patch("src.filtering._identify_doublets", side_effect=_fake_doublets):
        return data_filtering(adata, **kwargs)


# --------------------------------------------------------------------------
# Visium HD
# --------------------------------------------------------------------------


def test_visium_filtering_does_not_need_merfish_or_center_fields():
    adata = make_visium(with_center=False)
    assert "center_x" not in adata.obs and "pct_counts_blank" not in adata.obs

    result = run(adata, modality="visium_hd", min_genes=1, min_cells=1, max_mt_percentage=100)

    assert result.n_obs == adata.n_obs
    assert "is_blank" not in result.obs and "pct_counts_blank" not in result.obs
    assert "center_x" not in result.obs
    np.testing.assert_array_equal(result.obsm["spatial"], adata.obsm["spatial"])


def test_visium_requires_spatial_coordinates():
    adata = make_visium()
    del adata.obsm["spatial"]
    with pytest.raises(ValueError, match="obsm\\['spatial'\\]"):
        data_filtering(adata, modality="visium_hd")


def test_visium_requires_core_qc_metrics():
    adata = make_visium()
    del adata.obs["total_counts"]
    with pytest.raises(ValueError, match="total_counts"):
        data_filtering(adata, modality="visium_hd")


def test_visium_report_has_counts_for_every_step():
    counts = _counts()
    counts[0] = 0
    counts[1] = [100, 100] + [0] * 8
    adata = make_visium(counts=counts)
    mt_cutoff = 90.0

    result = run(adata, modality="visium_hd", min_genes=1, min_cells=1, max_mt_percentage=mt_cutoff)
    report = result.uns[REPORT_KEY].set_index("step")

    assert list(report.index) == ["min_genes", "min_cells", "max_mt_percentage", "doublets"]
    assert report.loc["min_genes", ["n_obs_before", "n_obs_after"]].tolist() == [40, 39]
    assert report.loc["max_mt_percentage", ["n_obs_before", "n_obs_after"]].tolist() == [39, 38]
    assert report.loc["min_genes", "n_vars_before"] == 10
    for step in report.index:
        assert report.loc[step, "n_obs_after"] <= report.loc[step, "n_obs_before"]
    assert report["n_obs_before"].iloc[0] == adata.n_obs
    assert report["n_obs_after"].iloc[-1] == result.n_obs
    assert report["n_vars_after"].iloc[-1] == result.n_vars


def test_visium_doublets_are_flagged_not_removed():
    adata = make_visium()
    result = run(adata, modality="visium_hd", min_genes=1, min_cells=1, max_mt_percentage=100)
    report = result.uns[REPORT_KEY].set_index("step")

    assert result.n_obs == adata.n_obs
    assert report.loc["doublets", "action"] == "flagged"
    assert report.loc["doublets", "n_obs_before"] == report.loc["doublets", "n_obs_after"]
    assert report.loc["doublets", "n_flagged"] == int(result.obs["predicted_doublet"].sum()) > 0
    assert report.loc["min_genes", "action"] == "removed"


def test_visium_missing_mitochondrial_column_is_skipped_explicitly():
    adata = make_visium()
    del adata.obs["pct_counts_mt"]

    result = run(adata, modality="visium_hd", min_genes=1, min_cells=1)
    row = result.uns[REPORT_KEY].set_index("step").loc["max_mt_percentage"]

    assert result.n_obs == adata.n_obs
    assert row["action"] == "skipped"
    assert "pct_counts_mt" in row["note"]


def test_visium_without_mitochondrial_genes_is_skipped_explicitly():
    adata = make_visium(genes=[f"G{i}" for i in range(10)])
    assert not adata.var["mt"].any()

    result = run(adata, modality="visium_hd", min_genes=1, min_cells=1)
    row = result.uns[REPORT_KEY].set_index("step").loc["max_mt_percentage"]

    assert row["action"] == "skipped"
    assert "mitochondrial genes" in row["note"]
    assert result.n_obs == adata.n_obs


def test_visium_filtering_is_sparse_and_preserves_identity():
    counts = _counts()
    counts[3] = 0
    adata = make_visium(counts=counts)
    original = adata.X.toarray().copy()

    result = run(adata, modality="visium_hd", min_genes=1, min_cells=1, max_mt_percentage=100)

    assert sparse.issparse(result.X)
    assert adata.obs_names[3] not in result.obs_names
    positions = adata.obs_names.get_indexer(result.obs_names)
    np.testing.assert_array_equal(result.obsm["spatial"], adata.obsm["spatial"][positions])
    np.testing.assert_array_equal(adata.X.toarray(), original)


def test_visium_empty_result_is_reported():
    result = run(make_visium(), modality="visium_hd", min_genes=10_000, min_cells=1)
    report = result.uns[REPORT_KEY].set_index("step")

    assert result.n_obs == 0
    assert report.loc["doublets", "action"] == "skipped"
    assert report.loc["min_genes", "n_obs_after"] == 0


def test_filtering_report_survives_h5ad_round_trip(tmp_path):
    result = run(make_visium(), modality="visium_hd", min_genes=1, min_cells=1)
    path = tmp_path / "filtered.h5ad"
    result.write_h5ad(path)

    restored = ad.read_h5ad(path).uns[REPORT_KEY]
    assert list(restored["step"]) == list(result.uns[REPORT_KEY]["step"])


def test_unknown_modality_is_rejected():
    with pytest.raises(ValueError, match="Unknown modality"):
        data_filtering(make_visium(), modality="visium")


def test_qc_metric_availability_reports_missing_and_unexpected():
    adata = make_visium()
    del adata.obs["pct_counts_hb"]
    adata.obs["pct_counts_blank"] = 0.0

    availability = qc_metric_availability(adata, "visium_hd")

    assert availability["missing_required"] == []
    assert availability["missing_optional"] == ["pct_counts_hb"]
    assert availability["unexpected_unavailable"] == ["pct_counts_blank"]
    assert "center_x" not in QC_METRIC_CONTRACT["visium_hd"]["required_obs"]


# --------------------------------------------------------------------------
# MERFISH (existing behaviour preserved)
# --------------------------------------------------------------------------


@pytest.mark.parametrize("kwargs", [{"is_spatial": True}, {"modality": "merfish"}])
def test_merfish_flags_blanks_without_removing_observations(kwargs):
    counts = _counts(n_vars=8)
    counts[5, -2:] = 500
    adata = make_merfish(counts)

    result = data_filtering(
        adata, min_genes=1, min_cells=1, blank_max_percentage=10.0, **kwargs
    )
    report = result.uns[REPORT_KEY].set_index("step")

    assert result.n_obs == adata.n_obs
    assert bool(result.obs["is_blank"].iloc[5])
    assert result.obs["is_blank"].sum() == report.loc["blank_probes", "n_flagged"]
    assert report.loc["blank_probes", "action"] == "flagged"
    assert list(report.index) == ["min_genes", "min_cells", "blank_probes"]


def test_merfish_does_not_use_single_cell_steps():
    adata = make_merfish()
    with patch("src.filtering._filter_mitochondrial_genes") as mito, patch(
        "src.filtering._identify_doublets"
    ) as doublets:
        result = data_filtering(adata, is_spatial=True, min_genes=1, min_cells=1)

    mito.assert_not_called()
    doublets.assert_not_called()
    assert "predicted_doublet" not in result.obs
    assert "pct_counts_mt" not in result.obs


def test_merfish_keeps_blank_columns_and_coordinates_aligned():
    counts = _counts(n_vars=8)
    counts[2] = 0
    adata = make_merfish(counts)

    result = data_filtering(adata, is_spatial=True, min_genes=1, min_cells=1)

    assert adata.obs_names[2] not in result.obs_names
    positions = adata.obs_names.get_indexer(result.obs_names)
    np.testing.assert_array_equal(result.obsm["spatial"], adata.obsm["spatial"][positions])
    np.testing.assert_allclose(
        result.obs["pct_counts_blank"], adata.obs["pct_counts_blank"].iloc[positions]
    )


def test_legacy_single_cell_call_still_filters_mito_and_flags_doublets():
    adata = make_visium()
    del adata.obsm["spatial"]

    result = run(adata, min_genes=1, min_cells=1, max_mt_percentage=100)

    assert list(result.uns[REPORT_KEY]["step"]) == [
        "min_genes", "min_cells", "max_mt_percentage", "doublets",
    ]
    assert "predicted_doublet" in result.obs
