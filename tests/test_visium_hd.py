import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from unittest.mock import patch

from src.filtering import data_filtering
from src.utils import load_data, validate_visium_hd_adata

N_SIDE = 6
N_GENES = 8
GENE_NAMES = ["MT-CO1", "MT-ND1", "RPS3", "RPL5", "HBB", "ACTB", "GAPDH", "EPCAM"]


def _counts(n_obs=N_SIDE * N_SIDE, seed=0):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(3.0, size=(n_obs, N_GENES)).astype(np.int32)
    counts[:, 5] += 1  # ACTB detected in every bin
    return counts


def make_visium_hd(*, counts=None, as_sparse=True, spatial=True, uns=True):
    counts = _counts() if counts is None else counts
    n_obs = counts.shape[0]
    X = sparse.csr_matrix(counts) if as_sparse else counts
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame(index=[f"s_008um_{i:05d}-1" for i in range(n_obs)]),
        var=pd.DataFrame(index=GENE_NAMES[: counts.shape[1]]),
    )
    if spatial:
        side = int(np.ceil(np.sqrt(n_obs)))
        adata.obsm["spatial"] = np.array(
            [[(i % side) * 8.0 + 100.0, (i // side) * 8.0 + 200.0] for i in range(n_obs)]
        )
    if uns:
        adata.uns["spatial"] = {
            "lib": {"images": {}, "scalefactors": {"spot_diameter_fullres": 8.0}}
        }
    return adata


def load(adata, tmp_path):
    path = tmp_path / "visium_hd.h5ad"
    adata.write_h5ad(path)
    return load_data(str(path), is_spatial=True, platform="Visium HD")["raw"]


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------


@pytest.mark.parametrize("as_sparse", [True, False])
def test_validate_accepts_valid_data(as_sparse):
    validate_visium_hd_adata(make_visium_hd(as_sparse=as_sparse))


def test_validate_accepts_missing_optional_spatial_metadata():
    validate_visium_hd_adata(make_visium_hd(uns=False))


def test_validate_requires_spatial_metadata_when_requested():
    with pytest.raises(ValueError, match="uns\\['spatial'\\]"):
        validate_visium_hd_adata(make_visium_hd(uns=False), require_spatial_metadata=True)


def test_validate_rejects_missing_coordinates():
    with pytest.raises(ValueError, match="obsm\\['spatial'\\]"):
        validate_visium_hd_adata(make_visium_hd(spatial=False))


def test_validate_rejects_wrong_coordinate_shape():
    adata = make_visium_hd()
    adata.obsm["spatial"] = np.zeros((adata.n_obs, 3))
    with pytest.raises(ValueError, match="shape"):
        validate_visium_hd_adata(adata)


def test_validate_rejects_non_finite_coordinates():
    adata = make_visium_hd()
    coords = adata.obsm["spatial"].copy()
    coords[0, 0] = np.nan
    adata.obsm["spatial"] = coords
    with pytest.raises(ValueError, match="non-finite"):
        validate_visium_hd_adata(adata)


def test_validate_rejects_non_numeric_coordinates():
    adata = make_visium_hd()
    adata.obsm["spatial"] = np.full((adata.n_obs, 2), "a", dtype=object)
    with pytest.raises(ValueError, match="numeric"):
        validate_visium_hd_adata(adata)


def test_validate_rejects_duplicate_identifiers():
    adata = make_visium_hd()
    adata.var_names = ["GENE"] * adata.n_vars
    with pytest.raises(ValueError, match="Gene identifiers"):
        validate_visium_hd_adata(adata)

    adata = make_visium_hd()
    adata.obs_names = ["spot"] * adata.n_obs
    with pytest.raises(ValueError, match="Observation identifiers"):
        validate_visium_hd_adata(adata)


def test_validate_rejects_empty_dataset():
    empty = ad.AnnData(np.zeros((0, N_GENES)))
    with pytest.raises(ValueError, match="no observations"):
        validate_visium_hd_adata(empty)


def test_validate_rejects_non_finite_expression():
    counts = _counts().astype(np.float64)
    counts[1, 1] = np.inf
    with pytest.raises(ValueError, match="non-finite"):
        validate_visium_hd_adata(make_visium_hd(counts=counts))


def test_validate_raw_count_checks_are_opt_in():
    fractional = make_visium_hd(counts=_counts().astype(np.float64) + 0.5)
    validate_visium_hd_adata(fractional)
    with pytest.raises(ValueError, match="integer-like"):
        validate_visium_hd_adata(fractional, require_raw_counts=True)

    negative = _counts().astype(np.float64)
    negative[0, 0] = -1.0
    with pytest.raises(ValueError, match="negative"):
        validate_visium_hd_adata(make_visium_hd(counts=negative), require_raw_counts=True)


def test_validate_rejects_malformed_spatial_metadata():
    adata = make_visium_hd()
    adata.uns["spatial"] = {"lib": "not-a-dict"}
    with pytest.raises(ValueError, match="malformed"):
        validate_visium_hd_adata(adata)

    adata.uns["spatial"] = {"lib": {"images": [1, 2]}}
    with pytest.raises(ValueError, match="Image metadata"):
        validate_visium_hd_adata(adata)


# --------------------------------------------------------------------------
# Loading and QC
# --------------------------------------------------------------------------


@pytest.mark.parametrize("as_sparse", [True, False])
def test_load_preserves_identity_and_coordinates(tmp_path, as_sparse):
    source = make_visium_hd(as_sparse=as_sparse)
    loaded = load(source, tmp_path)

    assert list(loaded.obs_names) == list(source.obs_names)
    assert list(loaded.var_names) == list(source.var_names)
    np.testing.assert_array_equal(loaded.obsm["spatial"], source.obsm["spatial"])
    assert loaded.shape == source.shape


def test_load_does_not_modify_expression_matrix(tmp_path):
    source = make_visium_hd()
    loaded = load(source, tmp_path)

    assert sparse.issparse(loaded.X)
    np.testing.assert_array_equal(loaded.X.toarray(), source.X.toarray())


def test_load_adds_single_cell_qc_and_center_columns(tmp_path):
    source = make_visium_hd()
    loaded = load(source, tmp_path)
    counts = source.X.toarray()

    for column in (
        "n_genes_by_counts",
        "total_counts",
        "log1p_total_counts",
        "log1p_n_genes_by_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
        "pct_counts_hb",
        "center_x",
        "center_y",
    ):
        assert column in loaded.obs.columns

    np.testing.assert_allclose(loaded.obs["total_counts"], counts.sum(axis=1))
    np.testing.assert_array_equal(loaded.obs["n_genes_by_counts"], (counts > 0).sum(axis=1))
    np.testing.assert_allclose(
        loaded.obs["total_counts_mt"], counts[:, :2].sum(axis=1)
    )
    np.testing.assert_allclose(
        loaded.obs["pct_counts_mt"],
        counts[:, :2].sum(axis=1) / counts.sum(axis=1) * 100,
    )
    np.testing.assert_allclose(loaded.obs["center_x"], source.obsm["spatial"][:, 0])
    np.testing.assert_allclose(loaded.obs["center_y"], source.obsm["spatial"][:, 1])


def test_load_does_not_create_merfish_blank_columns(tmp_path):
    loaded = load(make_visium_hd(), tmp_path)
    assert "blank_counts" not in loaded.obs.columns
    assert "pct_counts_blank" not in loaded.obs.columns


def test_load_without_optional_spatial_metadata(tmp_path):
    loaded = load(make_visium_hd(uns=False), tmp_path)
    assert "spatial" not in loaded.uns
    assert "total_counts" in loaded.obs.columns


def test_load_handles_zero_count_observations(tmp_path):
    counts = _counts()
    counts[[3, 10]] = 0
    loaded = load(make_visium_hd(counts=counts), tmp_path)

    zero_rows = [3, 10]
    assert (loaded.obs["total_counts"].iloc[zero_rows] == 0).all()
    assert (loaded.obs["n_genes_by_counts"].iloc[zero_rows] == 0).all()
    assert np.isfinite(loaded.obs["total_counts"]).all()
    assert loaded.n_obs == counts.shape[0]


def test_load_rejects_missing_coordinates(tmp_path):
    with pytest.raises(ValueError, match="obsm\\['spatial'\\]"):
        load(make_visium_hd(spatial=False), tmp_path)


def test_load_rejects_non_finite_coordinates(tmp_path):
    adata = make_visium_hd()
    coords = adata.obsm["spatial"].copy()
    coords[2, 1] = np.inf
    adata.obsm["spatial"] = coords
    with pytest.raises(ValueError, match="non-finite"):
        load(adata, tmp_path)


def test_load_platform_name_is_case_insensitive(tmp_path):
    path = tmp_path / "visium_hd.h5ad"
    make_visium_hd().write_h5ad(path)
    loaded = load_data(str(path), is_spatial=True, platform="visium hd")["raw"]
    assert "center_x" in loaded.obs.columns


# --------------------------------------------------------------------------
# Filtering
# --------------------------------------------------------------------------


def _filter(adata, **kwargs):
    # Scrublet needs a realistically sized matrix; its behaviour is not under test here.
    with patch("src.filtering._identify_doublets", side_effect=lambda a, rate: a) as doublets:
        result = data_filtering(adata, is_spatial=False, **kwargs)
    return result, doublets


def test_filtering_uses_single_cell_qc_path_for_visium_hd(tmp_path):
    loaded = load(make_visium_hd(), tmp_path)
    result, doublets = _filter(
        loaded, min_genes=1, min_cells=1, max_mt_percentage=100.0, expected_doublet_rate=0.05
    )

    doublets.assert_called_once()
    assert "is_blank" not in result.obs.columns


def test_filtering_keeps_obs_var_and_spatial_alignment(tmp_path):
    counts = _counts()
    counts[0, :] = 0
    counts[1, :] = 0
    counts[1, 5] = 1  # single detected gene
    counts[:, 7] = 0
    counts[4, 7] = 1  # gene seen in only one bin
    loaded = load(make_visium_hd(counts=counts), tmp_path)

    result, _ = _filter(
        loaded, min_genes=2, min_cells=2, max_mt_percentage=100.0, expected_doublet_rate=0.05
    )

    assert loaded.obs_names[0] not in result.obs_names
    assert loaded.obs_names[1] not in result.obs_names
    assert "EPCAM" not in result.var_names
    assert result.obsm["spatial"].shape == (result.n_obs, 2)

    position = loaded.obs_names.get_indexer(result.obs_names)
    np.testing.assert_array_equal(result.obsm["spatial"], loaded.obsm["spatial"][position])
    np.testing.assert_allclose(result.obs["center_x"], loaded.obs["center_x"].iloc[position])
    kept_vars = loaded.var_names.get_indexer(result.var_names)
    np.testing.assert_array_equal(
        result.X.toarray(), loaded.X.toarray()[np.ix_(position, kept_vars)]
    )
    assert sparse.issparse(result.X)


def test_filtering_removes_zero_count_observations(tmp_path):
    counts = _counts()
    counts[[3, 10]] = 0
    loaded = load(make_visium_hd(counts=counts), tmp_path)

    result, _ = _filter(
        loaded, min_genes=1, min_cells=1, max_mt_percentage=100.0, expected_doublet_rate=0.05
    )

    assert loaded.obs_names[3] not in result.obs_names
    assert loaded.obs_names[10] not in result.obs_names
    assert result.n_obs == loaded.n_obs - 2


def test_mitochondrial_filter_removes_high_mt_bins(tmp_path):
    counts = _counts()
    counts[5] = [50, 50, 0, 0, 0, 1, 0, 0]  # ~99% mitochondrial
    loaded = load(make_visium_hd(counts=counts), tmp_path)

    result, _ = _filter(
        loaded, min_genes=1, min_cells=1, max_mt_percentage=90.0, expected_doublet_rate=0.05
    )

    assert loaded.obs_names[5] not in result.obs_names
    assert (result.obs["pct_counts_mt"] < 90.0).all()


def test_filtering_does_not_mutate_source_matrix(tmp_path):
    loaded = load(make_visium_hd(), tmp_path)
    before = loaded.X.toarray().copy()
    before_names = list(loaded.obs_names)

    _filter(loaded, min_genes=3, min_cells=2, max_mt_percentage=50.0, expected_doublet_rate=0.05)

    np.testing.assert_array_equal(loaded.X.toarray(), before)
    assert list(loaded.obs_names) == before_names


def test_filtering_with_scrublet_marks_doublets_on_realistic_sized_data(tmp_path):
    rng = np.random.default_rng(1)
    counts = rng.poisson(1.0, size=(600, 250)).astype(np.int32)
    adata = ad.AnnData(
        sparse.csr_matrix(counts),
        obs=pd.DataFrame(index=[f"s_{i}-1" for i in range(600)]),
        var=pd.DataFrame(index=[f"G{i}" for i in range(250)]),
    )
    adata.obsm["spatial"] = rng.uniform(0, 1000, size=(600, 2))
    loaded = load(adata, tmp_path)

    result = data_filtering(
        loaded,
        min_genes=20,
        min_cells=3,
        max_mt_percentage=100.0,
        expected_doublet_rate=0.05,
        is_spatial=False,
    )

    assert "predicted_doublet" in result.obs.columns
    assert result.obs["predicted_doublet"].dtype == bool
    assert result.obsm["spatial"].shape == (result.n_obs, 2)


def test_filtering_that_removes_everything_returns_empty_dataset(tmp_path):
    loaded = load(make_visium_hd(), tmp_path)

    result, doublets = _filter(
        loaded, min_genes=10_000, min_cells=1, max_mt_percentage=100.0, expected_doublet_rate=0.05
    )

    assert result.n_obs == 0
    assert result.obsm["spatial"].shape == (0, 2)
    doublets.assert_not_called()


def test_gene_filter_that_removes_everything_leaves_no_genes(tmp_path):
    loaded = load(make_visium_hd(), tmp_path)

    result, _ = _filter(
        loaded, min_genes=1, min_cells=10_000, max_mt_percentage=100.0, expected_doublet_rate=0.05
    )

    assert result.n_vars == 0
    assert result.n_obs == loaded.n_obs
