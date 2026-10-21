import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from src.utils import load_data, validate_visium_hd_adata


def make(counts=None, sparse_format="csr", uns=True):
    if counts is None:
        counts = np.random.default_rng(0).poisson(2.0, (12, 5)).astype(np.int32) + 1
    X = counts
    if sparse_format:
        X = getattr(sparse, f"{sparse_format}_matrix")(counts)
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame(index=[f"s_{i}-1" for i in range(counts.shape[0])]),
        var=pd.DataFrame(index=[f"G{i}" for i in range(counts.shape[1])]),
    )
    adata.obsm["spatial"] = np.column_stack(
        [np.arange(counts.shape[0]) * 8.0, np.arange(counts.shape[0]) % 3 * 8.0]
    )
    if uns:
        adata.uns["spatial"] = {"lib": {"images": {}, "scalefactors": {"spot_diameter_fullres": 8.0}}}
    return adata


def write_and_load(adata, tmp_path):
    path = tmp_path / "in.h5ad"
    adata.write_h5ad(path)
    return load_data(str(path), is_spatial=True, platform="Visium HD")["raw"]


def force_x(adata, value):
    """Bypass AnnData's own checks so the validator's guards can be exercised."""
    forced = type("ForcedXAnnData", (ad.AnnData,), {"X": property(lambda self: value)})
    adata.__class__ = forced
    return adata


# --- valid inputs are accepted ------------------------------------------------


@pytest.mark.parametrize("fmt", ["csr", "csc", None])
def test_accepts_sparse_and_dense_matrices(fmt):
    adata = make(sparse_format=fmt)
    validate_visium_hd_adata(adata)


def test_accepts_sparse_matrix_with_stored_zeros_and_empty_rows():
    counts = np.zeros((6, 4), dtype=np.int32)
    counts[1, 2] = 3
    counts[4, 0] = 1
    adata = make(counts)
    adata.X.data[0] = 0  # explicit stored zero
    validate_visium_hd_adata(adata)


def test_accepts_float_and_unsigned_dtypes():
    for dtype in (np.float32, np.float64, np.uint16):
        counts = (np.random.default_rng(1).poisson(2.0, (8, 4)) + 1).astype(dtype)
        validate_visium_hd_adata(make(counts))


def test_missing_optional_metadata_is_accepted():
    adata = make(uns=False)
    assert "spatial" not in adata.uns
    validate_visium_hd_adata(adata)

    adata.uns["spatial"] = {"lib": {}}  # no images, no scalefactors
    validate_visium_hd_adata(adata)


def test_coordinate_dataframe_aligned_to_obs_is_accepted():
    adata = make()
    adata.obsm["spatial"] = pd.DataFrame(
        adata.obsm["spatial"], index=adata.obs_names, columns=["x", "y"]
    )
    validate_visium_hd_adata(adata)


def test_validation_does_not_modify_or_print(capsys):
    adata = make()
    before_x = adata.X.toarray().copy()
    before_coords = adata.obsm["spatial"].copy()
    before_names = list(adata.obs_names), list(adata.var_names)

    validate_visium_hd_adata(adata)

    assert capsys.readouterr().out == ""
    np.testing.assert_array_equal(adata.X.toarray(), before_x)
    np.testing.assert_array_equal(adata.obsm["spatial"], before_coords)
    assert (list(adata.obs_names), list(adata.var_names)) == before_names


# --- structure and dimensions -------------------------------------------------


def test_rejects_non_anndata_input():
    with pytest.raises(ValueError, match="Expected an AnnData"):
        validate_visium_hd_adata(np.zeros((3, 3)))


@pytest.mark.parametrize("shape", [(0, 5), (5, 0)])
def test_rejects_empty_datasets_and_reports_shape(shape):
    adata = ad.AnnData(np.zeros(shape))
    with pytest.raises(ValueError, match="no observations or no genes"):
        validate_visium_hd_adata(adata)


def test_rejects_missing_expression_matrix():
    adata = make()
    adata.X = None
    with pytest.raises(ValueError, match="no expression matrix"):
        validate_visium_hd_adata(adata)


def test_rejects_expression_matrix_with_wrong_dimensions():
    adata = force_x(make(), sparse.csr_matrix(np.ones((3, 3))))
    with pytest.raises(ValueError, match="does not match the dataset dimensions"):
        validate_visium_hd_adata(adata)


def test_rejects_unsupported_matrix_type():
    adata = force_x(make(), [[1, 2], [3, 4]])
    with pytest.raises(ValueError, match="NumPy array or SciPy sparse"):
        validate_visium_hd_adata(adata)


# --- identifiers --------------------------------------------------------------


def test_rejects_duplicate_observation_and_gene_identifiers():
    adata = make()
    adata.obs_names = ["spot"] * adata.n_obs
    with pytest.raises(ValueError, match="Observation identifiers must be unique"):
        validate_visium_hd_adata(adata)

    adata = make()
    adata.var_names = ["GENE"] * adata.n_vars
    with pytest.raises(ValueError, match="Gene identifiers must be unique"):
        validate_visium_hd_adata(adata)


def test_duplicate_error_names_examples():
    adata = make()
    names = list(adata.var_names)
    names[1] = names[0]
    adata.var_names = names
    with pytest.raises(ValueError, match="G0"):
        validate_visium_hd_adata(adata)


@pytest.mark.parametrize("bad", ["", "   ", "nan", "None"])
def test_rejects_empty_or_missing_identifiers(bad):
    adata = make()
    names = list(adata.obs_names)
    names[2] = bad
    adata.obs_names = names
    with pytest.raises(ValueError, match="Observation identifiers must be non-empty"):
        validate_visium_hd_adata(adata)

    adata = make()
    names = list(adata.var_names)
    names[1] = bad
    adata.var_names = names
    with pytest.raises(ValueError, match="Gene identifiers must be non-empty"):
        validate_visium_hd_adata(adata)


# --- spatial coordinates ------------------------------------------------------


def test_rejects_missing_coordinates():
    adata = make()
    del adata.obsm["spatial"]
    with pytest.raises(ValueError, match="obsm\\['spatial'\\]"):
        validate_visium_hd_adata(adata)


@pytest.mark.parametrize("columns", [1, 3])
def test_rejects_wrong_coordinate_width(columns):
    adata = make()
    adata.obsm["spatial"] = np.zeros((adata.n_obs, columns))
    with pytest.raises(ValueError, match="shape"):
        validate_visium_hd_adata(adata)


def test_rejects_non_numeric_and_boolean_coordinates():
    adata = make()
    adata.obsm["spatial"] = np.full((adata.n_obs, 2), "a", dtype=object)
    with pytest.raises(ValueError, match="numeric"):
        validate_visium_hd_adata(adata)

    adata.obsm["spatial"] = np.ones((adata.n_obs, 2), dtype=bool)
    with pytest.raises(ValueError, match="not boolean"):
        validate_visium_hd_adata(adata)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_rejects_non_finite_coordinates_with_count(bad):
    adata = make()
    coords = adata.obsm["spatial"].copy()
    coords[[0, 3], 1] = bad
    adata.obsm["spatial"] = coords
    with pytest.raises(ValueError, match="non-finite values in 2 observation"):
        validate_visium_hd_adata(adata)


def test_rejects_coordinate_dataframe_misaligned_with_observations():
    # AnnData itself refuses a coordinate DataFrame keyed to different observations.
    adata = make()
    with pytest.raises(ValueError, match="does not match"):
        adata.obsm["spatial"] = pd.DataFrame(
            adata.obsm["spatial"], index=adata.obs_names[::-1], columns=["x", "y"]
        )


def test_rejects_sparse_coordinates():
    adata = make()
    adata.obsm["spatial"] = sparse.csr_matrix(adata.obsm["spatial"])
    with pytest.raises(ValueError, match="dense numeric array"):
        validate_visium_hd_adata(adata)


def test_rejects_degenerate_identical_coordinates_but_not_single_observation():
    adata = make()
    adata.obsm["spatial"] = np.full((adata.n_obs, 2), 5.0)
    with pytest.raises(ValueError, match="no spatial extent"):
        validate_visium_hd_adata(adata)

    single = make(np.array([[1, 2, 3]], dtype=np.int32))
    validate_visium_hd_adata(single)


# --- expression values --------------------------------------------------------


def test_rejects_non_finite_expression_values():
    for bad in (np.nan, np.inf):
        counts = np.ones((6, 4), dtype=np.float64)
        counts[2, 2] = bad
        with pytest.raises(ValueError, match="non-finite"):
            validate_visium_hd_adata(make(counts))


def test_rejects_negative_expression_values():
    counts = np.ones((6, 4), dtype=np.float64)
    counts[1, 1] = -0.5
    with pytest.raises(ValueError, match="negative values"):
        validate_visium_hd_adata(make(counts))


def test_rejects_all_zero_matrix_but_accepts_zero_count_observations():
    with pytest.raises(ValueError, match="all zeros"):
        validate_visium_hd_adata(make(np.zeros((6, 4), dtype=np.int32)))

    counts = np.ones((6, 4), dtype=np.int32)
    counts[0] = 0
    validate_visium_hd_adata(make(counts))


def test_rejects_non_numeric_and_boolean_expression():
    for value in (np.full((12, 5), True), np.full((12, 5), "a", dtype=object)):
        adata = force_x(make(sparse_format=None), value)
        with pytest.raises(ValueError, match="numeric values"):
            validate_visium_hd_adata(adata)


def test_fractional_values_are_only_rejected_when_counts_are_required():
    counts = np.full((6, 4), 1.5)
    validate_visium_hd_adata(make(counts))
    with pytest.raises(ValueError, match="integer-like"):
        validate_visium_hd_adata(make(counts), require_raw_counts=True)


# --- optional metadata --------------------------------------------------------


def test_requires_spatial_metadata_only_when_requested():
    adata = make(uns=False)
    validate_visium_hd_adata(adata)
    with pytest.raises(ValueError, match="uns\\['spatial'\\]"):
        validate_visium_hd_adata(adata, require_spatial_metadata=True)


def test_rejects_malformed_spatial_metadata():
    adata = make()
    adata.uns["spatial"] = {"lib": "oops"}
    with pytest.raises(ValueError, match="malformed"):
        validate_visium_hd_adata(adata)

    adata.uns["spatial"] = {"lib": {"images": [1]}}
    with pytest.raises(ValueError, match="Image metadata"):
        validate_visium_hd_adata(adata)

    adata.uns["spatial"] = {"lib": {"scalefactors": [1]}}
    with pytest.raises(ValueError, match="Scale-factor metadata"):
        validate_visium_hd_adata(adata)


@pytest.mark.parametrize("bad", [0, -1.0, np.nan, "8", True])
def test_rejects_unusable_scale_factors(bad):
    adata = make()
    adata.uns["spatial"]["lib"]["scalefactors"]["spot_diameter_fullres"] = bad
    with pytest.raises(ValueError, match="spot_diameter_fullres"):
        validate_visium_hd_adata(adata)


def test_absent_scale_factors_are_accepted():
    adata = make()
    adata.uns["spatial"]["lib"]["scalefactors"] = {"tissue_hires_scalef": 0.05}
    validate_visium_hd_adata(adata)


# --- loader boundary: no silent overwrite ------------------------------------


def test_loader_validates_before_adding_columns(tmp_path):
    adata = make()
    del adata.obsm["spatial"]
    with pytest.raises(ValueError, match="obsm\\['spatial'\\]"):
        write_and_load(adata, tmp_path)


def test_loader_keeps_matching_center_columns(tmp_path):
    adata = make()
    adata.obs["center_x"] = adata.obsm["spatial"][:, 0]
    adata.obs["center_y"] = adata.obsm["spatial"][:, 1]

    loaded = write_and_load(adata, tmp_path)

    np.testing.assert_allclose(loaded.obs["center_x"], adata.obsm["spatial"][:, 0])


def test_loader_refuses_to_overwrite_conflicting_center_columns(tmp_path):
    adata = make()
    adata.obs["center_x"] = 999.0
    adata.obs["center_y"] = adata.obsm["spatial"][:, 1]

    with pytest.raises(ValueError, match="center_x.*differs from"):
        write_and_load(adata, tmp_path)


def test_loader_leaves_inputs_unchanged(tmp_path):
    adata = make()
    loaded = write_and_load(adata, tmp_path)

    assert list(loaded.obs_names) == list(adata.obs_names)
    assert list(loaded.var_names) == list(adata.var_names)
    np.testing.assert_array_equal(loaded.obsm["spatial"], adata.obsm["spatial"])
    np.testing.assert_array_equal(loaded.X.toarray(), adata.X.toarray())
