import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from unittest.mock import patch

from src.filtering import REPORT_KEY, data_filtering, preview_filter_impact
from src.utils import load_data, qc_metric_availability, validate_xenium_adata


def make_xenium(*, spatial_in_obsm=False, counts=None):
    if counts is None:
        counts = np.full((4, 3), 10, dtype=np.int32)
    obs = pd.DataFrame(
        {
            "control_probe_counts": [0, 10, 0, 0],
            "control_codeword_counts": [0, 0, 20, 0],
            "cell_area": [100.0, 100.0, 200.0, 50.0],
            "nucleus_area": [50.0, 60.0, 100.0, 25.0],
            "x_centroid": [10.0, 20.0, 30.0, 40.0],
            "y_centroid": [1.0, 2.0, 3.0, 4.0],
        },
        index=[f"cell-{i}" for i in range(counts.shape[0])],
    )
    adata = ad.AnnData(
        X=sparse.csr_matrix(counts),
        obs=obs,
        var=pd.DataFrame(index=[f"gene-{i}" for i in range(counts.shape[1])]),
    )
    if spatial_in_obsm:
        adata.obsm["spatial"] = obs[["x_centroid", "y_centroid"]].to_numpy()
    return adata


def load(adata, tmp_path):
    path = tmp_path / "xenium.h5ad"
    adata.write_h5ad(path)
    return load_data(str(path), is_spatial=True, platform="Xenium")["raw"]


def test_load_xenium_builds_qc_and_spatial_fields_from_centroids(tmp_path):
    loaded = load(make_xenium(), tmp_path)

    np.testing.assert_allclose(loaded.obsm["spatial"], [[10, 1], [20, 2], [30, 3], [40, 4]])
    np.testing.assert_allclose(loaded.obs["center_x"], [10, 20, 30, 40])
    np.testing.assert_allclose(loaded.obs["center_y"], [1, 2, 3, 4])
    np.testing.assert_allclose(loaded.obs["negative_probe_pct"], [0, 33.3333, 0, 0], rtol=1e-4)
    np.testing.assert_allclose(loaded.obs["negative_codeword_pct"], [0, 0, 66.6667, 0], rtol=1e-4)
    np.testing.assert_allclose(loaded.obs["nucleus_ratio"], [0.5, 0.6, 0.5, 0.5])


def test_load_xenium_keeps_existing_spatial_coordinates(tmp_path):
    adata = make_xenium(spatial_in_obsm=True)
    adata.obsm["spatial"] = np.array([[100, 1], [200, 2], [300, 3], [400, 4]])

    loaded = load(adata, tmp_path)

    np.testing.assert_array_equal(loaded.obsm["spatial"], adata.obsm["spatial"])
    np.testing.assert_array_equal(loaded.obs["center_x"], [100, 200, 300, 400])


def test_validate_xenium_requires_controls_areas_and_coordinates():
    adata = make_xenium()
    del adata.obs["control_probe_counts"]
    with pytest.raises(ValueError, match="control_probe_counts"):
        validate_xenium_adata(adata)

    adata = make_xenium()
    del adata.obs["x_centroid"]
    del adata.obs["y_centroid"]
    with pytest.raises(ValueError, match="requires adata.obsm"):
        validate_xenium_adata(adata)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -1.0])
def test_validate_xenium_rejects_invalid_control_counts(bad):
    adata = make_xenium()
    adata.obs["control_probe_counts"] = adata.obs[
        "control_probe_counts"
    ].astype(float)
    adata.obs.loc["cell-1", "control_probe_counts"] = bad
    with pytest.raises(ValueError, match="control_probe_counts"):
        validate_xenium_adata(adata)


def test_validate_xenium_rejects_invalid_coordinates_and_area():
    adata = make_xenium(spatial_in_obsm=True)
    coordinates = adata.obsm["spatial"].copy()
    coordinates[1, 0] = np.nan
    adata.obsm["spatial"] = coordinates
    with pytest.raises(ValueError, match="coordinates must be finite"):
        validate_xenium_adata(adata)

    adata = make_xenium()
    adata.obs.loc["cell-0", "cell_area"] = 0
    with pytest.raises(ValueError, match="cell_area.*greater than zero"):
        validate_xenium_adata(adata)


def test_xenium_qc_contract_matches_loader_output(tmp_path):
    loaded = load(make_xenium(), tmp_path)

    assert qc_metric_availability(loaded, "xenium")["missing_required"] == []
    assert "blank_counts" not in loaded.obs


def test_xenium_filtering_flags_either_control_without_removing_cells(tmp_path):
    loaded = load(make_xenium(), tmp_path)
    with patch("src.filtering._identify_doublets") as doublets:
        filtered = data_filtering(
            loaded,
            min_genes=1,
            min_cells=1,
            modality="xenium",
            negative_probe_max_percentage=5,
            negative_codeword_max_percentage=10,
        )

    doublets.assert_not_called()
    assert filtered.n_obs == loaded.n_obs
    assert filtered.obs["is_xenium_control"].tolist() == [False, True, True, False]
    report = filtered.uns[REPORT_KEY].set_index("step")
    assert list(report.index) == ["min_genes", "min_cells", "negative_controls"]
    assert report.loc["negative_controls", "n_flagged"] == 2
    assert report.loc["negative_controls", "action"] == "flagged"

    clean = filtered[~filtered.obs["is_xenium_control"]].copy()
    assert list(clean.obs_names) == ["cell-0", "cell-3"]
    assert clean.obsm["spatial"].shape == (2, 2)


def test_xenium_filtering_requires_loaded_control_metrics(tmp_path):
    loaded = load(make_xenium(), tmp_path)
    del loaded.obs["negative_probe_pct"]

    with pytest.raises(ValueError, match="Xenium filtering requires"):
        data_filtering(loaded, modality="xenium")


def test_preview_filter_impact_reports_xenium_control_flags(tmp_path):
    loaded = load(make_xenium(), tmp_path)

    impact = preview_filter_impact(
        loaded,
        min_genes=1,
        min_cells=1,
        negative_probe_max_percentage=5,
        negative_codeword_max_percentage=10,
    )

    assert impact["obs_flagged_controls"] == 2
    assert impact["obs_remaining"] == 4


@pytest.mark.parametrize(
    "kwargs",
    [
        {"negative_probe_max_percentage": -1},
        {"negative_codeword_max_percentage": 101},
    ],
)
def test_xenium_filtering_rejects_invalid_control_thresholds(tmp_path, kwargs):
    loaded = load(make_xenium(), tmp_path)

    with pytest.raises(ValueError, match="Invalid Xenium filtering thresholds"):
        data_filtering(loaded, modality="xenium", **kwargs)
