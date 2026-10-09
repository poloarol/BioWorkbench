import numpy as np
import pytest
import scanpy as sc
from unittest.mock import patch

from src.utils import (
    _calculate_qc_metrics,
    calculate_module_score,
    get_spatially_variable_genes,
    load_data,
    write_to_disk,
    validate_spatial_adata,
)


def _spatial_adata(**overrides):
    adata = sc.AnnData(
        X=np.array([[2, 1, 1], [4, 0, 2]], dtype=float),
        obs={"center_x": [1.0, 2.0], "center_y": [3.0, 4.0]},
    )
    adata.var_names = ["GeneA", "Blank-1", "GeneB"]
    for column, values in overrides.items():
        adata.obs[column] = values
    return adata


def test_validate_spatial_adata_accepts_valid_data():
    validate_spatial_adata(_spatial_adata())


def test_validate_spatial_adata_requires_coordinates():
    adata = _spatial_adata()
    del adata.obs["center_y"]
    with pytest.raises(ValueError, match="center_y"):
        validate_spatial_adata(adata)


def test_validate_spatial_adata_rejects_non_finite_coordinates():
    with pytest.raises(ValueError, match="non-finite"):
        validate_spatial_adata(_spatial_adata(center_x=[1.0, np.nan]))


def test_validate_spatial_adata_rejects_non_numeric_coordinates():
    with pytest.raises(ValueError, match="numeric"):
        validate_spatial_adata(_spatial_adata(center_x=["a", "b"]))


def test_validate_spatial_adata_requires_blank_genes():
    adata = _spatial_adata()
    adata.var_names = ["GeneA", "GeneC", "GeneB"]
    with pytest.raises(ValueError, match="Blank-"):
        validate_spatial_adata(adata)


def test_validate_spatial_adata_rejects_all_zero_counts():
    adata = _spatial_adata()
    adata.X = np.zeros((2, 3))
    with pytest.raises(ValueError, match="zero total counts"):
        validate_spatial_adata(adata)


def test_load_data_handles_zero_count_cells(tmp_path):
    adata = _spatial_adata()
    adata.X = np.array([[2, 1, 1], [0, 0, 0]], dtype=float)
    filepath = tmp_path / "zero.h5ad"
    adata.write_h5ad(filepath)

    loaded = load_data(str(filepath), is_spatial=True)["raw"]

    assert loaded.obs["pct_counts_blank"].tolist() == pytest.approx([25.0, 0.0])


def test_calculate_qc_metrics_for_single_cell_data(small_adata):
    _calculate_qc_metrics(small_adata)

    assert small_adata.var["mt"].tolist() == [False, False, False, False, True]
    assert small_adata.obs["total_counts"].tolist() == [15, 13, 18, 17, 12]
    assert "pct_counts_mt" in small_adata.obs


def test_calculate_qc_metrics_for_spatial_data(small_adata):
    _calculate_qc_metrics(small_adata, is_spatial=True)

    assert "total_counts" in small_adata.obs
    assert "n_genes_by_counts" in small_adata.obs
    assert "mt" not in small_adata.var


def test_load_data_prepares_single_cell_dataset(tmp_path, small_adata):
    filepath = tmp_path / "input.h5ad"
    small_adata.write_h5ad(filepath)

    result = load_data(str(filepath))

    assert result["raw"].n_obs == small_adata.n_obs
    assert result["subset"] is None
    assert np.array_equal(result["raw"].layers["raw"], small_adata.X)
    assert "pct_counts_mt" in result["raw"].obs
    assert result["raw"].var["mt"].tolist() == [
        False, False, False, False, True
    ]


def test_load_data_creates_requested_subset(tmp_path, small_adata):
    filepath = tmp_path / "input.h5ad"
    small_adata.write_h5ad(filepath)

    result = load_data(
        str(filepath),
        cell_columns=["Cell1", "Cell3"],
        genes=["GeneA", "GeneC"],
    )

    try:
        subset = result["subset"]
        assert list(subset.obs_names) == ["Cell1", "Cell3"]
        assert list(subset.var_names) == ["GeneA", "GeneC"]
        assert np.array_equal(subset.layers["raw"], subset.X)
        assert "pct_counts_mt" in subset.obs
    finally:
        result["raw"].file.close()


def test_load_data_prepares_spatial_metadata(tmp_path):
    adata = sc.AnnData(
        X=np.array([[2, 1, 1], [4, 0, 2]]),
        obs={
            "center_x": [10.0, 20.0],
            "center_y": [30.0, 40.0],
        },
    )
    adata.var_names = ["GeneA", "Blank-1", "GeneB"]
    filepath = tmp_path / "spatial.h5ad"
    adata.write_h5ad(filepath)

    result = load_data(str(filepath), is_spatial=True)
    loaded = result["raw"]

    assert result["subset"] is None
    assert loaded.obs["blank_counts"].tolist() == [1, 0]
    assert loaded.obs["pct_counts_blank"].tolist() == pytest.approx(
        [25.0, 0.0]
    )
    assert np.array_equal(
        loaded.obsm["spatial"],
        np.array([[10.0, 30.0], [20.0, 40.0]]),
    )
    assert loaded.obs["log1p_total_counts"].tolist() == pytest.approx(
        np.log1p([4, 6])
    )


@pytest.mark.parametrize(
    ("cell_columns", "genes", "expected_obs", "expected_vars"),
    [
        (["Cell1", "Cell4"], None, ["Cell1", "Cell4"], [
            "GeneA", "GeneB", "GeneC", "GeneD", "MT-GeneE"
        ]),
        (None, ["GeneA", "GeneD"], [
            "Cell1", "Cell2", "Cell3", "Cell4", "Cell5"
        ], ["GeneA", "GeneD"]),
    ],
)
def test_load_data_allows_subsetting_one_axis(
    tmp_path,
    small_adata,
    cell_columns,
    genes,
    expected_obs,
    expected_vars,
):
    filepath = tmp_path / "input.h5ad"
    small_adata.write_h5ad(filepath)

    result = load_data(
        str(filepath),
        cell_columns=cell_columns,
        genes=genes,
    )

    try:
        assert list(result["subset"].obs_names) == expected_obs
        assert list(result["subset"].var_names) == expected_vars
    finally:
        result["raw"].file.close()


def test_write_to_disk_round_trips_anndata(tmp_path, small_adata):
    filepath = tmp_path / "output.h5ad"

    write_to_disk(small_adata, str(filepath))

    loaded = sc.read_h5ad(filepath)
    assert list(loaded.obs_names) == list(small_adata.obs_names)
    assert list(loaded.var_names) == list(small_adata.var_names)
    assert np.array_equal(loaded.X, small_adata.X)


def test_calculate_module_score_uses_only_available_genes(small_adata):
    with patch("src.utils.sc.tl.score_genes") as mock_score_genes:
        result = calculate_module_score(
            small_adata,
            genes=["GeneA", "not-a-gene"],
            score_name="test_score",
        )

    mock_score_genes.assert_called_once_with(
        small_adata,
        gene_list=["GeneA"],
        score_name="test_score",
    )
    assert result is small_adata


def test_calculate_module_score_rejects_no_available_genes(small_adata):
    with pytest.raises(
        ValueError,
        match="None of the selected genes were found",
    ):
        calculate_module_score(small_adata, genes=["not-a-gene"])


def test_get_spatially_variable_genes_excludes_blank_probes(small_adata):
    small_adata.var_names = [
        "GeneA",
        "GeneB",
        "Blank-1",
        "GeneD",
        "MT-GeneE",
    ]

    def mark_selected_genes(adata, **kwargs):
        adata.var["highly_variable"] = [
            True, False, True, False, True
        ]

    with patch(
        "src.utils.sc.pp.highly_variable_genes",
        side_effect=mark_selected_genes,
    ) as mock_hvg:
        result = get_spatially_variable_genes(
            small_adata,
            n_top_genes=3,
        )

    mock_hvg.assert_called_once_with(
        small_adata,
        flavor="seurat_v3",
        n_top_genes=3,
        subset=False,
        inplace=True,
    )
    assert result == ["GeneA", "MT-GeneE"]
