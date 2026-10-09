import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from src import interactive
from src.clustering import run_spatially_variable_genes


@pytest.fixture
def spatial_adata():
    rng = np.random.default_rng(0)
    side = 15
    xs, ys = np.meshgrid(np.arange(side), np.arange(side))
    coords = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    n = len(coords)

    X = rng.poisson(3, (n, 30)).astype(float)
    X[:, 0] = rng.poisson(1 + 8 * (coords[:, 0] > side / 2))  # spatial pattern
    adata = ad.AnnData(X)
    adata.var_names = ["PATTERN"] + [f"G{i}" for i in range(1, 29)] + ["Blank-2"]
    adata.obsm["spatial"] = coords
    adata.obs["center_x"] = coords[:, 0]
    adata.obs["center_y"] = coords[:, 1]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=15)
    adata.var.loc["PATTERN", "highly_variable"] = True
    return adata


def test_hvg_table_and_plot(spatial_adata):
    table = interactive.hvg_table(spatial_adata)
    assert table["rank"].tolist() == list(range(1, len(table) + 1))
    assert set(table["gene"]) == set(
        spatial_adata.var_names[spatial_adata.var["highly_variable"]]
    )
    assert table["dispersions_norm"].is_monotonic_decreasing

    fig = interactive.hvg_plot(spatial_adata)
    assert sum(len(t.x) for t in fig.data[:2]) == spatial_adata.n_vars
    assert fig.layout.xaxis.type == "log"


def test_hvg_requires_annotation(spatial_adata):
    del spatial_adata.var["highly_variable"]
    assert interactive.hvg_table(spatial_adata).empty
    with pytest.raises(ValueError):
        interactive.hvg_plot(spatial_adata)


def test_spatially_variable_genes_end_to_end(spatial_adata):
    run_spatially_variable_genes(spatial_adata)

    tested = spatial_adata.uns["moranI"].index
    assert not any("Blank-" in g for g in tested)
    assert set(tested) <= set(
        spatial_adata.var_names[spatial_adata.var["highly_variable"]]
    )

    table = interactive.svg_table(spatial_adata)
    assert table.loc[0, "gene"] == "PATTERN"
    assert bool(table.loc[0, "significant"])

    fig = interactive.svg_plot(spatial_adata)
    assert sum(len(t.x) for t in fig.data[:2]) == len(table)


def test_svg_can_test_all_genes(spatial_adata):
    run_spatially_variable_genes(spatial_adata, only_highly_variable=False)
    assert len(spatial_adata.uns["moranI"]) == spatial_adata.n_vars - 1  # no blank


def test_svg_errors(spatial_adata):
    with pytest.raises(ValueError):
        interactive.svg_plot(spatial_adata)
    del spatial_adata.obsm["spatial"]
    with pytest.raises(ValueError):
        run_spatially_variable_genes(spatial_adata)
