import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from src import plotting


@pytest.fixture
def plotted_adata(small_adata):
    small_adata.obsm["X_pca"] = np.array(
        [[0, 0], [1, 0], [0, 1], [1, 1], [2, 1]],
        dtype=float,
    )
    small_adata.obsm["X_umap"] = np.array(
        [[0, 0], [1, 0], [0, 1], [1, 1], [2, 1]],
        dtype=float,
    )
    small_adata.obs["group"] = pd.Categorical(
        ["A", "A", "B", "B", "C"]
    )
    small_adata.obs["score"] = [0.1, 0.2, 0.3, 0.4, 0.5]
    small_adata.obs["center_x"] = [10, 20, 30, 40, 50]
    small_adata.obs["center_y"] = [5, 10, 15, 20, 25]
    small_adata.obs["volume"] = [1, 2, 3, 4, 5]
    small_adata.obs["n_genes_by_counts"] = [2, 3, 2, 4, 3]
    small_adata.obs["total_counts"] = [15, 13, 18, 17, 12]
    return small_adata


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def test_highest_expressed_genes_returns_figure(plotted_adata):
    with pytest.MonkeyPatch.context() as monkeypatch:
        called = {}

        def fake_plot(adata, **kwargs):
            called["adata"] = adata
            called.update(kwargs)

        monkeypatch.setattr(
            plotting.sc.pl,
            "highest_expr_genes",
            fake_plot,
        )

        fig = plotting.highest_expressed_genes(
            plotted_adata,
            n_top_genes=4,
        )

    assert called["adata"] is plotted_adata
    assert called["n_top"] == 4
    assert called["show"] is False
    assert isinstance(fig, plt.Figure)


def test_plot_qc_metrics_spatial_returns_three_panels(plotted_adata):
    fig = plotting.plot_qc_metrics(plotted_adata, is_spatial=True)

    assert len(fig.axes) == 3
    assert [ax.get_title() for ax in fig.axes] == [
        "Total transcripts per cell",
        "Unique transcripts per cell",
        "Volume of segmented cells",
    ]


def test_plot_qc_metrics_single_cell_returns_figure(plotted_adata):
    fig = plotting.plot_qc_metrics(plotted_adata)

    assert isinstance(fig, plt.Figure)
    assert len(fig.axes) == 3


def test_qc_summary_calculates_expected_statistics(plotted_adata):
    result = plotting.qc_summary(plotted_adata)

    assert result == {
        "n_cells": 5,
        "median_genes": 3.0,
        "mean_genes": 2.8,
        "median_counts": 15.0,
        "mean_counts": 15.0,
    }


@pytest.mark.parametrize(
    ("function", "coordinates", "color_by", "title"),
    [
        ("plot_pca", "X_pca", "group", "PCA — group"),
        ("plot_umap", "X_umap", "group", "UMAP — Group"),
        ("plot_umap", "X_umap", "score", "UMAP — Score"),
    ],
)
def test_embedding_plots_return_labeled_figures(
    plotted_adata,
    function,
    coordinates,
    color_by,
    title,
):
    fig = getattr(plotting, function)(plotted_adata, color_by=color_by)

    assert isinstance(fig, plt.Figure)
    assert fig.axes[0].get_title() == title
    assert coordinates in plotted_adata.obsm


def test_plot_spatial_handles_categorical_and_continuous_values(
    plotted_adata,
):
    fig = plotting.plot_spatial(
        plotted_adata,
        color_by=["group", "score"],
    )

    assert [ax.get_title() for ax in fig.axes[:2]] == [
        "Spatial — Group",
        "Spatial — Score",
    ]
    assert len(fig.axes) == 3


def test_plot_spatial_accepts_single_color_column(plotted_adata):
    fig = plotting.plot_spatial(plotted_adata, color_by="group")

    assert len(fig.axes) == 1
    assert fig.axes[0].get_title() == "Spatial — Group"


@pytest.mark.parametrize(
    ("genes", "titles"),
    [
        ("GeneA", ["UMAP — GeneA"]),
        (["GeneA", "GeneB"], ["UMAP — GeneA", "UMAP — GeneB"]),
    ],
)
def test_plot_umap_genes_accepts_one_or_multiple_genes(
    plotted_adata,
    genes,
    titles,
):
    fig = plotting.plot_umap_genes(plotted_adata, genes)

    assert [ax.get_title() for ax in fig.axes if ax.get_title()] == titles


@pytest.mark.parametrize(
    ("genes", "message"),
    [
        ([], "No genes were provided"),
        (["missing"], r"Gene\(s\) not found"),
    ],
)
def test_plot_umap_genes_rejects_invalid_genes(
    plotted_adata,
    genes,
    message,
):
    with pytest.raises(ValueError, match=message):
        plotting.plot_umap_genes(plotted_adata, genes)


def test_plot_umap_genes_requires_embedding(small_adata):
    with pytest.raises(ValueError, match="UMAP coordinates not found"):
        plotting.plot_umap_genes(small_adata, ["GeneA"])


def test_plot_spatial_genes_returns_figure(plotted_adata):
    fig = plotting.plot_spatial_genes(
        plotted_adata,
        ["GeneA", "GeneB"],
    )

    assert [ax.get_title() for ax in fig.axes if ax.get_title()] == [
        "Spatial — GeneA",
        "Spatial — GeneB",
    ]


@pytest.mark.parametrize(
    ("remove_columns", "genes", "message"),
    [
        (["center_x"], ["GeneA"], "Missing required spatial coordinate"),
        ([], ["missing"], r"Gene\(s\) not found"),
        ([], [], "No genes were provided"),
    ],
)
def test_plot_spatial_genes_rejects_invalid_input(
    plotted_adata,
    remove_columns,
    genes,
    message,
):
    for column in remove_columns:
        del plotted_adata.obs[column]

    with pytest.raises(ValueError, match=message):
        plotting.plot_spatial_genes(plotted_adata, genes)
