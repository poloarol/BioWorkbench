import numpy as np
import pandas as pd
import pytest

from src import interactive


@pytest.fixture
def adata(small_adata):
    small_adata.obsm["X_pca"] = np.array([[0, 0], [1, 0], [0, 1], [1, 1], [2, 1]], dtype=float)
    small_adata.obsm["X_umap"] = small_adata.obsm["X_pca"].copy()
    small_adata.obs["group"] = pd.Categorical(["A", "A", "B", "B", "C"])
    small_adata.obs["score"] = [0.1, 0.2, 0.3, 0.4, 0.5]
    small_adata.obs["center_x"] = [10, 20, 30, 40, 50]
    small_adata.obs["center_y"] = [5, 10, 15, 20, 25]
    small_adata.obs["total_counts"] = [15, 13, 18, 17, 12]
    return small_adata


@pytest.mark.parametrize("basis", ["umap", "pca", "spatial"])
def test_categorical_plot_has_one_trace_per_category(adata, basis):
    fig = interactive.interactive_plot(adata, basis, color_by="group")
    assert [t.name.split(" ")[0] for t in fig.data] == ["A", "B", "C"]
    assert sum(len(t.x) for t in fig.data) == adata.n_obs


def test_highlight_dims_other_groups(adata):
    fig = interactive.interactive_plot(adata, "umap", color_by="group", highlight=["B"])
    colors = {t.name.split(" ")[0]: t.marker.color for t in fig.data}
    assert colors["A"] == colors["C"] == interactive.DIM_COLOR
    assert colors["B"] != interactive.DIM_COLOR


def test_continuous_and_gene_plots_are_single_trace(adata):
    assert len(interactive.interactive_plot(adata, "umap", color_by="score").data) == 1
    gene = adata.var_names[0]
    assert len(interactive.interactive_plot(adata, "spatial", gene=gene).data) == 1


def test_plot_requires_color_and_valid_gene(adata):
    with pytest.raises(ValueError):
        interactive.interactive_plot(adata, "umap")
    with pytest.raises(ValueError):
        interactive.interactive_plot(adata, "umap", gene="not-a-gene")


def test_missing_embedding_raises(small_adata):
    with pytest.raises(ValueError):
        interactive.interactive_plot(small_adata, "umap", color_by=small_adata.obs.columns[0])


def test_group_statistics(adata):
    stats = interactive.group_statistics(adata, "group", groups=["A", "B"])
    assert list(stats["group"]) == ["A", "B"]
    assert list(stats["n_cells"]) == [2, 2]
    assert stats["pct_of_total"].tolist() == [40.0, 40.0]
    assert stats.loc[0, "total_counts (mean)"] == 14.0
    assert stats.loc[0, "centroid_x"] == 15.0


def test_group_statistics_all_groups_with_gene(adata):
    gene = adata.var_names[0]
    stats = interactive.group_statistics(adata, "group", gene=gene)
    assert len(stats) == 3
    assert f"{gene} % expressing" in stats.columns


def test_selection_statistics_and_event_parsing(adata):
    event = {"selection": {"points": [{"customdata": [0, "a"]}, {"customdata": [3, "b"]}]}}
    positions = interactive.selected_positions(event)
    assert positions.tolist() == [0, 3]

    stats = interactive.selection_statistics(adata, positions, group_col="group")
    assert stats["n_cells"].tolist() == [2, 1, 1]
    assert stats["selection"].tolist() == ["all selected", "group = A", "group = B"]
    assert interactive.selection_statistics(adata, []).empty
    assert interactive.selected_positions(None).size == 0
