import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from unittest.mock import patch
from src.clustering import (_run_pca, _run_umap, 
                            _run_clustering, _normalize_and_log_transform,
                            _identify_highly_variable_genes, identify_marker_genes,
                            identify_spatial_domains, run_clustering, run_pca)


# Normalization

def test_normalize_and_log_transform(small_adata):
    result = _normalize_and_log_transform(small_adata)

    expected_first_cell = np.log1p(
        np.array([10, 0, 0, 5, 0]) / 15 * 1e4
    )

    assert np.allclose(result.X[0], expected_first_cell)


def test_normalize_and_log_transform_returns_same_adata(small_adata):
    result = _normalize_and_log_transform(small_adata)

    assert result is small_adata


# Highly variable genes

# def test_identify_highly_variable_genes_adds_annotations(small_adata):
#     result = _identify_highly_variable_genes(
#         small_adata,
#         n_top_genes=2
#     )

#     assert "highly_variable" in result.var


# def test_identify_highly_variable_genes_selects_requested_number(small_adata):
#     result = _identify_highly_variable_genes(
#         small_adata,
#         n_top_genes=2
#     )

#     assert result.var["highly_variable"].sum() == 2


def test_identify_highly_variable_genes_passes_n_top_genes(small_adata):
    with patch("src.filtering.sc.pp.highly_variable_genes") as mock_hvg:
        _identify_highly_variable_genes(
            small_adata,
            n_top_genes=3
        )

        mock_hvg.assert_called_once_with(
            small_adata,
            n_top_genes=3
        )


# PCA

def test_run_pca_creates_embedding(small_adata):
    result = _run_pca(
        small_adata,
        n_comps=2
    )

    assert "X_pca" in result.obsm
    assert result.obsm["X_pca"].shape == (5, 2)


def test_run_pca_passes_n_comps(small_adata):
    with patch("src.filtering.sc.pp.pca") as mock_pca:
        result = _run_pca(
            small_adata,
            n_comps=3
        )

        mock_pca.assert_called_once_with(
            small_adata,
            n_comps=3
        )

    assert result is small_adata


# UMAP

def test_run_umap_creates_embedding(small_adata):
    result = _run_umap(
        small_adata,
        n_neighbors=3,
        min_dist=0.1
    )

    assert "X_umap" in result.obsm
    assert result.obsm["X_umap"].shape == (5, 2)


def test_run_umap_passes_parameters(small_adata):
    with patch("src.filtering.sc.pp.neighbors") as mock_neighbors, \
         patch("src.filtering.sc.tl.umap") as mock_umap:

        result = _run_umap(
            small_adata,
            n_neighbors=7,
            min_dist=0.2
        )

    mock_neighbors.assert_called_once_with(
        small_adata,
        n_neighbors=7
    )

    mock_umap.assert_called_once_with(
        small_adata,
        min_dist=0.2
    )

    assert result is small_adata


def test_run_umap_runs_neighbors_before_umap(small_adata):
    calls = []

    def fake_neighbors(*args, **kwargs):
        calls.append("neighbors")

    def fake_umap(*args, **kwargs):
        calls.append("umap")

    with patch(
        "src.filtering.sc.pp.neighbors",
        side_effect=fake_neighbors
    ), patch(
        "src.filtering.sc.tl.umap",
        side_effect=fake_umap
    ):
        _run_umap(
            small_adata,
            n_neighbors=3,
            min_dist=0.1
        )

    assert calls == ["neighbors", "umap"]


# Clustering

def test_run_clustering_creates_cluster_labels(small_adata):
    sc.pp.neighbors(small_adata, n_neighbors=2)

    result = _run_clustering(
        small_adata,
        resolution=1.0
    )

    assert "leiden" in result.obs
    assert result.obs["leiden"].notna().all()


def test_run_clustering_passes_resolution(small_adata):
    with patch("src.filtering.sc.tl.leiden") as mock_leiden:
        result = _run_clustering(
            small_adata,
            resolution=0.5
        )

    mock_leiden.assert_called_once_with(
        small_adata,
        resolution=0.5
    )

    assert result is small_adata

def test_run_pca_pipeline(small_adata):
    normalized = small_adata.copy()
    variable_genes = small_adata.copy()
    pca_result = small_adata.copy()

    with patch(
        "src.clustering._normalize_and_log_transform",
        return_value=normalized
    ) as mock_normalize, patch(
        "src.clustering._identify_highly_variable_genes",
        return_value=variable_genes
    ) as mock_hvg, patch(
        "src.clustering._run_pca",
        return_value=pca_result
    ) as mock_pca:

        result = run_pca(
            small_adata,
            n_top_genes=100,
            n_comps=10
        )

    mock_normalize.assert_called_once_with(small_adata)

    mock_hvg.assert_called_once_with(
        normalized,
        n_top_genes=100
    )

    mock_pca.assert_called_once_with(
        variable_genes,
        n_comps=10
    )

    assert result is pca_result


def test_run_pca_uses_default_parameters(small_adata):
    with patch(
        "src.clustering._normalize_and_log_transform",
        return_value=small_adata
    ) as mock_normalize, patch(
        "src.clustering._identify_highly_variable_genes",
        return_value=small_adata
    ) as mock_hvg, patch(
        "src.clustering._run_pca",
        return_value=small_adata
    ) as mock_pca:

        result = run_pca(small_adata)

    mock_normalize.assert_called_once_with(small_adata)

    mock_hvg.assert_called_once_with(
        small_adata,
        n_top_genes=2000
    )

    mock_pca.assert_called_once_with(
        small_adata,
        n_comps=50
    )

    assert result is small_adata


# Marker genes

def test_identify_marker_genes(small_adata):
    with patch("src.filtering.sc.tl.rank_genes_groups") as mock_rank_genes:
        result = identify_marker_genes(
            small_adata,
            groupby="cluster_label"
        )

    mock_rank_genes.assert_called_once_with(
        small_adata,
        groupby="cluster_label",
        method="t-test"
    )

    assert result is small_adata


# Spatial domains

def test_identify_spatial_domains(small_adata):
    gene_graph = np.eye(small_adata.n_obs)
    spatial_graph = np.ones(
        (small_adata.n_obs, small_adata.n_obs)
    )

    small_adata.obsp["connectivities"] = gene_graph
    small_adata.obsp["spatial_connectivities"] = spatial_graph

    with patch("src.clustering.sc.pp.neighbors") as mock_neighbors, \
         patch("src.clustering.sq.gr.spatial_neighbors") as mock_spatial, \
         patch("src.clustering.sc.tl.leiden") as mock_leiden:

        result = identify_spatial_domains(
            small_adata,
            alpha=0.2
        )

    mock_neighbors.assert_called_once_with(small_adata)
    mock_spatial.assert_called_once_with(small_adata)
    mock_leiden.assert_called_once()

    assert result is small_adata


def test_identify_spatial_domains_builds_joint_graph(small_adata):
    gene_graph = np.eye(small_adata.n_obs)
    spatial_graph = np.ones(
        (small_adata.n_obs, small_adata.n_obs)
    )

    small_adata.obsp["connectivities"] = gene_graph
    small_adata.obsp["spatial_connectivities"] = spatial_graph

    alpha = 0.2

    with patch("src.clustering.sc.pp.neighbors"), \
         patch("src.clustering.sq.gr.spatial_neighbors"), \
         patch("src.clustering.sc.tl.leiden") as mock_leiden:

        identify_spatial_domains(
            small_adata,
            alpha=alpha
        )

    _, kwargs = mock_leiden.call_args

    graph = kwargs["adjacency"]

    if hasattr(graph, "toarray"):
        graph = graph.toarray()

    expected_graph = (
        (1 - alpha) * gene_graph
        + alpha * spatial_graph
    )

    assert np.allclose(graph, expected_graph)


# Clustering pipeline

def test_run_clustering_pipeline(small_adata):
    umap_result = small_adata.copy()
    clustering_result = small_adata.copy()

    clustering_result.obs["leiden"] = pd.Categorical(
        ["0", "1", "0", "1", "0"]
    )

    with patch(
        "src.clustering._run_umap",
        return_value=umap_result
    ) as mock_umap, patch(
        "src.clustering._run_clustering",
        return_value=clustering_result
    ) as mock_clustering:

        result = run_clustering(
            small_adata,
            n_neighbors=10,
            min_dist=0.2,
            resolution=0.5
        )

    mock_umap.assert_called_once_with(
        small_adata,
        n_neighbors=10,
        min_dist=0.2
    )

    mock_clustering.assert_called_once_with(
        umap_result,
        resolution=0.5
    )

    assert result is clustering_result


def test_run_clustering_creates_cluster_labels(small_adata):
    small_adata.obs["leiden"] = pd.Categorical(
        ["0", "1", "0", "1", "0"]
    )

    with patch(
        "src.clustering._run_umap",
        return_value=small_adata
    ), patch(
        "src.clustering._run_clustering",
        return_value=small_adata
    ):

        result = run_clustering(small_adata)

    assert "cluster_label" in result.obs

    assert list(result.obs["cluster_label"]) == [
        "Cluster 0",
        "Cluster 1",
        "Cluster 0",
        "Cluster 1",
        "Cluster 0",
    ]

    assert str(result.obs["cluster_label"].dtype) == "category"
