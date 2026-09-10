
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt


def _normalize_and_log_transform(adata: sc.AnnData) -> sc.AnnData:
    """
    Normalize and log-transform the data in an AnnData object.

    Parameters:
    - adata: sc.AnnData
        The AnnData object containing gene expression data.
    
    Returns:
    - adata: sc.AnnData
        The normalized and log-transformed AnnData object.
    """

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    return adata

def _identify_highly_variable_genes(adata: sc.AnnData, n_top_genes: int = 2000) -> sc.AnnData:
    """
    Identify highly variable genes in an AnnData object.

    Parameters:
    - adata: sc.AnnData
        The AnnData object containing gene expression data.
    - n_top_genes: int, optional
        Number of top highly variable genes to identify. Default is 2000.
    
    Returns:
    - adata: sc.AnnData
        The AnnData object with highly variable genes annotated.
    """

    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)

    return adata

def _run_pca(adata: sc.AnnData, n_comps: int = 50) -> sc.AnnData:
    """
    Perform PCA dimensionality reduction on an AnnData object.

    Parameters:
    - adata: sc.AnnData
        The AnnData object containing gene expression data.
    - n_comps: int, optional
        Number of principal components to use for PCA. Default is 50.

    Returns:
    - adata: sc.AnnData
        The AnnData object with PCA embeddings.
    """

    sc.pp.pca(adata, n_comps=n_comps)

    return adata

def run_pca(adata: sc.AnnData, n_top_genes: int = 2000, n_comps: int = 50) -> sc.AnnData:
    """
    Run PCA on an AnnData object after normalizing, log-transforming, and identifying highly variable genes.

    Parameters:
    - adata: sc.AnnData
        The AnnData object containing gene expression data.
    - n_top_genes: int, optional
        Number of top highly variable genes to identify. Default is 2000.
    - n_comps: int, optional
        Number of principal components to use for PCA. Default is 50.

    Returns:
    - adata: sc.AnnData
        The AnnData object with PCA embeddings.
    """

    adata = _normalize_and_log_transform(adata)
    adata = _identify_highly_variable_genes(adata, n_top_genes=n_top_genes)
    adata = _run_pca(adata, n_comps=n_comps)

    return adata

def _run_umap(adata: sc.AnnData, n_neighbors: int = 15, min_dist: float = 0.1) -> sc.AnnData:
    """
    Perform UMAP dimensionality reduction on an AnnData object.

    Parameters:
    - adata: sc.AnnData
        The AnnData object containing gene expression data.
    - n_neighbors: int, optional
        Number of neighbors to use for UMAP. Default is 15.
    - min_dist: float, optional
        Minimum distance parameter for UMAP. Default is 0.1.

    Returns:
    - adata: sc.AnnData
        The AnnData object with UMAP embeddings.
    """

    sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    sc.tl.umap(adata, min_dist=min_dist)

    return adata

def _run_clustering(adata: sc.AnnData, resolution: float = 1.0) -> sc.AnnData:
    """
    Perform clustering on an AnnData object using the Leiden algorithm.

    Parameters:
    - adata: sc.AnnData
        The AnnData object containing gene expression data.
    - resolution: float, optional
        Resolution parameter for the Leiden algorithm. Default is 1.0.

    Returns:
    - adata: sc.AnnData
        The AnnData object with cluster labels.
    """

    sc.tl.leiden(adata, resolution=resolution)

    return adata

def run_clustering(adata: sc.AnnData, n_neighbors: int = 15, min_dist: float = 0.1, resolution: float = 1.0) -> sc.AnnData:
    """
    Run UMAP and clustering on an AnnData object.

    Parameters:
    - adata: sc.AnnData
        The AnnData object containing gene expression data.
    - n_neighbors: int, optional
        Number of neighbors to use for UMAP. Default is 15.
    - min_dist: float, optional
        Minimum distance parameter for UMAP. Default is 0.1.
    - resolution: float, optional
        Resolution parameter for the Leiden algorithm. Default is 1.0.

    Returns:
    - adata: sc.AnnData
        The AnnData object with UMAP embeddings and cluster labels.
    """

    adata = _run_umap(adata, n_neighbors=n_neighbors, min_dist=min_dist)
    adata = _run_clustering(adata, resolution=resolution)
    
    adata.obs["cluster_label"] = (
            "Cluster "
            + adata.obs["leiden"].astype(str)).astype("category")


    return adata

def identify_marker_genes(adata: sc.AnnData, groupby: str = "cluster_label") -> sc.AnnData:
    """
    Identify marker genes for each cluster in an AnnData object.

    Parameters:
    - adata: sc.AnnData
        The AnnData object containing gene expression data.
    - groupby: str, optional
        The column in adata.obs to group by for marker gene identification. Default is "cluster_label".

    Returns:
    - adata: sc.AnnData
        The AnnData object with marker genes identified.
    """

    sc.tl.rank_genes_groups(adata, groupby=groupby, method="t-test")

    return adata