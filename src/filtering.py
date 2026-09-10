import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


def _filter_cells(adata: sc.AnnData, min_genes: int = 200) -> sc.AnnData:
    """
    Filter cells in an AnnData object based on quality control metrics.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to filter.
    - min_genes: int, optional
        Minimum number of genes expressed in a cell to keep it. Default is 200.

    Returns:
    - filtered_adata: sc.AnnData
        The filtered AnnData object.
    """
    
    # Filter cells based on the specified criteria
    cell_subset, number_per_cell = sc.pp.filter_cells(
        adata,
        min_genes=min_genes,
        inplace=False
    )

    return adata[cell_subset, :].copy()

def _filter_genes(adata: sc.AnnData, min_cells: int = 3) -> sc.AnnData:
    """
    Filter genes in an AnnData object based on quality control metrics.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to filter.
    - min_cells: int, optional
        Minimum number of cells expressing a gene to keep it. Default is 3.

    Returns:
    - filtered_adata: sc.AnnData
        The filtered AnnData object.
    """


    gene_subset, number_per_gene = sc.pp.filter_genes(
        adata,
        min_cells=min_cells,
        inplace=False
    )

    return adata[:, gene_subset].copy()


def _filter_mitochondrial_genes(adata: sc.AnnData, max_mt_percentage: float = 20) -> sc.AnnData:
    """
    Filter cells in an AnnData object based on the fraction of mitochondrial genes.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to filter.
    - max_mt_percentage: float, optional
        Maximum percentage of mitochondrial genes allowed in a cell. Default is 20%.

    Returns:
    - filtered_adata: sc.AnnData
        The filtered AnnData object.
    """
    
    # Filter cells based on the specified criteria
    filtered_adata = adata[adata.obs["pct_counts_mt"] < max_mt_percentage].copy()
    
    return filtered_adata


def _identify_doublets(adata: sc.AnnData, expected_doublet_rate: float = 0.05) -> sc.AnnData:
    """
    Identify potential doublets in an AnnData object using Scrublet.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to analyze.
    - expected_doublet_rate: float, optional
        Expected doublet rate in the dataset. Default is 0.1 (10%).

    Returns:
    - adata_with_doublets: sc.AnnData
        The AnnData object with an additional column indicating potential doublets.
    """
    
    sc.pp.scrublet(adata, expected_doublet_rate=expected_doublet_rate)
    
    return adata

def data_filtering(adata, 
           min_genes: int = 200, 
           min_cells: int = 3, 
           max_mt_percentage: float = 20, 
           expected_doublet_rate: float = 0.05) -> sc.AnnData:
    """
    Filter an AnnData object based on quality control metrics.

    Parameters:
    - adata: sc.AnnData
        The AnnData object to filter.
    - min_genes: int, optional
        Minimum number of genes expressed in a cell to keep it. Default is 200.
    - min_cells: int, optional
        Minimum number of cells expressing a gene to keep it. Default is 3.
    - max_mt_percentage: float, optional
        Maximum percentage of mitochondrial genes allowed in a cell. Default is 20%.
    - expected_doublet_rate: float, optional
        Expected doublet rate in the dataset. Default is 0.05 (5%).

    Returns:
    - filtered_adata: sc.AnnData
        The filtered AnnData object.
    """
    
    # Filter cells based on the specified criteria
    adata = _filter_cells(adata, min_genes)
    
    # Filter genes based on the specified criteria
    adata = _filter_genes(adata, min_cells)
    
    # Filter cells based on mitochondrial gene percentage
    adata = _filter_mitochondrial_genes(adata, max_mt_percentage)
    
    # Identify potential doublets in the dataset
    adata = _identify_doublets(adata, expected_doublet_rate)
    
    return adata