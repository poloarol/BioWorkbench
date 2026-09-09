
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt


def highest_expressed_genes(adata: sc.AnnData, n_top_genes: int = 10) -> plt.Figure:
    """
    Plot the highest expressed genes in an AnnData object.

    Parameters:
    - adata: sc.AnnData
        The AnnData object containing gene expression data.
    - n_top_genes: int, optional
        Number of top expressed genes to display. Default is 10.
    Returns:
    - fig: plt.Figure
        The matplotlib figure object containing the plot.
    """

    return sc.pl.highest_expr_genes(adata, n_top=n_top_genes, show=False)

def plot_qc_metrics(adata: sc.AnnData) -> plt.Figure:
    """
    Plot quality control metrics for an AnnData object.

    Parameters:
    - adata: sc.AnnData
        The AnnData object containing gene expression data.
    
    Returns:
    - fig: plt.Figure
        The matplotlib figure object containing the QC metrics plot.
    """

    return sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            show=False
    )