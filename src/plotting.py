
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

    fig, ax = plt.subplots(figsize=(6, 4))

    sc.pl.highest_expr_genes(
        adata,
        n_top=n_top_genes,
        ax=ax,
        show=False,
    )

    return fig

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
            show=False,
    )

def qc_summary(adata):
    n_genes = adata.obs["n_genes_by_counts"]
    total_counts = adata.obs["total_counts"]

    return {
        "n_cells": adata.n_obs,
        "median_genes": n_genes.median(),
        "mean_genes": n_genes.mean(),
        "median_counts": total_counts.median(),
        "mean_counts": total_counts.mean(),
    }
