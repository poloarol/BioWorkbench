
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

# -----------------------------------------------------------------------------
# Visualization helpers
# -----------------------------------------------------------------------------

def plot_pca(
    adata: sc.AnnData,
    color_by: str,
):
    """Plot the first two principal components."""

    fig, ax = plt.subplots(figsize=(8, 6))

    values = adata.obs[color_by]

    # -------------------------------------------------------------------------
    # Categorical data
    # -------------------------------------------------------------------------

    if (
        values.dtype.name == "category"
        or values.dtype == object
    ):
        groups = values.astype(str)

        palette = sns.color_palette(
            "tab20",
            n_colors=groups.nunique(),
        )

        for color, group in zip(
            palette,
            sorted(groups.unique()),
        ):
            mask = groups == group

            ax.scatter(
                adata.obsm["X_pca"][mask, 0],
                adata.obsm["X_pca"][mask, 1],
                s=5,
                alpha=0.7,
                color=color,
                label=group,
            )

        ax.legend(
            title=color_by,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            markerscale=2,
        )

    # -------------------------------------------------------------------------
    # Numeric data
    # -------------------------------------------------------------------------

    else:
        scatter = ax.scatter(
            adata.obsm["X_pca"][:, 0],
            adata.obsm["X_pca"][:, 1],
            c=values,
            s=5,
            alpha=0.7,
            cmap="viridis",
        )

        fig.colorbar(
            scatter,
            ax=ax,
            label=color_by,
        )

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"PCA — {color_by}")

    return fig



def plot_umap(
    adata: sc.AnnData,
    color_by: str,
):
    """Plot UMAP colored by an obs column."""

    fig, ax = plt.subplots(figsize=(8, 6))

    values = adata.obs[color_by]

    if values.dtype.name == "category" or values.dtype == object:
        categories = values.astype(str)
        unique_categories = sorted(categories.unique())

        palette = sns.color_palette(
            "tab20",
            n_colors=len(unique_categories),
        )

        for color, category in zip(palette, unique_categories):
            mask = categories == category

            ax.scatter(
                adata.obsm["X_umap"][mask, 0],
                adata.obsm["X_umap"][mask, 1],
                s=5,
                alpha=0.7,
                color=color,
                label=category,
            )

        ax.legend(
            title=color_by,
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            markerscale=2,
        )

    else:
        scatter = ax.scatter(
            adata.obsm["X_umap"][:, 0],
            adata.obsm["X_umap"][:, 1],
            c=values,
            s=5,
            alpha=0.7,
            cmap="viridis",
        )

        fig.colorbar(
            scatter,
            ax=ax,
            label=color_by,
        )

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(f"UMAP — {color_by}")

    return fig
