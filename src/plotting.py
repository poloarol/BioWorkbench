import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import squidpy as sq
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


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

def plot_qc_metrics(
    adata: sc.AnnData,
    is_spatial: bool = False,
) -> plt.Figure:
    """
    Plot quality control metrics for an AnnData object.

    Parameters
    ----------
    adata : sc.AnnData
        The AnnData object containing gene expression data.
    is_spatial : bool
        Whether the data are spatial transcriptomics data.

    Returns
    -------
    plt.Figure
        The matplotlib figure containing the QC metrics plots.
    """

    if is_spatial:
        fig, axs = plt.subplots(1, 3, figsize=(10, 4))

        axs[0].set_title("Total transcripts per cell")
        sns.histplot(
            adata.obs["total_counts"],
            kde=False,
            ax=axs[0],
        )

        axs[1].set_title("Unique transcripts per cell")
        sns.histplot(
            adata.obs["n_genes_by_counts"],
            kde=False,
            ax=axs[1],
        )

        # axs[2].set_title("Transcripts per FOV")
        # transcripts_per_fov = (
        #     adata.obs.groupby("fov")["total_counts"].sum()
        # )
        # sns.histplot(
        #     transcripts_per_fov,
        #     kde=False,
        #     ax=axs[2],
        # )

        axs[2].set_title("Volume of segmented cells")
        sns.histplot(
            adata.obs["volume"],
            kde=False,
            ax=axs[2],
        )

        fig.tight_layout()
        return fig

    # Non-spatial QC
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )

    fig = plt.gcf()
    fig.tight_layout()
    return fig


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


def plot_spatial(
    adata,
    color_by,
    figsize=(18, 6),
):
    """Plot spatial representation of the data.

    Args:
        adata (sc.AnnData): Annotated data matrix.
        color_by (str or list of str): Column(s) in `adata.obs` to color the points by.
        figsize (tuple, optional): Figure size. Defaults to (18, 6).

    Returns:
        matplotlib.figure.Figure: The resulting figure.
    """
    if isinstance(color_by, str):
        color_by = [color_by]

    n_plots = len(color_by)

    fig, axes = plt.subplots(
        1,
        n_plots,
        figsize=figsize,
        squeeze=False,
    )

    axes = axes[0]

    for ax, variable in zip(axes, color_by):

        values = adata.obs[variable]

        # Categorical
        if (
            pd.api.types.is_categorical_dtype(values)
            or pd.api.types.is_object_dtype(values)
            or pd.api.types.is_bool_dtype(values)
        ):
            categories = values.astype("category")
            category_names = list(categories.cat.categories)
            n_categories = len(category_names)

            if n_categories == 2:
                color_map = {
                    category_names[0]: "red",
                    category_names[1]: "blue",
                }

            else:
                cmap = plt.get_cmap("tab20", n_categories)

                color_map = {
                    category: cmap(i)
                    for i, category in enumerate(category_names)
                }

            # Avoid pandas Series.map() MultiIndex issue
            colors = [
                color_map.get(category, "gray")
                for category in categories
            ]

            ax.scatter(
                adata.obs["center_x"],
                -adata.obs["center_y"],
                c=colors,
                s=1,
                alpha=0.7,
            )

            # Categorical legend
            legend_handles = [
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="none",
                    markerfacecolor=color_map[category],
                    markeredgecolor="none",
                    markersize=7,
                    label=str(category),
                )
                for category in category_names
            ]

            ax.legend(
                handles=legend_handles,
                title=variable.replace("_", " ").title(),
                loc="upper right",
                frameon=True,
            )

        # Continuous
        else:
            scatter = ax.scatter(
                adata.obs["center_x"],
                -adata.obs["center_y"],
                c=values,
                s=1,
                alpha=0.7,
                cmap="viridis",
            )

            fig.colorbar(
                scatter,
                ax=ax,
                fraction=0.046,
                pad=0.04,
                label=variable.replace("_", " ").title(),
            )
        
        ax.set_xlabel("Spatial X")
        ax.set_ylabel("Spatial Y")
        ax.set_title(f"Spatial — {variable.replace('_', ' ').title()}")
        
        
        ax.set_aspect("equal", adjustable="box")
        ax.margins(0)
        ax.axis("off")

    fig.tight_layout()

    return fig


def plot_umap(
    adata: sc.AnnData,
    color_by: str,
    figsize: tuple = (8, 6),
):
    """Plot UMAP colored by a categorical or continuous obs column.
    Args:
        adata (sc.AnnData): Annotated data matrix.
        color_by (str): Column in `adata.obs` to color the points by.
        figsize (tuple, optional): Figure size. Defaults to (8, 6).

    Returns:
        matplotlib.figure.Figure: The resulting figure.
    """

    fig, ax = plt.subplots(figsize=figsize)

    values = adata.obs[color_by]

    # Categorical variable
    if (
        pd.api.types.is_categorical_dtype(values)
        or pd.api.types.is_object_dtype(values)
        or pd.api.types.is_bool_dtype(values)
    ):
        categories = values.astype("category")
        category_names = categories.cat.categories
        n_categories = len(category_names)

        # Binary categorical variable
        if n_categories == 2:
            color_map = {
                category_names[0]: "red",
                category_names[1]: "blue",
            }

        # Multi-class categorical variable
        else:
            palette = sns.color_palette(
                "tab20",
                n_colors=n_categories,
            )

            color_map = {
                category: color
                for category, color in zip(
                    category_names,
                    palette,
                )
            }

        # Plot each category
        for category in category_names:
            mask = categories == category

            ax.scatter(
                adata.obsm["X_umap"][mask, 0],
                adata.obsm["X_umap"][mask, 1],
                s=5,
                alpha=0.7,
                color=color_map[category],
                label=str(category),
            )

        # Categorical legend
        ax.legend(
            title=color_by.replace("_", " ").title(),
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            markerscale=2,
        )

    # Continuous variable
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
            label=color_by.replace("_", " ").title(),
        )

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(f"UMAP — {color_by.replace('_', ' ').title()}")

    fig.tight_layout()

    return fig


def plot_umap_genes(
    adata,
    genes,
    figsize=(8, 6),
):
    """Plot UMAP expression of multiple genes.

    Args:
        adata (sc.AnnData): Annotated data matrix containing UMAP
            coordinates in `adata.obsm["X_umap"]`.
        genes (list[str]): List of gene names to plot.
        figsize (tuple, optional): Base figure size. Defaults to (8, 6).

    Returns:
        matplotlib.figure.Figure: The resulting figure.
    """

    if isinstance(genes, str):
        genes = [genes]

    if len(genes) == 0:
        raise ValueError("No genes were provided for UMAP plotting.")

    # Validate UMAP coordinates
    if "X_umap" not in adata.obsm:
        raise ValueError(
            "UMAP coordinates not found in adata.obsm['X_umap']."
        )

    # Validate genes
    missing_genes = [
        gene for gene in genes
        if gene not in adata.var_names
    ]

    if missing_genes:
        raise ValueError(
            f"Gene(s) not found in adata.var_names: {missing_genes}"
        )

    umap = adata.obsm["X_umap"]

    n_genes = len(genes)

    fig, axes = plt.subplots(
        1,
        n_genes,
        figsize=(figsize[0] * n_genes, figsize[1]),
        squeeze=False,
    )

    axes = axes[0]

    for ax, gene in zip(axes, genes):

        # Extract expression
        values = np.log1p(adata[:, gene].X)

        # Handle sparse matrices
        if hasattr(values, "toarray"):
            values = values.toarray()

        values = np.asarray(values).flatten()

        # Plot
        scatter = ax.scatter(
            umap[:, 0],
            umap[:, 1],
            c=values,
            s=1,
            alpha=0.7,
            cmap="viridis",
        )

        fig.colorbar(
            scatter,
            ax=ax,
            fraction=0.046,
            pad=0.04,
            label="Expression (log1p)",
        )

        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(f"UMAP — {gene}")

        ax.set_aspect("equal", adjustable="box")
        ax.margins(0.02)

    fig.tight_layout()

    return fig


def plot_spatial_genes(
    adata,
    genes,
    figsize=(8, 6),
):
    """Plot spatial expression of multiple genes.

    Args:
        adata (sc.AnnData): Annotated data matrix containing spatial
            coordinates in `adata.obs["center_x"]` and
            `adata.obs["center_y"]`.
        genes (list[str]): List of gene names to plot.
        figsize (tuple, optional): Base figure size. Defaults to (8, 6).

    Returns:
        matplotlib.figure.Figure: The resulting figure.
    """

    if isinstance(genes, str):
        genes = [genes]

    # Validate spatial coordinates
    required_coordinates = {"center_x", "center_y"}

    missing_coordinates = required_coordinates - set(adata.obs.columns)

    if missing_coordinates:
        raise ValueError(
            "Missing required spatial coordinate columns: "
            f"{sorted(missing_coordinates)}"
        )

    # Validate genes
    missing_genes = [gene for gene in genes if gene not in adata.var_names]

    if missing_genes:
        raise ValueError(
            f"Gene(s) not found in adata.var_names: {missing_genes}"
        )

    n_genes = len(genes)

    if n_genes == 0:
        raise ValueError("No genes were provided for spatial plotting.")

    fig, axes = plt.subplots(
        1,
        n_genes,
        figsize=(figsize[0] * n_genes, figsize[1]),
        squeeze=False,
    )

    axes = axes[0]

    for ax, gene in zip(axes, genes):

        # Extract expression values
        values = np.log1p(adata[:, gene].X)

        # Convert sparse matrix / 2D array to 1D numpy array
        if hasattr(values, "toarray"):
            values = values.toarray()

        values = np.asarray(values).flatten()

        scatter = ax.scatter(
            adata.obs["center_x"],
            -adata.obs["center_y"],
            c=values,
            s=1,
            alpha=0.7,
            cmap="viridis",
        )

        fig.colorbar(
            scatter,
            ax=ax,
            fraction=0.046,
            pad=0.04,
            label="Expression (log1p)",
        )

        ax.set_xlabel("Spatial X")
        ax.set_ylabel("Spatial Y")
        ax.set_title(
            f"Spatial — {gene}"
        )

        ax.set_aspect("equal", adjustable="box")
        ax.margins(0)
        ax.axis("off")

    fig.tight_layout()

    return fig