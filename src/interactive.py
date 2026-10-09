"""Interactive (Plotly) embeddings with group selection and statistics."""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

DIM_COLOR = "rgba(200,200,200,0.25)"
# Plot height follows the viewport; width already follows the column.
RESPONSIVE_CSS = """
<style>
div[data-testid="stPlotlyChart"],
div[data-testid="stPlotlyChart"] .js-plotly-plot,
div[data-testid="stPlotlyChart"] .plot-container {
    height: clamp(340px, 70vh, 900px) !important;
}
</style>
"""
QC_COLUMNS = ("total_counts", "n_genes_by_counts", "pct_counts_mt", "volume")


def is_categorical(values: pd.Series) -> bool:
    return (
        isinstance(values.dtype, pd.CategoricalDtype)
        or pd.api.types.is_object_dtype(values)
        or pd.api.types.is_string_dtype(values)
        or pd.api.types.is_bool_dtype(values)
    )


def get_coordinates(adata, basis: str) -> pd.DataFrame:
    """Return a DataFrame with x/y columns for 'umap', 'pca' or 'spatial'."""
    if basis in ("umap", "pca"):
        key = f"X_{basis}"
        if key not in adata.obsm:
            raise ValueError(f"{key} not found in adata.obsm.")
        coords = np.asarray(adata.obsm[key])[:, :2]
    elif basis == "spatial":
        missing = {"center_x", "center_y"} - set(adata.obs.columns)
        if missing:
            raise ValueError(
                f"Missing required spatial coordinate columns: {sorted(missing)}"
            )
        coords = np.column_stack(
            [adata.obs["center_x"].to_numpy(), -adata.obs["center_y"].to_numpy()]
        )
    else:
        raise ValueError(f"Unknown basis: {basis!r}")
    return pd.DataFrame(coords, columns=["x", "y"], index=adata.obs_names)


def gene_expression(adata, gene: str) -> np.ndarray:
    """log1p expression of one gene as a dense 1D array."""
    if gene not in adata.var_names:
        raise ValueError(f"Gene(s) not found in adata.var_names: {[gene]}")
    values = np.log1p(adata[:, gene].X)
    if hasattr(values, "toarray"):
        values = values.toarray()
    return np.asarray(values).ravel()


def _axis_labels(basis: str) -> tuple:
    return {
        "umap": ("UMAP 1", "UMAP 2"),
        "pca": ("PC1", "PC2"),
        "spatial": ("Spatial X", "Spatial Y"),
    }[basis]


def interactive_plot(
    adata,
    basis: str,
    color_by: str = None,
    gene: str = None,
    highlight=None,
    point_size: float = None,
    height: int = None,
) -> go.Figure:
    """Zoomable/pannable scatter of cells colored by an obs column or a gene.

    Args:
        adata: Annotated data matrix.
        basis: 'umap', 'pca' or 'spatial'.
        color_by: obs column to color by (categorical or numeric).
        gene: Gene to color by (log1p expression). Takes precedence over color_by.
        highlight: Categories of a categorical `color_by` to emphasise; all
            other cells are dimmed. None or empty shows everything.
        point_size: Marker size; chosen from the cell count when None.
        height: Fixed figure height in pixels; None (default) sizes to the viewport.

    Every point carries its integer row position in `customdata`, so selections
    made in the browser can be mapped back to cells.
    """
    if gene is None and color_by is None:
        raise ValueError("Provide either color_by or gene.")

    coords = get_coordinates(adata, basis)
    n = len(coords)
    size = point_size if point_size is not None else (2 if n > 50_000 else 4 if n > 5_000 else 6)
    positions = np.arange(n)
    xlabel, ylabel = _axis_labels(basis)

    fig = go.Figure()
    hover = "%{customdata[1]}<extra></extra>"
    names = adata.obs_names.to_numpy()

    if gene is not None:
        values = gene_expression(adata, gene)
        fig.add_trace(
            go.Scattergl(
                x=coords["x"],
                y=coords["y"],
                mode="markers",
                marker=dict(
                    size=size,
                    color=values,
                    colorscale="Viridis",
                    colorbar=dict(title="Expression (log1p)"),
                    opacity=0.8,
                ),
                customdata=np.column_stack([positions, names]),
                hovertemplate=f"%{{customdata[1]}}<br>{gene}: %{{marker.color:.2f}}<extra></extra>",
                showlegend=False,
            )
        )
        title = f"{basis.upper()} — {gene}"
    else:
        values = adata.obs[color_by]
        title = f"{basis.upper()} — {color_by.replace('_', ' ').title()}"
        if is_categorical(values):
            categories = values.astype("category")
            category_names = [str(c) for c in categories.cat.categories]
            labels = categories.astype(str).to_numpy()
            palette = px.colors.qualitative.Alphabet + px.colors.qualitative.Dark24
            highlight = {str(h) for h in (highlight or [])}
            for i, name in enumerate(category_names):
                mask = labels == name
                dimmed = bool(highlight) and name not in highlight
                fig.add_trace(
                    go.Scattergl(
                        x=coords["x"].to_numpy()[mask],
                        y=coords["y"].to_numpy()[mask],
                        mode="markers",
                        name=f"{name} ({mask.sum():,})",
                        marker=dict(
                            size=size,
                            color=DIM_COLOR if dimmed else palette[i % len(palette)],
                            opacity=1 if dimmed else 0.85,
                        ),
                        customdata=np.column_stack([positions[mask], names[mask]]),
                        hovertemplate=f"%{{customdata[1]}}<br>{name}<extra></extra>",
                    )
                )
            fig.update_layout(legend=dict(title=color_by.replace("_", " ").title(), itemsizing="constant"))
        else:
            fig.add_trace(
                go.Scattergl(
                    x=coords["x"],
                    y=coords["y"],
                    mode="markers",
                    marker=dict(
                        size=size,
                        color=values.to_numpy(dtype=float),
                        colorscale="Viridis",
                        colorbar=dict(title=color_by.replace("_", " ").title()),
                        opacity=0.8,
                    ),
                    customdata=np.column_stack([positions, names]),
                    hovertemplate=f"%{{customdata[1]}}<br>{color_by}: %{{marker.color:.3g}}<extra></extra>",
                    showlegend=False,
                )
            )

    fig.update_layout(
        title=title,
        height=height,
        margin=dict(l=10, r=10, t=50, b=10),
        dragmode="zoom",
        xaxis=dict(title=xlabel),
        yaxis=dict(title=ylabel, scaleanchor="x" if basis == "spatial" else None),
        template="plotly_white",
    )
    return fig


def selected_positions(event) -> np.ndarray:
    """Integer cell positions from a Streamlit plotly selection event."""
    if event is None:
        return np.array([], dtype=int)
    selection = event.get("selection", {}) if hasattr(event, "get") else {}
    points = selection.get("points", [])
    positions = []
    for point in points:
        custom = point.get("customdata")
        if custom:
            positions.append(int(custom[0]))
    return np.unique(np.asarray(positions, dtype=int))


def group_statistics(
    adata,
    group_col: str,
    groups=None,
    gene: str = None,
) -> pd.DataFrame:
    """Per-group summary: cell count, fraction, QC means and optional gene stats.

    Args:
        adata: Annotated data matrix.
        group_col: Categorical obs column defining the groups.
        groups: Groups to include; all groups when None/empty.
        gene: Optional gene to add mean log1p expression and % expressing.
    """
    labels = adata.obs[group_col].astype(str)
    if groups:
        wanted = [str(g) for g in groups]
    else:
        wanted = [str(c) for c in adata.obs[group_col].astype("category").cat.categories]

    qc_cols = [c for c in QC_COLUMNS if c in adata.obs.columns]
    expr = gene_expression(adata, gene) if gene is not None else None
    has_xy = {"center_x", "center_y"} <= set(adata.obs.columns)

    rows = []
    for group in wanted:
        mask = (labels == group).to_numpy()
        row = {
            group_col: group,
            "n_cells": int(mask.sum()),
            "pct_of_total": round(100 * mask.sum() / max(len(labels), 1), 2),
        }
        for col in qc_cols:
            col_values = adata.obs.loc[mask, col]
            row[f"{col} (mean)"] = round(float(col_values.mean()), 3) if mask.any() else np.nan
            row[f"{col} (median)"] = round(float(col_values.median()), 3) if mask.any() else np.nan
        if expr is not None and mask.any():
            row[f"{gene} mean (log1p)"] = round(float(expr[mask].mean()), 3)
            row[f"{gene} % expressing"] = round(100 * float((expr[mask] > 0).mean()), 2)
        if has_xy and mask.any():
            row["centroid_x"] = round(float(adata.obs.loc[mask, "center_x"].mean()), 2)
            row["centroid_y"] = round(float(adata.obs.loc[mask, "center_y"].mean()), 2)
        rows.append(row)

    return pd.DataFrame(rows)


def selection_statistics(adata, positions, group_col: str = None, gene: str = None) -> pd.DataFrame:
    """Summary of an arbitrary set of cells (e.g. a lasso selection).

    Returns one row for the whole selection, followed by one row per group of
    `group_col` present within the selection.
    """
    positions = np.asarray(positions, dtype=int)
    if positions.size == 0:
        return pd.DataFrame()

    subset = adata[positions]
    rows = [{"selection": "all selected", **_summary_row(subset, adata.n_obs, gene)}]
    if group_col is not None and group_col in subset.obs:
        labels = subset.obs[group_col].astype(str)
        for group in sorted(labels.unique()):
            rows.append(
                {
                    "selection": f"{group_col} = {group}",
                    **_summary_row(subset[(labels == group).to_numpy()], adata.n_obs, gene),
                }
            )
    return pd.DataFrame(rows)


def render_interactive(
    adata,
    basis: str,
    key: str,
    color_by: str = None,
    gene: str = None,
    group_by: str = None,
    height: int = None,
) -> None:
    """Render an interactive plot with a group picker and statistics in Streamlit.

    Zoom with the scroll wheel / box-zoom, pan by dragging, and use the box or
    lasso tool to select cells. Statistics are shown for the groups picked in
    the multiselect and, separately, for any box/lasso selection.

    When coloring by `color_by` (categorical), the picked groups are highlighted.
    When coloring by `gene`, `group_by` provides the groups; picking some
    restricts the plot to those cells.
    """
    import streamlit as st

    st.html(RESPONSIVE_CSS)

    if color_by is not None and is_categorical(adata.obs[color_by]):
        group_col = color_by
    else:
        group_col = group_by
    highlight = []
    if group_col is not None:
        options = [str(c) for c in adata.obs[group_col].astype("category").cat.categories]
        highlight = st.multiselect(
            (f"Highlight {group_col.replace('_', ' ')}" if group_col == color_by
             else f"Show only {group_col.replace('_', ' ')}"),
            options=options,
            key=f"{key}_highlight",
            placeholder="All (select one or more to focus)",
        )

    plotted = adata
    if gene is not None and highlight:
        plotted = adata[adata.obs[group_col].astype(str).isin(highlight).to_numpy()]

    fig = interactive_plot(
        plotted,
        basis=basis,
        color_by=color_by,
        gene=gene,
        highlight=highlight if group_col == color_by else None,
        height=height,
    )
    event = st.plotly_chart(
        fig,
        key=f"{key}_chart",
        on_select="rerun",
        selection_mode=("box", "lasso"),
        config={"scrollZoom": True, "displaylogo": False},
        width="stretch",
    )

    if group_col is not None:
        stats = group_statistics(adata, group_col, groups=highlight, gene=gene)
        with st.expander(
            "Group statistics" + (f" ({len(highlight)} selected)" if highlight else " (all)"),
            expanded=bool(highlight),
        ):
            st.dataframe(stats, hide_index=True, width="stretch")

    positions = selected_positions(event)
    if positions.size:
        st.caption(f"{positions.size:,} cells selected on the plot")
        st.dataframe(
            selection_statistics(plotted, positions, group_col=group_col, gene=gene),
            hide_index=True,
            width="stretch",
        )


def _summary_row(subset, total: int, gene: str = None) -> dict:
    row = {
        "n_cells": int(subset.n_obs),
        "pct_of_total": round(100 * subset.n_obs / max(total, 1), 2),
    }
    for col in QC_COLUMNS:
        if col in subset.obs.columns:
            row[f"{col} (mean)"] = round(float(subset.obs[col].mean()), 3)
    if gene is not None:
        expr = gene_expression(subset, gene)
        row[f"{gene} mean (log1p)"] = round(float(expr.mean()), 3)
        row[f"{gene} % expressing"] = round(100 * float((expr > 0).mean()), 2)
    return row
