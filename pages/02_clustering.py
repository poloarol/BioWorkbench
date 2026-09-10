import io
import os
import tempfile

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
import streamlit as st

from src.clustering import run_clustering, run_pca, identify_marker_genes
from src.plotting import plot_pca, plot_umap


# -----------------------------------------------------------------------------
# Sidebar
# -----------------------------------------------------------------------------

with st.sidebar:
    st.title("Processing parameters")

    st.subheader("PCA")

    n_top_genes = st.number_input(
        "Highly variable genes",
        min_value=100,
        max_value=10000,
        value=2000,
        step=100,
    )

    n_comps = st.number_input(
        "Number of PCs",
        min_value=2,
        max_value=100,
        value=50,
        step=1,
    )

    st.subheader("UMAP / Clustering")

    n_neighbors = st.number_input(
        "Neighbors",
        min_value=2,
        max_value=100,
        value=15,
        step=1,
    )

    min_dist = st.slider(
        "UMAP min_dist",
        min_value=0.0,
        max_value=1.0,
        value=0.1,
        step=0.05,
    )

    resolution = st.slider(
        "Leiden resolution",
        min_value=0.1,
        max_value=3.0,
        value=1.0,
        step=0.1,
    )

# -----------------------------------------------------------------------------
# Dataset overview
# -----------------------------------------------------------------------------

adata = st.session_state.adatas['filtered_wout_doublets']

# -----------------------------------------------------------------------------
# Pipeline controls
# -----------------------------------------------------------------------------

st.divider()

st.subheader("Analysis")

col1, col2 = st.columns(2)

with col1:
    if st.button(
        "▶ Run PCA",
        type="primary",
        use_container_width=True,
    ):
        with st.spinner("Running normalization, HVG selection and PCA..."):
            try:
                st.session_state.adata = run_pca(
                    adata=adata,
                    n_top_genes=n_top_genes,
                    n_comps=n_comps,
                )

                st.success("PCA completed.")

            except Exception as exc:
                st.error(f"PCA failed: {exc}")

with col2:
    if st.button(
        "▶ Run UMAP + Leiden",
        use_container_width=True,
    ):
        if "X_pca" not in st.session_state.adata.obsm:
            st.warning("Run PCA first.")
        else:
            with st.spinner("Running neighbors, UMAP and Leiden..."):
                try:
                    st.session_state.adata = run_clustering(
                        adata=adata,
                        n_neighbors=n_neighbors,
                        min_dist=min_dist,
                        resolution=resolution,
                    )

                    st.success("UMAP and clustering completed.")

                except Exception as exc:
                    st.error(f"Clustering failed: {exc}")

# -----------------------------------------------------------------------------
# Results
# -----------------------------------------------------------------------------

if "X_pca" in adata.obsm or "X_umap" in adata.obsm:

    st.divider()

    col_pca, col_umap = st.columns(2)

    # -------------------------------------------------------------------------
    # PCA
    # -------------------------------------------------------------------------

    with col_pca:

        if "X_pca" in adata.obsm:

            st.subheader("PCA")

            pca_color_by = st.selectbox(
                "Color cells by",
                options=list(adata.obs.columns),
                key="pca_color_by",
            )

            fig = plot_pca(
                adata,
                color_by=pca_color_by,
            )

            st.pyplot(
                fig,
                use_container_width=True,
            )

            plt.close(fig)

        else:
            st.info("PCA has not been run yet.")

    # -------------------------------------------------------------------------
    # UMAP
    # -------------------------------------------------------------------------

    with col_umap:

        if "X_umap" in adata.obsm:

            st.subheader("UMAP")

            umap_color_by = st.selectbox(
                "Color cells by",
                options=list(adata.obs.columns),
                key="umap_color_by",
            )

            fig = plot_umap(
                adata,
                color_by=umap_color_by,
            )

            st.pyplot(
                fig,
                use_container_width=True,
            )

            plt.close(fig)

        else:
            st.info("UMAP has not been run yet.")


# -----------------------------------------------------------------------------
# Highly variable genes
# -----------------------------------------------------------------------------

if "highly_variable" in adata.var:

    st.divider()

    st.subheader("Highly Variable Genes")

    hvg = adata.var[
        adata.var["highly_variable"]
    ].copy()

    st.write(
        f"Selected **{len(hvg):,}** highly variable genes."
    )

    st.dataframe(
        hvg,
        use_container_width=True,
    )


adata = identify_marker_genes(adata, groupby="cluster_label")

# -----------------------------------------------------------------------------
# Marker genes
# -----------------------------------------------------------------------------

if "rank_genes_groups" in adata.uns:

    st.divider()

    st.subheader("Top 10 Marker Genes")

    marker_df = sc.get.rank_genes_groups_df(
        adata,
        group=None,
    )

    top_markers = (
        marker_df
        .sort_values(
            ["group", "scores"],
            ascending=[True, False],
        )
        .groupby("group", sort=False)
        .head(10)
        .reset_index(drop=True)
    )

    # -------------------------------------------------------------------------
    # Format p-values
    # -------------------------------------------------------------------------

    def format_pvalue(value):
        if pd.isna(value):
            return "—"

        if value == 0:
            return "<1e-300"

        if value < 0.001:
            return f"{value:.2e}"

        return f"{value:.4f}"

    top_markers["pvals_adj"] = (
        top_markers["pvals_adj"]
        .apply(format_pvalue)
    )

    top_markers["pvals"] = (
        top_markers["pvals"]
        .apply(format_pvalue)
    )

    display_columns = [
        "group",
        "names",
        "scores",
        "logfoldchanges",
        "pvals_adj",
    ]

    st.dataframe(
        top_markers[display_columns],
        use_container_width=True,
        hide_index=True,
    )


# -----------------------------------------------------------------------------
# Cluster summary
# -----------------------------------------------------------------------------

if "cluster_label" in adata.obs:

    st.divider()

    st.subheader("Cluster Summary")

    cluster_counts = (
        adata.obs["cluster_label"]
        .value_counts()
        .sort_index()
        .rename("cells")
        .to_frame()
    )

    cluster_counts["percentage"] = (
        cluster_counts["cells"]
        / len(adata)
        * 100
    ).round(2)

    st.dataframe(
        cluster_counts,
        use_container_width=True,
    )


st.session_state.adatas['clustered'] = adata

# -----------------------------------------------------------------------------
# Download
# -----------------------------------------------------------------------------

st.divider()

st.subheader("Export")

download_col1, download_col2 = st.columns(2)


# -----------------------------------------------------------------------------
# Download processed AnnData
# -----------------------------------------------------------------------------

with download_col1:

    with tempfile.NamedTemporaryFile(
        suffix=".h5ad",
        delete=False,
    ) as tmp:
        tmp_path = tmp.name

    try:
        adata.write_h5ad(tmp_path)

        with open(tmp_path, "rb") as f:
            h5ad_data = f.read()

        st.download_button(
            label="⬇ Download processed AnnData",
            data=h5ad_data,
            file_name="processed.h5ad",
            mime="application/octet-stream",
            use_container_width=True,
        )

    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


# -----------------------------------------------------------------------------
# Download marker genes
# -----------------------------------------------------------------------------

with download_col2:

    if "rank_genes_groups" in adata.uns:

        marker_df = sc.get.rank_genes_groups_df(
            adata,
            group=None,
        )

        marker_csv = marker_df.to_csv(
            index=False
        ).encode("utf-8")

        st.download_button(
            label="⬇ Download marker genes",
            data=marker_csv,
            file_name="marker_genes.csv",
            mime="text/csv",
            use_container_width=True,
        )

    else:

        st.button(
            "⬇ Download marker genes",
            disabled=True,
            use_container_width=True,
        )

        st.caption(
            "Run marker-gene analysis first."
        )