
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import streamlit as st

if "is_spatial" not in st.session_state:
    st.session_state["is_spatial"] = False

from src.clustering import run_clustering, run_pca, identify_marker_genes, identify_spatial_domains, run_spatially_variable_genes, RANDOM_STATE
from src.interactive import render_interactive, render_hvg, render_svg
from src.utils import color_options


# -----------------------------------------------------------------------------
# Sidebar
# -----------------------------------------------------------------------------

alpha = 0.2


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

    if st.session_state['is_spatial']:
        st.subheader("Spatial Domains")

        alpha = st.slider(
            "Alpha (weight of spatial graph)",
            min_value=0.0,
            max_value=1.0,
            value=0.2,
            step=0.05,
        )

# -----------------------------------------------------------------------------
# Dataset overview
# -----------------------------------------------------------------------------

adata = None

filtered_key = (
    "filtered_wout_blanks"
    if st.session_state.get("is_spatial", False)
    else "filtered_wout_doublets"
)
filtered = st.session_state.get("adatas", {}).get(filtered_key)

if filtered is None:
    st.info("Apply filtering on the Cell and Gene Filtering page first.")
    st.stop()

# Downstream steps modify AnnData in place, so work on a copy that is rebuilt
# whenever the filtered object changes. Raw counts stay in layers["raw"].
work = st.session_state.get("clustering_work")
if work is None or st.session_state.get("clustering_source") is not filtered:
    restored = st.session_state.get("adatas", {}).get("clustered")
    if restored is not None and st.session_state.get("clustering_source") is None:
        work = restored
    else:
        work = filtered.copy()
    st.session_state["clustering_work"] = work
    st.session_state["clustering_source"] = filtered

adata = work

if st.session_state.get("is_spatial", False):
    col1, col2, col3 = st.columns(3)
else:
    col1, col2 = st.columns(2)

# -----------------------------------------------------------------------------
# Pipeline controls
# -----------------------------------------------------------------------------


with col1:
    if st.button(
        "▶ Run PCA",
        type="primary",
        use_container_width=True,
    ):
        with st.spinner("Running normalization, HVG selection and PCA..."):
            try:
                adata = run_pca(
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
        type="primary",
        use_container_width=True,
    ):
        if "X_pca" not in adata.obsm:
            st.warning("Run PCA first.")
        else:
            with st.spinner("Running neighbors, UMAP and Leiden..."):
                try:
                    adata = run_clustering(
                        adata=adata,
                        n_neighbors=n_neighbors,
                        min_dist=min_dist,
                        resolution=resolution,
                    )

                    st.success("UMAP and clustering completed.")

                except Exception as exc:
                    st.error(f"Clustering failed: {exc}")

if st.session_state['is_spatial']:
    with col3:
        if st.button(
            "▶ Run Identify Spatial Domains",
            type="primary",
            use_container_width=True,
        ):
            with st.spinner("Running spatial domain identification..."):
                try:
                    # Placeholder for the actual spatial domain identification function
                    adata = identify_spatial_domains(adata, alpha=alpha)
                    st.success("Spatial domain identification completed.")
                except Exception as exc:
                    st.error(f"Spatial domain identification failed: {exc}")


# -----------------------------------------------------------------------------
# Results
# -----------------------------------------------------------------------------

if "X_pca" in adata.obsm or "X_umap" in adata.obsm:

    st.divider()

    if st.session_state['is_spatial']:
        col_pca, col_umap, col_domains = st.columns(3)
    else:
        col_pca, col_umap = st.columns(2)

    # -------------------------------------------------------------------------
    # PCA
    # -------------------------------------------------------------------------

    with col_pca:

        if "X_pca" in adata.obsm:

            st.subheader("PCA")

            pca_color_by = st.selectbox(
                "Color cells by",
                options=color_options(adata),
                key="pca_color_by",
            )

            render_interactive(adata, "pca", key="pca", color_by=pca_color_by)

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
                options=color_options(adata),
                key="umap_color_by",
            )

            render_interactive(adata, "umap", key="umap", color_by=umap_color_by)

        else:
            st.info("UMAP has not been run yet.")

    if st.session_state['is_spatial']:
        with col_domains:
            
            if "squidpy_domains" not in adata.obs:
                st.info("Run UMAP Clustering and Spatial Domain analysis first.")
            else:
                st.subheader("Spatial")
                # st.caption(
                #     "Spatial domains are obtained from a joint expression/spatial graph, not validated biological compartments."
                # )
                color_by = st.selectbox(
                    "Color cells by",
                    options=color_options(adata),
                    placeholder="squidpy_domains",
                    key="spatial_color_by",
                )
            
                render_interactive(adata, "spatial", key="spatial", color_by=color_by)

            # st.dataframe(
            #     domain_counts,
            #     use_container_width=True,
            # )

# -----------------------------------------------------------------------------
# Highly variable genes
# -----------------------------------------------------------------------------

if "highly_variable" in adata.var:
    st.divider()
    st.subheader("Highly Variable Genes")
    render_hvg(adata, key="hvg")

if 'cluster_label' in adata.obs:
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


if 'squidpy_domains' in adata.obs:
    st.divider()

    st.subheader("Spatial Domains")

    domain_counts = (
        adata.obs["squidpy_domains"]
        .value_counts()
        .sort_index()
        .rename("cells")
        .to_frame()
    )

    domain_counts["percentage"] = (
        domain_counts["cells"]
        / len(adata)
        * 100
    ).round(2)

    st.dataframe(
        domain_counts,
        use_container_width=True,
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

if "params" not in st.session_state:
    st.session_state.params = {}

st.session_state.params['n_top_genes'] = n_top_genes
st.session_state.params['n_comps'] = n_comps
st.session_state.params['n_neighbors'] = n_neighbors
st.session_state.params['resolution'] = resolution
st.session_state.params['min_dist'] = min_dist
st.session_state.params['random_state'] = RANDOM_STATE
if st.session_state['is_spatial']:
    st.session_state.params['alpha'] = alpha

st.session_state["clustering_work"] = adata

# Annotation requires a complete embedding and clustering, not just a started run.
if (
    "X_pca" in adata.obsm
    and "X_umap" in adata.obsm
    and "cluster_label" in adata.obs
):
    st.session_state.adatas['clustered'] = adata
else:
    st.session_state.adatas.pop('clustered', None)