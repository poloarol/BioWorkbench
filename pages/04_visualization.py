import streamlit as st

from src.utils import calculate_module_score
from src.interactive import render_interactive, is_categorical

if "is_spatial" not in st.session_state:
    st.session_state["is_spatial"] = False


if 'module_score' not in st.session_state:
    st.session_state['module_score'] = []

adatas = st.session_state.get("adatas", {})
adata = adatas.get("annotated") if isinstance(adatas, dict) else None

if adata is None:
    st.info("Run cell annotation before opening the visualizer.")
    st.stop()

if "X_umap" not in adata.obsm:
    st.info("Run UMAP on the Clustering page before opening the visualizer.")
    st.stop()
    
st.sidebar.subheader("Genes")

genes = [
    gene for gene in adata.var_names
    if not gene.startswith("Blank-")
]

gene_search = st.sidebar.text_input(
    "Search genes",
     placeholder="e.g. GFAP, EGFR...",
)


if gene_search:
    genes_to_display = [
        gene for gene in genes
        if gene_search.lower() in gene.lower()
    ]
else:
    genes_to_display = genes

selected_genes = []

with st.sidebar.container(height=400):
    for gene in genes_to_display:
        if st.checkbox(gene, key=f"gene_{gene}"):
            selected_genes.append(gene)


module_name = st.sidebar.text_input(
    "Module name",
    placeholder="e.g. Astrocyte module",
    value="module_score"
)

if st.sidebar.button("Submit Gene Set"):
    module_name = module_name.strip()
    if selected_genes and module_name:
        score_name = f"module_{module_name}"
        adata = calculate_module_score(
            adata,
            genes=selected_genes,
            score_name=score_name,
        )

        if score_name not in st.session_state['module_score']:
            st.session_state['module_score'].append(score_name)
        
st.subheader("Gene Expression Visualizer")

if not genes:
    st.info("No genes are available to visualize.")
    st.stop()

color_by = st.selectbox(
    "Select Gene...",
    options=genes,
    key="gene_exp",
)

group_options = [
    column
    for column in adata.obs.columns
    if is_categorical(adata.obs[column]) and 1 < adata.obs[column].nunique() <= 100
]
group_by = st.selectbox(
    "Group cells by (for statistics / filtering)",
    options=group_options,
    index=None,
    key="gene_group_by",
    placeholder="Optional: cluster, domain, annotation...",
)

if st.session_state['is_spatial']:
    col1, col2 = st.columns(2)

    with col1:
        render_interactive(adata, "umap", key="gene_umap", gene=color_by, group_by=group_by)
    with col2:
        render_interactive(adata, "spatial", key="gene_spatial", gene=color_by, group_by=group_by)
else:
    render_interactive(adata, "umap", key="gene_umap", gene=color_by, group_by=group_by)

st.divider()

st.subheader("Module Scores Visualizer")

available_modules = [
    module
    for module in st.session_state["module_score"]
    if module in adata.obs.columns
]

if not available_modules:
    st.info("Submit a gene set to create a module score for visualization.")
else:
    color_by = st.selectbox(
        "Module Score",
        options=available_modules,
        key="module_color",
        placeholder="Select a module score",
    )

    if color_by is None:
        st.warning("Please select a module score to visualize.")
    elif st.session_state['is_spatial']:
        col1, col2 = st.columns(2)

        with col1:
            render_interactive(adata, "umap", key="module_umap", color_by=color_by)

        with col2:
            render_interactive(adata, "spatial", key="module_spatial", color_by=color_by)
    else:
        render_interactive(adata, "umap", key="module_umap", color_by=color_by)
