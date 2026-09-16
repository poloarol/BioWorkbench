import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

from src.utils import calculate_module_score
from src.plotting import plot_umap, plot_spatial

if st.session_state['is_spatial']:
    adata = st.session_state.adatas['annotated']
    
    
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

    if st.sidebar.button("Submit genes"):
        if selected_genes and module_name:
            adata = calculate_module_score(
                adata,
                genes=selected_genes,
                score_name=f"module_{module_name}",
            )
            # st.session_state.adatas['annotated'] = adata
            col1, col2 = st.columns(2)

            with col1:
                st.pyplot(
                    plot_umap(
                        adata,
                        color_by=f"module_{module_name}",
                        figsize=(12, 10),
                    )
                )

            with col2:
                st.pyplot(
                    plot_spatial(
                        adata,
                        color_by=f"module_{module_name}",
                        figsize=(12, 10),
                    )
                )
    else:
        col1, col2 = st.columns(2)

        with col1:
            st.pyplot(
                plot_umap(
                    adata,
                    color_by="celltypist_predicted_labels",
                    figsize=(12, 10),
                )
            )

        with col2:
            st.pyplot(
                plot_spatial(
                    adata,
                    color_by="celltypist_predicted_labels",
                    figsize=(12, 10),
                )
            )
