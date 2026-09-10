import tempfile

import pandas as pd
import scanpy as sc
import streamlit as st

import plotly.graph_objects as go


from src.utils import load_data
from src.filtering import data_filtering
from src.plotting import highest_expressed_genes, plot_qc_metrics, qc_summary

if "adatas" not in st.session_state:
    st.session_state.adatas = None
    st.session_state.filtered = False

col1, col2, col3 = st.columns(3)

with col1:
    adata_file = st.file_uploader(
        "Upload your .h5ad file",
        type=["h5ad"]
    )

with col2:
    cells_file = st.file_uploader(
        "Upload your cells csv file",
        type=["csv"]
    )

with col3:
    genes_file = st.file_uploader(
        "Upload your genes csv file",
        type=["csv"]
    )

try:
    if adata_file is not None:
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
            tmp.write(adata_file.getvalue())
            tmp_path = tmp.name
        st.session_state.adatas = load_data(tmp_path)
    if cells_file is not None:
        cells = pd.read_csv(cells_file)
    if genes_file is not None:
        genes = pd.read_csv(genes_file)
except Exception as e:
    print(f"Error reading uploaded file: {e}")



# --------------------------------------------------
# Analysis section
# --------------------------------------------------

st.divider()
st.header("Analysis")

filter_col, plot_col = st.columns([1, 3])

# --------------------------------------------------
# Filters
# --------------------------------------------------

with filter_col:
    st.subheader("Filters")

    with st.form("filter_form"):
        min_genes = st.number_input(
            "Minimum genes",
            min_value=0,
            max_value=1_000_000,
            value=200
        )

        min_cells = st.number_input(
            "Minimum cells",
            min_value=0,
            max_value=1_000_000,
            value=3
        )

        max_mt_percentage = st.number_input(
            "Maximum mitochondrial gene percentage",
            min_value=0.0,
            max_value=100.0,
            value=20.0
        )
        
        exp_doublet_rate = st.number_input(
            "Expected doublet rate",
            min_value=0.0,
            max_value=1.0,
            value=0.05
        )

        apply_filter = st.form_submit_button("Apply filters")

        if apply_filter:
            adata = (
                        st.session_state.adatas['subset']
                        if st.session_state.adatas['subset'] is not None
                        else st.session_state.adatas['raw']
                    )
            adata_filtered = data_filtering(
                adata,
                min_genes=float(min_genes),
                min_cells=float(min_cells),
                max_mt_percentage=float(max_mt_percentage),
                expected_doublet_rate=float(exp_doublet_rate)
            )
            
            st.session_state.adatas['filtered_w_doublets'] = adata_filtered
            st.session_state.adatas['filtered_wout_doublets'] = adata_filtered[~adata_filtered.obs['predicted_doublet']].copy()
            st.session_state.filtered = True

# --------------------------------------------------
# Plots
# --------------------------------------------------

with plot_col:
    
    selected = st.pills(
        "Analysis",
        ["Top Gene", "QC Metrics", "Statistics"],
        default="Top Gene",
    )

    if selected == "Top Gene":
        if st.session_state.filtered:
            adata = (
                        st.session_state.adatas['subset']
                        if st.session_state.adatas['subset'] is not None
                        else st.session_state.adatas['raw']
                    )
            
            
            fig_one = highest_expressed_genes(adata, n_top_genes=10)
            fig_two = highest_expressed_genes(st.session_state.adatas['filtered_w_doublets'], n_top_genes=10)
            fig_three = highest_expressed_genes(st.session_state.adatas['filtered_wout_doublets'], n_top_genes=10)

            col1, col2, col3 = st.columns(3)

            with col1:
                st.subheader("Subset / Raw")
                st.pyplot(fig_one, use_container_width=True)

            with col2:
                st.subheader("Filtered with Doublets")
                st.pyplot(fig_two, use_container_width=True)

            with col3:
                st.subheader("Filtered without Doublets")
                st.pyplot(fig_three, use_container_width=True)

        
    elif selected == "QC Metrics":
        if st.session_state.filtered:
            adata = (
                        st.session_state.adatas['subset']
                        if st.session_state.adatas['subset'] is not None
                        else st.session_state.adatas['raw']
                     )
                    
                    
            fig_one = plot_qc_metrics(adata)
            fig_two = plot_qc_metrics(st.session_state.adatas['filtered_w_doublets'])
            fig_three = plot_qc_metrics(st.session_state.adatas['filtered_wout_doublets'])

            col1, col2, col3 = st.columns(3)
        
            with col1:
                st.subheader("Subset / Raw")
                st.pyplot(fig_one, use_container_width=True)
        
            with col2:
                st.subheader("Filtered with Doublets")
                st.pyplot(fig_two, use_container_width=True)
        
            with col3:
                st.subheader("Filtered without Doublets")
                st.pyplot(fig_three, use_container_width=True)
    elif selected == "Statistics":
        if st.session_state.filtered:
            adata = (
                        st.session_state.adatas['subset']
                        if st.session_state.adatas['subset'] is not None
                        else st.session_state.adatas['raw']
                     )
            counts = st.session_state.adatas['filtered_w_doublets'].obs['predicted_doublet'].value_counts()

            singlets = counts.get(False, 0)
            doublets = counts.get(True, 0)

            total = singlets + doublets

            singlet_pct = singlets / total * 100
            doublet_pct = doublets / total * 100
            
            fig = go.Figure()

            fig.add_trace(
                go.Bar(
                    x=[singlet_pct],
                    y=["Cells"],
                    orientation="h",
                    name="Singlets",
                    width=0.35,
                    customdata=[[singlets]],
                    hovertemplate="Singlets: %{customdata[0]:,}<extra></extra>",
                    marker_color="#4CAF50",
                )
            )

            fig.add_trace(
                go.Bar(
                    x=[doublet_pct],
                    y=["Cells"],
                    orientation="h",
                    name="Doublets",
                    width=0.35,
                    customdata=[[doublets]],
                    hovertemplate="Doublets: %{customdata[0]:,}<extra></extra>",
                    marker_color="#E74C3C",
                )
            )

            fig.update_layout(
                barmode="stack",
                height=100,
                xaxis=dict(
                    range=[0, 100],
                    title="Percentage of cells (%)",
                ),
                yaxis=dict(
                    showticklabels=False,
                ),
                margin=dict(l=10, r=10, t=5, b=40),
                showlegend=False,
            )


            st.plotly_chart(fig, use_container_width=True)
            
            
            summaries = [
                ("Subset / Raw", qc_summary(adata)),
                ("Filtered + Doublets", qc_summary(st.session_state.adatas['filtered_w_doublets'])),
                ("Filtered − Doublets", qc_summary(st.session_state.adatas['filtered_wout_doublets'])),
            ]

            cols = st.columns(3)

            for col, (title, summary) in zip(cols, summaries):
                with col:
                    st.markdown(f"### {title}")

                    st.metric(
                        "Cells",
                        f"{summary['n_cells']:,}",
                    )

                    st.metric(
                        "Median genes / cell",
                        f"{summary['median_genes']:,.0f}",
                        help=(
                            "n_genes_by_counts is the number of genes detected in each "
                            "cell. A higher value generally indicates greater RNA "
                            "complexity, while very low values can indicate low-quality "
                            "or empty cells."
                        ),
                    )

                    st.caption(
                        f"Mean: {summary['mean_genes']:,.0f}"
                    )

                    st.metric(
                        "Median counts / cell",
                        f"{summary['median_counts']:,.0f}",
                        help=(
                            "total_counts is the total number of RNA counts (UMIs) "
                            "detected in each cell. It reflects the amount of RNA "
                            "captured and the sequencing depth of the cell."
                        ),
                    )

                    st.caption(
                        f"Mean: {summary['mean_counts']:,.0f}"
                    )

