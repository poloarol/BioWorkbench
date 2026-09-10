import tempfile

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from src.utils import load_data
from src.filtering import data_filtering
from src.plotting import (
    highest_expressed_genes,
    plot_qc_metrics,
    qc_summary,
)


# ============================================================
# Session state
# ============================================================

if "adatas" not in st.session_state:
    st.session_state.adatas = None

if "filtered" not in st.session_state:
    st.session_state.filtered = False

if "uploaded_adata_name" not in st.session_state:
    st.session_state.uploaded_adata_name = None

if "cells" not in st.session_state:
    st.session_state.cells = None

if "genes" not in st.session_state:
    st.session_state.genes = None


# ============================================================
# Helper functions
# ============================================================

def get_base_adata():
    """
    Return the subset AnnData if it exists, otherwise raw AnnData.

    Returns None if no usable AnnData object is loaded.
    """
    adatas = st.session_state.get("adatas")

    if not isinstance(adatas, dict):
        return None

    subset = adatas.get("subset")
    raw = adatas.get("raw")

    if subset is not None:
        return subset

    if raw is not None:
        return raw

    return None


def get_filtered_adata():
    """Return the filtered AnnData object, or None if unavailable."""
    adatas = st.session_state.get("adatas")

    if not isinstance(adatas, dict):
        return None

    return adatas.get("filtered_w_doublets")


def get_filtered_singlets():
    """Return the filtered singlet AnnData object, or None if unavailable."""
    adatas = st.session_state.get("adatas")

    if not isinstance(adatas, dict):
        return None

    return adatas.get("filtered_wout_doublets")


def clear_filtered_results():
    """Remove previous filtering results."""
    if not isinstance(st.session_state.adatas, dict):
        return

    st.session_state.adatas.pop("filtered_w_doublets", None)
    st.session_state.adatas.pop("filtered_wout_doublets", None)
    st.session_state.filtered = False


# ============================================================
# Sidebar
# ============================================================

with st.sidebar:
    st.title("Single-cell Analysis")

    # --------------------------------------------------------
    # Data upload
    # --------------------------------------------------------

    st.header("Data")

    adata_file = st.file_uploader(
        "Upload your .h5ad file",
        type=["h5ad"],
        key="adata_uploader",
    )

    cells_file = st.file_uploader(
        "Upload your cells CSV file",
        type=["csv"],
        key="cells_uploader",
    )

    genes_file = st.file_uploader(
        "Upload your genes CSV file",
        type=["csv"],
        key="genes_uploader",
    )

    # --------------------------------------------------------
    # Load uploaded files
    # --------------------------------------------------------

    if adata_file is not None:
        # Only reload when a different file is uploaded.
        if st.session_state.uploaded_adata_name != adata_file.name:
            try:
                with tempfile.NamedTemporaryFile(
                    suffix=".h5ad",
                    delete=False,
                ) as tmp:
                    tmp.write(adata_file.getvalue())
                    tmp_path = tmp.name

                loaded_adatas = load_data(tmp_path)

                if loaded_adatas is None:
                    st.error("The .h5ad file did not contain any data.")
                elif not isinstance(loaded_adatas, dict):
                    st.error(
                        "load_data() returned an unexpected object. "
                        "Expected a dictionary containing AnnData objects."
                    )
                else:
                    st.session_state.adatas = loaded_adatas
                    st.session_state.uploaded_adata_name = adata_file.name
                    st.session_state.filtered = False

            except Exception as e:
                st.error(f"Error loading .h5ad file: {e}")
                st.session_state.adatas = None
                st.session_state.uploaded_adata_name = None
                st.session_state.filtered = False

    if cells_file is not None:
        try:
            st.session_state.cells = pd.read_csv(cells_file)
        except Exception as e:
            st.error(f"Error reading cells CSV: {e}")
            st.session_state.cells = None

    if genes_file is not None:
        try:
            st.session_state.genes = pd.read_csv(genes_file)
        except Exception as e:
            st.error(f"Error reading genes CSV: {e}")
            st.session_state.genes = None

    # --------------------------------------------------------
    # Filter controls
    # --------------------------------------------------------

    st.divider()
    st.header("Filters")

    with st.form("filter_form"):
        min_genes = st.number_input(
            "Minimum genes",
            min_value=0,
            max_value=1_000_000,
            value=200,
            step=10,
        )

        min_cells = st.number_input(
            "Minimum cells",
            min_value=0,
            max_value=1_000_000,
            value=3,
            step=1,
        )

        max_mt_percentage = st.number_input(
            "Maximum mitochondrial gene percentage",
            min_value=0.0,
            max_value=100.0,
            value=20.0,
            step=1.0,
        )

        exp_doublet_rate = st.number_input(
            "Expected doublet rate",
            min_value=0.0,
            max_value=1.0,
            value=0.05,
            step=0.01,
        )

        apply_filter = st.form_submit_button(
            "Apply filters",
            use_container_width=True,
        )

    # --------------------------------------------------------
    # Apply filters
    # --------------------------------------------------------

    if apply_filter:
        adata = get_base_adata()

        if adata is None:
            st.warning(
                "Please upload a valid .h5ad file before applying filters."
            )
        elif adata.n_obs == 0:
            st.warning("The input dataset contains no cells.")
        elif adata.n_vars == 0:
            st.warning("The input dataset contains no genes.")
        else:
            try:
                with st.spinner("Applying filters..."):
                    adata_filtered = data_filtering(
                        adata,
                        min_genes=float(min_genes),
                        min_cells=float(min_cells),
                        max_mt_percentage=float(max_mt_percentage),
                        expected_doublet_rate=float(exp_doublet_rate),
                    )

                if adata_filtered is None:
                    st.error("Filtering returned no data.")
                    clear_filtered_results()

                elif adata_filtered.n_obs == 0:
                    st.warning(
                        "The filters removed all cells. "
                        "Try less restrictive filtering parameters."
                    )
                    st.session_state.adatas[
                        "filtered_w_doublets"
                    ] = adata_filtered

                    st.session_state.adatas[
                        "filtered_wout_doublets"
                    ] = adata_filtered.copy()

                    st.session_state.filtered = True

                else:
                    # Make sure predicted_doublet exists before using it.
                    if "predicted_doublet" not in adata_filtered.obs.columns:
                        st.warning(
                            "Filtering completed, but "
                            "'predicted_doublet' was not found. "
                            "Doublet-specific results may be unavailable."
                        )

                        adata_without_doublets = adata_filtered.copy()

                    else:
                        doublet_mask = (
                            adata_filtered.obs["predicted_doublet"]
                            .fillna(False)
                            .astype(bool)
                        )

                        adata_without_doublets = adata_filtered[
                            ~doublet_mask
                        ].copy()

                    st.session_state.adatas[
                        "filtered_w_doublets"
                    ] = adata_filtered

                    st.session_state.adatas[
                        "filtered_wout_doublets"
                    ] = adata_without_doublets

                    st.session_state.filtered = True

                    st.success(
                        f"Filtering complete: "
                        f"{adata_filtered.n_obs:,} cells remaining."
                    )

            except Exception as e:
                st.error(f"Error while filtering data: {e}")
                clear_filtered_results()


# ============================================================
# Main analysis area
# ============================================================

st.title("Analysis")

base_adata = get_base_adata()
filtered_adata = get_filtered_adata()
singlet_adata = get_filtered_singlets()


# ------------------------------------------------------------
# No data loaded
# ------------------------------------------------------------

if base_adata is None:
    st.info(
        "Upload an .h5ad file from the sidebar to begin the analysis."
    )
    st.stop()


# ------------------------------------------------------------
# Dataset overview
# ------------------------------------------------------------

st.caption(
    f"Current dataset: **{base_adata.n_obs:,} cells × "
    f"{base_adata.n_vars:,} genes**"
)

st.divider()


# ------------------------------------------------------------
# Analysis selector
# ------------------------------------------------------------

selected = st.pills(
    "Analysis",
    ["Top Gene", "QC Metrics", "Statistics"],
    default="Top Gene",
)


# ------------------------------------------------------------
# Require filtering for downstream analysis
# ------------------------------------------------------------

if not st.session_state.filtered:
    st.info(
        "Configure your filtering parameters in the sidebar "
        "and click **Apply filters** to view the analysis."
    )
    st.stop()


if filtered_adata is None:
    st.warning(
        "No filtered dataset is available. "
        "Please apply the filters again."
    )
    st.stop()


# ============================================================
# Top Gene
# ============================================================

if selected == "Top Gene":

    try:
        fig_one = highest_expressed_genes(
            base_adata,
            n_top_genes=10,
        )

        fig_two = highest_expressed_genes(
            filtered_adata,
            n_top_genes=10,
        )

        fig_three = highest_expressed_genes(
            singlet_adata,
            n_top_genes=10,
        )

        cols = st.columns(3)

        with cols[0]:
            st.subheader("Subset / Raw")
            st.pyplot(fig_one, use_container_width=True)

        with cols[1]:
            st.subheader("Filtered with Doublets")
            st.pyplot(fig_two, use_container_width=True)

        with cols[2]:
            st.subheader("Filtered without Doublets")
            st.pyplot(fig_three, use_container_width=True)

    except Exception as e:
        st.error(f"Unable to generate gene plots: {e}")


# ============================================================
# QC Metrics
# ============================================================

elif selected == "QC Metrics":

    try:
        fig_one = plot_qc_metrics(base_adata)
        fig_two = plot_qc_metrics(filtered_adata)
        fig_three = plot_qc_metrics(singlet_adata)

        cols = st.columns(3)

        with cols[0]:
            st.subheader("Subset / Raw")
            st.pyplot(fig_one, use_container_width=True)

        with cols[1]:
            st.subheader("Filtered with Doublets")
            st.pyplot(fig_two, use_container_width=True)

        with cols[2]:
            st.subheader("Filtered without Doublets")
            st.pyplot(fig_three, use_container_width=True)

    except Exception as e:
        st.error(f"Unable to generate QC plots: {e}")


# ============================================================
# Statistics
# ============================================================

elif selected == "Statistics":

    # --------------------------------------------------------
    # Doublet statistics
    # --------------------------------------------------------

    predicted_doublet = None

    if "predicted_doublet" in filtered_adata.obs.columns:
        predicted_doublet = (
            filtered_adata.obs["predicted_doublet"]
            .fillna(False)
            .astype(bool)
        )

    if predicted_doublet is not None:

        singlets = int((~predicted_doublet).sum())
        doublets = int(predicted_doublet.sum())
        total = singlets + doublets

        if total > 0:
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
                    hovertemplate=(
                        "Singlets: %{customdata[0]:,}"
                        "<extra></extra>"
                    ),
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
                    hovertemplate=(
                        "Doublets: %{customdata[0]:,}"
                        "<extra></extra>"
                    ),
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
                margin=dict(
                    l=10,
                    r=10,
                    t=5,
                    b=40,
                ),
                showlegend=False,
            )

            st.plotly_chart(
                fig,
                use_container_width=True,
            )

        else:
            st.info(
                "No cells remain after filtering, so doublet "
                "statistics cannot be calculated."
            )

    else:
        st.info(
            "The filtered dataset does not contain a "
            "'predicted_doublet' column."
        )

    # --------------------------------------------------------
    # QC summaries
    # --------------------------------------------------------

    datasets = [
        ("Subset / Raw", base_adata),
        ("Filtered + Doublets", filtered_adata),
        ("Filtered − Doublets", singlet_adata),
    ]

    cols = st.columns(3)

    for col, (title, adata) in zip(cols, datasets):

        with col:
            st.markdown(f"### {title}")

            if adata is None or adata.n_obs == 0:
                st.info("No cells available.")
                continue

            try:
                summary = qc_summary(adata)

                n_cells = summary.get("n_cells")
                median_genes = summary.get("median_genes")
                mean_genes = summary.get("mean_genes")
                median_counts = summary.get("median_counts")
                mean_counts = summary.get("mean_counts")

                if n_cells is not None:
                    st.metric(
                        "Cells",
                        f"{n_cells:,}",
                    )

                if median_genes is not None:
                    st.metric(
                        "Median genes / cell",
                        f"{median_genes:,.0f}",
                        help=(
                            "n_genes_by_counts is the number of genes "
                            "detected in each cell. A higher value "
                            "generally indicates greater RNA complexity, "
                            "while very low values can indicate "
                            "low-quality or empty cells."
                        ),
                    )

                if mean_genes is not None:
                    st.caption(
                        f"Mean: {mean_genes:,.0f}"
                    )

                if median_counts is not None:
                    st.metric(
                        "Median counts / cell",
                        f"{median_counts:,.0f}",
                        help=(
                            "total_counts is the total number of RNA "
                            "counts (UMIs) detected in each cell. "
                            "It reflects the amount of RNA captured "
                            "and the sequencing depth of the cell."
                        ),
                    )

                if mean_counts is not None:
                    st.caption(
                        f"Mean: {mean_counts:,.0f}"
                    )

            except Exception as e:
                st.error(
                    f"Unable to calculate summary: {e}"
                )
