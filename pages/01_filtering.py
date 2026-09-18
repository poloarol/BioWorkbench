import tempfile
import yaml

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
import streamlit as st

from plotly.subplots import make_subplots

from src.utils import load_data
from src.filtering import data_filtering
from src.plotting import (
    highest_expressed_genes,
    plot_qc_metrics,
    qc_summary,
    plot_spatial
)


# ============================================================
# Session state
# ============================================================

CONFIG_PATH = Path(__file__).parents[1] / "config" / "params.yaml"

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

if 'is_spatial' not in st.session_state:
    st.session_state.is_spatial = False


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

    # subset = adatas.get("subset")
    raw = adatas.get("raw")

    # if subset is not None:
    #     return subset

    if raw is not None:
        return raw

    return None


def get_adata(key: str):
    """Return the AnnData object corresponding to the given key, or None if unavailable."""
    adatas = st.session_state.get("adatas")

    if not isinstance(adatas, dict):
        return None

    return adatas.get(key)


def clear_filtered_results():
    """Remove previous filtering results."""
    if not isinstance(st.session_state.adatas, dict):
        return

    st.session_state.adatas.pop("filtered_w_doublets", None)
    st.session_state.adatas.pop("filtered_wout_doublets", None)
    st.session_state.adatas.pop("filtered_w_blanks", None)
    st.session_state.adatas.pop("filtered_wout_blanks", None)
    st.session_state.filtered = False


def load_filtering_config():
    with open(CONFIG_PATH, "r") as f:
        return yaml.safe_load(f)


def get_filtering_config(dataset_type, technology=None):
    config = load_filtering_config()

    if dataset_type == "single_cell":
        return config["single_cell"]

    if dataset_type == "spatial":
        return config["spatial"][technology]

    raise ValueError(f"Unsupported dataset type: {dataset_type}")

# ============================================================
# Sidebar
# ============================================================

with st.sidebar:
    st.title("Single-cell Analysis")
    
    spatial_technology = st.sidebar.selectbox(
        "Technology",
        [
            "None — single-cell dataset",
            "MERFISH",
            # "Xenium",
            # "CosMX",
            # "Visium HD",
        ],
    )
    
    is_spatial = spatial_technology != "None — single-cell dataset"
    st.session_state["spatial_technology"] = spatial_technology
    st.session_state["is_spatial"] = is_spatial
    
    config = get_filtering_config(
        "spatial" if is_spatial else "single_cell",
        technology=spatial_technology if is_spatial else None,
    )
    
    # --------------------------------------------------------
    # Data upload
    # --------------------------------------------------------

    st.subheader("Data")

    adata_file = st.file_uploader(
        "Upload your .h5ad file",
        type=["h5ad"],
        key="adata_uploader",
        max_upload_size = 500
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

                loaded_adatas = load_data(tmp_path, is_spatial=is_spatial)

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
                st.session_state.adatas = None
                st.session_state.uploaded_adata_name = None
                st.session_state.filtered = False
                
                import traceback

                st.error(
                    f"Error loading .h5ad file: "
                    f"{type(e).__name__}: {e}"
                )
                st.code(traceback.format_exc())

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
    
    min_genes: int = 0
    min_cells: int = 0
    blank_threshold: float = 0.0
    max_mt_percentage: float = 0.0
    exp_doublet_rate: float = 0.0

    with st.form("filter_form"):
        min_genes = st.number_input(
            "Minimum genes",
            min_value=0,
            max_value=1_000_000,
            value=config.get("min_genes", 200),
            step=10,
        )

        min_cells = st.number_input(
            "Minimum cells",
            min_value=0,
            max_value=1_000_000,
            value=config.get("min_cells", 3),
            step=1,
        )
        
        if is_spatial:
            
            blank_threshold = st.number_input(
                "Blank threshold",
                min_value=0.0,
                max_value=100.0,
                value=float(config.get('blank_threshold', 5.0)),
                step=1.0
            )
            
        else:

            max_mt_percentage = st.number_input(
                "Maximum mitochondrial gene percentage",
                min_value=0.0,
                max_value=100.0,
                value=float(config.get("max_mito_percent", 20.0)),
                step=1.0,
            )

            exp_doublet_rate = st.number_input(
                "Expected doublet rate",
                min_value=0.0,
                max_value=1.0,
                value=float(config.get("expected_doublet_rate", 0.05)),
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
                        is_spatial=is_spatial,
                        blank_max_percentage=float(blank_threshold)
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
                    if is_spatial:
                        if "pct_counts_blank" not in adata_filtered.obs.columns:
                            st.warning(
                                "Filtering completed, but "
                                "'pct_counts_blank' was not found. "
                                "Blank-specific results may be unavailable."
                            )
                            adata_without_blanks = adata_filtered.copy()
                        else:
                            adata_without_blanks = adata_filtered[
                                ~adata_filtered.obs['is_blank']
                            ].copy()
                            st.session_state.adatas[
                                "filtered_wout_blanks"
                            ] = adata_without_blanks
                            st.session_state.adatas[
                                "filtered_w_blanks"
                            ] = adata_filtered
                            
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

if is_spatial:
    filtered_adata = get_adata("filtered_w_blanks")
    singlet_adata = get_adata("filtered_wout_blanks")
else:
    filtered_adata = get_adata("filtered_w_doublets")
    singlet_adata = get_adata("filtered_wout_doublets")

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

if is_spatial:
    blank_percentage: int = (base_adata.obs["blank_counts"].to_numpy().sum() / base_adata.var["total_counts"].sum()) * 100
    st.caption(f"Blank gene percentage: {blank_percentage:.2f}%")

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


ttl: str = 'Blank' if is_spatial else 'Doublet'

print(filtered_adata)

if selected == "Top Gene":
    if is_spatial:
        x = base_adata.obs["pct_counts_blank"].dropna()

        fig = make_subplots(
            rows=2,
            cols=1,
            shared_xaxes=True,
            row_heights=[0.75, 0.25],
            vertical_spacing=0.03,
        )

        # Histogram
        fig.add_trace(
            go.Histogram(
                x=x,
                nbinsx=50,
                name="Cells",
                marker_color="steelblue",
                opacity=0.85,
                hovertemplate="pct_counts_blank: %{x}<br>Cells: %{y}<extra></extra>",
            ),
            row=1,
            col=1,
        )

        # Violin
        fig.add_trace(
            go.Violin(
                x=x,
                name="Distribution",
                box=dict(visible=True),
                meanline=dict(visible=True),
                points=False,
                fillcolor="lightblue",
                line_color="steelblue",
                hovertemplate="pct_counts_blank: %{x}<extra></extra>",
                orientation="h",
            ),
            row=2,
            col=1,
        )

        fig.update_layout(
            height=600,
            template="plotly_white",
            showlegend=False,
            title="Distribution of pct_counts_blank",
            dragmode="zoom",
        )

        fig.update_yaxes(title_text="Number of cells", row=1, col=1)
        fig.update_xaxes(title_text="pct_counts_blank (%)", row=2, col=1)
        
        st.plotly_chart(
            fig,
            use_container_width=True,
        )

    st.subheader("Highly Variable Genes")
    st.markdown(
        "The following plots show the top highly variable genes for each dataset subset."
    )

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
            st.subheader(f"Filtered with {ttl}s")
            st.pyplot(fig_two, use_container_width=True)

        with cols[2]:
            st.subheader(f"Filtered without {ttl}s")
            st.pyplot(fig_three, use_container_width=True)

    except Exception as e:
        st.error(f"Unable to generate gene plots: {e}")


# ============================================================
# QC Metrics
# ============================================================

elif selected == "QC Metrics":

    try:
        fig_one = plot_qc_metrics(base_adata, is_spatial=is_spatial)
        fig_two = plot_qc_metrics(filtered_adata, is_spatial=is_spatial)
        fig_three = plot_qc_metrics(singlet_adata, is_spatial=is_spatial)

        st.subheader("Subset / Raw")
        if is_spatial:
            col1, col2 = st.columns(2)
            with col1:
                st.pyplot(fig_one, use_container_width=True)
            with col2:
                spatial_qc_plt = plot_spatial(
                    filtered_adata,
                    color_by=[
                        "is_blank",
                        "log1p_n_genes_by_counts",
                        "log1p_total_counts",
                    ],
                )
                st.pyplot(spatial_qc_plt, use_container_width=True)     
        else:
            st.pyplot(fig_one, use_container_width=True)

        st.markdown("<br><br>", unsafe_allow_html=True)

        st.subheader(f"Filtered with {ttl}s")
        if is_spatial:
            col1, col2 = st.columns(2)
            with col1:
                st.pyplot(fig_two, use_container_width=True)
            with col2:
                spatial_qc_plt = plot_spatial(
                    filtered_adata,
                    color_by=[
                        "is_blank",
                        "log1p_n_genes_by_counts",
                        "log1p_total_counts",
                    ],
                )
                st.pyplot(spatial_qc_plt, use_container_width=True)
        else:
            st.pyplot(fig_two, use_container_width=True)

        st.markdown("<br><br>", unsafe_allow_html=True)

        st.subheader(f"Filtered without {ttl}s")
        if is_spatial:
            col1, col2 = st.columns(2)
            with col1:
                st.pyplot(fig_three, use_container_width=True)
            with col2:
                spatial_qc_plt = plot_spatial(
                    singlet_adata,
                    color_by=[
                        # "is_blank",
                        "log1p_n_genes_by_counts",
                        "log1p_total_counts",
                    ],
                )
                st.pyplot(spatial_qc_plt, use_container_width=True)
        else:
            st.pyplot(fig_three, use_container_width=True)

    except Exception as e:
        st.error(f"Unable to generate QC plots: {e}")

# ============================================================
# Statistics
# ============================================================

elif selected == "Statistics":

    # --------------------------------------------------------
    # Helpers
    # --------------------------------------------------------

    def get_bool_obs_column(adata, column):
        """
        Return an obs column as boolean.

        Supports:
        - actual boolean columns
        - numeric 0/1 columns

        Returns None if the column does not exist.

        Raises TypeError for unexpected data types so that
        invalid prediction columns are not silently converted.
        """
        if adata is None or column not in adata.obs.columns:
            return None

        values = adata.obs[column]

        # Already boolean
        if values.dtype == bool:
            return values

        # Numeric 0/1
        if values.dtype.kind in "biu":
            unique_values = values.dropna().unique()

            if all(value in (0, 1) for value in unique_values):
                return values.fillna(0).astype(bool)

        raise TypeError(
            f"{column!r} must contain boolean or numeric 0/1 values; "
            f"got dtype {values.dtype} with values "
            f"{values.dropna().unique()[:10]}"
        )


    def create_percentage_bar(
        values,
        labels,
        total_label,
        colors,
        counts,
    ):
        """
        Create a stacked horizontal percentage bar chart.

        `values` contains the counts used to calculate percentages.
        `counts` contains the counts displayed in the hover text.
        """

        if not (
            len(values)
            == len(labels)
            == len(colors)
            == len(counts)
        ):
            raise ValueError(
                "values, labels, colors, and counts must have "
                "the same length."
            )

        total = sum(values)

        if total <= 0:
            raise ValueError(
                "Cannot create a percentage bar with a total of zero."
            )

        percentages = [
            value / total * 100
            for value in values
        ]

        fig = go.Figure()

        for percentage, label, color, count in zip(
            percentages,
            labels,
            colors,
            counts,
        ):
            fig.add_trace(
                go.Bar(
                    x=[percentage],
                    y=[total_label],
                    orientation="h",
                    name=label,
                    width=0.35,
                    customdata=[[count]],
                    hovertemplate=(
                        f"{label}: "
                        "%{customdata[0]:,}"
                        "<extra></extra>"
                    ),
                    marker_color=color,
                )
            )

        fig.update_layout(
            barmode="stack",
            height=100,
            template="plotly_white",
            xaxis=dict(
                range=[0, 100],
                title=(
                    "Percentage of spots (%)"
                    if total_label == "Spots"
                    else "Percentage of cells (%)"
                ),
                ticksuffix="%",
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

        return fig


    # --------------------------------------------------------
    # Prediction statistics
    # --------------------------------------------------------

    if is_spatial:

        # ----------------------------------------------------
        # Predicted blanks
        # ----------------------------------------------------

        predicted_blanks = get_bool_obs_column(
            filtered_adata,
            "is_blank",
        )

        if predicted_blanks is None:

            st.info(
                "No predicted blank information is available."
            )

        else:

            total = len(predicted_blanks)

            if total == 0:

                st.info(
                    "No spots remain after filtering, so blank "
                    "statistics cannot be calculated."
                )

            else:
                
                is_blank = filtered_adata.obs["is_blank"].fillna(False).astype(bool)
                low_blank = (~is_blank).sum()
                high_blank = is_blank.sum()
                total = low_blank + high_blank
                
                perc_low_blank = low_blank / total * 100 if total > 0 else 0
                perc_high_blank = high_blank / total * 100 if total > 0 else 0
                
                fig = create_percentage_bar(
                    values=[
                        low_blank,
                        high_blank,
                    ],
                    labels=[
                        "Low Blank",
                        "High Blank",
                    ],
                    total_label="Spots",
                    colors=[
                        "#4CAF50",
                        "#E74C3C",
                    ],
                    counts=[
                        low_blank,
                        high_blank,
                    ],
                )

                st.subheader("Predicted Blank Statistics")

                st.plotly_chart(
                    fig,
                    use_container_width=True,
                )

                st.caption(
                    f"{high_blank:,} of {total:,} spots "
                    f"({perc_high_blank:.2f}%) "
                    "are predicted blanks."
                )
    else:

        # ----------------------------------------------------
        # Predicted doublets
        # ----------------------------------------------------

        predicted_doublet = get_bool_obs_column(
            filtered_adata,
            "predicted_doublet",
        )

        if predicted_doublet is None:

            st.info(
                "No predicted doublet information is available."
            )

        else:

            total = len(predicted_doublet)

            if total == 0:

                st.info(
                    "No cells remain after filtering, so doublet "
                    "statistics cannot be calculated."
                )

            else:

                doublets = int(predicted_doublet.sum())
                singlets = total - doublets

                fig = create_percentage_bar(
                    values=[
                        singlets,
                        doublets,
                    ],
                    labels=[
                        "Singlets",
                        "Doublets",
                    ],
                    total_label="Cells",
                    colors=[
                        "#4CAF50",
                        "#E74C3C",
                    ],
                    counts=[
                        singlets,
                        doublets,
                    ],
                )

                st.subheader("Predicted Doublet Statistics")

                st.plotly_chart(
                    fig,
                    use_container_width=True,
                )

                st.caption(
                    f"{doublets:,} of {total:,} cells "
                    f"({doublets / total * 100:.2f}%) "
                    "are predicted doublets."
                )


    # --------------------------------------------------------
    # QC summaries
    # --------------------------------------------------------

    datasets = [
        (
            "Subset / Raw",
            base_adata,
        ),
        (
            f"Filtered + {ttl}s",
            filtered_adata,
        ),
        (
            f"Filtered − {ttl}s",
            singlet_adata,
        ),
    ]

    cols = st.columns(3)

    for col, (title, adata) in zip(
        cols,
        datasets,
    ):

        with col:

            st.markdown(
                f"### {title}"
            )

            if adata is None or adata.n_obs == 0:

                st.info(
                    "No cells available."
                )

                continue

            try:

                summary = qc_summary(
                    adata
                )

                n_cells = summary.get(
                    "n_cells"
                )

                median_genes = summary.get(
                    "median_genes"
                )

                mean_genes = summary.get(
                    "mean_genes"
                )

                median_counts = summary.get(
                    "median_counts"
                )

                mean_counts = summary.get(
                    "mean_counts"
                )

                # --------------------------------------------
                # Number of cells / spots
                # --------------------------------------------

                if n_cells is not None:

                    st.metric(
                        "Cells",
                        f"{n_cells:,}",
                    )

                # --------------------------------------------
                # Median genes
                # --------------------------------------------

                if median_genes is not None:

                    st.metric(
                        "Median genes / cell",
                        f"{median_genes:,.0f}",
                        help=(
                            "n_genes_by_counts is the number "
                            "of genes detected in each cell. "
                            "A higher value generally indicates "
                            "greater RNA complexity, while very "
                            "low values can indicate low-quality "
                            "or empty cells."
                        ),
                    )

                # --------------------------------------------
                # Mean genes
                # --------------------------------------------

                if mean_genes is not None:

                    st.caption(
                        f"Mean: {mean_genes:,.0f}"
                    )

                # --------------------------------------------
                # Median counts
                # --------------------------------------------

                if median_counts is not None:

                    st.metric(
                        "Median counts / cell",
                        f"{median_counts:,.0f}",
                        help=(
                            "total_counts is the total number "
                            "of RNA counts (UMIs) detected in "
                            "each cell. It reflects the amount "
                            "of RNA captured and the sequencing "
                            "depth of the cell."
                        ),
                    )

                # --------------------------------------------
                # Mean counts
                # --------------------------------------------

                if mean_counts is not None:

                    st.caption(
                        f"Mean: {mean_counts:,.0f}"
                    )

            except Exception as e:

                st.error(
                    f"Unable to calculate summary: {e}"
                )
