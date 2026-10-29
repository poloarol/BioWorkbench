import tempfile
import hashlib
import yaml

from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from plotly.subplots import make_subplots

from src.utils import load_data, get_spatially_variable_genes
from src.filtering import (
    data_filtering,
    preview_filter_impact,
    validate_filter_thresholds,
)
from src.interactive import render_interactive
from src.plotting import (
    highest_expressed_genes,
    plot_qc_metrics,
    qc_summary,
    plot_spatial,
    plot_spatial_genes
)
from src.session_bundle import load_session_bundle


# ============================================================
# Session state
# ============================================================

CONFIG_PATH = Path(__file__).parents[1] / "config" / "params.yaml"

if "adatas" not in st.session_state:
    st.session_state.adatas = {}

if "loaded_adata_modality" not in st.session_state:
    st.session_state.loaded_adata_modality = None

if "loaded_adata_signature" not in st.session_state:
    st.session_state.loaded_adata_signature = None

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

if "spatial_technology" not in st.session_state:
    st.session_state.spatial_technology = "None — single-cell dataset"


# ============================================================
# Helper functions
# ============================================================

def get_base_adata():
    """ Return the raw AnnData object if available, otherwise None. """
    adatas = st.session_state.get("adatas")

    if not isinstance(adatas, dict):
        return None

    raw = adatas.get("raw")

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
    st.session_state.filtered = False
    if not isinstance(st.session_state.adatas, dict):
        return

    st.session_state.adatas.pop("filtered_w_doublets", None)
    st.session_state.adatas.pop("filtered_wout_doublets", None)
    st.session_state.adatas.pop("filtered_w_blanks", None)
    st.session_state.adatas.pop("filtered_wout_blanks", None)
    st.session_state.adatas.pop("filtered_w_xenium_controls", None)
    st.session_state.adatas.pop("filtered_wout_xenium_controls", None)


def clear_downstream_results():
    """Remove analysis outputs that depend on the currently loaded dataset."""
    clear_filtered_results()
    if not isinstance(st.session_state.adatas, dict):
        return
    st.session_state.adatas.pop("clustered", None)
    st.session_state.pop("clustering_work", None)
    st.session_state.adatas.pop("annotated", None)
    st.session_state["module_score"] = []


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


def get_saved_number(name, fallback, minimum, maximum, *, integer=False):
    """Return a bounded numeric parameter restored from an export bundle."""
    saved_params = st.session_state.get("params")
    if not isinstance(saved_params, dict):
        return fallback

    value = saved_params.get(name, fallback)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return fallback

    value = min(max(value, minimum), maximum)
    return int(value) if integer else float(value)


# ============================================================
# Sidebar
# ============================================================

with st.sidebar:
    st.title("Single-cell Analysis")

    st.subheader("Restore an export bundle")
    st.caption(
        "Upload a .wkb bundle to replace the current analysis objects "
        "and parameters."
    )
    bundle_file = st.file_uploader(
        "Upload a BioWorkbench bundle",
        type=["wkb"],
        key="session_bundle_upload",
        max_upload_size=2000,
    )
    if st.button(
        "Restore bundle",
        disabled=bundle_file is None,
        key="restore_session_bundle",
        use_container_width=True,
    ):
        try:
            imported_adatas, imported_params = load_session_bundle(
                bundle_file.getvalue()
            )

            spatial_keys = {
                "filtered_w_blanks",
                "filtered_wout_blanks",
                "filtered_w_xenium_controls",
                "filtered_wout_xenium_controls",
            }
            contains_spatial_data = bool(
                spatial_keys.intersection(imported_adatas)
            )
            contains_spatial_metadata = any(
                "spatial" in adata.obsm
                or {"center_x", "center_y"}.issubset(adata.obs.columns)
                or "pct_counts_blank" in adata.obs.columns
                for adata in imported_adatas.values()
            )
            is_spatial = (
                contains_spatial_data
                or contains_spatial_metadata
                or imported_params.get("is_spatial") is True
            )

            st.session_state.adatas = imported_adatas
            st.session_state.pop("clustering_work", None)
            st.session_state.pop("clustering_source", None)
            st.session_state.params = imported_params
            st.session_state.is_spatial = is_spatial
            saved_technology = imported_params.get("spatial_technology")
            st.session_state.pop("spatial_technology_widget", None)
            st.session_state["spatial_technology"] = (
                saved_technology
                if is_spatial and saved_technology in ("MERFISH", "Visium HD", "Xenium", "CosMx")
                else "None — single-cell dataset"
            )
            st.session_state.filtered = any(
                key in imported_adatas
                for key in (
                    "filtered_w_doublets",
                    "filtered_wout_doublets",
                    "filtered_w_blanks",
                    "filtered_wout_blanks",
                    "filtered_w_xenium_controls",
                    "filtered_wout_xenium_controls",
                )
            )
            st.session_state["module_score"] = sorted({
                column
                for adata in imported_adatas.values()
                for column in adata.obs.columns
                if column.startswith("module_")
            })
            st.session_state.loaded_adata_modality = is_spatial
            st.session_state.bundle_restored = True

            st.success(
                f"Restored {len(imported_adatas)} AnnData object(s) and "
                f"{len(imported_params)} parameter(s). The previous session "
                "state was replaced."
            )
        except (ValueError, OSError) as exc:
            st.error(
                f"Could not restore the bundle: {type(exc).__name__}: {exc}"
            )

    st.divider()
    
    technology_options = [
        "None — single-cell dataset",
        "MERFISH",
        "Visium HD",
        "Xenium",
        "CosMx",
    ]
    # Widget keys are dropped when the page is left, so persist under a plain key.
    spatial_technology = st.sidebar.selectbox(
        "Technology",
        technology_options,
        index=technology_options.index(st.session_state["spatial_technology"])
        if st.session_state["spatial_technology"] in technology_options
        else 0,
        key="spatial_technology_widget",
    )
    st.session_state["spatial_technology"] = spatial_technology

    is_spatial = st.session_state['spatial_technology'] != "None — single-cell dataset"
    st.session_state["is_spatial"] = is_spatial
    # Visium HD uses the single-cell QC (mitochondrial %, doublets); only
    # MERFISH is filtered on blank probes; Xenium uses negative controls.
    uses_blanks = st.session_state['spatial_technology'] == "MERFISH"
    uses_xenium = st.session_state['spatial_technology'] == "Xenium"
    
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

    # cells_file = st.file_uploader(
    #     "Upload your cells CSV file",
    #     type=["csv"],
    #     key="cells_uploader",
    # )

    # genes_file = st.file_uploader(
    #     "Upload your genes CSV file",
    #     type=["csv"],
    #     key="genes_uploader",
    # )

    # --------------------------------------------------------
    # Load uploaded files
    # --------------------------------------------------------

    if adata_file is None and st.session_state.get("bundle_restored"):
        st.session_state.bundle_restored = False
        st.session_state.uploaded_adata_name = None
        st.session_state.loaded_adata_signature = None
        st.session_state.loaded_adata_modality = None

    if adata_file is not None:
        uploaded_bytes = adata_file.getvalue()
        # The technology is part of the signature so switching it reloads the file.
        upload_signature = hashlib.sha256(
            uploaded_bytes + spatial_technology.encode()
        ).hexdigest()
        is_unchanged_upload_after_restore = (
            st.session_state.get("bundle_restored")
            and st.session_state.uploaded_adata_name == adata_file.name
            and st.session_state.loaded_adata_signature == upload_signature
        )
        if is_unchanged_upload_after_restore:
            pass
        elif (
            st.session_state.uploaded_adata_name != adata_file.name
            or st.session_state.loaded_adata_modality != is_spatial
            or st.session_state.loaded_adata_signature != upload_signature
        ):
            st.session_state.bundle_restored = False
            try:
                with tempfile.TemporaryDirectory() as temp_dir:
                    tmp_path = Path(temp_dir) / "upload.h5ad"
                    tmp_path.write_bytes(uploaded_bytes)
                    loaded_adatas = load_data(
                        str(tmp_path),
                        is_spatial=is_spatial,
                        platform=spatial_technology
                    )

                if loaded_adatas is None:
                    st.error("The .h5ad file did not contain any data.")
                    clear_downstream_results()
                    st.session_state.adatas["raw"] = None
                    st.session_state.uploaded_adata_name = None
                    st.session_state.loaded_adata_modality = None
                    st.session_state.loaded_adata_signature = None
                elif not isinstance(loaded_adatas, dict):
                    st.error(
                        "load_data() returned an unexpected object. "
                        "Expected a dictionary containing AnnData objects."
                    )
                    clear_downstream_results()
                    st.session_state.adatas["raw"] = None
                    st.session_state.uploaded_adata_name = None
                    st.session_state.loaded_adata_modality = None
                    st.session_state.loaded_adata_signature = None
                else:
                    clear_downstream_results()
                    st.session_state.adatas['raw'] = loaded_adatas['raw']
                    st.session_state.uploaded_adata_name = adata_file.name
                    st.session_state.loaded_adata_modality = is_spatial
                    st.session_state.loaded_adata_signature = upload_signature

            except Exception as e:
                clear_downstream_results()
                st.session_state.adatas['raw'] = None
                st.session_state.uploaded_adata_name = None
                st.session_state.loaded_adata_modality = None
                st.session_state.loaded_adata_signature = None
                
                if isinstance(e, ValueError):
                    st.error(f"Invalid dataset: {e}")
                else:
                    import traceback

                    st.error(
                        f"Error loading .h5ad file: "
                        f"{type(e).__name__}: {e}"
                    )
                    st.code(traceback.format_exc())

    # if cells_file is not None:
    #     try:
    #         st.session_state.cells = pd.read_csv(cells_file)
    #     except Exception as e:
    #         st.error(f"Error reading cells CSV: {e}")
    #         st.session_state.cells = None

    # if genes_file is not None:
    #     try:
    #         st.session_state.genes = pd.read_csv(genes_file)
    #     except Exception as e:
    #         st.error(f"Error reading genes CSV: {e}")
    #         st.session_state.genes = None

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
    cprobes_threshold: float = 5.0
    cwords_threshold: float = 5.0

    with st.form("filter_form"):
        min_genes_default = get_saved_number(
            "min_genes",
            config.get("min_genes", 200),
            0,
            1_000_000,
            integer=True,
        )
        min_genes = st.number_input(
            "Minimum genes",
            min_value=0,
            max_value=1_000_000,
            value=min_genes_default,
            step=10,
            help=(
                "Observations (cells, or bins/spots for spatial data) with "
                "fewer detected genes than this are REMOVED. A value equal "
                "to the threshold is kept. Defaults come from "
                "config/params.yaml and are not tuned to bin size or depth: "
                "check the QC plots first."
            ),
        )

        min_cells_default = get_saved_number(
            "min_cells",
            config.get("min_cells", 3),
            0,
            1_000_000,
            integer=True,
        )
        min_cells = st.number_input(
            "Minimum cells",
            min_value=0,
            max_value=1_000_000,
            value=min_cells_default,
            step=1,
            help=(
                "Genes detected in fewer observations than this (counted "
                "after the minimum-genes step) are REMOVED from the "
                "feature set. No observations are removed by this step."
            ),
        )
        
        if st.session_state['spatial_technology'].lower() == "merfish":
            
            blank_threshold_default = get_saved_number(
                "blank_threshold",
                float(config.get("blank_threshold", 5.0)),
                0.0,
                100.0,
            )
            blank_threshold = st.number_input(
                "Blank threshold",
                min_value=0.0,
                max_value=100.0,
                value=blank_threshold_default,
                step=1.0,
                help=(
                    "Observations with a blank-probe percentage above this "
                    "are FLAGGED in obs['is_blank'], and the page later "
                    "builds a blank-free copy. Percent, 0-100."
                ),
            )

        elif st.session_state['spatial_technology'].lower() == "xenium":
            cprobes_default = get_saved_number(
                "cprobes",
                float(config.get("cprobes", 5.0)),
                0.0,
                100.0,
            )

            cwords_default = get_saved_number(
                "cwords",
                float(config.get("cwords", 5.0)),
                0.0,
                100.0,
            )

            cprobes_threshold = st.number_input(
                "Negative DNA probes threshold",
                min_value=0.0,
                max_value=100.0,
                value=cprobes_default,
                step=1.0,
                help=(
                    "Flag cells with a negative-probe percentage above "
                    "this threshold. Percent, 0-100."
                ),
            )

            cwords_threshold = st.number_input(
                "Negative decoding threshold",
                min_value=0.0,
                max_value=100.0,
                value=cwords_default,
                step=1.0,
                help=(
                    "Flag cells with a negative-codeword percentage above "
                    "this threshold. Percent, 0-100."
                ),
            )
            
        else:

            max_mt_percentage_default = get_saved_number(
                "max_mt_percentage",
                float(config.get("max_mito_percent", 20.0)),
                0.0,
                100.0,
            )
            max_mt_percentage = st.number_input(
                "Maximum mitochondrial gene percentage",
                min_value=0.0,
                max_value=100.0,
                value=max_mt_percentage_default,
                step=1.0,
                help=(
                    "Observations whose mitochondrial percentage is at or "
                    "above this value are REMOVED (strict '<' is kept). "
                    "Percent, 0-100. Skipped, with a note in the filtering "
                    "report, when no human 'MT-' genes are annotated."
                ),
            )

            exp_doublet_rate_default = get_saved_number(
                "exp_doublet_rate",
                float(config.get("expected_doublet_rate", 0.05)),
                0.0,
                1.0,
            )
            exp_doublet_rate = st.number_input(
                "Expected doublet rate",
                min_value=0.0,
                max_value=1.0,
                value=exp_doublet_rate_default,
                step=0.01,
                help=(
                    "Fraction (0-1) passed to Scrublet. Predicted doublets "
                    "are FLAGGED here, and are removed only in the "
                    "'without doublets' copy used downstream."
                ),
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
            threshold_errors = validate_filter_thresholds(
                min_genes=min_genes,
                min_cells=min_cells,
                **(
                    {"blank_max_percentage": blank_threshold}
                    if uses_blanks
                    else {
                        "negative_probe_max_percentage": cprobes_threshold,
                        "negative_codeword_max_percentage": cwords_threshold,
                    }
                    if uses_xenium
                    else {
                        "max_mt_percentage": max_mt_percentage,
                        "expected_doublet_rate": exp_doublet_rate,
                    }
                ),
            )
            for message in threshold_errors:
                st.error(message)
            if threshold_errors:
                st.stop()

            impact = preview_filter_impact(
                adata,
                min_genes,
                min_cells,
                None if uses_blanks or uses_xenium else max_mt_percentage,
                negative_probe_max_percentage=(
                    cprobes_threshold if uses_xenium else None
                ),
                negative_codeword_max_percentage=(
                    cwords_threshold if uses_xenium else None
                ),
            )
            summary = (
                f"These thresholds remove {impact['obs_removed_min_genes']:,} "
                f"of {impact['n_obs']:,} observations (fewer than {min_genes} "
                f"genes) and {impact['vars_removed_min_cells']:,} of "
                f"{impact['n_vars']:,} genes (detected in fewer than "
                f"{min_cells} observations)"
            )
            if impact["obs_removed_mito"] is not None:
                summary += (
                    f", then {impact['obs_removed_mito']:,} observations "
                    f"with mitochondrial percentage >= {max_mt_percentage}"
                )
            remaining_after_controls = impact["obs_remaining"]
            if uses_xenium:
                flagged_controls = impact["obs_flagged_controls"] or 0
                remaining_after_controls -= flagged_controls
                summary += (
                    f", flagging {flagged_controls:,} cells with excessive "
                    f"negative-probe or negative-codeword percentages; "
                    f"{remaining_after_controls:,} pass both control thresholds"
                )
            summary += f". {impact['obs_remaining']:,} observations remain after filtering."
            if remaining_after_controls == 0 or impact["vars_remaining"] == 0:
                st.warning(summary + " Nothing would be left: relax the thresholds.")
            else:
                st.info(summary)
            try:
                with st.spinner("Applying filters..."):
                    adata_filtered = data_filtering(
                        adata,
                        min_genes=int(min_genes),
                        min_cells=int(min_cells),
                        max_mt_percentage=float(max_mt_percentage),
                        expected_doublet_rate=float(exp_doublet_rate),
                        is_spatial=uses_blanks,
                        modality=(
                            "visium_hd"
                            if spatial_technology == "Visium HD"
                            else "merfish"
                            if uses_blanks
                            else "xenium"
                            if uses_xenium
                            else "single_cell"
                        ),
                        blank_max_percentage=float(blank_threshold),
                        negative_probe_max_percentage=float(cprobes_threshold),
                        negative_codeword_max_percentage=float(cwords_threshold),
                    )

                # New filtering results invalidate any analysis built on the old ones.
                clear_downstream_results()
                st.session_state["filtering_report"] = (
                    adata_filtered.uns.get("filtering_report")
                    if adata_filtered is not None
                    else None
                )

                if adata_filtered is None:
                    st.error("Filtering returned no data.")
                    clear_downstream_results()

                elif adata_filtered.n_obs == 0:
                    st.warning(
                        "The filters removed all cells. "
                        "Try less restrictive filtering parameters."
                    )
                    if uses_blanks:
                        st.session_state.adatas[
                            "filtered_w_blanks"
                        ] = adata_filtered
                        st.session_state.adatas[
                            "filtered_wout_blanks"
                        ] = adata_filtered.copy()
                    elif uses_xenium:
                        st.session_state.adatas[
                            "filtered_w_xenium_controls"
                        ] = adata_filtered
                        st.session_state.adatas[
                            "filtered_wout_xenium_controls"
                        ] = adata_filtered.copy()
                    else:
                        st.session_state.adatas[
                            "filtered_w_doublets"
                        ] = adata_filtered
                        st.session_state.adatas[
                            "filtered_wout_doublets"
                        ] = adata_filtered.copy()
                    st.session_state.filtered = True

                else:
                    if uses_blanks:
                        if "pct_counts_blank" not in adata_filtered.obs.columns:
                            st.warning(
                                "Filtering completed, but "
                                "'pct_counts_blank' was not found. "
                                "Blank-specific results may be unavailable."
                            )
                            adata_without_blanks = adata_filtered.copy()
                            st.session_state.adatas[
                                "filtered_wout_blanks"
                            ] = adata_without_blanks
                            st.session_state.adatas[
                                "filtered_w_blanks"
                            ] = adata_filtered
                            st.session_state.filtered = True
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
                    elif uses_xenium:
                        control_mask = (
                            adata_filtered.obs["is_xenium_control"]
                            .fillna(False)
                            .astype(bool)
                        )
                        st.session_state.adatas[
                            "filtered_w_xenium_controls"
                        ] = adata_filtered
                        st.session_state.adatas[
                            "filtered_wout_xenium_controls"
                        ] = adata_filtered[~control_mask].copy()
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
                clear_downstream_results()


# ============================================================
# Main analysis area
# ============================================================

st.title("Analysis")

base_adata = get_base_adata()

if st.session_state['spatial_technology'] == "MERFISH":
    filtered_adata = get_adata("filtered_w_blanks")
    singlet_adata = get_adata("filtered_wout_blanks")
elif st.session_state['spatial_technology'] == "Xenium":
    filtered_adata = get_adata("filtered_w_xenium_controls")
    singlet_adata = get_adata("filtered_wout_xenium_controls")
elif st.session_state['spatial_technology'] == "None — single-cell dataset" or st.session_state['spatial_technology'] == "Visium HD":
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

if st.session_state['spatial_technology'] == "MERFISH":
    total_counts = base_adata.obs["total_counts"].sum()
    blank_percentage = (
        base_adata.obs["blank_counts"].sum() / total_counts * 100
        if total_counts
        else 0
    )
    st.caption(f"Blank gene percentage: {blank_percentage:.2f}%")

st.divider()

filtering_report = st.session_state.get("filtering_report")
if filtering_report is not None and st.session_state.filtered:
    with st.expander("Filtering steps (counts before / after)"):
        st.dataframe(filtering_report, hide_index=True, use_container_width=True)

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


ttl: str = (
    "Blank"
    if uses_blanks
    else "Xenium-control"
    if uses_xenium
    else "Doublet"
)

if selected == "Top Gene":
    if st.session_state['spatial_technology'] == "MERFISH":
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
    
    if st.session_state['is_spatial']:
        st.subheader("Spatially Variable Genes")
        st.markdown(
            "The following lists show the top spatially variable genes for each dataset subset."
        )
        
        gene_set1 = (
            get_spatially_variable_genes(base_adata, n_top_genes=10)
            if base_adata.n_obs and base_adata.n_vars
            else []
        )
        gene_set2 = (
            get_spatially_variable_genes(filtered_adata, n_top_genes=10)
            if filtered_adata.n_obs and filtered_adata.n_vars
            else []
        )
        gene_set3 = (
            get_spatially_variable_genes(singlet_adata, n_top_genes=10)
            if singlet_adata.n_obs and singlet_adata.n_vars
            else []
        )

        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.subheader("Subset / Raw")
            color_by = st.selectbox(
                "Module Score",
                options=gene_set1,
                key="svg1",
                placeholder="Select a Spatially Variable Gene"
            )
            
            if color_by is not None:
                render_interactive(base_adata, "spatial", key="svg1_map", gene=color_by, height=450)
            else:
                st.info("No spatially variable genes are available.")

        with col2:
            st.subheader(f"Filtered with {ttl}s")
            color_by = st.selectbox(
                "Module Score",
                options=gene_set2,
                key="svg2",
                placeholder="Select a Spatially Variable Gene"
            )
            if color_by is not None:
                render_interactive(filtered_adata, "spatial", key="svg2_map", gene=color_by, height=450)
            else:
                st.info("No spatially variable genes are available.")

        with col3:
            st.subheader(f"Filtered without {ttl}s")
            color_by = st.selectbox(
                "Module Score",
                options=gene_set3,
                key="svg3",
                placeholder="Select a Spatially Variable Gene"
            )
            if color_by is not None:
                render_interactive(singlet_adata, "spatial", key="svg3_map", gene=color_by, height=450)
            else:
                st.info("No spatially variable genes are available.")

    else:
        
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
        fig_one = plot_qc_metrics(base_adata, is_spatial=is_spatial, platform=st.session_state['spatial_technology'])
        fig_two = plot_qc_metrics(filtered_adata, is_spatial=is_spatial, platform=st.session_state['spatial_technology'])
        fig_three = plot_qc_metrics(singlet_adata, is_spatial=is_spatial, platform=st.session_state['spatial_technology'])
        
        if st.session_state['spatial_technology'] in ("MERFISH", "CosMx"):
            color_by = ["is_blank", "log1p_n_genes_by_counts", "log1p_total_counts"]
        elif uses_xenium:
            color_by = [
                "negative_probe_pct",
                "negative_codeword_pct",
                "log1p_n_genes_by_counts",
            ]
        else:
            color_by = ["log1p_n_genes_by_counts", "log1p_total_counts", "pct_counts_mt"]

        st.subheader("Subset / Raw")
        if is_spatial:
            col1, col2 = st.columns(2)
            with col1:
                st.pyplot(fig_one, use_container_width=True)
            with col2:
                spatial_qc_plt = plot_spatial(
                    base_adata,
                    color_by=color_by
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
                    color_by=color_by,
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
                    color_by=color_by,
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

    if uses_blanks:

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
    elif uses_xenium:
        predicted_controls = get_bool_obs_column(
            filtered_adata,
            "is_xenium_control",
        )

        if predicted_controls is None:
            st.info("No Xenium negative-control flag information is available.")
        elif len(predicted_controls) == 0:
            st.info(
                "No cells remain after filtering, so control statistics "
                "cannot be calculated."
            )
        else:
            high_controls = int(predicted_controls.sum())
            passing_controls = len(predicted_controls) - high_controls
            high_percentage = high_controls / len(predicted_controls) * 100

            fig = create_percentage_bar(
                values=[passing_controls, high_controls],
                labels=["Within thresholds", "Above threshold"],
                total_label="Cells",
                colors=["#4CAF50", "#E74C3C"],
                counts=[passing_controls, high_controls],
            )
            st.subheader("Xenium Negative-Control Statistics")
            st.plotly_chart(fig, use_container_width=True)
            st.caption(
                f"{high_controls:,} of {len(predicted_controls):,} cells "
                f"({high_percentage:.2f}%) exceed at least one control threshold."
            )
    else:

        # ----------------------------------------------------
        # Predicted doublets
        # ----------------------------------------------------
        
        filtered_adata.obs['is_doublet'] = filtered_adata.obs['predicted_doublet'] > np.median(filtered_adata.obs['predicted_doublet'])

        predicted_doublet = get_bool_obs_column(
            filtered_adata,
            "is_doublet",
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

if "params" not in st.session_state or not isinstance(
    st.session_state.params,
    dict,
):
    st.session_state.params = {}

st.session_state.params.update({
    "min_genes": min_genes,
    "min_cells": min_cells,
    "blank_threshold": blank_threshold,
    "cprobes": cprobes_threshold,
    "cwords": cwords_threshold,
    "max_mt_percentage": max_mt_percentage,
    "exp_doublet_rate": exp_doublet_rate,
    "is_spatial": is_spatial,
    "spatial_technology": spatial_technology,
})