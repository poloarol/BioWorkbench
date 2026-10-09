import json

import celltypist
import streamlit as st

from celltypist import models
from matplotlib import pyplot as plt

from src.interactive import render_interactive

if "is_spatial" not in st.session_state:
    st.session_state["is_spatial"] = False


adatas = st.session_state.get("adatas", {})
adata = adatas.get("clustered") if isinstance(adatas, dict) else None

if adata is None:
    st.info("Run PCA and clustering on the Clustering page before annotation.")
    st.stop()


@st.cache_data
def get_celltypist_models():
    """Retrieve the available CellTypist model catalogue."""

    return models.models_description()

selected_model = None
model_info = get_celltypist_models()

with st.sidebar:

    st.subheader("Celltypist Reference Models")

    selected_model = st.selectbox(
        "CellTypist model",
        options=model_info["model"].tolist(),
        key="celltypist_model",
    )

    selected_model_info = model_info[
        model_info["model"] == selected_model
    ].iloc[0]

    with st.expander("Model description", expanded=True):
        st.write(
            selected_model_info["description"]
        )
    # st.session_state.params['model'] = selected_model
        
    st.subheader("Custom Annotations")

    annotation_file = st.file_uploader(
        "Upload annotations",
        type=["json"],
        key="annotation_upload",
        help="JSON file mapping cell IDs to annotation labels.",
    )

    st.divider()

    run_annotation = st.button(
        "Run Annotation",
        type="primary",
        use_container_width=True,
    )

if run_annotation:

    with st.spinner(
        f"Annotating cells with {selected_model}..."
    ):

        try:

            # Load the selected model
            model = models.Model.load(
                model=selected_model,
            )

            # Run annotation
            predictions = celltypist.annotate(
                adata,
                model=model,
                majority_voting=True,
            )
            
            tmp = predictions.to_adata()

            # Store predicted labels
            adata.obs["celltypist_predicted_labels"] = (
                tmp.obs[
                    "predicted_labels"
                ].values
            )

            adata.obs["celltypist_majority_voting"] = (
                tmp.obs[
                    "majority_voting"
                ].values
            )
            
            adata.obs["celltypist_conf_score"] = (
                            tmp.obs[
                                "conf_score"
                            ].values
                        )

            st.success(
                f"Annotation completed using {selected_model}."
            )

        except Exception as e:
            st.error(
                f"Annotation failed: {e}"
            )


if annotation_file is not None:

    try:
        # Load custom annotations from the uploaded JSON file
        custom_annotations = json.load(annotation_file)
        if not isinstance(custom_annotations, dict):
            raise ValueError("Expected a JSON object mapping cell IDs to labels.")
        adata.obs["custom_annotation"] = adata.obs_names.map(
            custom_annotations
        )

    except Exception as e:
        st.error(
            f"Failed to load custom annotations: {e}"
        )


# -----------------------------------------------------------------------------
# Annotation UMAPs
# -----------------------------------------------------------------------------

annotation_columns = [
    "custom_annotation",
    "celltypist_predicted_labels",
    "celltypist_majority_voting",
    "celltypist_conf_score",
]

available_annotations = [
    column
    for column in annotation_columns
    if column in adata.obs.columns
]

if "X_umap" in adata.obsm and available_annotations:

    st.divider()
    st.subheader("Annotation Visualizations")

    for color_by in available_annotations:

        st.markdown(f"### {color_by}")

        if st.session_state["is_spatial"]:

            plot_cols = st.columns(2)

            with plot_cols[0]:
                render_interactive(
                    adata, "umap", key=f"ann_umap_{color_by}", color_by=color_by
                )

            with plot_cols[1]:
                render_interactive(
                    adata, "spatial", key=f"ann_spatial_{color_by}", color_by=color_by
                )

        else:
            render_interactive(
                adata, "umap", key=f"ann_umap_{color_by}", color_by=color_by
            )

st.session_state.adatas['annotated'] = adata
