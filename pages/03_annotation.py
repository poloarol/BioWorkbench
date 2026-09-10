import json

import celltypist
import streamlit as st

from celltypist import models
from matplotlib import pyplot as plt

from src.plotting import plot_umap

@st.cache_data
def get_celltypist_models():
    """Retrieve the available CellTypist model catalogue."""

    return models.models_description()

selected_model = None
model_info = get_celltypist_models()
adata = st.session_state.adatas['clustered']

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
        adata.obs["custom_annotation"] = adata.obs["annotation"].map(custom_annotations)

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

    st.subheader("Annotation UMAPs")

    # Display two plots per row
    for i in range(0, len(available_annotations), 2):

        columns = available_annotations[i:i + 2]

        plot_cols = st.columns(2)

        for plot_col, color_by in zip(plot_cols, columns):

            with plot_col:

                st.caption(color_by)

                fig = plot_umap(
                    adata,
                    color_by=color_by,
                )

                st.pyplot(
                    fig,
                    use_container_width=True,
                )

                plt.close(fig)
