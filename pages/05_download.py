import os
import streamlit as st
import tempfile

import pandas as pd
import scanpy as sc


# -----------------------------------------------------------------------------
# Downloads
# -----------------------------------------------------------------------------

st.divider()
st.subheader("Download")

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

st.markdown("### Parameters")

params = st.session_state.get("params")

if params:
    params_df = pd.DataFrame(
        list(params.items()),
        columns=["parameter", "value"],
    )

    csv_data = params_df.to_csv(index=False).encode("utf-8")

    st.download_button(
        label="Download parameters",
        data=csv_data,
        file_name="parameters.csv",
        mime="text/csv",
        use_container_width=True,
    )
else:
    st.info("No parameters are available to download.")


# -----------------------------------------------------------------------------
# AnnData objects
# -----------------------------------------------------------------------------

st.markdown("### AnnData")

adatas = st.session_state.get("adatas")

if isinstance(adatas, dict):

    available_adatas = {
        name: adata
        for name, adata in adatas.items()
        if adata is not None
    }

    if available_adatas:

        for name, adata in available_adatas.items():

            col1, col2 = st.columns([3, 1])

            with col1:
                st.write(f"**{name}**")

                if hasattr(adata, "n_obs") and hasattr(adata, "n_vars"):
                    st.caption(
                        f"{adata.n_obs:,} cells × {adata.n_vars:,} genes"
                    )

            with col2:

                # AnnData.write_h5ad() needs a filename, so use a
                # temporary in-memory BytesIO buffer.
                buffer = io.BytesIO()

                adata.write_h5ad(buffer)

                buffer.seek(0)

                st.download_button(
                    label="Download",
                    data=buffer,
                    file_name=f"{name}.h5ad",
                    mime="application/octet-stream",
                    key=f"download_adata_{name}",
                    use_container_width=True,
                )

    else:
        st.info("No AnnData objects are available.")

else:
    st.info("No AnnData objects are available.")
