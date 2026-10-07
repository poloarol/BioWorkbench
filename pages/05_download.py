import streamlit as st

# -----------------------------------------------------------------------------
# Downloads
# -----------------------------------------------------------------------------

st.divider()
st.subheader("Download")

params = st.session_state.get("params")
adatas = st.session_state.get("adatas")


# -----------------------------------------------------------------------------
# Selection
# -----------------------------------------------------------------------------

selected_params = []
selected_adatas = []

col_params, col_adatas = st.columns(2)


# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

with col_params:
    st.markdown("### Parameters")

    if isinstance(params, dict) and params:
        select_all_params = st.checkbox(
            "Select all parameters",
            value=True,
            key="select_all_params",
        )

        key_col, value_col = st.columns([1, 2])

        with key_col:
            st.markdown("**Key**")

        with value_col:
            st.markdown("**Value**")

        for parameter, value in params.items():
            key_col, value_col = st.columns([1, 2])

            with key_col:
                selected = st.checkbox(
                    str(parameter),
                    value=select_all_params,
                    key=f"download_param_{parameter}",
                )

            with value_col:
                st.write(str(value))

            if selected:
                selected_params.append(parameter)

    else:
        st.info("No parameters are available.")


# -----------------------------------------------------------------------------
# AnnData
# -----------------------------------------------------------------------------

with col_adatas:
    st.markdown("### AnnData")

    if isinstance(adatas, dict):
        available_adatas = {
            name: adata
            for name, adata in adatas.items()
            if adata is not None
        }

        if available_adatas:
            select_all_adatas = st.checkbox(
                "Select all AnnData objects",
                value=True,
                key="select_all_adatas",
            )

            name_col, size_col = st.columns([1, 1])

            with name_col:
                st.markdown("**Object**")

            with size_col:
                st.markdown("**Size**")

            for name, adata in available_adatas.items():
                name_col, size_col = st.columns([1, 1])

                with name_col:
                    selected = st.checkbox(
                        str(name),
                        value=select_all_adatas,
                        key=f"download_adata_{name}",
                    )

                with size_col:
                    if hasattr(adata, "n_obs") and hasattr(adata, "n_vars"):
                        st.write(f"{adata.n_obs:,} × {adata.n_vars:,}")
                    else:
                        st.write("AnnData")

                if selected:
                    selected_adatas.append(name)

        else:
            st.info("No AnnData objects are available.")

    else:
        st.info("No AnnData objects are available.")


# -----------------------------------------------------------------------------
# Download
# -----------------------------------------------------------------------------

st.divider()

n_params = len(selected_params)
n_adatas = len(selected_adatas)
total_selected = n_params + n_adatas

st.write(
    f"Selected: **{n_params} parameter(s)** and "
    f"**{n_adatas} AnnData object(s)**"
)

if total_selected == 0:
    st.info("Select at least one item to download.")

else:
    selected_params_data = {
        str(key): params[key]
        for key in selected_params
    }
    selected_adata_data = {
        str(name): adatas[name]
        for name in selected_adatas
    }
    zip_data = create_session_bundle(
        selected_adata_data,
        selected_params_data,
    )

    # -------------------------------------------------------------------------
    # Download
    # -------------------------------------------------------------------------

    st.download_button(
        label="Download selected data (.wkb)",
        data=zip_data,
        file_name="bioworkbench.wkb",
        mime="application/zip",
        use_container_width=True,
    )