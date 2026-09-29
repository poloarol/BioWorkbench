import base64

import streamlit as st
from pathlib import Path


st.set_page_config(
    page_title="BioWorkbench",
    page_icon="docs/images/bioworkbench-logo-small.png",
    # layout="wide"
)

def app_header():
    logo_path = Path(__file__).parent / "docs" / "images" / "bioworkbench-logo.png"

    with open(logo_path, "rb") as f:
        logo_base64 = base64.b64encode(f.read()).decode()

    st.markdown(
        f"""
        <div style="
            width: 100%;
            height: 220px;
            display: flex;
            align-items: center;
            justify-content: center;
            overflow: hidden;
            margin-bottom: 2rem;
        ">
            <img
                src="data:image/png;base64,{logo_base64}"
                style="
                    width: 90%;
                    max-height: 210px;
                    object-fit: contain;
                "
            >
        </div>
        """,
        unsafe_allow_html=True,
    )

app_header()

pages = [
    st.Page("pages/01_filtering.py", title="01 - Cell and Gene Filtering"),
    st.Page("pages/02_clustering.py", title="02 - Exploratory Analysis & Clustering"),
    st.Page("pages/03_annotation.py", title="03 - Celltype Annotation"),
    st.Page("pages/04_visualization.py", title="04 - Gene & Module Visualizer"),
    st.Page("pages/05_download.py", title="05 - Download")
]

pg = st.navigation(pages)
pg.run()
