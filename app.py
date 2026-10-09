import base64
from pathlib import Path
import streamlit as st

st.set_page_config(
    page_title="BioWorkbench",
    page_icon="docs/images/bioworkbench-logo-small.png",
    layout="wide",
)

st.logo(
    "docs/images/bioworkbench-logo.png",
    icon_image="docs/images/bioworkbench-logo-small.png",
    size="large",
)

# Remove Streamlit's default content width/padding
st.markdown(
    """
    <style>
        .stAppViewContainer .main .block-container {
            max-width: 100%;
            padding: 0;
        }

        /* Sidebar logo: lower, and ~90% of the sidebar width */
        [data-testid="stSidebarHeader"] {
            height: auto;
            display: block;
            width: 100%;
            padding: 2.5rem 0 1rem 0;
        }

        [data-testid="stSidebarHeader"] > div,
        [data-testid="stSidebarHeader"] a {
            display: block;
            width: 100% !important;
            max-width: 100%;
            height: auto;
        }

        [data-testid="stSidebar"] {
            container-type: inline-size;
        }

        [data-testid="stSidebarHeader"] img {
            display: block;
            width: 90cqw !important;
            height: auto !important;
            max-width: none !important;
            max-height: none !important;
            margin: 0 0 0 calc((100% - 90cqw) / 2);
        }

        /* Make the header edge-to-edge */
        .bioworkbench-header {
            width: 100vw;
            margin-left: calc(50% - 50vw);
            margin-right: calc(50% - 50vw);
            height: 220px;

            display: flex;
            align-items: center;
            justify-content: center;

            overflow: hidden;
            margin-bottom: 2rem;
        }

        .bioworkbench-header img {
            width: 100%;
            height: 100%;
            object-fit: contain;
        }
    </style>
    """,
    unsafe_allow_html=True,
)


def app_header():
    logo_path = Path(__file__).parent / "docs" / "images" / "bioworkbench-logo.png"

    with open(logo_path, "rb") as f:
        logo_base64 = base64.b64encode(f.read()).decode()

    st.markdown(
        f"""
        <div class="bioworkbench-header">
            <img
                src="data:image/png;base64,{logo_base64}"
                alt="BioWorkbench"
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
