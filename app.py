import streamlit as st

st.set_page_config(
    page_title="BioWorkbench",
    page_icon="🧬",
    layout="wide"
)

def app_header():
    st.markdown("""
        <h1 style="font-size: 2.5rem; margin-bottom: 0.5rem;">
            BioWorkbench - A lightweight single-cell genomics analysis platform
           </h1>
        """,
        unsafe_allow_html=True,)

app_header()

pages = [
    st.Page("pages/01_filtering.py", title="01 - Cell and Gene Filtering"),
    st.Page("pages/02_clustering.py", title="02 - Clustering & Exploratory Analysis"),
    st.Page("pages/03_annotation.py", title="03 - Annotation"),
    st.Page("pages/04_visualization.py", title="04 - Visualizer"),
    st.Page("pages/05_download.py", title="05 - Download")
]

pg = st.navigation(pages)
pg.run()
