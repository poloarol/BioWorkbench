from pathlib import Path

from streamlit.testing.v1 import AppTest


def test_application_starts_without_streamlit_exceptions():
    app_path = Path(__file__).parents[1] / "app.py"

    app = AppTest.from_file(str(app_path)).run(timeout=60)

    assert not app.exception


def test_analysis_pages_handle_missing_prior_workflow_data():
    app_path = Path(__file__).parents[1] / "app.py"
    app = AppTest.from_file(str(app_path)).run(timeout=60)

    for page in (
        "pages/02_clustering.py",
        "pages/03_annotation.py",
        "pages/04_visualization.py",
        "pages/05_download.py",
    ):
        app.switch_page(page).run(timeout=60)
        assert not app.exception, page


def test_clustering_page_prompts_when_filtering_has_not_run():
    app_path = Path(__file__).parents[1] / "app.py"
    app = AppTest.from_file(str(app_path)).run(timeout=60)

    app.switch_page("pages/02_clustering.py").run(timeout=60)

    assert not app.exception
    assert any(
        "Apply filtering" in element.value
        for element in app.info
    )
