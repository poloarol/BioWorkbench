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


def test_upload_page_offers_bundle_restore_in_sidebar():
    app_path = Path(__file__).parents[1] / "app.py"
    app = AppTest.from_file(str(app_path)).run(timeout=60)

    app.switch_page("pages/01_filtering.py").run(timeout=60)

    assert not app.exception
    assert any(button.label == "Restore bundle" for button in app.button)


def test_filtering_page_uses_restored_numeric_parameters():
    app_path = Path(__file__).parents[1] / "app.py"
    app = AppTest.from_file(str(app_path)).run(timeout=60)
    app.session_state["params"] = {
        "min_genes": 123,
        "min_cells": 7,
        "max_mt_percentage": 12.0,
        "exp_doublet_rate": 0.08,
    }

    app.switch_page("pages/01_filtering.py").run(timeout=60)

    assert not app.exception
    assert next(
        element
        for element in app.number_input
        if element.label == "Minimum genes"
    ).value == 123
    assert next(
        element
        for element in app.number_input
        if element.label == "Minimum cells"
    ).value == 7
    assert (
        next(
            element
            for element in app.number_input
            if element.label == "Maximum mitochondrial gene percentage"
        ).value
        == 12.0
    )
