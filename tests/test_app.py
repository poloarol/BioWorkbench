from pathlib import Path

from streamlit.testing.v1 import AppTest


def test_application_starts_without_streamlit_exceptions():
    app_path = Path(__file__).parents[1] / "app.py"

    app = AppTest.from_file(str(app_path)).run(timeout=60)

    assert not app.exception
