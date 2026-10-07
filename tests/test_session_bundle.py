import io
import json
import zipfile

import numpy as np
import pytest

from src.session_bundle import create_session_bundle, load_session_bundle


def test_session_bundle_round_trips_data_and_parameters(small_adata):
    objects = {
        "raw/processed": small_adata,
        "raw\\processed": small_adata.copy(),
    }
    params = {"min_genes": 5, "is_spatial": True}

    bundle = create_session_bundle(objects, params)
    loaded_adatas, loaded_params = load_session_bundle(bundle)

    assert loaded_params == params
    assert list(loaded_adatas) == list(objects)
    for name, original in objects.items():
        np.testing.assert_array_equal(loaded_adatas[name].X, original.X)
        assert loaded_adatas[name].obs_names.equals(original.obs_names)

    with zipfile.ZipFile(io.BytesIO(bundle)) as archive:
        manifest = json.loads(archive.read("manifest.json"))
    assert len(set(manifest["adata_files"].values())) == len(objects)


def test_session_bundle_reads_legacy_manifest(small_adata, tmp_path):
    adata_path = tmp_path / "raw.h5ad"
    small_adata.write_h5ad(adata_path)
    bundle_buffer = io.BytesIO()
    with zipfile.ZipFile(bundle_buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(
            "manifest.json",
            json.dumps({
                "format": "BioWorkbench",
                "version": 1,
                "params": ["min_genes"],
                "adatas": ["raw"],
            }),
        )
        archive.writestr("params.json", json.dumps({"min_genes": 10}))
        archive.write(adata_path, "adatas/raw.h5ad")

    loaded_adatas, loaded_params = load_session_bundle(bundle_buffer.getvalue())

    assert loaded_params == {"min_genes": 10}
    np.testing.assert_array_equal(loaded_adatas["raw"].X, small_adata.X)


def test_session_bundle_rejects_non_zip_data():
    with pytest.raises(ValueError, match="valid ZIP"):
        load_session_bundle(b"not a ZIP archive")


def test_session_bundle_rejects_manifest_path_traversal():
    bundle_buffer = io.BytesIO()
    with zipfile.ZipFile(bundle_buffer, "w") as archive:
        archive.writestr(
            "manifest.json",
            json.dumps({
                "format": "BioWorkbench",
                "version": 1,
                "params": [],
                "adatas": ["raw"],
                "adata_files": {"raw": "../raw.h5ad"},
            }),
        )
        archive.writestr("params.json", "{}")

    with pytest.raises(ValueError, match="Invalid AnnData archive path"):
        load_session_bundle(bundle_buffer.getvalue())


def test_session_bundle_rejects_unreadable_anndata_member():
    bundle_buffer = io.BytesIO()
    with zipfile.ZipFile(bundle_buffer, "w") as archive:
        archive.writestr(
            "manifest.json",
            json.dumps({
                "format": "BioWorkbench",
                "version": 1,
                "params": [],
                "adatas": ["raw"],
                "adata_files": {"raw": "adatas/raw.h5ad"},
            }),
        )
        archive.writestr("params.json", "{}")
        archive.writestr("adatas/raw.h5ad", "not an H5AD file")

    with pytest.raises(ValueError, match="could not be read"):
        load_session_bundle(bundle_buffer.getvalue())


def test_session_bundle_requires_anndata_objects(small_adata):
    with pytest.raises(ValueError, match="not an AnnData"):
        create_session_bundle({"raw": object()}, {})
