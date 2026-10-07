import io
import json
import re
import tempfile
import zipfile
import zlib
from pathlib import Path, PurePosixPath
from typing import BinaryIO

import anndata as ad
import scanpy as sc


_FORMAT_NAME = "BioWorkbench"
_FORMAT_VERSION = 1
_MAX_ARCHIVE_ENTRIES = 1024
_MAX_UNCOMPRESSED_SIZE = 2 * 1024**3
_MAX_PARAMS_SIZE = 10 * 1024**2
_CHUNK_SIZE = 1024**2


def _safe_filename(name: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._")
    return sanitized or "adata"


def create_session_bundle(
    adatas: dict[str, ad.AnnData],
    params: dict[str, object],
) -> bytes:
    """Serialize selected AnnData objects and parameters as a .wkb archive."""
    if not isinstance(adatas, dict) or not isinstance(params, dict):
        raise ValueError("AnnData objects and parameters must be dictionaries.")

    file_mapping: dict[str, str] = {}
    for index, (name, value) in enumerate(adatas.items()):
        if not isinstance(name, str) or not name:
            raise ValueError("AnnData object names must be non-empty strings.")
        if not isinstance(value, ad.AnnData):
            raise ValueError(f"Selected object {name!r} is not an AnnData object.")
        file_mapping[name] = f"adatas/{index}_{_safe_filename(name)}.h5ad"

    manifest = {
        "format": _FORMAT_NAME,
        "version": _FORMAT_VERSION,
        "params": list(params),
        "adatas": list(adatas),
        "adata_files": file_mapping,
    }

    archive_buffer = io.BytesIO()
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        with zipfile.ZipFile(
            archive_buffer,
            mode="w",
            compression=zipfile.ZIP_DEFLATED,
        ) as archive:
            archive.writestr(
                "manifest.json",
                json.dumps(manifest, indent=2, default=str),
            )
            archive.writestr(
                "params.json",
                json.dumps(params, indent=2, default=str),
            )

            for index, (name, adata) in enumerate(adatas.items()):
                adata_path = temp_path / f"{index}.h5ad"
                adata.write_h5ad(adata_path)
                archive.write(adata_path, arcname=file_mapping[name])

    return archive_buffer.getvalue()


def _read_json_member(
    archive: zipfile.ZipFile,
    info: zipfile.ZipInfo,
    *,
    max_size: int,
) -> object:
    if info.file_size > max_size:
        raise ValueError(f"Archive member {info.filename!r} is too large.")
    try:
        with archive.open(info) as member:
            data = member.read(max_size + 1)
        if len(data) > max_size:
            raise ValueError(f"Archive member {info.filename!r} is too large.")
        return json.loads(data)
    except (
        EOFError,
        UnicodeDecodeError,
        json.JSONDecodeError,
        NotImplementedError,
        RuntimeError,
        zipfile.BadZipFile,
        zlib.error,
    ) as exc:
        raise ValueError(f"Archive member {info.filename!r} is not valid JSON.") from exc


def _validate_bundle(
    archive: zipfile.ZipFile,
) -> tuple[dict[str, object], dict[str, zipfile.ZipInfo]]:
    infos = archive.infolist()
    if len(infos) > _MAX_ARCHIVE_ENTRIES:
        raise ValueError("The bundle contains too many archive entries.")
    if len({info.filename for info in infos}) != len(infos):
        raise ValueError("The bundle contains duplicate archive entries.")

    total_size = sum(info.file_size for info in infos)
    if total_size > _MAX_UNCOMPRESSED_SIZE:
        raise ValueError("The uncompressed bundle exceeds the 2 GiB size limit.")
    if any(info.flag_bits & 0x1 for info in infos):
        raise ValueError("Encrypted bundle entries are not supported.")

    info_by_name = {info.filename: info for info in infos}
    manifest_info = info_by_name.get("manifest.json")
    params_info = info_by_name.get("params.json")
    if manifest_info is None or params_info is None:
        raise ValueError("The bundle must contain manifest.json and params.json.")

    manifest = _read_json_member(
        archive,
        manifest_info,
        max_size=_MAX_PARAMS_SIZE,
    )
    if not isinstance(manifest, dict):
        raise ValueError("The bundle manifest must be a JSON object.")
    if manifest.get("format") != _FORMAT_NAME:
        raise ValueError("This archive is not a BioWorkbench export bundle.")
    version = manifest.get("version")
    if type(version) is not int or version != _FORMAT_VERSION:
        raise ValueError(f"Unsupported BioWorkbench bundle version: {version!r}.")

    names = manifest.get("adatas")
    parameter_names = manifest.get("params")
    if (
        not isinstance(names, list)
        or any(not isinstance(name, str) or not name for name in names)
        or len(names) != len(set(names))
    ):
        raise ValueError("The manifest AnnData list is invalid.")
    if (
        not isinstance(parameter_names, list)
        or any(not isinstance(name, str) or not name for name in parameter_names)
        or len(parameter_names) != len(set(parameter_names))
    ):
        raise ValueError("The manifest parameter list is invalid.")

    params = _read_json_member(
        archive,
        params_info,
        max_size=_MAX_PARAMS_SIZE,
    )
    if not isinstance(params, dict) or set(params) != set(parameter_names):
        raise ValueError("params.json does not match the manifest.")

    file_mapping = manifest.get("adata_files")
    if file_mapping is None:
        file_mapping = {
            name: f"adatas/{name.replace('/', '_').replace(chr(92), '_')}.h5ad"
            for name in names
        }
    if (
        not isinstance(file_mapping, dict)
        or set(file_mapping) != set(names)
        or any(not isinstance(path, str) for path in file_mapping.values())
    ):
        raise ValueError("The manifest AnnData file mapping is invalid.")

    expected_members = {"manifest.json", "params.json"}
    adata_infos: dict[str, zipfile.ZipInfo] = {}
    for name, member_name in file_mapping.items():
        member_path = PurePosixPath(member_name)
        if (
            member_path.is_absolute()
            or ".." in member_path.parts
            or member_path.parts[0:1] != ("adatas",)
            or len(member_path.parts) != 2
            or member_path.suffix != ".h5ad"
            or "\\" in member_name
            or member_path.as_posix() != member_name
        ):
            raise ValueError(f"Invalid AnnData archive path: {member_name!r}.")
        info = info_by_name.get(member_name)
        if info is None or info.is_dir():
            raise ValueError(f"AnnData object {name!r} is missing from the bundle.")
        if info.compress_type not in {
            zipfile.ZIP_STORED,
            zipfile.ZIP_DEFLATED,
            zipfile.ZIP_BZIP2,
            zipfile.ZIP_LZMA,
        }:
            raise ValueError(f"AnnData object {name!r} uses unsupported compression.")
        expected_members.add(member_name)
        adata_infos[name] = info

    if set(info_by_name) != expected_members:
        raise ValueError("The bundle contains unexpected archive entries.")

    return params, adata_infos


def load_session_bundle(
    source: bytes | bytearray | BinaryIO,
) -> tuple[dict[str, ad.AnnData], dict[str, object]]:
    """Validate and load a BioWorkbench .wkb bundle without extracting paths."""
    if isinstance(source, (bytes, bytearray)):
        source = io.BytesIO(source)

    try:
        archive = zipfile.ZipFile(source)
    except (zipfile.BadZipFile, OSError) as exc:
        raise ValueError("The uploaded file is not a valid ZIP bundle.") from exc

    with archive:
        params, adata_infos = _validate_bundle(archive)
        loaded_adatas: dict[str, ad.AnnData] = {}
        with tempfile.TemporaryDirectory() as temp_dir:
            for index, (name, info) in enumerate(adata_infos.items()):
                adata_path = Path(temp_dir) / f"{index}.h5ad"
                copied_size = 0
                try:
                    with (
                        archive.open(info) as source_member,
                        adata_path.open("wb") as target,
                    ):
                        while chunk := source_member.read(_CHUNK_SIZE):
                            copied_size += len(chunk)
                            if copied_size > info.file_size:
                                raise ValueError(
                                    f"AnnData object {name!r} exceeds its "
                                    "declared archive size."
                                )
                            target.write(chunk)
                except (
                    EOFError,
                    NotImplementedError,
                    OSError,
                    RuntimeError,
                    zipfile.BadZipFile,
                    zlib.error,
                ) as exc:
                    raise ValueError(
                        f"AnnData object {name!r} could not be read from the bundle."
                    ) from exc
                if copied_size != info.file_size:
                    raise ValueError(
                        f"AnnData object {name!r} is truncated in the bundle."
                    )
                try:
                    loaded_adatas[name] = sc.read_h5ad(adata_path)
                except (OSError, ValueError, KeyError) as exc:
                    raise ValueError(
                        f"AnnData object {name!r} could not be read from the bundle."
                    ) from exc

        return loaded_adatas, params
