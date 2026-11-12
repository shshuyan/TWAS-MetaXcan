"""Frozen public-dataset registry and checksum-verifiable downloads."""
from __future__ import annotations

import json
from pathlib import Path, PurePosixPath
import shutil
from urllib.parse import urlparse
from urllib.request import Request, urlopen

import yaml

from .screen_io import sha256


REQUIRED = {"id", "url", "filename", "role"}


def load_registry(path):
    """Load and validate a frozen registry without accessing the network."""
    path = Path(path).resolve()
    registry = yaml.safe_load(path.read_text())
    if not isinstance(registry, dict) or registry.get("schema_version") != 1:
        raise ValueError("Dataset registry requires schema_version: 1")
    if not registry.get("frozen_on") or not isinstance(registry.get("datasets"), list):
        raise ValueError("Dataset registry requires frozen_on and datasets")
    seen = set()
    for item in registry["datasets"]:
        if not isinstance(item, dict):
            raise ValueError("Each dataset entry must be a mapping")
        missing = sorted(REQUIRED - set(item))
        if missing:
            raise ValueError("Dataset entry missing " + ", ".join(missing))
        if item["id"] in seen:
            raise ValueError("Duplicate dataset id: " + str(item["id"]))
        seen.add(item["id"])
        rel = PurePosixPath(str(item["filename"]))
        if rel.is_absolute() or ".." in rel.parts or not rel.name:
            raise ValueError("Dataset filename must be a safe relative path: " + str(rel))
        if item.get("bytes") is not None and (not isinstance(item["bytes"], int) or item["bytes"] < 0):
            raise ValueError("Dataset bytes must be a nonnegative integer")
        checksum = item.get("sha256")
        if checksum is not None and (len(str(checksum)) != 64 or any(c not in "0123456789abcdef" for c in str(checksum).lower())):
            raise ValueError("Dataset sha256 must contain 64 hexadecimal characters")
    registry["_registry_path"] = str(path)
    return registry


def _matches(path, item):
    if item.get("bytes") is not None and path.stat().st_size != item["bytes"]:
        return False
    return item.get("sha256") is None or sha256(path) == item["sha256"].lower()


def _download(url, partial):
    """Resume local/HTTP downloads; restart when a server ignores Range."""
    offset = partial.stat().st_size if partial.exists() else 0
    parsed = urlparse(url)
    if parsed.scheme == "file":
        source = Path(parsed.path)
        if offset > source.stat().st_size:
            partial.unlink()
            offset = 0
        with source.open("rb") as src, partial.open("ab" if offset else "wb") as dst:
            src.seek(offset)
            shutil.copyfileobj(src, dst, length=1024 * 1024)
        return
    request = Request(url, headers={"Range": f"bytes={offset}-"} if offset else {})
    with urlopen(request) as response:
        append = offset > 0 and getattr(response, "status", None) == 206
        with partial.open("ab" if append else "wb") as dst:
            shutil.copyfileobj(response, dst, length=1024 * 1024)


def download_snapshot(registry, destination, selected_ids=None):
    """Download a selected immutable snapshot and write its audit manifest."""
    destination = Path(destination).resolve()
    selected = set(selected_ids) if selected_ids is not None else None
    known = {item["id"] for item in registry["datasets"]}
    if selected is not None and selected - known:
        raise ValueError("Unknown dataset ids: " + ", ".join(sorted(selected - known)))
    records = []
    for item in registry["datasets"]:
        if selected is not None and item["id"] not in selected:
            continue
        target = destination.joinpath(*PurePosixPath(item["filename"]).parts)
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if not _matches(target, item):
                raise ValueError(f"Existing file does not match registry; refusing to overwrite: {target}")
            status = "verified_existing"
        else:
            partial = target.with_name(target.name + ".part")
            _download(item["url"], partial)
            if not _matches(partial, item):
                raise ValueError(f"Downloaded file does not match registry; partial retained: {partial}")
            partial.replace(target)
            status = "downloaded"
        records.append({
            "id": item["id"], "role": item["role"], "disease": item.get("disease"),
            "url": item["url"], "path": str(target.relative_to(destination)),
            "bytes": target.stat().st_size, "sha256": sha256(target), "status": status,
        })
    destination.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema_version": 1,
        "registry_frozen_on": str(registry["frozen_on"]),
        "registry_path": registry.get("_registry_path"),
        "files": records,
    }
    (destination / "snapshot_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    )
    return destination
