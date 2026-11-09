"""Local configuration, table readers and provenance helpers for TWAS drug discovery."""
import hashlib
import json
from pathlib import Path

import pandas as pd
import yaml


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def read_table(path):
    return pd.read_csv(path, sep="," if str(path).endswith((".csv", ".csv.gz")) else "\t")


def symbol(value):
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def load_settings(path, variables=None):
    path = Path(path).resolve()
    settings = yaml.safe_load(path.read_text())
    if not isinstance(settings, dict):
        raise ValueError("Drug-screen config must be a mapping")
    env = {"root": str(Path(__file__).resolve().parents[1]), "config_dir": str(path.parent), **(variables or {})}
    def expand(value):
        if isinstance(value, str):
            for k, v in env.items():
                value = value.replace("${" + k + "}", str(v))
            if "${" in value:
                raise ValueError("Unresolved config variable: " + value)
        elif isinstance(value, dict):
            return {k: expand(v) for k, v in value.items()}
        elif isinstance(value, list):
            return [expand(v) for v in value]
        return value
    return expand(settings), path


def input_path(value, base):
    p = Path(value).expanduser()
    p = p if p.is_absolute() else base / p
    p = p.resolve()
    if not p.is_file():
        raise ValueError(f"Missing input file: {p}")
    return p


def json_write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
