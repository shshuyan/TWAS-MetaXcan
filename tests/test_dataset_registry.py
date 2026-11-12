"""Behavioral tests for frozen public-input snapshots."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import yaml

from twaslib.dataset_registry import download_snapshot, load_registry


def digest(data):
    return hashlib.sha256(data).hexdigest()


class DatasetRegistryTests(unittest.TestCase):
    def registry_file(self, root, datasets):
        path = root / "registry.yaml"
        path.write_text(yaml.safe_dump({
            "schema_version": 1,
            "frozen_on": "2026-11-12",
            "datasets": datasets,
        }))
        return path

    def test_invalid_registry_is_rejected_before_network_access(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = self.registry_file(root, [{"id": "missing-fields"}])
            with self.assertRaisesRegex(ValueError, "missing.*filename.*role.*url"):
                load_registry(path)

    def test_local_url_download_writes_verified_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = b"small public fixture\n"
            source = root / "source.tsv.gz"
            source.write_bytes(payload)
            path = self.registry_file(root, [{
                "id": "AD_PRIMARY",
                "disease": "AD",
                "role": "primary_gwas",
                "url": source.as_uri(),
                "filename": "AD/source.tsv.gz",
                "sha256": digest(payload),
                "bytes": len(payload),
                "build": "GRCh38",
                "ancestry": "European",
            }])

            registry = load_registry(path)
            snapshot = download_snapshot(registry, root / "snapshot")

            self.assertEqual((snapshot / "AD/source.tsv.gz").read_bytes(), payload)
            manifest = json.loads((snapshot / "snapshot_manifest.json").read_text())
            self.assertEqual(manifest["registry_frozen_on"], "2026-11-12")
            self.assertEqual(manifest["files"][0]["sha256"], digest(payload))
            self.assertEqual(manifest["files"][0]["bytes"], len(payload))
            self.assertEqual(manifest["files"][0]["status"], "downloaded")

    def test_completed_mismatch_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.bin"
            source.write_bytes(b"expected")
            path = self.registry_file(root, [{
                "id": "CMAP",
                "role": "perturbation_reference",
                "url": source.as_uri(),
                "filename": "cmap/data.bin",
                "sha256": digest(b"expected"),
                "bytes": 8,
            }])
            destination = root / "snapshot" / "cmap" / "data.bin"
            destination.parent.mkdir(parents=True)
            destination.write_bytes(b"user-data")

            with self.assertRaisesRegex(ValueError, "does not match.*refusing to overwrite"):
                download_snapshot(load_registry(path), root / "snapshot")
            self.assertEqual(destination.read_bytes(), b"user-data")

    def test_partial_local_download_resumes_without_duplicating_prefix(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = b"0123456789"
            source = root / "source.bin"
            source.write_bytes(payload)
            path = self.registry_file(root, [{
                "id": "RESUME",
                "role": "test_fixture",
                "url": source.as_uri(),
                "filename": "data/source.bin",
                "sha256": digest(payload),
                "bytes": len(payload),
            }])
            partial = root / "snapshot" / "data" / "source.bin.part"
            partial.parent.mkdir(parents=True)
            partial.write_bytes(payload[:4])

            download_snapshot(load_registry(path), root / "snapshot")

            self.assertEqual((root / "snapshot" / "data" / "source.bin").read_bytes(), payload)
            self.assertFalse(partial.exists())

    def test_cli_lists_frozen_entries_without_downloading(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.bin"
            source.write_bytes(b"x")
            path = self.registry_file(root, [{
                "id": "LSC_PRIMARY",
                "disease": "LSC",
                "role": "primary_gwas",
                "url": source.as_uri(),
                "filename": "LSC/source.bin",
            }])
            cli = Path(__file__).resolve().parents[1] / "twas.py"
            result = subprocess.run(
                [sys.executable, str(cli), "datasets", str(path), "--list"],
                capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("LSC_PRIMARY", result.stdout)
            self.assertIn("primary_gwas", result.stdout)
            self.assertFalse((root / "snapshot").exists())


if __name__ == "__main__":
    unittest.main()
