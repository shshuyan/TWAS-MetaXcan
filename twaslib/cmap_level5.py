"""Bounded-memory reader for LINCS/CMap Level 5 GCTX matrices."""
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from .screen_io import read_table, symbol


def _decode(values):
    return [value.decode() if isinstance(value, bytes) else str(value) for value in np.asarray(values).ravel()]


class Level5Reference:
    """Index a Level 5 matrix and stream selected signatures in bounded blocks."""

    def __init__(self, gctx, sig_info, gene_info, pert_info, cells, compound_only=True):
        self.gctx = Path(gctx).resolve()
        self.source_paths = [self.gctx, Path(sig_info).resolve(), Path(gene_info).resolve(), Path(pert_info).resolve()]
        if self.gctx.suffix == ".gz":
            raise ValueError("GCTX must be decompressed before random-access screening; keep the .gz as the archived source")
        if not isinstance(cells, list) or not cells or len(cells) != len(set(cells)):
            raise ValueError("cells must be a unique nonempty list")
        with h5py.File(self.gctx) as handle:
            try:
                self.row_ids = _decode(handle["/0/META/ROW/id"][:])
                self.column_ids = _decode(handle["/0/META/COL/id"][:])
                shape = handle["/0/DATA/0/matrix"].shape
            except KeyError as exc:
                raise ValueError("Unsupported GCTX: require /0/DATA/0/matrix and row/column id metadata") from exc
        if shape == (len(self.row_ids), len(self.column_ids)):
            self.orientation = "genes_by_signatures"
        elif shape == (len(self.column_ids), len(self.row_ids)):
            self.orientation = "signatures_by_genes"
        else:
            raise ValueError(f"GCTX dimensions {shape} do not match metadata ({len(self.row_ids)}, {len(self.column_ids)})")

        genes = read_table(gene_info)
        gene_id = "pr_gene_id" if "pr_gene_id" in genes else "gene_id"
        gene_symbol = "pr_gene_symbol" if "pr_gene_symbol" in genes else "gene_symbol"
        if {gene_id, gene_symbol} - set(genes):
            raise ValueError("gene_info requires gene ID and symbol columns")
        alias = dict(zip(genes[gene_id].astype(str), genes[gene_symbol].map(symbol)))
        mapped = [alias.get(gene_id_value, "") for gene_id_value in self.row_ids]
        valid = [g for g in mapped if g]
        if len(valid) != len(set(valid)):
            raise ValueError("CMap gene_info contains a duplicate gene symbol")
        self.genes = mapped
        self.gene_index = {gene: i for i, gene in enumerate(mapped) if gene}

        sig = read_table(sig_info).copy()
        required = {"sig_id", "pert_id", "pert_iname", "pert_type", "cell_id"}
        if required - set(sig):
            raise ValueError("sig_info missing columns: " + ", ".join(sorted(required - set(sig))))
        if sig.sig_id.duplicated().any():
            raise ValueError("sig_info contains duplicate sig_id")
        sig = sig[sig.cell_id.isin(cells)]
        if compound_only:
            sig = sig[sig.pert_type == "trt_cp"]
        by_id = sig.set_index("sig_id", drop=False)
        keep = [sid for sid in self.column_ids if sid in by_id.index]
        if not keep:
            raise ValueError("Configured cells have no compound signatures in GCTX")
        self.signature_index = {sid: i for i, sid in enumerate(self.column_ids)}
        self.signatures = by_id.loc[keep].reset_index(drop=True)

        pert = read_table(pert_info)
        if {"pert_id", "pert_iname"}.issubset(pert):
            names = dict(zip(pert.pert_id.astype(str), pert.pert_iname.astype(str)))
            self.signatures["pert_iname"] = self.signatures.apply(
                lambda row: names.get(str(row.pert_id), row.pert_iname), axis=1
            )

    def iter_signature_blocks(self, gene_symbols, block_size=256):
        """Yield `(signature_metadata, signature × gene MODZ)` blocks."""
        if not isinstance(block_size, int) or block_size < 1:
            raise ValueError("block_size must be a positive integer")
        missing = [gene for gene in gene_symbols if gene not in self.gene_index]
        if missing:
            raise ValueError("Genes absent from CMap reference: " + ", ".join(missing[:5]))
        gene_ix = np.array([self.gene_index[g] for g in gene_symbols], dtype=int)
        signature_ix = np.array([self.signature_index[s] for s in self.signatures.sig_id], dtype=int)
        with h5py.File(self.gctx) as handle:
            matrix = handle["/0/DATA/0/matrix"]
            for start in range(0, len(signature_ix), block_size):
                chosen = signature_ix[start:start + block_size]
                low, high = int(chosen.min()), int(chosen.max()) + 1
                local = chosen - low
                if self.orientation == "genes_by_signatures":
                    dense = np.asarray(matrix[:, low:high], dtype=np.float32)[gene_ix][:, local].T
                else:
                    dense = np.asarray(matrix[low:high, :], dtype=np.float32)[local][:, gene_ix]
                if not np.isfinite(dense).all():
                    raise ValueError("Level 5 matrix contains nonfinite values in selected block")
                yield self.signatures.iloc[start:start + len(chosen)].copy(), dense
