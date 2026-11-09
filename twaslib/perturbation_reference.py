"""Validated legacy CMap rank matrix; compound/cell grouping is done by the screen."""
from collections import Counter

import h5py
import numpy as np
from scipy.stats import rankdata

from .screen_io import read_table, symbol


class CMapReference:
    def __init__(self, matrix, mapping):
        with h5py.File(matrix) as f:
            raw = np.asarray(f["assay"], dtype=np.float32)
            ids = [x.decode() if isinstance(x, bytes) else str(x) for x in f["rownames"][:].ravel()]
            self.names = [x.decode() if isinstance(x, bytes) else str(x) for x in f["colnames"][:].ravel()]
        if raw.shape != (len(self.names), len(ids)):
            if raw.shape == (len(ids), len(self.names)):
                raw = raw.T
            else:
                raise ValueError("CMap dimensions do not match row/column names")
        if not np.isfinite(raw).all():
            raise ValueError("Nonfinite CMap ranks")
        # Input contract: legacy build02, lower ranks denote higher expression.
        self.ranks = rankdata(raw, axis=1, method="average").astype(np.float32)
        d = read_table(mapping).fillna("")
        alias = dict(zip(d.cmap_gene_id.astype(str), d.symbol.map(symbol)))
        genes = [alias.get(g, "") for g in ids]
        counts = Counter(genes)
        self.gi = {g: i for i, g in enumerate(genes) if g and counts[g] == 1}
        self.genes = genes
        if any(len(name.split("__")) < 2 for name in self.names):
            raise ValueError("Malformed CMap compound/cell signature ID")
