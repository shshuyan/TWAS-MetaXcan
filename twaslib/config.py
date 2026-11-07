"""Load the pipeline config and a trait config into one object."""
import os
import re
from dataclasses import dataclass, field

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

LIFTOVER_CHAINS = {"hg17": "hg17ToHg38.over.chain.gz",
                   "hg18": "hg18ToHg38.over.chain.gz",
                   "hg19": "hg19ToHg38.over.chain.gz"}
BUILD_ALIASES = {"grch37": "hg19", "b37": "hg19", "37": "hg19", "hg19": "hg19",
                 "grch38": "hg38", "b38": "hg38", "38": "hg38", "hg38": "hg38",
                 "hg18": "hg18", "ncbi36": "hg18", "hg17": "hg17"}


def _subst(obj, env):
    """Expand ${var} using env (and environment variables) recursively."""
    if isinstance(obj, str):
        def rep(m):
            k = m.group(1)
            return str(env.get(k, os.environ.get(k, m.group(0))))
        prev = None
        while prev != obj:
            prev, obj = obj, re.sub(r"\$\{(\w+)\}", rep, obj)
        return obj
    if isinstance(obj, dict):
        return {k: _subst(v, env) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_subst(v, env) for v in obj]
    return obj


def _deep_merge(a, b):
    out = dict(a)
    for k, v in (b or {}).items():
        out[k] = _deep_merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


@dataclass
class Config:
    trait_id: str
    trait: dict
    pipe: dict
    ref: dict
    work: str
    results: str
    extra: dict = field(default_factory=dict)

    # ---- convenience paths ----
    def p(self, *parts):
        return os.path.join(self.work, *parts)

    @property
    def build(self):
        b = str(self.trait["input"].get("build", "hg38")).lower()
        if b not in BUILD_ALIASES:
            raise SystemExit(f"Unknown genome build '{b}'")
        return BUILD_ALIASES[b]

    @property
    def raw_vcf(self):
        return self.p("00_gwasvcf", f"{self.trait_id}.{self.build}.vcf.gz")

    @property
    def export_tsv(self):
        return self.p("01_export", f"{self.trait_id}.{self.build}.tsv.gz")

    @property
    def harmonized(self):
        return self.p("02_harmonized", f"{self.trait_id}.txt.gz")

    @property
    def harmonized_vcf(self):
        return self.p("02_harmonized", f"{self.trait_id}.hg38.vcf.gz")

    @property
    def imputed(self):
        return self.p("04_imputed", f"imputed_{self.trait_id}.txt.gz")

    @property
    def twas_input(self):
        """GWAS used by S-PrediXcan: imputed unless imputation is turned off."""
        return self.harmonized if self.trait.get("skip_imputation") else self.imputed

    def spredixcan_dir(self, mtype):
        return self.p("05_spredixcan", mtype)

    def smultixcan_out(self, mtype):
        return self.p("06_smultixcan", mtype, f"{self.trait_id}_smultixcan.txt")


def load(trait_path, pipeline_path=None):
    pipeline_path = pipeline_path or os.path.join(ROOT, "config", "pipeline.yaml")
    with open(pipeline_path) as fh:
        pipe = yaml.safe_load(fh)
    with open(trait_path) as fh:
        trait = yaml.safe_load(fh)

    env = {"root": ROOT}
    env.update({k: v for k, v in pipe.get("paths", {}).items() if isinstance(v, str)})
    env = _subst(env, env)
    pipe = _subst(pipe, env)

    # Per-trait override of the reference resources (e.g. an EAS panel later)
    ref = _deep_merge(pipe["reference"], trait.get("reference"))
    env.update({k: v for k, v in ref.items() if isinstance(v, str)})
    ref = _subst(ref, env)
    trait = _subst(trait, _subst(env, env))

    tid = trait["id"]
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", tid):
        raise SystemExit(f"Trait id '{tid}' may only contain letters, digits, '_', '.', '-'")

    work = os.path.join(pipe["paths"]["work_dir"], tid)
    results = os.path.join(pipe["paths"]["results_dir"], tid)
    cfg = Config(tid, trait, pipe, ref, work, results)

    anc_t = str(trait.get("ancestry", "EUR")).upper()
    anc_r = str(ref.get("ancestry", "EUR")).upper()
    if anc_t != anc_r and not trait.get("allow_ancestry_mismatch"):
        raise SystemExit(
            f"[{tid}] GWAS ancestry {anc_t} != reference ancestry {anc_r}. LD and expression models "
            f"would not match. Give the trait a matching `reference:` block, or set "
            f"allow_ancestry_mismatch: true to run anyway.")
    return cfg
