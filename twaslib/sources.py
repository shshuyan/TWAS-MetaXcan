"""Search public GWAS repositories and turn a hit into a trait config.

Sources:
  finngen       FinnGen R12 public summary statistics (hg38, Finnish/EUR)
  pan_ukbb      Pan-UK Biobank per-ancestry summary statistics (hg19)
  gwas_catalog  NHGRI-EBI GWAS Catalog studies with full summary statistics
                (harmonised hg38 file when available, otherwise the submitted file)

Catalog tables are downloaded once into raw_gwas/_catalogs and reused.
"""
import gzip
import io
import logging
import os
import re
import subprocess
import threading
import time
import urllib.request

import pandas as pd
import yaml

from .presets import PRESETS, autodetect

FINNGEN_MANIFEST = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv"
PANUKB_MANIFEST = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz"
GWASCAT_STUDIES = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-download-studies-v1.0.3.1.tsv"
GWASCAT_ANCESTRY = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-download-ancestries-v1.0.3.1.tsv"
PANUKB_POPS = ["EUR", "EAS", "CSA", "AFR", "AMR", "MID"]
ANCESTRY_MAP = {"European": "EUR", "East Asian": "EAS", "South Asian": "CSA", "African": "AFR",
                "African American or Afro-Caribbean": "AFR", "Hispanic or Latin American": "AMR",
                "Greater Middle Eastern": "MID", "Asian unspecified": "EAS"}

_lock = threading.Lock()
_cache = {}


def _get(url, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": "twas-pipeline/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def _get_range(url, nbytes=200_000):
    req = urllib.request.Request(url, headers={"User-Agent": "twas-pipeline/1.0", "Range": f"bytes=0-{nbytes}"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def _cached_file(catalog_dir, name, url, max_age_days=30, refresh=False):
    os.makedirs(catalog_dir, exist_ok=True)
    path = os.path.join(catalog_dir, name)
    stale = not os.path.exists(path) or (time.time() - os.path.getmtime(path)) > max_age_days * 86400
    if refresh or stale:
        logging.info("Downloading catalog %s", url)
        tmp = path + ".part"
        subprocess.run(["curl", "-fsSL", "--retry", "3", "-o", tmp, url], check=True)
        os.replace(tmp, path)
    return path


def _first_n_lines(url, n=2):
    raw = _get_range(url)
    if url.endswith((".gz", ".bgz")):
        d = gzip.GzipFile(fileobj=io.BytesIO(raw))
        text = b""
        try:
            while text.count(b"\n") < n:
                chunk = d.read(4096)
                if not chunk:
                    break
                text += chunk
        except (EOFError, OSError):
            pass
    else:
        text = raw
    return text.decode("utf-8", "replace").splitlines()[:n]


def _sniff_sep(line):
    return "\t" if "\t" in line else ("," if line.count(",") > line.count(" ") else r"\s+")


# --------------------------------------------------------------------------- catalogs

def load_catalogs(catalog_dir, refresh=False):
    with _lock:
        if _cache and not refresh:
            return _cache
        fg = pd.read_csv(_cached_file(catalog_dir, "finngen_R12_manifest.tsv", FINNGEN_MANIFEST, refresh=refresh),
                         sep="\t")
        pu = pd.read_csv(_cached_file(catalog_dir, "panukb_manifest.tsv.bgz", PANUKB_MANIFEST, refresh=refresh),
                         sep="\t", compression="gzip", low_memory=False)
        st = pd.read_csv(_cached_file(catalog_dir, "gwascat_studies.tsv", GWASCAT_STUDIES, refresh=refresh),
                         sep="\t", dtype=str, low_memory=False)
        st = st[st["FULL SUMMARY STATISTICS"].str.lower() == "yes"].copy()
        an = pd.read_csv(_cached_file(catalog_dir, "gwascat_ancestry.tsv", GWASCAT_ANCESTRY, refresh=refresh),
                         sep="\t", dtype=str, low_memory=False)
        an = an[an["STAGE"].str.lower() == "initial"]
        anc = (an.groupby("STUDY ACCESSION")["BROAD ANCESTRAL CATEGORY"]
                 .apply(lambda x: ",".join(sorted(set(", ".join(x.dropna()).split(", "))))))
        st["ancestry"] = st["STUDY ACCESSION"].map(anc).fillna("")
        st = st.drop_duplicates("STUDY ACCESSION")
        _cache.update(finngen=fg, pan_ukbb=pu, gwas_catalog=st)
        return _cache


def _anc_codes(text):
    codes = {ANCESTRY_MAP.get(a.strip(), "OTH") for a in re.split(r"[,;]", text) if a.strip()}
    return ",".join(sorted(codes)) if codes else ""


def search(catalog_dir, query, ancestry=None, sources=None, limit=200):
    """Search all sources. Returns a list of dicts with a common set of keys."""
    cat = load_catalogs(catalog_dir)
    sources = sources or ["finngen", "pan_ukbb", "gwas_catalog"]
    terms = [t for t in re.split(r"\s+", query.strip().lower()) if t]

    def match(series):
        m = pd.Series(True, index=series.index)
        s = series.fillna("").str.lower()
        for t in terms:
            m &= s.str.contains(re.escape(t))
        return m

    out = []
    if "finngen" in sources:
        fg = cat["finngen"]
        for r in fg[match(fg.phenotype + " " + fg.phenocode)].itertuples():
            out.append(dict(source="finngen", id=r.phenocode, trait=r.phenotype, ancestry="EUR",
                            n_cases=int(r.num_cases), n_controls=int(r.num_controls), build="hg38",
                            detail=r.category, pops=["EUR"]))
    if "pan_ukbb" in sources:
        pu = cat["pan_ukbb"]
        for r in pu[match(pu.description.astype(str) + " " + pu.phenocode.astype(str))].itertuples():
            pops = [p for p in str(r.pops).split(",") if p]
            if ancestry and ancestry not in pops:
                continue
            cases = {p: getattr(r, f"n_cases_{p}") for p in pops}
            pid = f"{r.trait_type}-{r.phenocode}" + (f"-{r.coding}" if pd.notna(r.coding) and r.coding != "" else "") \
                  + (f"-{r.modifier}" if pd.notna(r.modifier) and r.modifier != "" else "")
            out.append(dict(source="pan_ukbb", id=pid, trait=str(r.description), ancestry=",".join(pops),
                            n_cases=(int(cases.get(ancestry or "EUR")) if pd.notna(cases.get(ancestry or "EUR", float("nan"))) else None),
                            n_controls=None, build="hg19",
                            detail=f"{r.trait_type}; cases by ancestry: " + ", ".join(
                                f"{p} {int(c)}" for p, c in cases.items() if pd.notna(c)), pops=pops,
                            row=int(r.Index)))
    if "gwas_catalog" in sources:
        st = cat["gwas_catalog"]
        m = match(st["MAPPED_TRAIT"].fillna("") + " " + st["DISEASE/TRAIT"].fillna("") + " " + st["STUDY ACCESSION"])
        for r in st[m].itertuples():
            a = _anc_codes(r.ancestry)
            out.append(dict(source="gwas_catalog", id=r[st.columns.get_loc("STUDY ACCESSION") + 1],
                            trait=r[st.columns.get_loc("DISEASE/TRAIT") + 1], ancestry=a,
                            n_cases=None, n_controls=None, build="?",
                            detail=f"{r[st.columns.get_loc('FIRST AUTHOR') + 1]} {str(r.DATE)[:4]}, "
                                   f"{r[st.columns.get_loc('INITIAL SAMPLE SIZE') + 1]}",
                            pops=a.split(",") if a else []))
    if ancestry:
        out = [o for o in out if ancestry in o["ancestry"].split(",")]
    order = {"finngen": 0, "pan_ukbb": 1, "gwas_catalog": 2}
    # Single-ancestry matches first when filtering by ancestry, then by source and size
    out.sort(key=lambda o: (bool(ancestry) and o["ancestry"] != ancestry, order[o["source"]], -(o["n_cases"] or 0)))
    return out[:limit]


# --------------------------------------------------------------------------- resolve

def _safe_id(*parts):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", "_".join(str(p) for p in parts if p)).strip("_")[:80]


def resolve(catalog_dir, source, ident, pop=None):
    """Work out URL, format, build, column mapping and sample sizes for one dataset."""
    cat = load_catalogs(catalog_dir)
    if source == "finngen":
        r = cat["finngen"].set_index("phenocode").loc[ident]
        url = r.path_https
        header = _first_n_lines(url, 1)[0].split("\t")
        return dict(source=source, source_id=ident, trait_id=_safe_id(ident, "FinnGen_R12"),
                    description=f"{r.phenotype}, FinnGen R12 endpoint {ident}", ancestry="EUR",
                    url=url, format="finngen", build="hg38", separator="\t",
                    columns={k: v for k, v in PRESETS["finngen"].items() if v in header},
                    n_cases=int(r.num_cases), n_controls=int(r.num_controls),
                    sample_size=int(r.num_cases + r.num_controls), header=header, warnings=[],
                    source_label="FinnGen R12")
    if source == "pan_ukbb":
        pu = cat["pan_ukbb"]
        cand = [r for r in pu.itertuples() if _pan_id(r) == ident]
        if not cand:
            raise KeyError(ident)
        r = cand[0]
        pops = [p for p in str(r.pops).split(",") if p]
        pop = pop or ("EUR" if "EUR" in pops else pops[0])
        if pop not in pops:
            raise ValueError(f"{ident} has no {pop} analysis (available: {pops})")
        url = r.aws_path.replace("s3://pan-ukb-us-east-1/", "https://pan-ukb-us-east-1.s3.amazonaws.com/")
        header = _first_n_lines(url, 1)[0].split("\t")
        af = next((c for c in (f"af_controls_{pop}", f"af_{pop}") if c in header), None)
        cols = {"chrom": "chr", "pos": "pos", "ea": "alt", "nea": "ref", "eaf": af,
                "beta": f"beta_{pop}", "se": f"se_{pop}", "mlogp": f"neglog10_pval_{pop}"}
        if cols["mlogp"] not in header and f"pval_{pop}" in header:
            cols.pop("mlogp")
            cols["pval"] = f"pval_{pop}"  # some releases: natural-log p; check the warning
        warn = [f"column {v} not in file" for v in cols.values() if v and v not in header]
        nc = getattr(r, f"n_cases_{pop}")
        nctl = getattr(r, f"n_controls_{pop}")
        return dict(source=source, source_id=ident, trait_id=_safe_id(ident.replace("-both_sexes", ""), "PanUKBB", pop),
                    description=f"{r.description}, Pan-UKBB {ident}, {pop} analysis", ancestry=pop,
                    url=url, format="custom", build="hg19", separator="\t",
                    columns={k: v for k, v in cols.items() if v},
                    exclude_flag=f"low_confidence_{pop}" if f"low_confidence_{pop}" in header else None,
                    n_cases=int(nc) if pd.notna(nc) else None, n_controls=int(nctl) if pd.notna(nctl) else None,
                    sample_size=int(nc + nctl) if pd.notna(nc) and pd.notna(nctl) else None,
                    header=header, warnings=warn, source_label="Pan-UK Biobank")
    if source == "gwas_catalog":
        st = cat["gwas_catalog"].set_index("STUDY ACCESSION").loc[ident]
        base = str(st["SUMMARY STATS LOCATION"]).replace("http://", "https://").rstrip("/")
        files = _list_dir(base + "/")
        harm = [f for f in _list_dir(base + "/harmonised/") if f.endswith(".h.tsv.gz")] if "harmonised/" in files else []
        warn = []
        if harm:
            url, build = f"{base}/harmonised/{harm[0]}", "hg38"
        else:
            data = [f for f in files if re.search(r"\.(tsv|txt|csv)(\.gz|\.bgz)?$", f) and "meta" not in f
                    and not f.endswith((".log", ".tbi"))]
            if not data:
                raise ValueError(f"No summary-statistics file found in {base}")
            url = f"{base}/{data[0]}"
            build = "?"
            meta = [f for f in files if f.endswith("-meta.yaml")]
            if meta:
                try:
                    m = yaml.safe_load(_get(f"{base}/{meta[0]}"))
                    build = {"GRCh37": "hg19", "GRCh38": "hg38", "NCBI36": "hg18"}.get(str(m.get("genome_assembly")), "?")
                except Exception as e:  # noqa: BLE001
                    warn.append(f"could not read metadata: {e}")
            if build == "?":
                warn.append("genome build unknown: set it before running")
        lines = _first_n_lines(url, 2)
        sep = _sniff_sep(lines[0])
        header = re.split(sep, lines[0]) if sep == r"\s+" else lines[0].split(sep)
        cols, w = autodetect(header)
        warn += w
        size = str(st["INITIAL SAMPLE SIZE"])
        toint = lambda xs: [int(x.replace(",", "")) for x in xs if x.replace(",", "").isdigit()]
        nums = toint(re.findall(r"(\d[\d,]*) [^,]*?cases", size))
        ctl = toint(re.findall(r"(\d[\d,]*) [^,]*?controls", size))
        allnum = toint(re.findall(r"(\d[\d,]*)", size))
        anc = _anc_codes(str(st["ancestry"]))
        alist = anc.split(",") if anc else ["EUR"]
        primary = "EUR" if "EUR" in alist else alist[0]
        if len(alist) > 1:
            warn.append(f"multi-ancestry meta-analysis ({anc}); treated as {primary}. The LD reference "
                        f"panel is single-ancestry, so imputation/TWAS are approximate")
        return dict(source=source, source_id=ident, trait_id=_safe_id(ident, re.sub(r"\W+", "_", str(st["MAPPED_TRAIT"]))[:40]),
                    description=f"{st['DISEASE/TRAIT']} ({st['FIRST AUTHOR']} {str(st['DATE'])[:4]}, {ident})",
                    ancestry=primary, url=url, format="custom", build=build,
                    separator=sep, columns=cols, n_cases=sum(nums) or None, n_controls=sum(ctl) or None,
                    sample_size=(sum(nums) + sum(ctl)) or (sum(allnum) or None), header=header,
                    warnings=warn, source_label="GWAS Catalog")
    raise ValueError(f"unknown source {source}")


def _pan_id(r):
    pid = f"{r.trait_type}-{r.phenocode}"
    if pd.notna(r.coding) and r.coding != "":
        pid += f"-{r.coding}"
    if pd.notna(r.modifier) and r.modifier != "":
        pid += f"-{r.modifier}"
    return pid


def _list_dir(url):
    html = _get(url).decode("utf-8", "replace")
    return [h for h in re.findall(r'href="([^"?/][^"]*)"', html)]


# --------------------------------------------------------------------------- config writer

def write_trait_config(path, spec, raw_dir):
    """spec: output of resolve() after the user's edits, plus run options."""
    ext = ".gz" if spec["url"].endswith((".gz", ".bgz")) else ""
    fname = re.sub(r"[^A-Za-z0-9_.-]", "_", os.path.basename(spec["url"]).split("?")[0])
    cfg = {
        "id": spec["trait_id"],
        "description": spec.get("description", ""),
        "ancestry": spec.get("ancestry", "EUR"),
        "input": {
            "path": os.path.join("${raw_dir}", spec["trait_id"], fname if fname.endswith((".gz", ".bgz")) else fname + ext),
            "url": spec["url"],
            "format": spec.get("format", "custom"),
            "build": spec["build"],
            "separator": spec.get("separator", "\t"),
            "columns": spec.get("columns") or None,
            "sample_size": spec.get("sample_size"),
            "n_cases": spec.get("n_cases"),
            "n_controls": spec.get("n_controls"),
            "source": f"{spec.get('source_label', spec.get('source'))} {spec.get('source_id', '')}".strip(),
        },
        "model_types": spec.get("model_types", ["eqtl"]),
        "tissues": spec.get("tissues", "all"),
        "skip_imputation": bool(spec.get("skip_imputation", False)),
    }
    if spec.get("exclude_flag"):
        cfg["input"]["exclude_flag"] = spec["exclude_flag"]
    if spec.get("allow_ancestry_mismatch"):
        cfg["allow_ancestry_mismatch"] = True
    if cfg["input"]["format"] != "custom":
        cfg["input"].pop("columns")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(f"# Generated by the TWAS web app from {cfg['input']['source']}\n")
        yaml.safe_dump(cfg, fh, sort_keys=False, allow_unicode=True)
    return cfg
