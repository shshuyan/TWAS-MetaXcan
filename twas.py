#!/usr/bin/env python
"""Config-driven TWAS pipeline.

  python twas.py run config/traits/AD_Bellenguez2022.yaml                # all steps
  python twas.py run config/traits/X.yaml --steps standardize export     # some steps
  python twas.py run config/traits/*.yaml                                # several traits in turn
  python twas.py compare config/traits/AD_Bellenguez2022.yaml config/traits/AD_Kunkle2019.yaml  # cross-dataset
  python twas.py presets                                                  # list input formats
"""
import argparse
import logging
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from twaslib import config, steps  # noqa: E402
from twaslib.presets import PRESETS  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pipeline", help="pipeline config (default config/pipeline.yaml)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run pipeline steps for one or more traits")
    r.add_argument("traits", nargs="+")
    r.add_argument("--steps", nargs="+", choices=steps.STEPS, default=steps.STEPS)
    r.add_argument("--from", dest="from_step", choices=steps.STEPS, help="run this step and all later ones")
    r.add_argument("--force", action="store_true", help="redo steps even if output exists")
    c = sub.add_parser("compare", help="compare S-MultiXcan results across traits")
    c.add_argument("traits", nargs="+")
    c.add_argument("--model_type", default="eqtl")
    c.add_argument("--out", default=os.path.join(config.ROOT, "results", "cross_trait"))
    sub.add_parser("presets", help="list input format presets")
    d = sub.add_parser("drug-screen", help="discover compounds from signed TWAS results across LINCS/CMap; no shortlist/GEO required")
    d.add_argument("config", help="TWAS drug-screen YAML; paths relative to this file")
    d.add_argument("--output", help="new output directory; must not already exist")
    a = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    if a.cmd == "drug-screen":
        from twaslib.drug_screen import run_drug_screen
        result = run_drug_screen(a.config, a.output)
        print(f"Completed TWAS drug screen: {result}")
        return
    if a.cmd == "presets":
        for k, v in PRESETS.items():
            print(f"{k}: {v}")
        print("gwasvcf: input is already GWAS-VCF (e.g. MRC-IEU OpenGWAS)")
        return
    cfgs = [config.load(t, a.pipeline) for t in a.traits]
    if a.cmd == "run":
        todo = steps.STEPS[steps.STEPS.index(a.from_step):] if a.from_step else a.steps
        for cfg in cfgs:
            steps.run_steps(cfg, todo, a.force)
    else:
        from twaslib.summarize import compare
        compare(cfgs, a.out, a.model_type)


if __name__ == "__main__":
    main()
