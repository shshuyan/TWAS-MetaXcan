"""Orchestrate four independent TWAS × CMap Level 5 disease screens."""
from pathlib import Path
import json
import platform
import sys

import h5py
import numpy as np
import pandas as pd
import yaml

from .cmap_level5 import Level5Reference
from .continuous_screen import screen_disease
from .screen_io import input_path, load_settings, sha256


DISEASES = ["AD", "LSC", "atopic_dermatitis", "IPF"]


def _manifest_record(path, role, root=None):
    path = Path(path).resolve()
    shown = str(path.relative_to(root)) if root is not None else str(path)
    return {"role": role, "file": shown, "bytes": path.stat().st_size, "sha256": sha256(path)}


def _write_report(path, disease, summary):
    demo = summary["input_kind"] == "synthetic_test"
    lines = [
        f"# {disease}：TWAS × CMap Level 5 连续反转筛选", "",
        "**DEMO / synthetic：仅验证软件闭环，不代表真实候选药物。**" if demo else "本报告来自配置声明的真实输入，仍不构成临床建议。", "",
        f"- 状态：`{summary['status']}`",
        f"- 主组织：`{summary['tissue']}`；主细胞：{', '.join(summary['primary_cells'])}",
        f"- 共享基因：{summary['shared_genes']}；|TWAS Z| 截断：{summary['z_clip']}",
        f"- 置换次数：{summary['permutations']}；随机种子：{summary['seed']}",
        f"- 主分析假设数：{summary.get('primary_hypotheses', 0)}；候选数：{summary.get('candidate_count', 0)}", "",
        "主分析采用带权 cosine similarity；负值表示药物 MODZ 与疾病 TWAS Z 方向相反。",
        "每个疾病单独进行 BH-FDR，secondary cells 不能把未通过主分析的化合物升级为候选。",
        "TWAS 是遗传预测表达关联，不等于患者差异表达或因果治疗效应。", "",
    ]
    path.write_text("\n".join(lines))


def run_four_disease_screen(config_path, output=None):
    cfg, config_path = load_settings(config_path)
    if cfg.get("schema_version") != 1 or not isinstance(cfg.get("diseases"), list):
        raise ValueError("Four-disease config requires schema_version: 1 and diseases")
    by_disease = {item.get("disease"): item for item in cfg["diseases"]}
    if set(by_disease) != set(DISEASES) or len(cfg["diseases"]) != len(DISEASES):
        raise ValueError("diseases must contain exactly: " + ", ".join(DISEASES))
    analysis = cfg.get("analysis", {})
    required_analysis = {"z_clip", "minimum_genes", "permutations", "seed", "alpha"}
    if required_analysis - set(analysis):
        raise ValueError("analysis missing: " + ", ".join(sorted(required_analysis - set(analysis))))
    cmap = cfg.get("cmap", {})
    required_cmap = {"gctx", "sig_info", "gene_info", "pert_info", "cells"}
    if required_cmap - set(cmap):
        raise ValueError("cmap missing: " + ", ".join(sorted(required_cmap - set(cmap))))
    base = config_path.parent
    reference_paths = {key: input_path(cmap[key], base) for key in ["gctx", "sig_info", "gene_info", "pert_info"]}
    reference = Level5Reference(
        reference_paths["gctx"], reference_paths["sig_info"], reference_paths["gene_info"],
        reference_paths["pert_info"], cells=cmap["cells"],
    )
    output_value = output or cfg.get("output")
    if not output_value:
        raise ValueError("Set output in config or pass output")
    out = Path(output_value).expanduser()
    out = (out if out.is_absolute() else base / out).resolve()
    out.mkdir(parents=True, exist_ok=False)

    resolved = {**cfg, "output": str(out), "cmap": {**cmap, **{key: str(value) for key, value in reference_paths.items()}}}
    input_manifest = [_manifest_record(config_path, "batch_config")]
    input_manifest += [_manifest_record(path, "cmap_" + key) for key, path in reference_paths.items()]
    implementation_paths = [Path(__file__).resolve(), Path(screen_disease.__code__.co_filename).resolve(),
                            Path(Level5Reference.__init__.__code__.co_filename).resolve()]
    input_manifest += [_manifest_record(path, "implementation:" + path.name) for path in implementation_paths]
    cross_frames, summaries = [], []
    for disease in DISEASES:
        disease_cfg = {**analysis, **by_disease[disease], "input_kind": cfg.get("input_kind", "unspecified")}
        disease_cfg["twas"] = str(input_path(disease_cfg["twas"], base))
        disease_out = screen_disease(disease_cfg, reference, out / disease)
        summary = json.loads((disease_out / "summary.json").read_text())
        summaries.append(summary)
        _write_report(disease_out / "report_zh.md", disease, summary)
        (disease_out / "resolved_config.yaml").write_text(yaml.safe_dump(disease_cfg, sort_keys=False, allow_unicode=True))
        disease_inputs = [_manifest_record(config_path, "batch_config"),
                          _manifest_record(disease_cfg["twas"], "spredixcan_full")]
        disease_inputs += [_manifest_record(path, "cmap_" + key) for key, path in reference_paths.items()]
        disease_inputs += [_manifest_record(path, "implementation:" + path.name) for path in implementation_paths]
        (disease_out / "input_manifest.json").write_text(json.dumps(disease_inputs, indent=2, sort_keys=True) + "\n")
        disease_outputs = [_manifest_record(path, "output", disease_out) for path in sorted(disease_out.iterdir())
                           if path.is_file() and path.name != "output_manifest.json"]
        (disease_out / "output_manifest.json").write_text(json.dumps(disease_outputs, indent=2, sort_keys=True) + "\n")
        input_manifest.append(_manifest_record(disease_cfg["twas"], disease + "_spredixcan_full"))
        primary = pd.read_csv(disease_out / "all_primary_tests.tsv", sep="\t")
        if len(primary):
            cross_frames.append(primary)

    (out / "resolved_config.yaml").write_text(yaml.safe_dump(resolved, sort_keys=False, allow_unicode=True))
    (out / "input_manifest.json").write_text(json.dumps(input_manifest, indent=2, sort_keys=True) + "\n")
    runtime = {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
               "pandas": pd.__version__, "h5py": h5py.__version__, "pyyaml": yaml.__version__}
    (out / "runtime_versions.json").write_text(json.dumps(runtime, indent=2, sort_keys=True) + "\n")
    cross = pd.concat(cross_frames, ignore_index=True) if cross_frames else pd.DataFrame()
    cross.to_csv(out / "cross_disease_descriptive.tsv", sep="\t", index=False)
    summary = {
        "status": "complete", "input_kind": cfg.get("input_kind", "unspecified"), "diseases": DISEASES,
        "correction_families": len(DISEASES), "pooled_candidate_ranking": False,
        "candidate_counts": {item["disease"]: item.get("candidate_count", 0) for item in summaries},
        "limitations": [
            "TWAS association is not causal mediation or patient differential expression",
            "CMap cell models are proxies and IPF lacks an ideal primary lung fibroblast in GSE92742",
            "four disease FDR families and candidate lists are independent and must not be pooled post hoc",
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    demo = cfg.get("input_kind") == "synthetic_test"
    (out / "report_zh.md").write_text("\n".join([
        "# 四疾病 TWAS × CMap Level 5 筛选", "",
        "**DEMO / synthetic：仅验证软件闭环，不代表真实候选药物。**" if demo else "四个疾病独立分析；结果不构成临床建议。", "",
        "AD、LSC、atopic dermatitis 与 IPF 各自拥有独立候选表和全局 FDR family；本目录不生成跨疾病 pooled rank。", "",
    ]))
    output_manifest = [_manifest_record(path, "output", out) for path in sorted(out.rglob("*"))
                       if path.is_file() and path != out / "output_manifest.json"]
    (out / "output_manifest.json").write_text(json.dumps(output_manifest, indent=2, sort_keys=True) + "\n")
    return out
