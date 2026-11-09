"""Explain the TWAS-led screen without treating associations as treatment effects."""


def write_screen_report(out, summary, candidates, queries):
    primary = summary["primary"]
    lines = ["# AD：从 TWAS 出发的药物筛选", "",
             f"输入类型（运行配置声明）：`{summary['input_kind']}`。`synthetic_test` 是软件测试，不能解释为真实 AD 药物发现结果。", "",
             "GWAS → 标准化与对齐 → S-PrediXcan → S-MultiXcan → 有方向的 TWAS 签名 → LINCS/CMap → 候选药物。", "",
             "本步骤不需要预先提供候选药物或 GEO；筛选范围为输入参考库覆盖的全部配置细胞中的化合物。", "",
             f"- 主分析：{primary['tissue']} / {primary['source']} / {primary['cell']}，在运行前配置。",
             f"- 实际检验 {summary['tested_groups']} 组，共 {summary['total_groups']} 组；主分析候选 {len(candidates)} 个。",
             f"- {summary['permutations']} 次随机基因查询；所有来源×组织×药物×细胞统一 BH 校正。", "",
             "## TWAS 签名如何产生", "",
             "方向仅来自指定组织的 S-PrediXcan Z-score。正号表示遗传预测表达升高与 AD 风险正相关，负号表示负相关；不是患者实际表达上调/下调。",
             "S-MultiXcan 提供跨组织关联证据，不使用 z_mean 构造统一方向；单组织和联合结果的阈值见配置及完整注释表。",
             f"跨组织方向冲突基因 {len(summary['twas']['conflicting_genes'])} 个，按配置保留标记或排除，不平均冲突方向。", "",
             "| 组织 | 参考 | 正相关查询基因 | 负相关查询基因 | 状态 |",
             "|---|---|---:|---:|---|"]
    for q in queries:
        lines.append(f"| {q['tissue']} | {q['source']} | {len(q['risk_positive'])} | {len(q['risk_negative'])} | {q['status']} |")
    if summary["isolated_hit_resolution_limited"]:
        lines += ["", f"随机查询精度限制：最小经验 P 为 {summary['minimum_empirical_p']:.3g}。"
                  f"若仅一个药物条件组呈极端信号、其他组无信号，其 BH q 约为 {summary['isolated_hit_bh_resolution']:.3g}，"
                  "当前次数不足以使这样的孤立信号通过阈值。零候选不能排除该类信号；需预先增加随机次数后另开一轮分析。"]
    lines += ["", "## 主分析候选", "", "候选要求：主分析中分数为负且全筛选范围校正后 q 小于配置阈值。其他组织/细胞的结果同时保留，不代替主分析。", ""]
    if candidates:
        lines += ["| 药物 | 来源 ID | 分数 | 全局 q | 扰动条件数 |", "|---|---|---:|---:|---:|"]
        for r in candidates:
            lines.append(f"| {r['drug']} | {r['compound_id']} | {r['score']:.4g} | {r['q_bh_global']:.4g} | {r['signatures']} |")
    else:
        lines.append("没有候选通过主分析。检查查询状态区分证据不足与统计未通过；不会自动放宽阈值。")
    lines += ["", "## 方法与边界", "",
              "LINCS：风险正相关基因被药物压低、风险负相关基因被提高计为反向；使用上下调尾部重叠分数。",
              "旧版 CMap：完整基因排名的双向有符号 KS；非神经细胞仅作为补充参考。每个药物×细胞分别跨条件取中位数，保留全部条件记录。",
              "随机基因查询固定正负集合大小，背景为该组织全部合格 TWAS 受检基因与参考测量背景的交集；默认无每方向基因数截断。",
              "BH 覆盖本次全部有效检验；零分或同向药物不能因 p 值小而进入候选。没有参考覆盖与查询不足不能解释为药物无效。",
              "这是候选发现，不是独立药效验证；遗传表达关联反转不是因果治疗证据。后续需要共定位、靶点作用方向、独立神经细胞实验、暴露与安全性证据。",
              "随机基因零模型未校正 LD、共表达和模型覆盖偏差；本模块没有执行 LD 校正的遗传力富集。", "",
              "## 输出", "",
              "- `candidate_drugs.tsv/json`：主分析候选及其他组织/细胞结果。",
              "- `drug_screen.tsv/json`：全库药物×细胞×组织结果，包含未通过记录。",
              "- `spredixcan_all_annotated.tsv.gz`、`smultixcan_all_annotated.tsv.gz`：完整输入与注释；`twas_signed_genes.tsv`：指定组织的方向与选择标记。",
              "- `twas_queries.json`：每个组织/参考实际使用的正负查询基因和背景。",
              "- `signature_evidence.jsonl`：逐扰动条件分数；LINCS 另保存同向/反向基因与剂量时间。",
              "- `virtual_experiments.jsonl`：问题、假设、实验、结果、解释与后续问题。",
              "- `input_manifest.json`、`output_manifest.json`：输入/实现与输出哈希；`summary.json`：成功结束标记。", ""]
    (out / "report_zh.md").write_text("\n".join(lines))
