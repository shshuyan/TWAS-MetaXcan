"""Plain-language report that keeps statistical evidence and gaps separate."""
from collections import Counter
import json


LABELS = {"supported": "支持本项统计标准", "not_supported": "未获支持", "not_covered": "参考库未覆盖",
          "not_testable": "背景不足，无法检验", "insufficient_query": "疾病查询不足",
          "both_signals_exploratory": "两项探索性信号同时满足", "both_association_signals": "两项关联信号同时满足",
          "reversal_only_genetics_not_supported": "仅反转获支持", "genetics_only_reversal_not_reproduced": "仅遗传富集获支持",
          "incomplete_reversal_evidence": "反转证据不完整", "incomplete_genetic_evidence": "遗传证据不完整",
          "neither_check_supported": "两项均未获支持"}


def write_report(out, candidates, genetic, reversal, metadata):
    m = metadata["methods"]; cohorts = [s["id"] for s in m["cohorts"]]
    def f(value):
        return "缺失" if value is None else f"{value:.4g}"
    headers = ["候选", "遗传富集q", "遗传结论", *[m['primary_cell'] + "反转 " + c for c in cohorts], "联合记录"]
    lines = ["# AD 候选药物：GWAS/TWAS 与 CMap/LINCS 复核", "",
             f"固定候选数：{len(candidates)}。本模块复核已有药筛名单，不新增或替换候选。", "",
             "**这些结果不构成因果靶点、独立药效验证或进入湿实验的充分条件。**", "",
             "## 数据和遗传证据层级", "",
             f"遗传输入模式：`{m['genetics']['mode']}`。{m['genetics']['mapping_strength']}。", ""]
    if m["genetics"]["mode"] == "gwas_mapped":
        lines += ["本次只使用指定GWAS研究的显著位点映射基因；没有运行全量GWAS→TWAS。富集是药物响应基因集合对位点映射基因的过度代表检验，不是LD校正的遗传力富集。没有遗传因果方向。", ""]
    else:
        lines += ["输入必须包含完整受检基因与P值；按全部受检基因计算Bonferroni阈值，再取各扰动参考库可测基因作为背景。S-MultiXcan联合P值只用于关联筛选，不用来指定药物反转方向。", ""]
    for r in m["genetic_backgrounds"]:
        lines += [f"- {r['source']}：遗传检验背景{r['background_genes']}基因，其中{r['genetic_hits_in_background']}个遗传命中。"]
    lines += [""]
    lines += ["疾病上/下调查询独立来自每个指定GEO队列的完整差异表；查询用于计算反转，遗传命中集合用于计算富集，两者不混为一个方向签名。", "",
              "## 固定候选复核", "", "| " + " | ".join(headers) + " |",
              "| " + " | ".join(["---"] * len(headers)) + " |"]
    for c in candidates:
        row = [c["drug"], f(c["genetic_q"]), LABELS.get(c["genetic_status"], c["genetic_status"])]
        for cohort in cohorts:
            prefix = "lincs_" + m["primary_cell"] + "_" + cohort
            row.append(f"{LABELS.get(c[prefix+'_status'], c[prefix+'_status'])}；分数={f(c[prefix+'_score'])}，q={f(c[prefix+'_q'])}")
        row.append(LABELS[c["decision"]]); lines.append("| " + " | ".join(row) + " |")
    lines += ["", "## 方法与解释范围", "",
              f"- 遗传富集：单侧Fisher；对本次全部可检验的固定候选×NEU/legacy-CMap基因集合统一BH。保留2×2计数、重合基因、期望重合、富集倍数和比值比。基因之间LD/共调控及扰动数据库的收录偏差没有被校正。",
              "- 药物响应基因不是药物直接结合靶点。`target_genetic_context.tsv`单独列出上游已提供靶点与遗传命中的交集，不能以响应基因富集替代直接靶点证据。",
              f"- LINCS：精确pert_id匹配；{m['primary_cell']}为主、其他配置细胞为补充，保留细胞、剂量、时间与签名ID。上调疾病基因被药物下调、下调疾病基因被药物上调计为反向；跨条件取中位数。",
              "- CMap：旧版build02，非神经细胞参考；完整排名做有符号KS，跨细胞取中位数。本地分数不是官方CLUE tau。只接受精确名称或配置中明确列出的别名；未核验结构一致性。",
              f"- 反转校准：每队列{m['permutations']}次随机基因查询，固定上/下集合大小，种子{m['seed']}。LINCS按候选×细胞组BH，CMap按全部参考化合物BH。负分数且q<{m['alpha']}才标为本项获支持。",
              "- 随机基因零模型不保留共表达结构，经验q值不是患者标签置换结果，也不是临床有效概率。",
              f"- LINCS是否已用于上游药筛：{m['upstream_reference_reused']}。复用时只属于一致性检查，不能算新增独立证据；即使跨GEO队列仍复用同一药物扰动库。供体跨队列独立性未核实。",
              "- 未测得、未匹配、查询不足与统计未获支持分别记录；均不直接等于药物无效。联合结果不合并P值，也不设药效总分。",
              "- 多候选/多条件的本次校正不能消除上游药筛使用相同LINCS数据造成的选择偏差；应在冻结方法后使用新参考/新实验验证。", "",
              "## CMap覆盖与每队列复核", "",
              "| 候选 | 队列 | 匹配名称 | 状态 | 分数 | q |", "| --- | --- | --- | --- | --- | --- |"]
    for r in reversal:
        if r["source"] == "cmap":
            lines.append(f"| {r['drug']} | {r['cohort']} | {r.get('matched_compound', '无')} | {LABELS.get(r['status'], r['status'])} | {f(r.get('score'))} | {f(r.get('q_bh'))} |")
    lines += ["", "## 输出和追溯", "",
              "`candidate_recheck.tsv/json`为最终证据矩阵；`genetic_enrichment.tsv/json`为集合检验；`reversal_by_cohort.tsv/json`为分队列反转；`signature_evidence.tsv/json`保留具体条件及同向/反向基因；`disease_queries.tsv/json`保存实际查询；`input_manifest.json`记录输入路径和SHA-256；`summary.json`最后写出，标志成功结束。", "",
              "`virtual_experiments.jsonl`把每项复核按Question→Hypothesis→Experiment→Result→Interpretation→Next Question保存；记录由统计结果和明确规则生成，不是LLM编造的科研推理。", "",
              "完整DE输入和全基因背景保留在清单中，仅在构建扰动查询时做预设FDR/大小筛选。缺少完整TWAS、共定位、结构身份、脑暴露和安全性等证据时，模块不会自动批准湿实验。", ""]
    (out / "report_zh.md").write_text("\n".join(lines))
    experiments = []
    for r in genetic:
        experiments.append({"candidate_id": r["candidate_id"], "drug": r["drug"], "source": r["source"],
                            "question": "该固定候选的药物响应基因是否富集AD遗传关联基因？",
                            "hypothesis": "响应集合与指定遗传命中集合的重合超过声明背景下的随机期望。",
                            "experiment": {"method": m["genetic_test"], "genetics": m["genetics"]},
                            "result": r, "interpretation": "exploratory_support" if r["support"] else r["status"] + "_without_support",
                            "next_question": "对重合位点进行共定位/精细定位，并区分直接靶点与下游响应基因；缺失覆盖时先补数据。"})
    for r in reversal:
        experiments.append({"candidate_id": r["candidate_id"], "drug": r["drug"], "source": r["source"], "cohort": r["cohort"],
                            "question": "该固定候选是否反向影响本队列的AD表达查询？",
                            "hypothesis": "跨参考条件的中位连接分数为负，且比匹配大小随机基因查询更负。",
                            "experiment": {"method": m[r["source"]], "permutations": m["permutations"], "seed": m["seed"], "alpha": m["alpha"]},
                            "result": r, "interpretation": r["status"],
                            "next_question": "用未参与筛选的神经扰动数据检验方向，并核实细胞状态、剂量与毒性；缺失覆盖时先补参考。"})
    (out / "virtual_experiments.jsonl").write_text("".join(json.dumps(r, ensure_ascii=False, allow_nan=False) + "\n" for r in experiments))
