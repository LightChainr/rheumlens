#!/usr/bin/env python3
"""Build the final, manuscript-neutral RheumLens information package."""

from pathlib import Path
import csv
import hashlib
import importlib.metadata
import platform
import subprocess

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "evidence_package"


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    # A compact DE table for direct information retrieval.
    de = pd.read_csv(OUT / "tcell_de_pydeseq2/batch_age_gender_adjusted_all_genes.csv")
    de[["gene", "baseMean", "log2FoldChange", "pvalue", "padj"]].head(50).to_csv(
        OUT / "tcell_de_pydeseq2/batch_age_gender_adjusted_top50.csv", index=False)

    env = OUT / "environment"
    env.mkdir(exist_ok=True)
    rows = sorted((d.metadata["Name"], d.version) for d in importlib.metadata.distributions())
    with (env / "python_host_packages.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["package", "version"]); w.writerows(rows)
    (env / "system.txt").write_text(
        f"python={platform.python_version()}\nplatform={platform.platform()}\n", encoding="utf-8")
    for name in [".venv-stats", ".venv-de"]:
        pip = ROOT / name / "bin/pip"
        if pip.exists():
            text = subprocess.run([str(pip), "freeze"], check=True, capture_output=True, text=True).stdout
            (env / f"{name[1:]}_pip_freeze.txt").write_text(text, encoding="utf-8")

    report = r'''# RheumLens 最终信息包（不含成稿）

日期：2026-06-19。优先级：`results/fold_contained/` 是三队列主 AUC 与 paired-bootstrap CI 的唯一权威来源；本信息包中的新增结果均为补充审计或敏感性分析。

## 1. 已完成、可直接使用

### 主模型与统计稳健性

- GSE135779（44 donors）：scGPT 0.854（0.733–0.950）；PCA 0.948（0.871–1.000）；HVG-pseudobulk 0.898（0.755–0.994）。PCA−scGPT +0.094（+0.028–+0.182）。
- GSE285773（26 donors）：scGPT 0.950（0.850–1.000）；PCA 0.988（0.944–1.000）；HVG-pseudobulk 0.975（0.913–1.000）。
- GSE174188（261 donors）：scGPT 0.978（0.963–0.990）；PCA 0.981（0.967–0.993）；HVG-pseudobulk 0.984（0.970–0.995）。
- 三队列、三模型均完成 1,000 次 donor-label permutation；经验 P 均 ≤0.004。
- 30 个重复分层 CV seed：GSE135779 中 PCA 与 HVG 分别在 30/30 次高于 scGPT；GSE285773 分别为 23/30、25/30；GSE174188 中 PCA 无稳定优势，HVG 为 30/30。

### IFN 残差

- 固定 15-gene ISG panel 在三个队列中均完整存在；旧报告“后两队列不可用”已被原始矩阵核验推翻。
- GSE135779：scGPT 0.854→残差 0.760，衰减 0.094；主审计置换 P=0.01998（固定权威值）。
- GSE285773：0.950→0.844，衰减 0.106，P=0.00799。
- GSE174188：0.978→0.910，衰减 0.068，P=0.000999。
- 结论只能表述为 IFN 信号解释一部分、但不是全部判别能力；不能表述为 IFN-independent biology 已被因果证明。

### 组成、lineage 与 DE

- 原作者补充表 20-cluster donor composition 基线：AUC 0.915（0.821–0.983），P=0.000999；scGPT−composition −0.061（−0.174–+0.041）。旧 composition AUC 0.796 应弃用。
- 项目固定 marker 规则（不是原作者 cell labels）：T-lineage scGPT 0.826（0.675–0.948），P=0.00599；Monocyte 0.882（0.774–0.967），P=0.000999。旧 T=0.904、Monocyte=0.871 应弃用。
- 同一 marker-defined T-lineage 原始计数 donor pseudobulk，PyDESeq2：疾病单变量 1,188 个 FDR<0.05；Batch+Age+Gender 调整后 305 个。调整后最强信号包括 IFI27、ODF3B、IFIT5、HERC5、HELZ2；完整表与 top50 均已保存。该结果不能称为“原作者注释 T-cell DE”。

### 混杂敏感性（关键限制）

- GSE135779：可用技术/人口学协变量自身 AUC 0.857；协变量残差后的 scGPT/PCA/HVG AUC 为 0.551/0.631/0.766。
- GSE285773：仅 QC 协变量可用；其自身 AUC 0.781；残差后的 scGPT/PCA/HVG 为 0.819/0.869/0.775。
- GSE174188：QC、年龄、样本/文库/悬液数量、sex、ethnicity、processing cohort 自身 AUC 0.928；残差后的 scGPT/PCA/HVG 为 0.704/0.652/0.672。
- 因此三队列结果应称为“within-cohort disease-associated discrimination”，不能据此声称模型学习到无混杂的疾病机制或已实现临床泛化。
- GSE174188 managed/control-only：scGPT/PCA/HVG 0.979/0.975/0.983；single-sample-donor-only：0.960/0.958/0.977。processing-cohort 1 只有 HC，无法独立计算 AUC；cohort 2–4 的 scGPT 为 0.905/0.877/0.795。

### 跨年龄迁移

- GSE135779 pediatric→adult（成人仅 12 donors）：严格 source-scaler scGPT AUC 0.714（0.371–1.000）；source-fitted PCA 同为 0.714（0.371–1.000）。
- 不缩放 scGPT 为 0.829（0.486–1.000），显示结果对预处理敏感。旧 0.950 transfer 是硬编码/身份混淆，应弃用。

### 正则化与供体聚合敏感性

- 固定 folds 下扫描 logistic `C=0.01/0.1/1/10/100`。scGPT AUC 范围：GSE135779 0.846–0.909，GSE285773 0.931–0.969，GSE174188 0.965–0.978。小队列对正则化较敏感，因此不能根据 OOF AUC 事后挑选 C；主结果继续固定 C=1。
- 将 scGPT 细胞均值改为逐维精确中位数：GSE135779 AUC 0.821（0.680–0.939；median−mean −0.033，配对 CI −0.088–0.003）；GSE285773 0.950（0.850–1.000；差 0.000，−0.025–0.025）；GSE174188 0.980（0.965–0.991；差 +0.002，−0.002–0.005）。聚合结论总体稳定，但 GSE135779 的均值略优。

### coloc

- 使用 R `coloc` 5.2.3 标准 `coloc.abf` 重跑。默认 p12=1e-5：STAT4 0.063、IFIH1 0.006、PRDM1 0.051、TET2 0.009；均不支持强共定位。
- p12=1e-4 时 PP4 为 0.403/0.057/0.349/0.087，显示 prior sensitivity。
- BCL2、FCGR2A 是 0 个 shared variants，必须报告为“不可估计”，不能报告 PP4=0。
- 历史自定义 Python coloc 与标准实现不等价，旧“所有 PP4=0”弃用。

### 有效 Geneformer V2 重跑（GPU，matched sensitivity）

- 使用官方 Geneformer revision `04c2b2e…`、V2-104M、官方 rank-value tokenizer、原始 counts、CLS embedding、倒数第二层；固定抽取每供体 500 cells（44 donors，共 22,000 cells）。16,088/20,633 输入 gene symbols 映射到 Ensembl ID；token 序列长度中位数 944，95th percentile 1,782。
- Geneformer donor AUC 0.920（0.824–0.986，置换 P=0.000999）；同一 22,000 cells 的 matched scGPT AUC 0.873（0.758–0.961，P=0.000999）。Geneformer−scGPT +0.047（配对 CI −0.017–+0.127），未证明二者存在性能差异。
- 对项目 marker-lineage（非原作者标签）的 K-means ARI：Geneformer raw/standardized 0.255/0.302；scGPT 0.749/0.760。Geneformer 保留中等 lineage 结构但明显弱于 scGPT；其 768 个维度无常数维。
- 这是固定抽样 matched sensitivity，不能与使用全部细胞的三队列主 AUC直接混作同一 estimand；但它推翻了“Geneformer 完全退化/技术失败”的旧叙述。

## 2. 已获取的公共元数据

- GSE135779：GEO SOFT、论文补充工作簿、donor 临床表（Batch、Age、Gender、Race、Ethnicity、SLEDAI、治疗和实验室变量）、测序 QC、官方 cluster/subcluster donor counts 均已结构化。
- GSE285773：26 个独立供体的 GEO/BioSample 信息已归档；sorted CD4 T，10x 3′ v3，每样本目标上样 50,000 cells、NovaSeq 6000 单 lane、目标 50,000 reads/cell。公开记录未提供 age、sex、ancestry、treatment、site 或 batch。
- GSE174188：GEO 元数据及项目 H5AD/CELLxGENE 中可用 donor/sample 状态已结构化；重复样本 donor 映射已保存。

## 3. 明确弃用

- 全部旧 GPU/pre-fold HVG/PCA 主 AUC 与旧 IFN residual 结果。
- 历史 Geneformer 比较：旧脚本把表达矩阵列索引直接作为 token ID 输入 `BertForMaskedLM`，未使用 rank-value tokenizer/token dictionary；旧 AUC/ARI 从构造上无效。仅保留上述官方 V2 重跑。
- 旧 T-cell、Monocyte、composition AUC；旧 322 DE gene 数；旧 0.950 transfer；旧 custom-Python coloc PP4。
- unmatched PBMC→sorted-CD4 transfer 不作为验证结论。

## 4. 无法从现有公开来源可靠取得，按用户指示不再追逐

- GSE285773 的 age、sex、ancestry、治疗、采样时间、site/batch；公开 GEO、GSM、BioSample 均无这些字段。
- GSE135779 与 GSE285773 的 author-provided barcode-level cell labels；公开补充表只有 donor-level cluster counts。
- BCL2/FCGR2A 当前 coloc 输入的 shared variants，以及全部位点的完整 allele-level harmonization（历史 merged tables 缺 eQTL alleles）。
- 受控访问原始数据（dbGaP）与任何需要项目方账号授权的材料。
- 作者、单位、ORCID、funding、最终发布账号：由用户最后提供，本包不推断。

## 5. 文件入口

- 权威主结果：`results/fold_contained/all_cohorts_auc_bootstrap_ci.csv`
- 主模型置换：`results/evidence_package/all_primary_models_permutation_tests.csv`
- repeated CV：`results/evidence_package/repeated_cv/all_cohorts_repeated_cv_summary.csv`
- IFN：`results/evidence_package/ifn_residual_all_cohorts/all_cohorts_ifn_residual_summary.csv`
- 混杂：`results/evidence_package/confounding_complete/confounding_sensitivity_summary.csv`
- 组成/lineage：`results/evidence_package/secondary_foldcontained/`
- T-cell DE：`results/evidence_package/tcell_de_pydeseq2/`
- 跨年龄迁移：`results/evidence_package/cross_age_transfer/`
- coloc：`results/evidence_package/coloc_standard_rerun/`
- 公共元数据：`results/evidence_package/public_metadata/structured/`
- 环境：`results/evidence_package/environment/`
- 模型敏感性：`results/evidence_package/model_sensitivity/`
- 有效 Geneformer：`results/evidence_package/geneformer_valid/`
- 文件哈希清单：`results/evidence_package/FINAL_MANIFEST_SHA256.csv`
'''
    (OUT / "FINAL_INFORMATION_PACKAGE_20260619.md").write_text(report, encoding="utf-8")

    files = []
    for base in [ROOT / "results/fold_contained", OUT, ROOT / "analysis"]:
        for p in base.rglob("*"):
            if p.is_file() and p.name != "FINAL_MANIFEST_SHA256.csv":
                files.append(p)
    with (OUT / "FINAL_MANIFEST_SHA256.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["relative_path", "size_bytes", "sha256"])
        for p in sorted(set(files)):
            w.writerow([p.relative_to(ROOT), p.stat().st_size, sha256(p)])


if __name__ == "__main__":
    main()
