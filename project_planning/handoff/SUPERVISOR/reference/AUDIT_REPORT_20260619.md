# RheumLens 论文审计报告

**日期**: 2026-06-19  
**权威主结果**: `/Users/lc/Documents/RheumLens/results/fold_contained/`  
**审计来源**: fold-contained statistical audit (donor 单位, training-fold-only HVG/PCA/scaling, 10,000 bootstrap)  

---

## A. 执行摘要

### ✅ 已确认

| 项目 | 状态 |
|------|------|
| 主结果 (0.854/0.950/0.978) 来自 fold-contained 审计 | CONFIRMED |
| 残差 AUC=0.760, P=0.020 (旧值 0.609/P=0.184 作废) | CONFIRMED |
| 15-gene ISG panel 已确认 | CONFIRMED |
| 手稿 v0.3 已使用所有修正值 | CONFIRMED |
| Legacy transfer 0.950 已被标记为 unsupported + 已删除 | CONFIRMED |
| GSE174188 cell_type 可用; GSE135779/GSE285773 无 | CONFIRMED |
| 3 个 folds 文件存在, GSE285773/GSE174188 由 StratifiedKFold(seed=42) 生成 | CONFIRMED |

### ❌ 被否定

| 项目 | 旧值 | 新值 | 来源 |
|------|------|------|------|
| GSE135779 scGPT AUC | 0.890 | **0.854** | fold_contained |
| GSE285773 scGPT AUC | 0.894 → | **0.950** (within-cohort, 非 transfer) | fold_contained |
| GSE174188 scGPT AUC | 0.852 | **0.978** | fold_contained |
| IFN 残差 AUC | 0.609 | **0.760** | fold_contained |
| IFN 残差 P | 0.184 (NS) | **0.020** (显著) | fold_contained |
| 旧结论 "IFN removal eliminates discrimination" | — | **作废** — 残差仍显著 | AUDIT_AND_RESULTS.md:54 |

### ⚠️ 仍缺失

| 项目 |
|------|
| GSE285773 的 batch/site/ancestry/age/sex/treatment 元数据 (archived H5AD 无这些字段) |
| 匹配 CD4→CD4 transfer (需要细胞类型注释, GSE135779 无 cell_type) |
| Secondary results (0.904/0.871/0.796/322 DE) 的 fold-contained 重跑 |
| 环境锁文件 (conda-lock / pip freeze) |
| 公共代码仓库 URL / DOI |
| 最终出版级图表 (从 corrected OOF predictions 重新生成) |
| 作者信息、通讯作者、资金来源 |
| 独立引用审计 |

### 🔴 可能改变结论

1. **IFN 残差结论已改变**: 旧手稿声称 IFN 移除消除了判别 (P=0.184, NS)。权威审计显示残差 AUC=0.760, P=0.020 — **非 IFN 信号在统计学上显著保留**。
2. **旧 transfer 0.950 为硬编码值**: 实际 PBMC→CD4 transfer AUC=0.63 (scaled) 或 0.94 (unscaled), 非 0.950。0.950 是 GSE285773 within-cohort OOF AUC。
3. **T-cell-only AUC=0.904、monocyte AUC=0.871、322 DE genes**: 这些值在手稿 v0.3 中, 但来自旧 GPU 服务器运行, **未在 fold-contained 审计中重跑**。标记为 DEPRECATED, 建议重跑或删除。

---

## B. audit_matrix.csv

已生成至 `/Users/lc/Documents/RheumLens/results/audit_matrix.csv`

---

## C. metadata_dictionary.csv

已生成至 `/Users/lc/Documents/RheumLens/results/metadata_dictionary.csv`

---

## D. secondary_results_provenance.csv

已生成至 `/Users/lc/Documents/RheumLens/results/secondary_results_provenance.csv`

---

## E. exact_text_hits.csv

已生成至 `/Users/lc/Documents/RheumLens/results/exact_text_hits.csv`
