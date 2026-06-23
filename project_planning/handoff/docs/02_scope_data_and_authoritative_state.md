# 数据范围、当前权威状态与项目边界

## 1. 数据版图判断

RheumLens 使用 SLECA 2026 数据冻结作为公开 SLE 单细胞版图锚点。按“公开可得、血液/PBMC、donor-level SLE-vs-HC、样本量可用于模型评估”的纳入标准，当前主要数据已经基本见底。

主 Benchmark 队列固定为：

- `GSE135779`：儿童 PBMC，33 SLE / 11 HC；
- `GSE285773`：儿童 sorted CD4，16 SLE / 10 HC；
- `GSE174188`：成人 CD4 subset，项目纳入 261 donors。

辅助或审计数据：

- `GSE137029`：与 GSE174188 高度 overlap-sensitive，完成 donor 去重前不得称为独立验证；
- `GSE250024`：3 名 SLE donor 的疫苗纵向数据，仅作方法灵敏度或纵向辅助任务；
- 小样本、皮肤、空间或非 SLE 数据不进入主 Benchmark。

## 2. 当前权威结果摘要

### 三队列主要 OOF AUC

| Cohort | scGPT | Expression PCA | Donor-mean HVG |
|---|---:|---:|---:|
| GSE135779 | 0.854 | 0.948 | 0.898 |
| GSE285773 | 0.950 | 0.988 | 0.975 |
| GSE174188 | 0.978 | 0.981 | 0.984 |

主结论：frozen scGPT 捕获可重复的 within-cohort SLE 相关信号，但未显示稳定超过强表达基线。

### 有效 Geneformer matched sensitivity

在 GSE135779 中固定每 donor 500 cells：

- Geneformer V2-104M：AUC 0.920；
- matched scGPT：AUC 0.873；
- 差值 +0.047，配对 CI 跨零。

Geneformer 已被有效运行，但仅属于一个队列、固定抽样的 matched sensitivity，不能与三队列全细胞主 estimand 混合解释。

### 关键限制

- IFN 解释部分、但不是全部判别信号；
- 可用协变量本身具有较高判别能力；
- GSE174188 processing cohort 结构强；
- pediatric-to-adult transfer 样本小、区间宽、对预处理敏感；
- 当前结果不能声称无混杂机制或临床部署价值。

## 3. 权威与新 Benchmark 的关系

### Authoritative manuscript estimand

用于复现和修订现有论文，必须使用：

- 现有 donor inclusion；
- 现有 authoritative folds；
- 现有 fold-contained preprocessing；
- 已归档 OOF、bootstrap 和 permutation 资产。

### V2/V3 benchmark estimands

用于方法扩展，可包括：

- matched-500 cell sandbox；
- repeated-CV estimand；
- CD4 common-compartment cross-cohort estimand；
- fixed-cell-count distribution estimand；
- leave-one-cohort-out estimand。

任何新 estimand 必须单独命名，不得覆盖原论文权威主结果。
