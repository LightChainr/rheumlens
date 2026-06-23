# RheumLens 数据补齐与分析完成方案

日期：2026-06-19

原则：本地 `results/fold_contained/` 为统计结果权威来源。新信息报告用于补充公开身份和历史线索；与本地文件冲突时采用本地可复核证据。

## 0. 新信息报告中需要纠正的项目

- GSE135779 当前分析矩阵为 **20,633 genes**，不是 24,069。
- GSE285773 三个主模型的 1,000 次置换已经完成，不再是缺失项。
- 三个队列的 donor-feature matrices 均包含全部 15 个 ISG；GSE285773/GSE174188 的 `run_metadata` 空列表只表示尚未运行 IFN residual，不表示基因不可用。
- SSD 上的 `p3_transfer.py` 没有删除；第 60 行仍是硬编码 transfer AUC=0.950。
- coloc JSON 和脚本不在当前 workspace，但存在于 SSD 的 `rheumlens_sub_20260618_0152.tar.gz`，可以恢复。
- v0.3 已删除旧主 AUC 和旧 IFN 数值，但仍保留 unsupported transfer 0.950、cross-age 0.829、旧 composition/T/monocyte/DE 结果。
- HiSeq 4000 与 NovaSeq 6000 证明使用了不同测序仪器；两者仍是相关的 10x/droplet 工作流，不能把它直接写成广义“cross-platform generalization”。建议写成 **two sequencing instruments / independently generated cohorts**。
- 成人 GSE135779 归档对象实际为 **7 SLE/5 HD**；公开样本总表与最终项目纳入集应分别描述。

## A. 可以在 Mac 上完成：已有数据即可

### A1. 完成主模型统计闭环（最高优先级）

所需数据：三个 donor-feature NPZ、固定 folds、现有 OOF predictions。

方案：

1. 完成 GSE174188 的 scGPT、expression-PCA、donor-mean-HVG 1,000 次 donor-label permutation。
2. 汇总三个队列九个主模型的 observed AUC、exceedances、empirical P、null mean/SD。
3. 更新 SHA256 manifest 和 run metadata。

预期产物：统一 permutation table、三个 null 文件、哈希和审计说明。

### A2. Repeated-CV sensitivity

所需数据：三个 donor-feature NPZ。

方案：用 30–50 个预设随机种子完整重跑 5-fold donor CV；每次在训练 fold 内重选 HVG、重拟合 PCA/scaler/classifier。报告各模型 AUC 与 paired delta 的中位数、IQR、范围和种子级明细。

目的：量化固定 folds 未覆盖的不确定性，尤其是 26-donor GSE285773。

### A3. 三队列 IFN residual

所需数据：现有 donor expression/scGPT matrices；15 个 ISG 在三个队列均全部存在。

方案：对 GSE285773 和 GSE174188复用 GSE135779 的 training-fold-only ISG score residualization；生成 OOF predictions、10,000 paired bootstrap 和 1,000 permutations。另做 15-gene score 的病例/对照效应量及与 OOF score 的相关。

注意：这是项目内部 curated panel 的跨队列敏感性，不是原研究 IFN module 的严格复现。

### A4. GSE174188 混杂敏感性闭环

所需数据：官方 12.2 GB CELLxGENE H5AD、已生成 donor covariates、donor features。

方案：

1. 对 scGPT、PCA、donor-mean HVG 全部做 training-fold-only available-covariate residualization。
2. leave-processing-cohort-out；仅在测试 cohort 同时包含病例和对照时计算 AUC，否则报告不可识别。
3. processing-cohort 分层或匹配分析。
4. 排除 flare/treated donors；仅 managed SLE vs HC 重跑。
5. 对 11 个多 sample donors 做 first/predefined-sample-only、单 suspension 和排除多样本 donor 的敏感性。
6. 输出 disease×sex/ethnicity/processing cohort/age 列联及 PC1–PC5 关联的多重校正结果。

当前证据：technical-only AUC=0.926，必须在主文中明确队列结构可能解释相当部分判别。

### A5. GSE285773 已有 QC 混杂分析

所需数据：归档 H5AD。

方案：完善 cells/donor、UMI、genes/cell、mitochondrial fraction 的病例/对照比较；对 technical-only AUC=0.781 做 bootstrap CI/permutation；完成三主模型对这些 QC 协变量的 residual sensitivity。

局限：不能代替缺失的 batch/sex/age/treatment 审计。

### A6. 恢复并复核 coloc

所需数据：SSD 备份 `rheumlens_sub_20260618_0152.tar.gz`。

可恢复资产：`a1_coloc_results.json`、六个 merged TSV、Python/R 脚本、debug 脚本。

方案：

1. 解压到隔离的 provenance 目录并计算原文件 SHA256。
2. 核对 genome build、等位基因方向、palindromic SNP、MAF、case fraction 和样本量。
3. 使用原 R `coloc.abf` 重跑，并与 JSON 对比；区分 PP4=0、数值下溢和 too-few-shared。
4. 对 p12 等 priors 做敏感性分析。

已知：STAT4/IFIH1/PRDM1/TET2 有结果；BCL2/FCGR2A 是 `too_few_shared`，不能写成 PP4=0。

### A7. 恢复并复核 Geneformer extraction

所需数据：同一 SSD 备份中的 GF embedding、parquet、日志和诊断脚本，以及当前 scGPT H5AD。

方案：

1. 固定 matched cells 和 marker-derived labels。
2. 查清 checkpoint、版本、token dictionary、pooling、embedding shape/variance。
3. 对两种 embedding 使用相同 PCA/scaling/KMeans k、n_init 和多个 seeds。
4. 输出 ARI 分布、contingency 和 exact cell-index hash。

必要性：现有报告称 GF ARI=-0.007，而备份日志中至少一处为 0.0011，必须先统一定义和运行版本。

### A8. Cross-age transfer 重跑

所需数据：现有 GSE135779 pediatric/adult H5AD（成人归档为 7 SLE/5 HD）。

方案：从 donor matrices 重新生成 scGPT transfer predictions；明确是否 scaling；重跑 expression baseline，并输出 prediction、coefficient/scaler hash、bootstrap CI。旧 0.829 在无预测文件前不保留。

### A9. 环境、图表和发布资产

方案：

- 生成 `requirements.lock`/`pip freeze`、R `sessionInfo()`、输入输出 SHA256。
- 建立 figure provenance manifest，将 corrected、legacy、deprecated 分开。
- 所有保留图从最终 CSV/NPZ 重绘。
- 构建一键运行入口和 README；准备公开仓库包及 Zenodo 上传清单。
- 全文统一将当前 log-expression mean 改称 donor-mean/donor-aggregated expression。

## B. 可以在 Mac 上尝试：需要联网下载或恢复公开数据

### B1. GSE285773 完整公开元数据

目标数据：26 个 donor/sample 的年龄、性别、疾病状态、治疗、batch、site、FACS 信息、10x chemistry、library、lane、SRA run、processing date。

方案：

1. 下载 GEO series matrix、family SOFT/XML、全部 GSM 页面和 SRA RunInfo。
2. 下载并只读检查 `GSE285773_RAW.tar`（约 1.7 GB）中的文件名、barcodes 和可能的 metadata。
3. 查询 BioProject PRJNA1205970、BioSample 和 SRA 元数据并建立 donor↔sample↔run 映射。
4. 若 supplementary 中含 Seurat/metadata 表，提取作者 CD4 subtype 和 donor covariates。

若公开文件仍无这些字段，则转入“不能仅靠 Mac 完成”。

### B2. GSE135779 原作者 cell-type annotation

目标数据：cell barcode 到原作者 cell type/subtype 的映射，以及 donor/sample metadata。

方案：

1. 检查 GEO supplementary、论文 supplementary、公开 cell atlas/HCA/CELLxGENE 镜像和作者代码仓库。
2. 按 barcode/sample 对齐当前 H5AD；验证覆盖率和一对一性。
3. 若找到原作者标签，重跑 composition、T-cell、monocyte 和 matched CD4 transfer。
4. 若没有原作者标签，可在 Mac 用 marker scoring、CellTypist/Azimuth 等重新注释，但必须标为 project-derived annotation，并做 marker/一致性验证。

### B3. GSE174188 visit/treatment 映射

目标数据：多 sample donor 的 timepoint、flare/pre/post-treatment、sample accession 和 treatment。

方案：从 Perez 2022 supplementary tables、GEO GSM、CELLxGENE `sample_uuid/suspension_uuid/disease_state` 和公开作者元数据建立映射；优先采用论文定义的 baseline/flare/post-flare 规则。

### B4. Secondary analyses 重跑

依赖：B2 的可信 cell labels。

方案：

- composition：固定 broad-type harmonization，在相同 donor folds 中重跑并 bootstrap。
- T/monocyte scGPT：对同一 cell labels 子集 mean-pool，再做完整 fold-contained CV。
- donor-level DE：从 counts 层按 donor×cell type 求和，使用 edgeR/DESeq2/limma-voom；预先定义过滤、design、可用 covariates 和 FDR；输出完整 gene table，不沿用“322”常数。
- matched CD4 transfer：训练端与目标端使用同一 broad CD4 定义，输出 donor predictions 和 preprocessing provenance。

### B5. 文献、方法和 ISG 来源审计

方案：联网核查 GEO/PubMed/PMC、下载 supplementary methods；整理 Nehar-Belaid IFN modules M1.2/M3.4/M5.12 的原始来源。15-gene panel应描述为项目 curated panel；若找不到精确来源，不再追求虚构的单一文献出处。

## C. 不能保证仅靠 Mac 完成

### C1. 公开记录中不存在的临床/技术元数据

若 GEO/BioSample/SRA/supplementary 均未提供，则以下信息只能向数据提交者索取：

- GSE285773 的年龄、性别、ancestry、治疗、site、batch、chemistry、lane、processing date 和 sample/visit 对应关系。
- GSE135779 的 sex、ancestry、治疗、site/batch，以及受控 dbGaP 中的个体级临床信息。
- GSE174188 未公开的具体药物、精确治疗时间和 sequencing lane/date/site。

### C2. 受控访问原始数据

GSE135779 和 GSE174188 的 dbGaP 原始/个体级资料需要机构批准、数据使用协议和授权账号。Mac 可以下载已获批数据，但不能替代访问审批。

### C3. 权威 matched CD4 标签

若 GSE135779 没有公开原作者 cell labels，则 Mac 上的重新注释只能提供探索性 matched transfer，不能声称是原作者定义的 ground-truth CD4 subset。

### C4. 已丢失的历史决策证据

若无 git/log/聊天记录，无法证明：

- 15-gene panel 是否在看结果前预设；
- 旧“24-gene”完整名单及删除原因；
- GSE135779 folds 是否曾被用于调参或从多个 fold sets 中选择。

可采用的诚实处理：称 15-gene panel 为 project-defined；将固定-fold结果配合 repeated-CV；不声称完全独立预注册。

### C5. 投稿者身份和外部发布动作

作者、单位、通讯作者、ORCID、贡献、funding/conflict 需要用户提供。GitHub/Zenodo 最终发布和 DOI 生成需要相应账号授权；Mac 可以准备完整发布包，但不能替用户决定作者身份或法律声明。

## 建议执行顺序

1. A1 主置换闭环。
2. A4 GSE174188 混杂与多样本敏感性。
3. A2 repeated-CV。
4. A3 三队列 IFN residual。
5. B1/B2/B3 并行下载元数据和 cell labels。
6. A6/A7 恢复 coloc/Geneformer provenance。
7. B4 决定 secondary results 是重跑还是删除。
8. A8 cross-age transfer；unsupported cross-cohort transfer 从主文移除。
9. A9 全文、图、环境和发布资产统一冻结。
