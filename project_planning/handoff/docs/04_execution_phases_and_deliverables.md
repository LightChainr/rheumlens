# 执行阶段、里程碑与交付物

## Phase 0：启动和资产冻结

目标：服务器具备可复现计算条件，旧项目权威资产完整迁移。

任务：

- 建立 `/autodl-fs/data/rheumlens` 目录；
- 创建 core、scGPT、Geneformer 和 R 环境；
- 复制不可变资产；
- 生成 SHA256；
- 核查 raw counts、donor labels、folds 和 cell IDs；
- 保存 GPU、CUDA、OS 和软件清单。

验收：

- 所有必须材料状态为 present/verified；
- authoritative result inputs 的 hash 与历史清单一致；
- 1000-cell smoke test 成功。

## Phase 1：原稿复现与 matched benchmark

任务：

- 复现三队列 scGPT/PCA/HVG 主结果；
- 复现 permutation、repeated CV、IFN、confounding 和 transfer；
- 复现有效 Geneformer matched-500；
- 补 matched PCA、HVG、ISG 和 pseudobulk；
- 补 Geneformer ISG/confounding/repeated-CV 审计。

交付：

- authoritative reproduction report；
- matched-500 method table；
- 所有 OOF predictions 和 paired differences；
- 结果差异审计日志。

## Phase 2：统一方法 Benchmark

### Family A：经典基线

- ISG-only；
- donor-mean HVG；
- donor-expression PCA；
- raw-count pseudobulk；
- cell-type pseudobulk；
- composition；
- scFeatures-style multi-view；
- scGPT/Geneformer mean。

### Family B：固定分布方法

- variance/MAD/skewness；
- quantile/tail fraction；
- KME/RFF/multi-scale KME；
- robust KME；
- Sliced-Wasserstein；
- Bures/低维高斯距离；
- prototype histogram 和 distance-to-reference。

### Family C：已有患者级模型

- Deep Sets；
- Set Transformer；
- attention MIL；
- MixMIL；
- ProtoCell4P；
- hierarchical MIL/scMILD 风格方法。

交付：

- 统一 method registry；
- 模型卡和运行日志；
- within-cohort 和 matched-cell 结果矩阵；
- 计算成本和失败记录。

## Phase 3：原创模型

优先顺序：

1. FOCUS-Lite；
2. ccKME-U；
3. RED；
4. U-DER；
5. GDS；
6. IDE；
7. DonorCLR；
8. Language-compiled FOCUS。

每个原创模型必须包含：

- 数学定义；
- 明确输入输出；
- fold-contained fit/transform；
- 主超参数固定或 inner CV；
- 简单嵌套基线；
- 消融；
- 随机或 surrogate null；
- 与 cell count、batch、ISG 和 covariates 的关系；
- 跨队列冻结测试。

## Phase 4：泛化和临床表型

任务：

- pediatric→adult transfer；
- CD4 common-compartment transfer；
- leave-one-cohort-out；
- calibration 和 decision-curve（仅在适当 estimand 下）；
- witness/FOCUS cell-state localization；
- composition vs within-type decomposition；
- 临床与生物学解释。

## Phase 5：论文与发布

交付物：

- 主文和补充材料；
- 统一代码仓库；
- 环境锁定和容器；
- OOF prediction/summary results；
- model cards；
- data provenance 和 hashes；
- figure source tables；
- 方法选择指南。
