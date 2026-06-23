# RheumLens 项目章程：统一 donor-level 单细胞 Benchmark 与原创聚合模型开发

版本：2026-06-20

## 1. 项目定位

RheumLens 研究单细胞 RNA-seq 基础模型和患者级建模方法在系统性红斑狼疮（SLE）中的真实价值。项目的统计单位始终是 donor，而不是 cell。

项目不以“证明某个模型一定优于其他模型”为目标，也不以追求阳性结果为前提。核心目标是：

> 建立一个严格、可复现、可审计的 donor-level 单细胞方法 Benchmark；系统比较传统表达方法、单细胞基础模型、分布方法、患者级集合模型和原创 donor 聚合方法；明确性能来自细胞编码、供体聚合、细胞组成、稀有状态、队列结构还是技术混杂。

RheumLens 的长期方法学命题是：

> 单细胞基础模型在临床表型预测中的瓶颈，可能不只是细胞编码器本身，而是如何将同一位患者的大量细胞压缩为稳定、可解释、可迁移的 donor representation。

## 2. 总体目标

### 2.1 复现与修订现有论文

1. 以最新权威结果为唯一基准，复现三队列主分析；
2. 彻底弃用历史错误 Geneformer、旧 AUC、旧 transfer、旧 coloc 和其他已否定结果；
3. 将有效 Geneformer V2 matched-cell sensitivity 纳入修订稿；
4. 统一术语：donor-mean HVG expression 与 raw-count pseudobulk 严格区分；
5. 保持 donor-level splitting、fold-contained preprocessing、paired bootstrap、permutation、repeated CV 和 confounding sensitivity。

### 2.2 建立统一 Benchmark

将所有方法拆成可交换的三部分：

```text
cell representation
      ↓
donor aggregation
      ↓
donor estimator
```

在相同 donor、相同 folds、相同细胞纳入和相同评估预算下，比较：

- 传统表达基线；
- scGPT 和 Geneformer；
- 均值、方差、分位数、KME、OT 和 prototype 方法；
- Deep Sets、Set Transformer、MIL 等患者级方法；
- RheumLens 原创方法。

Benchmark 不只输出排行榜，还要提炼方法学规律：

- 什么时候 mean/PCA 足够；
- 什么时候分布信息有增量；
- 什么时候稀有状态和 cell composition 重要；
- 什么时候复杂模型因小 donor 数而不稳定；
- within-cohort 与 cross-cohort 排名是否一致；
- foundation-model utility 是否主要受 donor pooling 限制。

### 2.3 开发原创 donor-level 方法

原创路线以“低容量、可解释、可审计、跨模型通用”为原则，重点包括：

1. **FOCUS-Lite / Language-compiled FOCUS**：用固定机制 query 在海量细胞中定位疾病相关证据，分别表示强度、丰度、尾部和异质性；
2. **ccKME-U**：细胞数校准、带不确定性的核均值 donor embedding；
3. **RED**：相对训练健康参考的状态偏离；
4. **GDS**：连续细胞状态图上的 donor signature；
5. **U-DER**：表示均值与估计不确定性共同进入患者级预测；
6. **IDE**：以跨 cohort 稳定性而非单队列最高 AUC 为目标的证据选择；
7. **DonorCLR**：同一 donor 不同细胞子样本间保持表示一致的无监督 donor encoder。

这些模型可使用 HVG、PCA、scGPT 或 Geneformer 作为细胞输入，从而区分“新编码器”与“新聚合算子”的贡献。

## 3. 核心论文主线

### 主线 A：基础模型审计

- frozen scGPT 能捕获可重复的 within-cohort SLE 相关信号；
- 但 PCA 和 donor-mean HVG 是强基线，且当前主分析未显示 scGPT 的稳定增量价值；
- 正确运行的 Geneformer V2 同样具备 donor-level 判别能力；
- cell-lineage recovery 和 donor-level disease discrimination 不是同一评价目标。

### 主线 B：供体聚合是关键瓶颈

- mean pooling 只保留平均状态；
- variance、quantile、KME、OT、prototype 和 FOCUS 测试均值之外的异质性、尾部、稀有状态和分布几何；
- 通过相同细胞、相同编码器和相同分类器的析因比较，识别增益来源。

### 主线 C：从分类到患者级免疫表型

最终输出不只是一项 AUC，而是患者级可解释表型：

- 全局 IFN 轴；
- 少量 IFN-high 细胞负担；
- 细胞类型组成；
- 细胞类型内状态偏移；
- 异质性和多峰；
- 与健康参考的偏离；
- 对技术因素、细胞数和队列变化的稳定性。

## 4. 项目原则

1. donor 是唯一独立监督样本；
2. 细胞数量不能被当成独立临床样本量；
3. 所有数据适应步骤必须在训练 donor 内拟合；
4. 外部验证队列标签不得参与特征或超参数选择；
5. 简单基线必须获得与复杂模型公平的计算和调参预算；
6. 所有结果均保留，包括阴性结果、失败模型和不稳定结果；
7. 预先区分 confirmatory、secondary 和 exploratory analyses；
8. 不因新方法未获胜而停止实现完整 Benchmark；
9. 不以高 AUC 自动解释为疾病机制、临床泛化或可部署诊断；
10. 所有环境、模型、输入、fold、OOF prediction 和结果均可追溯。

## 5. 预期成果

### 科学成果

- 一套 donor-level 单细胞方法分类体系；
- 一套 foundation-model 与传统方法公平比较的 Benchmark；
- 对 SLE 疾病信号结构的系统结论；
- 2–4 个具有清晰原创性的 donor aggregation 模型；
- 关于模型复杂度、聚合方式、细胞数和跨队列泛化的通用方法学规律。

### 工程成果

- 可复现代码库；
- 固定 folds 和统一 method registry；
- scGPT/Geneformer embedding cache；
- 标准结果 schema；
- OOF prediction、bootstrap、permutation 和 surrogate-null 资产；
- 可复用 CLI、配置和测试；
- 容器/环境 lockfile 和 SHA256 manifest。

### 论文成果

理想定位：

> A systematic donor-level benchmark reveals that aggregation, rather than cell encoding alone, governs the utility of single-cell foundation models for disease prediction.

并进一步提出：

> Calibrated distributional and mechanism-focused donor representations recover patient-level information discarded by conventional mean pooling.

## 6. 项目阶段

- **Phase 0**：资产、环境、数据范围和权威结果冻结；
- **Phase 1**：原稿复现与有效 Geneformer matched benchmark；
- **Phase 2**：经典、分布和患者级方法完整 Benchmark；
- **Phase 3**：原创模型开发、消融和 surrogate-null 审计；
- **Phase 4**：跨队列迁移、临床表型解释与论文整合；
- **Phase 5**：多疾病外部扩展（若后续数据和时间允许）。

