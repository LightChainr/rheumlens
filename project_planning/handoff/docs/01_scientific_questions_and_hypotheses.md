# RheumLens 核心科学问题与假设

## 1. 主要科学问题

### Q1. Foundation model 是否真正增加 donor-level 信息？

比较 scGPT、Geneformer、PCA、donor-mean HVG 和 pseudobulk，判断基础模型 embedding 是否在强表达基线之外提供稳定增量。

### Q2. Cell encoder 与 donor aggregator 谁更重要？

在 HVG、PCA、scGPT 和 Geneformer 上应用相同聚合算子；同时在相同编码器上替换 mean、variance、quantile、KME、OT、prototype、MIL 和 FOCUS，进行析因比较。

### Q3. SLE 信号主要属于哪类结构？

区分：

- 全体细胞的均值移动；
- 细胞组成改变；
- 同一细胞类型内状态改变；
- 供体内异质性增加；
- 少量极端异常细胞；
- 多峰或连续状态轨迹改变；
- IFN 与非 IFN 相关信号；
- 技术和人口学队列结构。

### Q4. Within-cohort 高性能能否跨队列迁移？

比较固定队列内 OOF 性能、跨年龄 transfer、CD4 共同空间 transfer、leave-one-cohort-out 和 source-only fitting。

### Q5. 如何将模型输出转为可解释的患者级免疫表型？

使用 cell-type stratification、FOCUS query、reference deviation、witness function 和消融，定位哪些细胞状态驱动预测。

## 2. 预先声明的核心假设

### H1. 聚合瓶颈假设

单细胞基础模型包含疾病相关信息，但普通 mean pooling 丢失供体内分布结构，因此其 donor-level 增量价值受到限制。

### H2. 分布增量假设

当病例与对照具有相似均值、但在方差、尾部、稀有状态或多峰结构上不同，distribution-aware aggregation 将优于 mean pooling。

### H3. SLE 适配假设

SLE 同时存在强 IFN 均值轴、细胞组成变化和少量高激活状态，因此适合检验从一阶均值到高阶分布的完整信息阶梯。

### H4. 机制聚焦假设

固定、可审计的机制 query 能从海量细胞中识别对 donor phenotype 更有意义的局部证据，且优于随机 query 和无语义分布统计。

### H5. 复杂度边界假设

在 donor 数仅 26–261 的场景下，高容量 MIL/attention 模型可能不稳定；低容量、冻结编码器和显式统计聚合可能更可靠。

### H6. 泛化排序假设

队列内最高 AUC 的模型不一定具有最佳跨队列性能；跨队列稳定性、混杂敏感性和校准应与 AUC 同等重要。

## 3. 可证伪结果

下列结果都被视为有效科学结论：

- PCA/HVG 始终优于基础模型；
- mean pooling 已经接近充分统计量；
- 分布方法只在模拟中有效、真实队列无增量；
- FOCUS 与随机 query 相同；
- MIL 在小 n 下不稳定；
- Geneformer 与 scGPT 性能无显著差异；
- 方法队列内性能高，但跨队列迁移失败；
- 所有模型主要捕获 IFN 或技术结构。

项目成功不依赖某一原创模型获胜，而依赖结论是否经公平、完整、可复现的评估获得。
