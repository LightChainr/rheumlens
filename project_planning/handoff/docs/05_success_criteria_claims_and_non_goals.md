# 成功标准、结论边界与非目标

## 1. 工程成功标准

- 三队列和 matched-500 数据可从 manifest 唯一重建；
- 所有方法实现统一 `fit/transform/predict` 或 end-to-end 接口；
- 训练和测试 donor 无信息泄漏；
- 所有 OOF predictions、fold、seed、配置和 hash 可追溯；
- 长任务支持断点续跑；
- 单元测试、smoke test 和结果 schema 验证通过；
- 失败模型和异常均生成结构化状态记录。

## 2. 统计成功标准

- 相同 estimand 使用相同 donor、fold 和评价流程；
- 主比较使用 outcome-stratified paired donor bootstrap；
- 关键模型进行 full-pipeline donor-label permutation；
- 新分布方法进行 cell-count matched sensitivity；
- 高阶增量主张需接受 mean-sufficient 或 moment-2 surrogate null；
- 复杂模型的调参预算公开且不超过预注册范围；
- cross-cohort 测试中所有拟合步骤只使用 source donor。

## 3. 科学成功标准

项目可因以下任一结果成功：

- 证明简单表达基线是最可靠方案；
- 证明 cell encoder 的差异小于 donor aggregator 的差异；
- 证明分布方法能恢复均值遗漏的信息；
- 证明 FOCUS 能稳定定位稀有机制证据；
- 证明复杂 MIL 在小 n 下无优势；
- 证明 within-cohort 排名不能预测 cross-cohort 排名；
- 建立一套可推广的方法选择规则。

不要求原创模型一定获得最高 AUC。

## 4. 可主张结论的层级

### 当前可主张

- 基础模型包含 SLE 相关 donor-level 信息；
- scGPT 未稳定优于 PCA/HVG；
- Geneformer 有效运行后表现具有竞争力；
- IFN、composition 和 covariates 均与判别相关；
- donor-level 和 fold-contained 评估是必要条件。

### 完成 Benchmark 后可主张

- 不同模型族在均值、分布、稀有状态和泛化上的系统差异；
- 聚合方式相对编码器的重要性；
- 新模型是否提供 mean 之外的可重复信息；
- 某些患者级免疫表型是否跨队列稳定。

### 需要未来数据才能主张

- 临床诊断价值；
- 对疾病活动、预后、复发或治疗反应的预测；
- 个体用药指导；
- 生物学因果机制；
- 前瞻性医疗决策收益。

## 5. 明确非目标

- 不从头训练通用大型语言模型；
- 不将百万细胞当成百万个监督样本；
- 不以单个队列的最高 AUC 选择最终方法；
- 不无限搜索更多模型变体；
- 不重新追逐已确认无法公开取得的临床字段；
- 不将 GSE137029 在未去重前作为独立验证；
- 不将皮肤和血液队列混为相同 estimand；
- 不因结果阴性而删除方法；
- 不把 UMAP、ARI 或 attention 权重直接当成机制证明。
