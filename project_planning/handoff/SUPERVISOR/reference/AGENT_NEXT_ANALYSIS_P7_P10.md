# RheumLens 第三台服务器：下一步 P7–P10 自主分析指令

当前恢复状态已经确认：

- P4 Geneformer：`ACCEPTED`；
- P5 GSE135779 matched-500：`ACCEPTED`，27/27；
- P6 GSE285773：`ACCEPTED`，11/11；
- P6 GSE174188：`COMPLETED_TECHNICAL`，17/17，261 OOF rows/method，0 NaN；
- 独立 P7 目录不存在；P6 内执行过 6 个 originals，其中第六个是 `cckme_u_weighted`；
- P8–P10 尚未开始；
- 当前无 tmux、无活跃任务，A800 空闲。

不要重跑 P4、P5、GSE285773 或 GSE174188 的17方法固定 benchmark。接下来按本文自主推进，一次只运行一个重任务。

## 1. 先纠正阶段状态

“originals 在 P6 中运行过”不等于整个 P7 已完成。按原任务定义重新标记：

```text
P7A 可学习患者级模型：NOT_STARTED
P7B 原创方法：PARTIAL（仅以实际存在的6个 method_id 为准）
P8：NOT_STARTED
P9：NOT_STARTED
P10：NOT_STARTED
```

从实际 P6 summary/config 恢复 6 个 originals 的精确 method_id，不要再使用简称推测。将阶段映射写入：

```text
$RL_ROOT/results/STAGE_MAPPING_P6_P7.md
```

说明这些方法虽然计算产物位于 P6 目录，但在最终方法分类中属于 P7B；不得复制一份结果制造第二套 OOF。

## 2. Gate A：GSE174188 scGPT provenance 恢复

当前 scGPT NPZ 位于 `QUARANTINED_INTERRUPTED/`。禁止直接因为“结果看起来很好”就把它作为正式资产，也禁止立刻重提取。先做轻量 provenance 审计：

1. 找到 NPZ 的准确路径、大小、mtime、SHA256、keys、shape、dtype。
2. 检查 donor_id、cell_id/source-row、label、embedding 行之间的对齐信息。
3. 找到生成该 NPZ 的脚本、命令、日志、输入 H5AD、模型/checkpoint/vocab、参数和中断发生位置。
4. 确认 NPZ 是否在中断前已经完整关闭并可全量读取，而不是只写了一部分。
5. 核验：
   - 261 donors 与 authoritative folds 精确集合匹配；
   - 每 donor cell 数与源 H5AD/生成日志一致；
   - embedding 维度、NaN/Inf、常数维、重复行；
   - labels 在 donor 内一致；
   - NPZ 的 donor-level mean 能复现现有 `scgpt_mean` OOF 预测或指标，允许浮点误差但不允许重新训练口径漂移；
   - 现有 P6 config 是否直接引用了隔离路径。
6. 如有源 H5AD row identifiers，随机抽取固定种子样本核对 source-row mapping；不要全量重读两遍大文件。

判定：

- 全部通过：状态 `PROVENANCE_RECOVERED`。将 NPZ 以 `.partial` 复制/硬链接策略（按文件系统能力）放入正式 canonical 目录，校验 SHA256 后原子改名；原隔离文件保留，不删除。更新 config/manifest 的正式路径，但不得改变数据内容。
- 无法证明完整性或 row/donor 对齐：状态 `PROVENANCE_UNRESOLVED`。保留 P6 表达式方法验收，所有 GSE174188 scGPT 派生方法标为 provisional；停止依赖这些结果的 P7/P8/P9 分支，报告最小缺口。不要擅自重提取，除非现有明确脚本和输入足以确定性复现且磁盘预算允许。

输出：

```text
$RL_ROOT/results/GSE174188_SCGPT_PROVENANCE_AUDIT.md
$RL_ROOT/results/GSE174188_SCGPT_PROVENANCE_MANIFEST.json
```

## 3. Gate B：P6 GSE174188 科学验收

在 provenance 通过后，将 `COMPLETED_TECHNICAL` 升级为 `ACCEPTED` 前核验：

- 17 个唯一 method_id；
- 每方法 261 OOF donors，每 donor 恰好一次；
- train/test donor 互斥；
- 所有方法使用相同 authoritative folds；
- 标准化、HVG、PCA、特征选择和任何监督步骤仅在 outer-train 拟合；
- OOF AUC、PR-AUC、Brier、fold metrics 可由 OOF 表重算一致；
- labels、predictions、folds 无重复或缺失；
- config 不引用 `QUARANTINED`；
- REPORT、manifest、summary、OOF hashes 一致。

由于所有17个 AUC 均较高，执行一个小型泄漏 sanity gate：

- 对 `raw_pseudobulk`、`scgpt_mean`、表现最低的 `donorclr` 各做 100 次 donor-label permutation smoke；
- 必须重跑完整 fold-contained pipeline，不能只打乱最终 prediction；
- 固定 seed 并保留 null AUC；
- 这100次只用于 sanity，不得冒充 P8 的正式1,000次 permutation。

若 observed AUC 明显位于 null 尾部、无泄漏证据，则标记 P6 GSE174188 `ACCEPTED_PENDING_FORMAL_STATS`，并继续。若 null 异常偏高，停止并审计管线。

## 4. P6 三队列统一总结

建立一个不修改原始结果的统一层：

- GSE135779；
- GSE285773；
- GSE174188_CD4。

输出：

```text
P6_ALL_COHORTS_METHOD_SUMMARY.csv
P6_MODALITY_AVAILABILITY.csv
P6_FAILED_SKIPPED_REGISTRY.csv
P6_INPUT_AND_FOLD_MANIFEST.csv
P6_CONSOLIDATED_REPORT.md
```

规则：

- 不可用模态记录 `UNAVAILABLE`，不计为模型失败；
- amendment folds 与 authoritative folds 必须明显区分；
- 不把不同队列不存在的 method 当作零分；
- 保留 cohort、compartment、modality、method_family、fold_source、status 字段；
- 所有数字从现有 OOF 重算，不手工抄排行榜。

## 5. P7A：可学习患者级模型

从现有代码和冻结配置中寻找并运行预注册实现：

- Deep Sets；
- attention/top-k MIL；
- Set Transformer；
- MixMIL approximation。

执行规则：

1. 不得临时编写新的替代模型；若仓库没有实现，记录 `SKIPPED_NOT_IMPLEMENTED`。
2. 不得修改冻结模型定义或根据 test AUC 调参。
3. 超参数选择、early stopping、scaler 和 representation adaptation 只能在 outer-train donors 内完成。
4. 先选择最小且材料完整的队列做单 fold smoke；通过后再完成该队列 OOF。
5. 每次只跑一个模型、一个队列；初始 DataLoader workers ≤2。
6. 保存训练曲线、seed、outer/inner splits、checkpoint hash、OOF、失败原因和峰值资源。
7. 小样本模型若不稳定，按冻结 seed 报告，不得反复挑选“最好一次”。

完成每个实现真正支持的队列；不要求让所有模型强行适配所有模态。

## 6. P7B：原创方法完整性与缺失项

先将 P6 内已有的6个 originals 注册为 P7B 已完成产物，引用原 OOF 和 hash，不复制数据。

原始任务候选包括：

- ccKME-U；
- U-DER；
- FOCUS-Lite；
- Language-FOCUS；
- RED；
- GDS；
- DonorCLR；
- IDE（仅训练 cohorts 选择）。

对未出现在6个既有 method_id 中的方法：

1. 检查仓库是否有正式实现、冻结配置和所需材料；
2. 三者齐全才运行；
3. 缺实现记 `SKIPPED_NOT_IMPLEMENTED`；
4. 缺材料记 `SKIPPED_MISSING_ASSET`；
5. IDE 的任何选择只允许使用 training cohorts/outer-train，不能查看 target/test labels；
6. 不得临时发明 Language-FOCUS 输入或方法替代品。

每个方法生成 method card：定义、输入、训练范围、复杂度、失败模式、ablation 状态和结果路径。

## 7. P8：正式统计与稳健性审计

只有 Gate A、Gate B 和 P7 状态表完成后开始。不得把100次 sanity permutation 混入正式统计。

### P8.1 Paired stratified bootstrap

- 10,000 次；
- donor-level、outcome-stratified；
- 同一 cohort 内所有方法使用完全相同的 resampled donor indices；
- 从冻结 OOF predictions 计算，不重新训练；
- 输出各方法 AUC 95% CI 和预注册 paired method differences；
- 保存 bootstrap indices/seed 或可验证的紧凑表示。

这是低成本步骤，优先完成三队列所有通过验收的 OOF 方法。

### P8.2 Label permutation

- 正式规模 1,000 次；
- donor label permutation；
- 保持 fold 结构；
- 每次必须重跑完整 fold-contained pipeline；
- empirical P=`(1 + null_auc >= observed_auc 的次数)/(1 + N)`；
- 按冻结 primary-method 清单执行，不要未经定义地对全部方法扩张计算。

先从现有 config/task documents 恢复 primary-method 清单。若多个文档冲突，停止该小阶段并报告冲突，不自行挑选显著方法。

### P8.3–P8.8

按冻结定义依次执行：

1. repeated CV：30 seeds；
2. cell-count sensitivity：仅预注册下采样水平；
3. surrogate nulls：mean/moment2；
4. covariate sensitivity：只使用实际存在变量，缺失保持 NA；
5. kernel diagnostics：冻结网格；
6. random-query / focused-cell removal controls。

每个小阶段先做输出结构 smoke，再跑完整规模。不得根据中间显著性改变 seed、次数、方法或下采样水平。每个小阶段单独写 GO/NO-GO 和 manifest，可断点续跑。

## 8. P9：跨队列分析

仅在 P8 primary analyses 通过后执行：

1. 只选择生物学舱室可比、模态和特征定义兼容的 source/target；
2. 所有特征选择、缩放、模型拟合和阈值选择只使用 source；
3. target labels 只用于一次最终评价；
4. 分开报告 discrimination、calibration 和 failure modes；
5. 不把 within-cohort 高 AUC 表述为跨队列泛化；
6. age/sex/treatment/ancestry/site/batch 缺失保持 NA，不推断；
7. GSE285773 amendment folds 只用于 within-cohort，不把它错误用于 transfer fitting；
8. 输出 source-target 合法性矩阵，对不合法组合给出明确原因。

先完成 source-target eligibility table，再启动任何训练。

## 9. P10：汇总交付

最终生成：

- OOF master table；
- method summary；
- paired method differences；
- within-cohort/cross-cohort matrices；
- robustness dashboard；
- method-family figures；
- fold diagnostics；
- failed/skipped registry；
- paper-ready tables；
- code/data/model/environment manifests；
- `RESULTS_MANIFEST_SHA256.csv`；
- 最终执行报告。

旧 manuscript 或历史结果只能作为对照，不能覆盖本次冻结 OOF 和正式统计。任何与旧数字不一致之处必须列出来源、管线差异和权威结论。

## 10. 运行纪律与监控

- 使用 tmux；一次只运行一个重任务；
- 每60秒低开销记录 PID、D-state、PSI、RAM、GPU温度/显存、本地盘和共享盘；
- 环境、cache、tmp、checkpoint 优先放本地盘；最终小型结果再原子同步共享盘；
- 不同时下载、解压、hash 和 benchmark；
- 本地可用空间低于10GB、GPU>85°C、OOM、Xid、持续D-state或共享盘持续超时时优雅暂停；
- 日志暂时不增长不等于卡死，结合 PID state、CPU time、I/O 和输出 mtime 判断；
- 每个阶段记录命令、config hash、input hash、seed、PID、开始/结束时间、退出码和峰值资源。

## 11. 立即执行顺序

```text
恢复精确6个 original method_id 与阶段映射
→ Gate A scGPT provenance
→ Gate B P6科学验收 + 100次 permutation sanity
→ P6三队列统一总结
→ P7A冻结实现
→ P7B缺失项/skip registry
→ P8.1 bootstrap
→ P8.2 formal permutation
→ P8.3–P8.8
→ P9 eligibility + transfer
→ P10汇总
```

首次回复只需报告：当前执行步骤、主 PID、预计产物、是否发现 provenance/配置阻塞。没有阻塞时不要等待用户确认，继续自主推进。
