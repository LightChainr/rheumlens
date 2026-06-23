# P8–P10 执行计划

## P8.2 正式 permutation 恢复

### 第一步：恢复持久状态

- 核验 `scgpt_mean` 1000 reps，不因第三台实例失联自动重跑；
- 隔离失败48-worker运行；
- 恢复准确 primary cohort×method 矩阵；
- 冻结 code/config/input/folds/seed hashes。

### 第二步：实现合格 runner

- 同一核心函数支持不同 method_id，不复制易漂移脚本；
- method目录隔离；
- 只加载方法需要的输入；
- deterministic rep queue；
- 本地rep checkpoint + 持久快照；
- errors.jsonl和严格最终断言；
- serial/parallel、kill/resume和错误config resume测试。

### 第三步：执行

1. `donor_expression_pca`；
2. `kme_multiscale@scgpt`；
3. 合并三 primary methods；
4. 输出 observed、null、empirical P、运行资源和SHA256。

## P8.3–P8.8

按冻结任务定义依次推进，不得直接跳P9：

| 阶段 | 任务 | 正式规模/原则 |
|---|---|---|
| P8.3 | repeated CV | 30 seeds；冻结定义 |
| P8.4 | cell-count sensitivity | 仅预注册下采样水平 |
| P8.5 | surrogate nulls | mean/moment2 |
| P8.6 | covariate sensitivity | 只用实际变量，缺失为NA |
| P8.7 | kernel diagnostics | 冻结网格 |
| P8.8 | query/removal controls | random/removal；防事后挑选 |

每项先做输出结构 smoke，随后正式规模。若不适用或缺材料，必须写 `SKIPPED_*` 和证据，不能静默略过。

## P9 跨队列

- 先生成 source-target eligibility matrix；
- 仅比较生物学舱室和模态兼容组合；
- 所有选择、缩放、拟合和阈值只使用source；
- target label仅用于最终评价；
- discrimination、calibration和failure mode分开报告；
- within-cohort高AUC不得表述为临床泛化。

## P10 交付

- OOF master table；
- method summary；
- paired differences；
- within/cross-cohort matrices；
- robustness dashboard；
- figures和paper tables；
- failed/skipped registry；
- code/data/model/environment manifest；
- `RESULTS_MANIFEST_SHA256.csv`；
- 最终执行报告与旧结果差异说明。
