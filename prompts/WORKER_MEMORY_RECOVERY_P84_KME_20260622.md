# RheumLens Worker 记忆恢复与 P8.4 KME 修复指令

日期：2026-06-22

你现在是 RheumLens Worker Agent。你刚刚丢失了部分上下文。本文是当前任务的恢复依据。不要恢复 full scheduler，直到本文的 PRE-FULL gates 全部通过并收到 Supervisor 新授权。

## 1. 当前阶段与范围

当前执行 P8.4 cell-count sensitivity。P8.3 repeated CV 已完成；P8.5–P10 尚未开始。

P8.4 方法：

```text
scgpt_mean
donor_expression_pca
kme_multiscale@scgpt
```

levels：25、50、100、200、500；sampling repeats：10。

两个 estimands：

- `common_support`：固定237 donors，141 SLE / 96 controls，primary；
- `available_donor`：各level最大可用donor，composition-changing secondary。

## 2. 已确认事故

### 2.1 旧 dense PCA 已排除

旧runner把整个cell×gene矩阵转为dense后做cell-level PCA。level200运行1小时以上、RSS约43–46GB，已终止并标记：

```text
LEGACY_DENSE_SMOKE_TERMINATED_NOT_FORMAL
```

权威 `donor_expression_pca` 必须使用：

```text
sparse cells×genes
→ donor mean
→ outer-train选择2000 HVG
→ outer-train StandardScaler
→ outer-train PCA(25)
→ LogisticL2Estimator
```

对应类：`rheumlens.providers.expression.DonorExpressionPCAProvider`。

修复后的level50验证：AUC=0.9766、OOF=261、约187秒、峰值RSS约3GB。

### 2.2 KME identity collapse 已确认

错误V4验证中，`scgpt_mean` 与 `kme_multiscale@scgpt` 的261名donor预测逐行完全相同，`abs_diff=0.0`。

根因：V4 runner调用了 `build_method(method_id)`，但没有使用返回的method对象。它对所有方法都手工执行 `X.mean(0) + StandardScaler + LogisticRegression`，所以KME退化成mean pooling。

所有受影响输出必须标记并隔离：

```text
PRE_FULL_UNVERIFIED_NOT_FORMAL
KME_IDENTITY_COLLAPSE_EXCLUDED
```

full scheduler已暂停；GPU空闲；不得自动resume。

## 3. 你缺失的精确 API 信息

### 3.1 build_method 返回结构

`build_method()` 返回：

```python
@dataclass(frozen=True)
class RegisteredMethod:
    method: FixedMethod | BagMethod
    data_key: str
    expression_key: str | None = None
    learned: bool = False
```

正确层级：

```text
registered.method
registered.method.provider
registered.method.aggregator
registered.method.estimator
registered.data_key
```

不是 `registered.aggregator`，也不是 `registered.method.method.aggregator`。

三个方法应解析为：

```text
scgpt_mean:
  data_key=scgpt
  provider=FrozenCellProvider
  aggregator=MeanAggregator
  estimator=LogisticL2Estimator

donor_expression_pca:
  data_key=lognorm
  provider=DonorExpressionPCAProvider
  aggregator=None
  estimator=LogisticL2Estimator

kme_multiscale@scgpt:
  data_key=scgpt
  provider=CellPCAProvider
  aggregator=MultiScaleRFFKMEAggregator
  estimator=LogisticL2Estimator
```

冻结参数：

```yaml
estimator_C: 1.0
n_hvg: 2000
donor_pca_dim: 25
cell_pca_dim: 32
cell_pca_max_fit_cells: 22000
kme_rff_dim: 256
kme_scales: [0.5, 1.0, 2.0]
kme_max_bandwidth_points: 4000
kme_max_diagnostic_points: 512
```

registry对KME设置 `include_linear_mean=True`。

### 3.2 KME 精确签名

```python
MultiScaleRFFKMEAggregator.fit(
    data: CellDataset,
    train_donors: list[str],
    context: FitContext,
) -> MultiScaleRFFKMEAggregator

MultiScaleRFFKMEAggregator.transform(
    data: CellDataset,
    donors: list[str],
    context: FitContext,
) -> DonorFeatures
```

不要在runner中手工调用它；由 `run_fixed_oof()` 管理。

### 3.3 run_fixed_oof 是正确入口

是。精确签名：

```python
run_fixed_oof(
    data: CellDataset,
    folds: list[OuterFold],
    method: FixedMethod,
    cohort: str,
    estimand: str,
    expression_data: CellDataset | None = None,
    donor_covariates: pandas.DataFrame | None = None,
    random_state: int = 0,
) -> tuple[pandas.DataFrame, list[dict]]
```

必须传 `method=registered.method`，不是 `registered`。

它会在每个fold中deepcopy并依次执行provider、aggregator、estimator，保证所有拟合限定在outer-train donors。

结论：用 `run_fixed_oof()` 完整替换V4 runner中的手工mean/scaler/logistic逻辑。

## 4. 首先恢复实际环境

使用计划运行任务的同一个Python环境：

```bash
which python
python -V
python - <<'PY'
import inspect, rheumlens
print(rheumlens.__file__)
print(inspect.getfile(rheumlens))
PY
python -m pip show rheumlens || true
conda info --envs || true
```

若有旧PID，检查其真实解释器：

```bash
readlink -f /proc/<PID>/exe
tr '\0' '\n' </proc/<PID>/environ | grep -E '^(PATH|CONDA|VIRTUAL_ENV|PYTHONPATH)='
```

若import失败，定位源码：

```bash
rg --files /autodl-fs/data /root/autodl-tmp 2>/dev/null \
  | rg '/src/rheumlens/(registry|engine|kme|expression)\.py$'
```

不要直接安装未知版本；优先恢复项目既有editable install。

运行并保存API审计：

```python
import inspect
from rheumlens.registry import build_method
from rheumlens.evaluation.engine import run_fixed_oof
from rheumlens.aggregators.kme import MultiScaleRFFKMEAggregator

print("run_fixed_oof", inspect.signature(run_fixed_oof))
print("KME.fit", inspect.signature(MultiScaleRFFKMEAggregator.fit))
print("KME.transform", inspect.signature(MultiScaleRFFKMEAggregator.transform))

defaults = {
    "estimator_C": 1.0,
    "n_hvg": 2000,
    "donor_pca_dim": 25,
    "cell_pca_dim": 32,
    "cell_pca_max_fit_cells": 22000,
    "kme_rff_dim": 256,
    "kme_scales": [0.5, 1.0, 2.0],
    "kme_max_bandwidth_points": 4000,
    "kme_max_diagnostic_points": 512,
}

for method_id in ("scgpt_mean", "donor_expression_pca", "kme_multiscale@scgpt"):
    reg = build_method(method_id, defaults=defaults, seed=20260619)
    print({
        "method_id": method_id,
        "method_type": type(reg.method).__name__,
        "provider": type(reg.method.provider).__name__,
        "aggregator": None if reg.method.aggregator is None else type(reg.method.aggregator).__name__,
        "estimator": type(reg.method.estimator).__name__,
        "data_key": reg.data_key,
        "learned": reg.learned,
    })
```

输出保存到 `results/P8_4_cell_count_v4/recovery/API_AND_ENVIRONMENT_AUDIT.txt`。

## 5. 恢复scheduler和manifests

只读定位：

```bash
ps -eo pid,ppid,stat,etimes,%cpu,%mem,rss,args | grep -E 'p84|P8_4|kme' | grep -v grep || true
rg --files /autodl-fs/data /root/autodl-tmp 2>/dev/null \
  | rg 'p84_full_staged\.py|P8_4|sampling|heartbeat|status\.|errors\.jsonl' \
  | sort
```

报告：

- scheduler path/hash；
- 活跃PID与auto-resume；
- job总数、final/partial/error数量；
- corrected manifests的path、columns、levels、repeats和hash；
- 哪些输出由错误V4 runner生成。

错误输出移入或标记到：

```text
results/P8_4_cell_count_v4/QUARANTINED_KME_IDENTITY_COLLAPSE/
```

不删除，不resume。

如果旧 `KME_IDENTITY_AUDIT.md` 为空，不要阻塞；用本轮实际证据重新生成。

## 6. 检查 sampling 是否 nested

同一donor、同一repeat必须满足：

```text
L25 ⊂ L50 ⊂ L100 ⊂ L200 ⊂ L500
```

如果现有50个manifests不是nested：

1. 原文件标记 `NON_NESTED_PRE_FULL_EXCLUDED`；
2. 每个 `donor × repeat` 只生成一个seed和一个cell permutation；
3. 各level取同一permutation前N个；
4. seed不得包含level；
5. 新输出写入 `sampling_nested_v1/`；
6. 保存nested validator与完整SHA256。

sampling必须label-blind；clinical label只能来自 `raw_counts.npz:y`；同一level/repeat的三个方法必须使用完全相同cell IDs。

## 7. 新建V5 runner，不再修补V4

新路径：

```text
results/P8_4_cell_count_v5/code/p84_v5_runner.py
```

旧V4代码与输出保持不变。

核心实现：

```python
from __future__ import annotations

import numpy as np
from rheumlens.evaluation.engine import run_fixed_oof
from rheumlens.registry import build_method


def run_one_method(method_id, *, datasets, folds, defaults, seed, estimand):
    registered = build_method(method_id, defaults=defaults, seed=int(seed))

    if registered.learned:
        raise ValueError(f"learned method not permitted: {method_id}")
    if registered.data_key not in datasets:
        raise ValueError(f"missing data_key={registered.data_key}")

    data = datasets[registered.data_key]
    expression_data = (
        datasets.get(registered.expression_key)
        if registered.expression_key is not None
        else None
    )

    oof, diagnostics = run_fixed_oof(
        data=data,
        folds=folds,
        method=registered.method,
        cohort="GSE174188_CD4",
        estimand=estimand,
        expression_data=expression_data,
        random_state=int(seed),
    )

    expected = [str(d) for fold in folds for d in fold.test_donors]
    if len(expected) != len(set(expected)):
        raise AssertionError("test donors are not exactly-once")
    if len(oof) != len(expected):
        raise AssertionError(f"OOF rows={len(oof)} expected={len(expected)}")
    if oof.donor_id.astype(str).duplicated().any():
        raise AssertionError("duplicate donor OOF")
    if set(oof.donor_id.astype(str)) != set(expected):
        raise AssertionError("OOF donor set mismatch")
    if not oof.status.eq("SUCCESS").all():
        raise RuntimeError(oof.loc[~oof.status.eq("SUCCESS")].to_string(index=False))
    if oof.score.isna().any() or not np.isfinite(oof.score).all():
        raise AssertionError("non-finite OOF score")
    return oof, diagnostics, registered
```

runner自身不得出现统一处理所有方法的：

```text
X.mean(axis=0)
StandardScaler
LogisticRegression
```

这些只能存在于项目已有provider/aggregator/estimator内部。

## 8. folds

fold CSV是donor-fold long table：261 donors × 5 folds；每donor四次train、一次test；fold sizes 53/52/52/52/52。

优先使用：

```python
from rheumlens.evaluation.splits import load_folds
folds = load_folds(folds_path)
```

再按estimand eligible donors限制每个fold，并断言train/test互斥、覆盖完整、每donor恰好test一次。

禁止 `drop_duplicates("donor_id")` 和 `dict.get(donor, 0)`。

## 9. 只允许一个 PRE-FULL smoke

仅运行：

```text
estimand=common_support
level=50
sampling_repeat=0
methods=3
```

Gate A：类型必须为：

```text
KME provider=CellPCAProvider
KME aggregator=MultiScaleRFFKMEAggregator
KME estimator=LogisticL2Estimator
```

Gate B：KME diagnostics至少包含：

```text
bandwidth
effective_rank
offdiag_q01_q50_q99
feature_norm_cellcount_corr
```

Gate C：按donor比较KME和mean：

```python
delta = np.abs(score_scgpt - score_kme)
print({
    "max_abs_diff": float(delta.max()),
    "mean_abs_diff": float(delta.mean()),
    "n_exact_equal": int((delta == 0).sum()),
    "correlation": float(np.corrcoef(score_scgpt, score_kme)[0, 1]),
})
```

若全部donor仍完全相同，立即BLOCK。不能只凭AUC不同验收。

Gate D：PCA runner不得对完整cell×gene调用toarray，RSS必须<8GB。

Gate E：OOF完整、唯一、无NaN/Inf、无fold leakage、labels正确。

Gate F：`.partial → validator → atomic rename`，kill/resume只补缺失job；V4错误输出不能被V5识别为完成项。

## 10. 最终统计定义

每个 `estimand × level × repeat × method` 单独计算AUC、PR-AUC和Brier。

不得把10个repeat拼成10倍donor数计算一个AUC。

同一repeat内计算：

```text
attenuation = AUC(level, repeat) - AUC(500, same repeat)
method difference = AUC(methodA, repeat) - AUC(methodB, repeat)
```

再对10 repeats汇总mean、SD、median、min、max、q025、q975。

common_support与available_donor分开。available_donor不能解释为纯cell-count效应。两个estimand的level500完全相同，只保留一个validated计算结果并引用，不能重复产生两个版本。

## 11. full执行资源规则

PRE-FULL通过后仍需Supervisor重新授权。

建议：workers=2起步，最多4；max PCA=1；max KME=1。

```bash
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
```

停止条件：PCA RSS>8GB、持续iowait>15%、D-state、SSH lag、label/fold错误、NaN/Inf、method/cache collision、KME再次identity collapse。

## 12. 现在的执行顺序

1. 恢复Python环境和rheumlens路径；
2. 保存API_AND_ENVIRONMENT_AUDIT；
3. 确认scheduler未运行且不会auto-resume；
4. 隔离错误V4输出；
5. 定位corrected manifests并验证nested；
6. 必要时生成全新nested manifests；
7. 新建V5 runner并调用 `run_fixed_oof(registered.method)`；
8. 完成单元测试和resume测试；
9. 只跑common_support level50 repeat0三方法；
10. 输出KME identity、diagnostics、资源、OOF和SHA256报告；
11. 汇报并等待Supervisor；
12. 不自行启动full。

## 13. 汇报格式

```markdown
[WORKER REPORT]

Stage: P8.4 memory recovery + KME repair
Status: COMPLETED_TECHNICAL / BLOCKED

Environment:
- python:
- rheumlens path:
- code SHA256:

Method identity:
- scgpt provider/aggregator/estimator:
- PCA provider/aggregator/estimator:
- KME provider/aggregator/estimator:

Scheduler:
- active PID:
- auto-resume:
- wrong outputs quarantined:

Sampling:
- manifest root:
- nested:
- levels/repeats:
- SHA256:

V5 smoke:
- scgpt metrics:
- PCA metrics:
- KME metrics:
- KME vs mean max_abs_diff:
- n_exact_equal:
- KME diagnostics:
- PCA peak RSS:

Validation:
- OOF coverage:
- donor uniqueness:
- fold leakage:
- labels:
- NaN/Inf:
- resume test:

Outputs:
- runner/config/reports/manifest:

Recommendation:
- START_FULL_RECOMMENDED: yes/no
- blockers:
```
