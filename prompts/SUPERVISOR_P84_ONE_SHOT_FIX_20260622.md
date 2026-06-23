# RheumLens P8.4 一次性修复与优化方案

日期：2026-06-22

用途：直接发给 Supervisor，再由 Supervisor 原样下发给 Worker。目标是停止远端临时补丁，回到已经验收的 RheumLens canonical pipeline。

## 1. 结论与禁止事项

当前 smoke 可以自然完成，但只能标记为 `ENGINEERING_SMOKE_LEGACY_DENSE`。

正式 P8.4 不得使用以下流程：

```text
selected cells × 61,497 genes
→ 全矩阵 .toarray()
→ cell-level PCA
```

原因不是单纯“太慢”，而是该流程改变了 `donor_expression_pca` 的方法定义。

权威实现位于：

```text
src/rheumlens/providers/expression.py
class DonorExpressionPCAProvider
```

其冻结定义是：

```text
sparse cell expression
→ 每名 donor 的 gene mean
→ outer-train donors 上选择2000个高方差基因
→ outer-train StandardScaler
→ outer-train PCA(25)
→ LogisticL2Estimator
```

正式 runner 必须调用现有 `build_method()` 和 `run_fixed_oof()`，不得复制或重写 PCA 数学过程。

## 2. 立即停止临时调试

当前 smoke 结束后：

1. 保存已有6项结果、日志、峰值RSS和代码hash；
2. 标记整个目录为 `LEGACY_DENSE_SMOKE_NOT_FORMAL`；
3. 禁止继续在 `/tmp/*.py` 上打补丁；
4. 新建一个正式分支目录，不覆盖旧输出：

```text
code/p8/p84_v4/
results/P8_4_cell_count_v4/
```

## 3. 唯一允许的数据流

### 3.1 输入只加载一次

加载：

- `lognorm.npz`：供 `donor_expression_pca`；
- `GSE174188_CD4_scgpt.npz`：供 `scgpt_mean` 和 `kme_multiscale@scgpt`；
- corrected sampling manifests；
- authoritative folds long table。

所有数据使用项目已有的：

```python
from rheumlens.data.io import load_npz_dataset
```

不得自行重新解析 sparse NPZ。

### 3.2 同一manifest跨模态选择完全相同的cell IDs

每个 corrected manifest 应至少包含：

```text
cell_id, donor_id, level, sampling_repeat
```

选中集合必须同时存在于 lognorm 和 scGPT：

```python
def subset_exact_cells(data, selected_cell_ids):
    selected = {str(x).strip() for x in selected_cell_ids}
    if len(selected) != len(selected_cell_ids):
        raise ValueError("sampling manifest contains duplicate cell IDs")

    available = {str(x).strip() for x in data.cell_ids}
    missing = selected - available
    if missing:
        raise ValueError(f"{len(missing)} selected cells absent from {data.name}")

    mask = np.fromiter(
        (str(x).strip() in selected for x in data.cell_ids),
        dtype=bool,
        count=len(data.cell_ids),
    )
    subset = data.subset_cells(mask)
    if subset.X.shape[0] != len(selected):
        raise AssertionError("selected cell count mismatch")
    return subset
```

不要假设两个NPZ row order相同；按 `cell_id` 选择。

## 4. folds 的唯一正确解析方式

该CSV是 donor-fold assignment long table，不是细胞级文件：

- 每名 donor 5行；
- 每个fold中该donor为train或test；
- 每名 donor 恰好一次test；
- `split_id` 可恒为0；
- 真正fold编号在 `fold`。

使用以下实现：

```python
from __future__ import annotations

import pandas as pd
from rheumlens.types import OuterFold


def canonical_id(value) -> str:
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    value = str(value).strip()
    if not value:
        raise ValueError("empty donor ID")
    return value


def load_authoritative_folds(path, eligible_donors) -> list[OuterFold]:
    frame = pd.read_csv(path)
    required = {"split_id", "fold", "role", "donor_id"}
    if not required.issubset(frame.columns):
        raise ValueError(f"missing fold columns: {required - set(frame.columns)}")

    frame = frame.copy()
    frame["donor_id"] = frame["donor_id"].map(canonical_id)
    frame["role"] = frame["role"].astype(str).str.strip().str.lower()
    frame["fold"] = frame["fold"].astype(int)
    eligible = {canonical_id(x) for x in eligible_donors}

    if set(frame.role.unique()) != {"train", "test"}:
        raise ValueError(f"unexpected role values: {sorted(frame.role.unique())}")

    counts = frame.groupby("donor_id").size()
    if not (counts == 5).all():
        raise ValueError("each donor must have exactly 5 fold rows")

    test_counts = frame.assign(is_test=frame.role.eq("test")).groupby("donor_id")["is_test"].sum()
    if not (test_counts == 1).all():
        raise ValueError("each donor must be test exactly once")

    source_donors = set(frame.donor_id)
    missing = eligible - source_donors
    if missing:
        raise ValueError(f"eligible donors absent from fold file: {sorted(missing)[:10]}")

    folds = []
    for fold_id in sorted(frame.fold.unique()):
        block = frame.loc[frame.fold.eq(fold_id)]
        train = tuple(sorted(set(block.loc[block.role.eq("train"), "donor_id"]) & eligible))
        test = tuple(sorted(set(block.loc[block.role.eq("test"), "donor_id"]) & eligible))
        if set(train) & set(test):
            raise ValueError(f"donor overlap in fold {fold_id}")
        if set(train) | set(test) != eligible:
            raise ValueError(f"fold {fold_id} does not cover eligible donor set")
        folds.append(OuterFold(split_id="0", fold=int(fold_id), train_donors=train, test_donors=test))

    test_union = [d for fold in folds for d in fold.test_donors]
    if len(test_union) != len(eligible) or set(test_union) != eligible:
        raise ValueError("OOF test assignment is not exactly-once")
    return folds
```

禁止 `drop_duplicates("donor_id")`，因为这通常保留第一条train行，不代表test fold。

禁止 `dict.get(donor, 0)`。

## 5. 正式 runner 必须复用项目引擎

核心运行函数如下。Worker只需要接入路径、manifest解析与checkpoint，不要重写方法逻辑。

```python
from __future__ import annotations

import numpy as np
import pandas as pd

from rheumlens.evaluation.engine import run_fixed_oof
from rheumlens.registry import build_method


PRIMARY_METHODS = (
    "scgpt_mean",
    "donor_expression_pca",
    "kme_multiscale@scgpt",
)


def run_one_method(
    method_id,
    *,
    lognorm_subset,
    scgpt_subset,
    folds,
    defaults,
    seed,
    estimand,
):
    registered = build_method(
        method_id,
        defaults=defaults,
        seed=int(seed),
    )
    if registered.learned:
        raise ValueError(f"P8.4 primary runner does not accept learned method {method_id}")

    datasets = {
        "lognorm": lognorm_subset,
        "scgpt": scgpt_subset,
    }
    if registered.data_key not in datasets:
        raise ValueError(f"unsupported data key {registered.data_key}")

    data = datasets[registered.data_key]
    oof, diagnostics = run_fixed_oof(
        data=data,
        folds=folds,
        method=registered.method,
        cohort="GSE174188_CD4",
        estimand=estimand,
        random_state=int(seed),
    )

    expected = {str(x) for fold in folds for x in fold.test_donors}
    observed = set(oof.loc[oof.status.eq("SUCCESS"), "donor_id"].astype(str))
    if observed != expected:
        raise ValueError(
            f"OOF coverage mismatch: missing={len(expected-observed)}, extra={len(observed-expected)}"
        )
    if len(oof) != len(expected):
        raise ValueError("OOF must contain exactly one row per donor")
    if oof.score.isna().any() or not np.isfinite(oof.score).all():
        raise ValueError("OOF contains non-finite scores")
    if oof.donor_id.duplicated().any():
        raise ValueError("duplicate donor OOF rows")
    return oof, diagnostics
```

这段代码会自动保证：

- `donor_expression_pca` 使用 `DonorExpressionPCAProvider`；
- 先donor mean，后HVG/scaler/PCA；
- 所有拟合限定在outer-train；
- `scgpt_mean` 使用冻结embedding均值；
- `kme_multiscale@scgpt` 使用现有KME与冻结参数；
- LogisticL2定义与P6/P8.3一致。

## 6. 为什么内存会立即下降

canonical `_donor_means()` 对每名 donor 的 sparse block调用：

```python
np.asarray(block.mean(axis=0)).ravel()
```

它只生成每名donor的一行均值，不会将全部cells转成dense。

最大 donor matrix 约：

```text
261 × 61,497 × 4 bytes ≈ 61 MB
```

选择2000 HVG后仅约：

```text
261 × 2,000 × 8 bytes ≈ 4 MB
```

因此正式 `donor_expression_pca` 的合理峰值通常应远低于当前43–46GB。若单方法仍超过约8GB，立即视为实现错误，不继续正式运行。

## 7. Primary 与 Secondary estimand

### 7.1 Primary common-support

固定500-cell level可用的237名donor：

```text
141 SLE / 96 controls
```

对所有 levels 都只使用这237人。这样不同cell-count level之间严格配对。

### 7.2 Secondary available-donor

允许各level使用最大可用集合：

```text
25: 261
50: 261
100: 260
200: 257
500: 237
```

必须标记为 `composition-changing secondary estimand`，不得直接将AUC变化解释成cell-count效应。

每个输出必须含：

```text
estimand_id
level
sampling_repeat
method_id
donor_id
y_true
score
fold
n_cells
```

## 8. 一次性测试矩阵

不得再靠“跑一下看看”调试。正式前完成以下自动测试。

### T1 folds解析

- 261 donors；
- 每donor 5行；
- 每donor 4 train + 1 test；
- folds大小53/52/52/52/52；
- OOF test union恰好261；
- 删除任意一行必须失败。

### T2 ID规范化

- `int`、`np.int64`、`str`、`np.str_`、`bytes`；
- 规范化后无碰撞；
- 空字符串、缺失ID、重复ID必须失败。

### T3 cell manifest跨模态对齐

- 同一manifest在lognorm和scGPT中选中完全相同cell IDs；
- 每donor cell数等于level；
- clinical label为162/99来源；
- y不参与sampling seed或排序。

### T4 sparse安全

对实际level=50运行 `donor_expression_pca`：

- 峰值RSS建议 <8GB；
- 不得出现对 `subset.X.toarray()` 的调用；
- OOF=eligible donors；
- NaN/Inf=0。

### T5 canonical一致性

使用P6 accepted config和未下采样/既有reference输入，确认：

```text
donor_expression_pca AUC ≈ 0.9856590597
```

容差先按项目既有复现规则；如果没有冻结规则，至少要求AUC绝对差 `<1e-8` 且OOF probability最大差 `<1e-6`，否则报告差异，不修改容差。

注意：不要要求新canonical runner匹配错误dense smoke的PCA概率；错误dense smoke不是reference。

### T6 resume

- 完成一个 `estimand × level × repeat × method` 后写临时文件；
- validator通过后原子rename；
- kill后重启只补缺失job；
- 已完成文件hash不改变。

## 9. checkpoint 与任务队列

一个job定义为：

```text
estimand × level × sampling_repeat × method
```

输出示例：

```text
oof/common_support/L0050/R000/scgpt_mean.parquet
oof/common_support/L0050/R000/donor_expression_pca.parquet
oof/common_support/L0050/R000/kme_multiscale_at_scgpt.parquet
```

每个job状态：

```text
PENDING → RUNNING → VALIDATED → COMPLETE
                    ↘ FAILED
```

只有validator通过才能从 `.partial` 原子改名为正式文件。

不要将多项方法写入同一个不可恢复的大文件。

## 10. 资源参数

正式v4先单进程跑完 `level=50, repeat=0, 3 methods`。

通过后：

- worker processes：2；
- `OMP_NUM_THREADS=1`；
- `MKL_NUM_THREADS=1`；
- `OPENBLAS_NUM_THREADS=1`；
- `NUMEXPR_NUM_THREADS=1`；
- 同一时间最多一个expression PCA job；
- 同一时间最多一个KME job；
- 数据先复制到本地SSD；
- 不启用48 workers；
- 不在共享盘反复解压NPZ。

即使服务器有1TB RAM，也不允许用内存掩盖错误算法。

## 11. Supervisor 的启动门

Worker必须先提交：

```text
P8_4_V4_IMPLEMENTATION_REPORT.md
P8_4_V4_TEST_REPORT.md
P8_4_V4_JOB_MATRIX.tsv
config.frozen.yaml
MANIFEST_SHA256.pre_run.tsv
```

验收条件：

- T1–T6全部通过；
- level50三方法完成；
- PCA峰值RSS降至合理范围；
- canonical reference复现；
- common-support与available-donor明确分离；
- 没有读取旧错误manifest或legacy dense结果。

满足后 Supervisor 可授权 `START_P8_4_FULL_V4`。

## 12. 给 Worker 的直接指令

```text
不要继续修补当前dense smoke。让它结束并归档为LEGACY_DENSE_SMOKE_NOT_FORMAL。

按P8.4 v4方案新建独立runner：
1. 使用load_npz_dataset加载数据；
2. 使用corrected manifest按cell_id跨模态选择；
3. 使用role train/test long table构建OuterFold；
4. 使用build_method和run_fixed_oof运行三个方法；
5. 禁止全cell矩阵toarray；
6. 完成T1–T6；
7. 只运行level50 repeat0三方法作为v4 validation；
8. 提交实现、测试、资源与hash报告；
9. 不自动启动full。
```
