# RheumLens Supervisor 下一步计划：P8.4 与 Attention 增量支线

日期：2026-06-22

## 0. 当前决策

P8.3 不需要重跑。先完成归档签署，再按冻结阶段定义进入 P8.4。

本地权威执行计划中的阶段定义为：

- P8.3：repeated CV，30 seeds；
- P8.4：cell-count sensitivity，donor 内下采样；
- P8.5：surrogate nulls，mean/moment2；
- P8.6：covariate sensitivity；
- P8.7：kernel diagnostics；
- P8.8：random-query / focused-cell removal controls。

不要将 P8.4 改名为 ISG residualization。ISG residualization 是已有的独立稳健性证据模块，不应覆盖冻结阶段编号。

## 1. P8.3 最终签署

Worker 报告足以支持 `ACCEPT_RECOMMENDED`，但签署前补两项小修正：

1. `142/143` 必须明确唯一未纳入 manifest 的文件确实是持续变化的 `run.log`，并在 manifest policy 中解释；
2. completion report 出现两个截断/格式异常的 SHA256，重新输出完整的 `path<TAB>64位SHA256`，不得只给前缀。

若其余142项重新验证通过，Supervisor 可签署 P8.3 为 `ACCEPTED`。不要重跑计算。

## 2. 主线任务：P8.4 cell-count sensitivity

### 2.1 科学问题

判断 scGPT mean、expression PCA、KME 等 donor 表示的性能是否依赖每名 donor 可用的细胞数量，以及复杂分布方法是否只是从更多细胞中获益。

### 2.2 不得擅自确定下采样水平

Worker 必须先搜索 frozen config、stage ledger、旧 task spec 和代码中的预注册 cell-count levels。

输出：

```text
results/P8_4_cell_count_v1/P8_4_DEFINITION_RECOVERY.md
results/P8_4_cell_count_v1/P8_4_CONFIG_CANDIDATES.tsv
```

如果找到了唯一冻结定义，按其执行。如果不同文档冲突，列出路径、hash 和定义，暂停正式运行。若完全没有定义，只能生成设计提案，不能事后根据结果挑选水平。

### 2.3 执行原则

- donor 是抽样单位；下采样发生在 donor 内；
- 同一 `cohort × repeat × cell_count` 下所有方法使用相同 cell IDs；
- 每个 cell-count level 使用预先生成的 deterministic sampling manifest；
- sampling manifest 不依赖疾病标签；
- train/test donor folds 不变；
- PCA、scaler、kernel参数、原型和所有可学习步骤只能在 outer-train 内拟合；
- 对细胞数不足的 donor 按冻结规则处理，不得临时改为有放回采样；
- 每个 level 报告保留 donor 数和病例/对照数；
- 输出 AUC、PR-AUC、Brier、paired AUC difference 和相对 full-cell/matched reference 的衰减；
- 保存 donor-level OOF，不只保存汇总表。

### 2.4 推荐输出结构

```text
results/P8_4_cell_count_v1/
  config.frozen.yaml
  definition_recovery.md
  sampling/
    sample_manifest_<cohort>_<level>_<repeat>.parquet
  oof/
    <cohort>/<method>/<level>/<repeat>.parquet
  summaries/
    cell_count_method_summary.csv
    paired_attenuation.csv
    donor_retention.csv
  diagnostics/
    leakage_checks.json
    sampling_balance.csv
  logs/
  MANIFEST_SHA256.tsv
  COMPLETION_REPORT.md
  ACCEPTANCE_PACKET.md
```

### 2.5 最小采样代码

仅在恢复冻结 levels 后使用。`replace=False` 是强制默认值；若冻结定义另有规定，以冻结定义为准。

```python
from __future__ import annotations

import hashlib
import numpy as np


def stable_seed(base_seed: int, *parts: object) -> int:
    payload = "|".join(map(str, (base_seed, *parts))).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def sample_donor_cells(
    cell_ids: np.ndarray,
    n_cells: int,
    *,
    base_seed: int,
    cohort: str,
    donor_id: str,
    repeat: int,
) -> np.ndarray:
    cell_ids = np.asarray(cell_ids).astype(str)
    if len(cell_ids) < n_cells:
        raise ValueError(
            f"donor {donor_id} has {len(cell_ids)} cells, requires {n_cells}"
        )
    rng = np.random.default_rng(
        stable_seed(base_seed, cohort, donor_id, n_cells, repeat)
    )
    selected = rng.choice(len(cell_ids), size=n_cells, replace=False)
    return np.sort(cell_ids[selected])
```

不要使用 Python 内置 `hash()` 生成 seed，因为跨进程/会话不保证一致。

## 3. 5090 创新支线：P7A-N neural aggregation audit

该支线与 P8.4 分开注册：

```text
stage_id: P7A-N
status: EXPLORATORY_PREDEFINED
```

不得把它追加成新的 primary method，也不得依据测试 AUC 改架构。目标是回答：在冻结 scGPT cell embeddings 上，跨细胞 attention 是否稳定优于 mean pooling。

### 3.1 先审计现有实现

现有代码：

```text
src/rheumlens/bag_models/torch_models.py
```

已有：

- Deep Sets；
- gated attention MIL；
- top-k attention；
- Set Transformer；
- cells_per_bag=512；
- donor-level internal validation。

正式扩展前必须记录现有 P7A 结果、config、folds、seed 和 checkpoint。不得把旧结果覆盖掉。

现有实现有三个需要修复后再做重复验证的问题：

1. `predict_score()` 每个 donor 只随机抽取一次细胞，测试预测含较大 Monte Carlo 噪声；
2. 当前 Set Transformer 只有 attention residual，没有 LayerNorm 和 feed-forward block；
3. attention 权重没有导出，无法做稳定性和混杂审计。

### 3.2 冻结比较矩阵

建议只比较以下低容量模型：

| method_id | 作用 |
|---|---|
| `scgpt_mean` | 不可学习基准 |
| `deepsets@scgpt` | 非线性但无跨细胞 attention |
| `gated_attention_mil@scgpt` | 单池化注意力 |
| `set_transformer_isab@scgpt` | 低复杂度跨细胞 attention |
| `kme_multiscale@scgpt` | 非神经分布基准 |

top-k attention 只作为 ablation；MixMIL 保持原定义，不混入 attention 主比较。

主队列优先 GSE174188 CD4（261 donors）。26和44 donor队列仅作受限外部敏感性，不允许独立大规模调参。

### 3.3 固定探索性超参数

以下是新 P7A-N 的候选冻结配置，不是原 P7A 的追溯定义。Supervisor 必须在首次查看新结果前保存 config 和 SHA256：

```yaml
input: frozen_scgpt_cell_embeddings
hidden_dim: 64
n_heads: 4
n_inducing: 32
dropout: 0.10
cells_per_bag: 512
batch_size: 8
epochs_max: 150
patience: 15
learning_rate: 0.001
weight_decay: 0.0001
mc_test_passes: 8
outer_splits: reuse_P8_3
outer_repeats: 30
inner_validation: donor_stratified_train_only
selection_metric: inner_validation_log_loss
```

不做超参数网格。若显存不足，只能减 batch size，不能根据外层测试表现改 hidden/head/inducing points。

### 3.4 参考实现：Induced Set Attention MIL

该实现将复杂度从完整 self-attention 的 `O(N²)` 降为约 `O(NM)`，其中 `M` 是 inducing points 数量。

```python
from __future__ import annotations

import torch
from torch import nn


class MAB(nn.Module):
    def __init__(self, dim: int, heads: int, dropout: float = 0.1):
        super().__init__()
        self.attn = nn.MultiheadAttention(
            dim, heads, dropout=dropout, batch_first=True
        )
        self.norm1 = nn.LayerNorm(dim)
        self.ff = nn.Sequential(
            nn.Linear(dim, 2 * dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(2 * dim, dim),
        )
        self.norm2 = nn.LayerNorm(dim)

    def forward(self, q, kv, kv_padding_mask=None, need_weights=False):
        z, weights = self.attn(
            q, kv, kv,
            key_padding_mask=kv_padding_mask,
            need_weights=need_weights,
            average_attn_weights=False,
        )
        h = self.norm1(q + z)
        h = self.norm2(h + self.ff(h))
        return h, weights


class ISAB(nn.Module):
    def __init__(self, dim: int, heads: int, n_inducing: int, dropout: float):
        super().__init__()
        self.inducing = nn.Parameter(torch.randn(1, n_inducing, dim) * 0.02)
        self.to_inducing = MAB(dim, heads, dropout)
        self.to_cells = MAB(dim, heads, dropout)

    def forward(self, x, padding_mask):
        inducing = self.inducing.expand(x.shape[0], -1, -1)
        h, _ = self.to_inducing(
            inducing, x, kv_padding_mask=padding_mask
        )
        out, _ = self.to_cells(x, h, kv_padding_mask=None)
        return out


class PMA(nn.Module):
    def __init__(self, dim: int, heads: int, dropout: float):
        super().__init__()
        self.seed = nn.Parameter(torch.randn(1, 1, dim) * 0.02)
        self.mab = MAB(dim, heads, dropout)

    def forward(self, x, padding_mask, need_weights=False):
        seed = self.seed.expand(x.shape[0], -1, -1)
        return self.mab(
            seed,
            x,
            kv_padding_mask=padding_mask,
            need_weights=need_weights,
        )


class InducedSetAttentionMIL(nn.Module):
    def __init__(
        self,
        input_dim: int,
        hidden_dim: int = 64,
        heads: int = 4,
        n_inducing: int = 32,
        dropout: float = 0.1,
    ):
        super().__init__()
        if hidden_dim % heads:
            raise ValueError("hidden_dim must be divisible by heads")
        self.project = nn.Sequential(
            nn.LayerNorm(input_dim),
            nn.Linear(input_dim, hidden_dim),
            nn.GELU(),
        )
        self.isab1 = ISAB(hidden_dim, heads, n_inducing, dropout)
        self.isab2 = ISAB(hidden_dim, heads, n_inducing, dropout)
        self.pma = PMA(hidden_dim, heads, dropout)
        self.head = nn.Sequential(
            nn.LayerNorm(hidden_dim),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, x, mask, return_attention=False):
        # x: [batch, cells, input_dim]; mask=True means valid cell.
        padding_mask = ~mask.bool()
        h = self.project(x)
        h = self.isab1(h, padding_mask)
        h = self.isab2(h, padding_mask)
        pooled, weights = self.pma(
            h, padding_mask, need_weights=return_attention
        )
        logits = self.head(pooled[:, 0]).squeeze(-1)
        if return_attention:
            # weights: [batch, heads, 1, cells]
            return logits, weights
        return logits
```

### 3.5 多次确定性测试抽样

不能继续使用单次随机 test bag。每个 test donor 固定做8次独立、可重现的无放回抽样，平均 probability。

```python
import numpy as np
import torch


@torch.no_grad()
def predict_mc(model, make_batch, bags, donor_ids, device, base_seed, passes=8):
    model.eval()
    predictions = []
    attention = []
    for mc in range(passes):
        xb, mask = make_batch(
            bags,
            donor_ids,
            seed=stable_seed(base_seed, "test", mc),
            replace=False,
        )
        logits, weights = model(
            xb.to(device), mask.to(device), return_attention=True
        )
        predictions.append(torch.sigmoid(logits).cpu().numpy())
        attention.append(weights.cpu().numpy())
    return np.mean(predictions, axis=0), attention
```

保存每个 MC pass 的预测。最终 AUC 使用8次 probability均值；同时报告 MC 标准差，量化 cell sampling noise。

### 3.6 训练与验证约束

- outer-test donor 在模型选择、early stopping、标准化和采样策略选择中完全不可见；
- inner validation 必须按 donor 和 label 分层；
- class imbalance 使用 outer-train 计算的 `pos_weight`，不能用全队列；
- checkpoint 选择依据 inner validation log loss，不按 outer AUC；
- 每个 outer fold 保存 train/validation/test donor IDs；
- deterministic algorithms 若不支持某算子，记录具体非确定性来源，不得静默关闭审计；
- 训练失败保留日志，不反复更换 seed 直到成功。

### 3.7 Attention 审计

attention 权重只能作为模型诊断，不能直接写成机制解释。

至少输出：

- 不同 folds/repeats 中同一 donor 的 attention entropy；
- 各 lineage 获得的总 attention mass；
- attention mass 与 ISG score、UMI、genes/cell、mitochondrial fraction 的相关；
- 高权重细胞在不同 seed 间的 Jaccard 稳定性；
- 打乱 cell order 后 prediction 数值不变性；
- 随机初始化/label permutation 下 attention 模式作为负对照。

若 attention 主要追踪 QC 或 ISG，不得称为新细胞状态发现。

### 3.8 P7A-N 停止规则

满足任一条件即停止扩展，不再增加模型容量：

- Set Transformer/attention MIL 的 repeated-CV paired difference 未稳定超过 `scgpt_mean`；
- 95% paired interval跨0且效应极小；
- outer-repeat 方差明显大于 mean/KME；
- attention 对 seed 或 cell sampling 极不稳定；
- 主要权重由 cell count、UMI、mitochondrial fraction或ISG解释；
- 只能通过查看测试结果后调参获得优势。

只有在预定义配置下出现稳定、配对、可复现增益，才允许提出第二版层级 lineage-aware attention。

## 4. 并行与资源计划

P8.4 的经典方法主要是 CPU 工作；P7A-N 使用5090。允许两条任务并行，但共享盘读取必须错峰。

- 大型 embeddings 先串行复制到本地 SSD；
- GPU任务最多 `DataLoader num_workers=2`；
- 不启动48进程；
- CPU正式并发初始8，稳定后最多16；
- GPU训练按一个 `cohort × method × outer fold` 串行队列；
- checkpoint和日志持续写本地，按fold原子同步到持久盘；
- 任何阶段不得直接从 `.incoming_*` 读取。

## 5. 给 Worker 的第一条指令

Supervisor 首先只下发以下任务：

```text
1. 修复P8.3两个归档元数据问题：完整completion-report hashes，以及142/143 manifest policy说明；不重跑。
2. 只读恢复P8.4预注册cell-count levels、methods、cohorts、repeats和seed定义。
3. 审计现有P7A神经模型结果与torch_models.py，不启动新训练。
4. 输出P8_4_DEFINITION_RECOVERY.md和P7A_N_IMPLEMENTATION_AUDIT.md。
5. 报告后等待Supervisor分别签署P8.4 frozen config和P7A-N exploratory config。
```

这两份定义未冻结前，不允许启动正式 P8.4 或新的 attention run。
