# RheumLens A800 项目 Agent 任务说明与阶段任务表

版本：2026-06-20  
适用环境：Ubuntu 22.04.4 / NVIDIA A800-SXM4-80GB / 144 CPU cores / 1 TiB RAM  
项目根目录：`/autodl-fs/data/rheumlens`

---

## 1. 文档目的

本文件是 RheumLens 项目的服务器端执行任务书，用于约束 Agent 按固定研究设计完成：

1. 环境和代码验收；
2. 真实材料恢复与真实性审计；
3. 权威结果复现；
4. Geneformer 全链路修复与复现；
5. matched-500 沙盒验证；
6. 三队列完整 benchmark；
7. 统计、稳健性、跨队列和临床相关性分析；
8. 最终结果、图表、论文材料和可复现 manifest 交付。

本项目不以获得阳性结果为目标。任何胜、平、负结果均须保留，不得为了提高 AUC 修改预注册分析、数据纳入、fold、标签或方法定义。

---

## 2. 项目核心问题

项目要回答的核心问题是：在严格 donor-level、防泄漏、强基线和跨队列审计框架下，单细胞基础模型的疾病预测价值来自细胞编码器本身，还是主要取决于供体级聚合方式？

重点包括：

- 比较 raw expression、pseudobulk、PCA、ISG、composition、scGPT、Geneformer 等表示；
- 比较均值、矩、分位数、核均值嵌入、距离方法、MIL 和原创 donor aggregation 方法；
- 检验细胞数、混杂变量、模型各向异性、随机 query 和细胞移除对结论的影响；
- 区分队列内预测、跨队列泛化和潜在临床意义，避免过度主张。

---

## 3. 不可违反的研究规则

1. 独立统计单位必须是 donor，不能把 cell 作为 train/test split 单位。
2. 每个 donor 在每个 OOF 设计中必须恰好获得一次测试预测。
3. PCA、HVG、scaler、threshold、prototype、query selection、feature selection 等任何数据适应步骤，只能在 outer-train donors 内拟合。
4. target cohort 标签不得用于跨队列特征选择、调参或模型选择。
5. raw counts、log-normalized expression、scGPT embedding 和 Geneformer embedding 必须使用不同字段和文件，禁止混用。
6. 不得覆盖 `immutable/`、权威 folds、原始 manifest 或已完成运行。
7. 不得因为结果不显著而改变方法、超参数、纳入标准或评价指标。
8. 缺失材料不得用模板、随机数据、重新划分的 folds 或推断标签代替。
9. 未通过阶段验收门时不得进入下一阶段。
10. 所有运行必须记录配置、输入 SHA256、代码版本、环境版本、日志和输出 manifest。

---

## 4. 固定服务器路径

| 用途 | 路径 |
|---|---|
| 项目根目录 | `/autodl-fs/data/rheumlens` |
| 运行代码 | `/autodl-fs/data/rheumlens/code/RheumLens_A800_full_code` |
| v5 原始执行包 | `/autodl-fs/data/rheumlens/packages/RheumLens_A800_v5_original` |
| Agent 起始提示词 | `/autodl-fs/data/rheumlens/AGENT_START_PROMPT.md` |
| 环境交接说明 | `/autodl-fs/data/rheumlens/V5_ENV_HANDOFF.md` |
| 主配置 | `/autodl-fs/data/rheumlens/configs/project.a800.yaml` |
| matched-500 配置 | `/autodl-fs/data/rheumlens/configs/project.matched500.yaml` |
| 材料清单 | `/autodl-fs/data/rheumlens/manifests/materials_manifest.tsv` |
| 权威资产 | `/autodl-fs/data/rheumlens/immutable` |
| metadata | `/autodl-fs/data/rheumlens/metadata` |
| folds | `/autodl-fs/data/rheumlens/splits` |
| 数据 | `/autodl-fs/data/rheumlens/data` |
| embeddings | `/autodl-fs/data/rheumlens/embeddings` |
| 模型 | `/autodl-fs/data/rheumlens/models` |
| 结果 | `/autodl-fs/data/rheumlens/results` |
| 日志 | `/autodl-fs/data/rheumlens/logs` |
| 初始材料审计 | `/autodl-fs/data/rheumlens/results/preflight/materials_audit_v5_initial.csv` |

Micromamba 环境实际位于：

```text
/root/autodl-tmp/rheumlens/micromamba_root
```

共享盘 inode 较紧张，禁止将 Python/R 环境迁回共享盘。数据、模型和结果仍保存在共享盘。

---

## 5. 当前已知状态

### 5.1 已完成

- v5 原始包 SHA256 校验通过；
- `rheumlens 4.0.0` 已 editable 安装；
- Python 静态编译通过；
- 11 项单元/集成测试通过；
- synthetic smoke 通过；
- PyTorch 2.5.1+cu124 可识别 A800；
- 四个隔离环境已建立：core、scGPT、Geneformer、R；
- GSE135779 matched-500 已跑完 27 个方法；
- folds 结构已解释为 44 donors × 5 folds；
- 已验证每个 donor 恰好作为 test donor 一次。

### 5.2 当前主要阻塞

Geneformer matched-500 当前运行未正确复用已验证的权威 embedding：

- 当前 `geneformer_mean` AUC 约 0.48；
- 无 scaler 的临时诊断约 0.63，但不能据此改变权威管线；
- 已确认另一台服务器曾由另一 Agent 官方重提取并验证 `GSE135779_v2_104M_all.csv`；该文件不假定存在于当前服务器；
- 该 CSV 的 SHA256 为 `745dac02033ac60b59f33345c82e09a00fd39524f2775631a7a122b56d7ef2b7`；
- 该权威 CSV 对应 Geneformer donor AUC 0.920、permutation P=0.000999、lineage ARI 0.255–0.302，且 768 个维度无常数维；
- 当前异常运行使用的 `GSE135779_v2_104m_matched500.npz` 并非该 Agent 直接生成，本地脚本也未引用它；
- NPZ 可能是其他 Worker 对 CSV 的二次封装，其生成脚本、source-row 对齐和 SHA256 尚未确认；
- 当前不能断言根因是 scaler、embedding 维度或 embedding 无信号；
- 当前服务器应重新执行官方提取链路，并以外部已验证 CSV 的哈希和指标作为复现参照；同时审计当前来源不明 NPZ 的派生链路。

### 5.3 材料状态注意事项

Agent 曾报告完整项目“0 项必需缺失”，但该结论必须重新验证。现有材料检查脚本在某些路径为空时可能沿用 manifest 原状态，因此不能仅凭 `status=ready` 判断资产真实存在。

所有必需资产必须同时满足：

- 有明确绝对路径；
- 路径实际存在；
- 文件大小合理；
- 可成功读取；
- 实算 SHA256；
- donor/cell/gene 维度符合预期；
- provenance 可追溯。

### 5.4 Geneformer 提取管道最新状态

截至当前，Agent 报告以下事实：

| 审计项 | 当前记录 |
|---|---|
| Geneformer revision | `04c2b2e`，需在最终 manifest 中保存完整 commit SHA |
| 模型 | V2-104M，文件约 418 MB |
| checkpoint SHA256 | `fff5cba29ddd8...`，最终报告必须保存完整 64 位哈希 |
| 模型结构 | hidden size 768，12 层，12 heads |
| token dictionary | `token_dictionary_gc104M.pkl`，约 426 KB |
| gene median | `gene_median_dictionary_gc104M.pkl`，约 1.5 MB |
| gene mapping | GRCh38，20,634 genes |
| 当前 extraction 设置 | `emb_mode=cls`、`emb_layer=-1` |
| tokenizer | rank-value、`collapse_gene_ids`、chunk size 512 |
| transformers | 从不兼容版本调整为 4.46.3 |
| Geneformer 安装 | 本地源码 editable install |
| 提取状态 | NPZ → AnnData → tokenize → extract 正在运行 |

另一台服务器上的已验证参考资产为：

```text
文件：GSE135779_v2_104M_all.csv
SHA256：745dac02033ac60b59f33345c82e09a00fd39524f2775631a7a122b56d7ef2b7
donor AUC：0.920
permutation P：0.000999
lineage ARI：0.255–0.302
维度：768，无常数维
```

该 CSV 不假定存在于当前服务器。当前任务是用相同模型、代码 revision、tokenizer、gene median、mapping、输入细胞和提取定义重新生成当前服务器版本，并比较输出文件哈希、结构和下游指标。外部 SHA256 是参考标识；只有字节级输入和输出完全相同时才要求哈希一致。

现有 `embeddings/geneformer/GSE135779_v2_104m_matched500.npz` 来源不明，应立即标记为 `QUARANTINED_UNKNOWN_PROVENANCE`，不得覆盖权威 CSV 或继续作为正式输入。任何重新提取结果也必须写入带 run ID 的新路径，不能替换权威 CSV。

当前存在两条并行任务：第一，在当前服务器重新完成官方提取；第二，调查未知 Worker 的 CSV/embedding → NPZ 二次封装是否发生 source-row、cell ID、donor ID、label 或行顺序错位。只有新提取与派生链路审计完成后，才能更新根因结论。

`emb_mode=cls` 和 `emb_layer=-1` 也必须有官方文档、历史有效 manifest 或冻结方案支持。不得因为某一层 AUC 更高而事后选择提取层。

---

## 6. 每次会话的固定启动步骤

```bash
export RL_ROOT=/autodl-fs/data/rheumlens
source "$RL_ROOT/env.sh"
cd "$RL_ROOT/code/RheumLens_A800_full_code"
micromamba env list
nvidia-smi
```

Python 和 CLI 默认使用：

```bash
micromamba run -n rheumlens-core <command>
```

长任务必须使用 `tmux`，并将 stdout/stderr 写入独立日志目录：

```text
$RL_ROOT/logs/v4/<UTC_RUN_ID>/
```

---

## 7. 大阶段总览

| 大阶段 | 目标 | 进入条件 | 退出条件 | 当前状态 |
|---|---|---|---|---|
| P0 | 环境与代码冻结 | 服务器可连接 | 测试、smoke、GPU 全通过 | 已完成 |
| P1 | 材料真实性审计 | P0 通过 | 所有必需资产有路径、SHA256 和结构报告 | 待重新验收 |
| P2 | 数据契约与 folds 审计 | P1 通过 | donor、label、cell、fold 无泄漏且可追溯 | 部分完成 |
| P3 | 权威结果复现 | P2 通过 | 主要历史结果在容差内复现 | 未完成 |
| P4 | Geneformer 全链路审计 | P2 通过 | 权威 embedding 或重提取结果可复现 | 阻塞中 |
| P5 | matched-500 正式验收 | P3、P4 通过 | 所有预注册方法和模态通过验收 | 未通过 |
| P6 | 三队列固定 benchmark | P5 通过 | 固定方法 OOF 与汇总完成 | 未开始 |
| P7 | 可学习及原创方法 | P6 通过 | 模型与原创方法完成 | 未开始 |
| P8 | 统计与稳健性审计 | P6/P7 通过 | bootstrap、permutation、nulls 完成 | 未开始 |
| P9 | 跨队列与临床分析 | P8 通过 | 合法 source-target 分析完成 | 未开始 |
| P10 | 汇总、图表与交付 | P9 通过 | 结果、图表、manifest、报告齐全 | 未开始 |

---

## 8. 详尽任务表

### P0：环境与代码冻结

| 小阶段 | 任务 | 检查方法 | 交付物 | 通过标准 |
|---|---|---|---|---|
| P0.1 | 加载固定环境 | `micromamba env list` | 环境列表 | 四个环境存在 |
| P0.2 | 验证 GPU | `nvidia-smi`、PyTorch CUDA | GPU 报告 | A800 可见，CUDA=True |
| P0.3 | 静态编译 | `python -m compileall -q src scripts` | 日志 | 无错误 |
| P0.4 | 单元测试 | `pytest -q tests` | `unit_tests_latest.txt` | 11/11 通过 |
| P0.5 | synthetic smoke | `rheumlens --stage smoke` | smoke 输出 | 全部 synthetic 方法成功 |
| P0.6 | 冻结代码来源 | 对比 v5 原始包 | diff 报告 | 核心实现无未授权修改 |

Agent 不得重新实现 benchmark、aggregator、split 或 evaluation engine。新增诊断代码必须放在独立 `diagnostics/` 目录，不能修改正式实现。

### P1：材料真实性与 provenance 审计

| 小阶段 | 任务 | 必须检查 | 交付物 | 通过标准 |
|---|---|---|---|---|
| P1.1 | 枚举所有资产 | 绝对路径、类型、大小 | `asset_inventory.csv` | 无空路径 |
| P1.2 | 计算哈希 | SHA256 实算 | `asset_sha256.csv` | 所有必需文件有哈希 |
| P1.3 | 验证可读性 | H5AD/NPZ/CSV/checkpoint 实际打开 | `asset_readability.json` | 无损坏文件 |
| P1.4 | 审计来源 | 原始目录、压缩包、下载 URL、日期 | `asset_provenance.csv` | 每项来源可追溯 |
| P1.5 | 冻结不可变资产 | donor、labels、folds、cell IDs | `immutable_manifest.csv` | 权威资产只读 |
| P1.6 | 区分项目范围 | matched-500 与 full project 分开 | 两份审计表 | 不混称“0 缺失” |

必须重点核验：

- authoritative donor inclusion 和 disease labels；
- primary folds 和 repeated-CV splits；
- GSE135779、GSE285773、GSE174188 数据；
- raw counts 与 lognorm；
- scGPT embeddings；
- Geneformer embeddings、checkpoint、token dictionary、gene median；
- GSE135779 matched-500 cell IDs；
- covariates、composition、query bank；
- `FINAL_MANIFEST_SHA256.csv` 或同等权威记录。

禁止仅根据 manifest 中已有 `status` 字段判定资产就绪。

### P2：数据契约、对齐与 folds 审计

| 小阶段 | 任务 | 检查项 | 交付物 | 通过标准 |
|---|---|---|---|---|
| P2.1 | donor 表审计 | donor 唯一性、标签一致性 | `donor_audit.csv` | 无重复或冲突 |
| P2.2 | cell 对齐 | cell ID 集合和顺序 | `cell_alignment.json` | matched 模态完全一致 |
| P2.3 | gene 对齐 | symbol/Ensembl、重复、缺失 | `gene_alignment.json` | 映射可解释 |
| P2.4 | folds 审计 | train/test 互斥 | `fold_audit.csv` | 无 donor 泄漏 |
| P2.5 | OOF 唯一性 | 每 donor 测试次数 | `oof_uniqueness.csv` | 每 donor 恰好一次 |
| P2.6 | 类别分布 | 各 fold SLE/HC | `fold_class_balance.csv` | 无单类测试 fold |
| P2.7 | modality contract | raw/lognorm/embedding 类型 | `modality_contract.json` | 无混用 |

GSE135779 当前 folds 的 220 行必须保留以下已验证结构：44 donors × 5 folds；每个 fold 35/36 train donors 和 8/9 test donors；每个 donor 恰好作为 test donor 一次。

### P3：权威结果复现

| 小阶段 | 方法 | 目标 | 交付物 | 停止条件 |
|---|---|---|---|---|
| P3.1 | donor-mean HVG | 复现权威基线 | OOF、AUC、配置 | 超出容差 |
| P3.2 | donor-expression PCA | 复现权威基线 | 同上 | 超出容差 |
| P3.3 | raw pseudobulk | 复现权威基线 | 同上 | 超出容差 |
| P3.4 | ISG scalar | 复现权威基线 | 同上 | 超出容差 |
| P3.5 | composition | 复现权威基线 | 同上 | 超出容差 |
| P3.6 | scGPT mean | 复现权威结果 | 同上 | 超出容差 |
| P3.7 | Geneformer mean | 复现权威结果 | 同上 | 超出容差并进入 P4 |

每项必须报告：

- donor 数和类别数；
- fold-wise AUC；
- pooled OOF AUC；
- OOF 行数与 donor 唯一性；
- 使用的 scaler、classifier、`C` 和 class weight；
- 与权威结果的绝对差值；
- 输入和配置 SHA256。

### P4：Geneformer 全链路专项审计

这是当前最高优先级阻塞阶段。

#### P4.1 历史权威结果定位

当前已获得另一台服务器约 0.92 AUC 的参考文件名、完整 SHA256 和验证指标。当前服务器不要求先找到该文件，而是要复现其提取定义和验证结果。

需要交付：

- 原始 OOF predictions；
- 原始运行配置；
- 原始 folds；
- embedding manifest；
- checkpoint/tokenizer/input SHA256；
- 生成日期、代码 revision 和命令。

如果当前服务器使用完全相同的输入、版本和序列化方式，输出 SHA256 应与参考哈希比较；若序列化或输入字节不同，则不能仅凭哈希不一致判定失败，必须继续比较 cell/source-row 对齐、数值误差和固定下游指标。

#### P4.1A CSV → NPZ 派生链路审计

对 `GSE135779_v2_104m_matched500.npz` 必须回答：

1. 是谁、何时、用哪个脚本生成；
2. 输入是否为权威 `GSE135779_v2_104M_all.csv`；
3. 输入 CSV 的 SHA256 是否与权威哈希一致；
4. CSV 每行代表 cell、donor 还是其他实体；
5. NPZ 中 embedding、cell ID、donor ID、label 的字段名和 shape；
6. source-row 到 NPZ row 的一一映射如何保存；
7. matched-500 cell ID 筛选发生在封装前还是封装后；
8. 是否排序、去重、merge、groupby 或 reset index；
9. 是否发生字符串类型、前导零、大小写或 donor alias 变化；
10. NPZ 的完整 SHA256 及生成日志。

必须用 cell ID 或明确的 source-row key 做 join，禁止按当前行位置假设 CSV 与 NPZ 对齐。

#### P4.2 checkpoint 与软件 provenance

检查并记录：

- Geneformer repository commit；
- V2-104M checkpoint 精确路径和 SHA256；
- config 中 hidden size、层数、模型类别；
- token dictionary SHA256；
- gene median SHA256；
- transformers、datasets、torch 版本；
- 是否使用官方 tokenizer 和官方 embedding API。

#### P4.3 输入 H5AD 审计

检查：

- 500 cells/donor 是否来自权威 matched cell IDs；
- 总细胞数是否为 22,000；
- donor、label、cell ID 是否完整；
- gene ID 是 symbol 还是 Ensembl；
- Ensembl version suffix 是否处理一致；
- duplicate genes 如何处理；
- raw counts 是否为非负整数；
- tokenizer 所需字段是否正确；
- 输入细胞顺序和输出 embedding 顺序是否可回溯。

#### P4.4 提取定义审计

不得自行采用“layer -4”或“cosine gap > 0.05”。必须以历史有效 manifest 或官方固定方案为依据。

需要明确：

- 提取层；
- token-level 到 cell-level 的 pooling；
- 是否排除特殊 token；
- 长序列截断；
- batch size、dtype 和 device；
- 输出是否为 cell embedding；
- 输出维度及其来源；
- 是否发生全零、重复行、NaN/Inf 或顺序错位。

当前正在运行的提取任务还必须记录：

- 完整启动命令；
- UTC run ID；
- tmux session 和日志路径；
- 输入 NPZ/AnnData 的 SHA256；
- tokenized dataset 路径和 SHA256/manifest；
- 完整 commit SHA，而不是短 SHA；
- checkpoint、config、token dictionary、gene median 和 mapping 的完整 SHA256；
- transformers、torch、datasets、anndata 和 Geneformer 安装版本；
- batch size、CUDA device、峰值显存、开始/结束时间和退出码。

#### P4.5 embedding 数值质量审计

生成：

- shape、dtype；
- NaN/Inf/zero row 数；
- 每维方差分布；
- 每 cell 范数分布；
- 重复行比例；
- PCA explained variance；
- donor mean 距离；
- 按 cell type 分层后的 within/between donor 距离；
- 随机抽样可视化。

原始 cosine 高或 within/between cosine 差小，不能单独证明 embedding 无效。Transformer embedding 可能存在各向异性，必须结合 centered、PCA 和下游结果判断。

#### P4.6 最小受控复现

只运行 `geneformer_mean`：

- 使用完全相同的 44 donors；
- 使用完全相同的 22,000 cells；
- 使用 authoritative folds；
- 使用权威 classifier、scaler、`C` 和 class weight；
- 禁止根据结果选择有无 scaler；
- 输出每折和 pooled OOF AUC。

如果要比较有无 scaler，只能作为标记清楚的诊断 sensitivity，不能替换主分析。

新 embedding 完成后，必须先通过以下技术门，再允许运行 `geneformer_mean`：

1. shape 与预期一致，目标应为 22,000 cells × 768 features；
2. 22,000 个 cell IDs 完整、唯一且与 matched-500 权威列表一致；
3. 每 donor 恰好 500 cells，共 44 donors；
4. 输出顺序可通过 cell ID manifest 回溯，不能假设 tokenizer 保持输入顺序；
5. donor labels 与 folds 100% 对齐；
6. 无 NaN、Inf、全零行或异常高重复行比例；
7. 文件写入带 run ID 的新路径，并生成 SHA256；
8. 未知来源旧 embedding 被隔离但不删除；
9. 主分析使用冻结的 scaler/classifier/`C`/folds；
10. 输出每折 AUC、pooled OOF AUC 和每 donor 唯一 OOF 检查。

#### P4.7 P4 验收标准

满足以下任一条件才可退出 P4：

1. 当前服务器官方重提取完成，并用固定 donor/folds 管线复现接近参考值的 donor AUC；且
2. CSV → matched-500 输入的 cell/source-row 对齐有可验证 manifest；且
3. 新派生 NPZ 的输出与 CSV 子集在数值和行级 provenance 上一致；或
4. 若派生 NPZ 被证明错误，则隔离 NPZ，恢复以权威 CSV 为源的正确转换流程。

### P5：matched-500 正式验收

| 小阶段 | 任务 | 交付物 | 通过标准 |
|---|---|---|---|
| P5.1 | 固定输入 | matched cell manifest | 44×500 cells，无错位 |
| P5.2 | 固定基线 | raw、PCA、ISG、composition | 全部完成 |
| P5.3 | scGPT 方法 | mean、moments、KME、FOCUS、RED 等 | OOF 完整 |
| P5.4 | Geneformer 方法 | mean、moments、KME、FOCUS 等 | P4 已验收 |
| P5.5 | 完整 27 方法审计 | success/fail/skip 表 | 每方法状态明确 |
| P5.6 | 汇总 | method summary、fold AUC | 文件齐全 |

不得将 “27 个命令退出码为 0” 等同于 “27 个方法科学验收通过”。还要检查：

- OOF 行数；
- donor 唯一性；
- NaN/Inf；
- 单折失败；
- 常数预测；
- AUC 方向；
- 方法输入模态是否正确。

### P6：三队列固定 benchmark

主队列：

- GSE135779；
- GSE285773；
- GSE174188。

执行顺序：

1. 每队列单独 validate；
2. 每队列强基线；
3. frozen foundation-model mean；
4. moment ladder；
5. quantile/tail；
6. KME/robust KME/MMD；
7. Bures/Sliced-Wasserstein；
8. prototypes；
9. cell-type pseudobulk/scFeatures；
10. 汇总固定方法结果。

每个队列必须单独生成：配置快照、输入 manifest、OOF、fold metrics、diagnostics 和错误清单。

### P7：可学习患者级模型与原创方法

#### P7A 可学习模型

- Deep Sets；
- attention/top-k MIL；
- Set Transformer；
- MixMIL approximation；
- 其他预注册低容量患者级模型。

要求：超参数只能在 outer-train 内选择；小样本下禁止无限扩大模型容量。

#### P7B 原创方法

- ccKME-U；
- U-DER；
- FOCUS-Lite；
- Language-FOCUS；
- RED；
- GDS；
- DonorCLR；
- IDE（仅训练 cohorts 选择）。

每个原创方法必须有：方法定义、训练范围、输入依赖、复杂度、失败模式和 ablation。

### P8：统计与稳健性审计

| 小阶段 | 任务 | 默认规模 | 说明 |
|---|---|---:|---|
| P8.1 | paired stratified bootstrap | 10,000 | donor-level 配对 |
| P8.2 | label permutation | 1,000 | fold 结构保持 |
| P8.3 | repeated CV | 30 seeds | 使用冻结定义 |
| P8.4 | cell-count sensitivity | 预注册水平 | donor 内下采样 |
| P8.5 | surrogate nulls | mean/moment2 | 检验额外分布信息 |
| P8.6 | covariate sensitivity | 固定方案 | 不推断缺失变量 |
| P8.7 | kernel diagnostics | 固定网格 | 检查退化 |
| P8.8 | query controls | random/removal | 防止事后挑选 |

只有主 OOF 管线通过后才能运行本阶段。探索性 Geneformer 诊断 permutation 不得混入正式统计结果。

### P9：跨队列与临床分析

1. 只选择生物学舱室可比的 source/target；
2. 所有拟合和特征选择只使用 source；
3. target 标签仅用于最终评价；
4. 分别报告 discrimination、calibration 和 failure modes；
5. 队列内高 AUC 不得表述为临床泛化；
6. 对 treatment、age、sex、ancestry、site、batch 缺失保持 NA，不推断。

### P10：汇总、图表和最终交付

最终至少需要：

- OOF master table；
- method summary；
- paired method differences；
- within-cohort 和 cross-cohort matrices；
- robustness dashboard；
- method-family figures；
- fold-level diagnostics；
- failed/skipped run registry；
- paper tables；
- code manifest；
- data manifest；
- model manifest；
- environment manifest；
- `RESULTS_MANIFEST_SHA256.csv`；
- 最终执行报告。

---

## 9. 运行状态定义

每个阶段只能使用以下状态：

- `NOT_STARTED`：尚未开始；
- `RUNNING`：正在执行；
- `BLOCKED_MATERIALS`：缺数据、模型、fold 或 provenance；
- `BLOCKED_VALIDATION`：输出存在科学或技术异常；
- `FAILED`：执行失败；
- `COMPLETED_TECHNICAL`：命令完成但尚未科学验收；
- `ACCEPTED`：技术和科学验收均通过。

禁止仅用“成功”描述存在关键异常的阶段。例如当前 matched-500 应标为 `BLOCKED_VALIDATION`，而不是成功。

---

## 10. Agent 每次汇报模板

```text
大阶段：
小阶段：
状态：
开始/结束时间：
执行命令：
代码版本：
配置路径与 SHA256：
输入路径与 SHA256：
donor/cell/gene 数：
fold 结构：
输出路径：
关键指标：
验证结果：
异常与失败：
日志路径：
是否满足退出条件：
下一步：
```

涉及模型 embedding 时必须额外报告：

```text
模型名称：
checkpoint SHA256：
代码 revision：
token dictionary SHA256：
gene median SHA256：
输入 H5AD SHA256：
提取层：
pooling：
shape/dtype：
NaN/Inf/zero/duplicate rows：
cell ID alignment：
```

---

## 11. 当前 Agent 的立即任务

当前不得启动 bootstrap、正式 permutation 或全 benchmark。立即任务按顺序为：

1. 冻结现有 matched-500 输出和日志，生成 SHA256；
2. 生成完整项目和 matched-500 两套真实资产清单；
3. 为每项必需资产提供绝对路径、大小、可读性和 SHA256；
4. 保存外部参考信息：文件名、完整 SHA256、AUC、P 值、ARI 和维度检查；
5. 保持当前官方重提取任务写入独立 run-ID 路径，不覆盖旧 NPZ；
6. 记录当前提取的完整输入、版本、命令、tokenized manifest 和 cell/source-row 映射；
7. 审计未知 NPZ 的生成脚本、输入哈希和 source-row/cell ID 对齐；
8. 提取完成后执行 shape、cell ID、donor、数值质量和顺序对齐检查；
9. 用新提取结果运行固定 `geneformer_mean`，禁止事后改 scaler/层/超参数；
10. 提交 fold-wise AUC、pooled OOF、完整 provenance，并与外部参考指标比较；
11. 等待 P4 验收决定，再恢复 P5。

不得让未知来源 NPZ 覆盖当前新提取结果。外部参考 CSV 不在当前服务器并不构成材料缺失；它是复现基准。不得不断尝试不同层、pooling、scaler 或超参数来追逐 0.92。

---

## 12. 完成标准

项目只有在以下条件全部满足时才能声明完成：

- 环境、代码和包版本可重建；
- 所有真实输入和权威资产有 SHA256；
- donor-level 防泄漏审计通过；
- 权威结果已复现，或未复现原因有可验证证据；
- 所有方法有完整 OOF 和运行状态；
- 统计与稳健性分析使用冻结方案；
- 跨队列分析无 target leakage；
- 所有失败和负结果被保留；
- 图表和论文结论不超出证据边界；
- 最终结果目录有完整 manifest。
