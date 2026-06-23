# RheumLens A800 环境、存储与材料准备清单

版本：2026-06-20  
目标服务器：Ubuntu 22.04；A800-SXM4 80GB ×1；144 CPU cores；1 TiB RAM；`/autodl-fs/data` 可用约 929GB。

## 0. 数据版图结论与项目边界

RheumLens 将 SLECA 2026 预印本作为公开数据版图锚点。SLECA 在其数据冻结时声明整合所有公开可得的 SLE scRNA-seq 数据，共 366 个样本、超过 400 万细胞、8 项研究。

项目据此冻结以下判断：

1. 在 SLECA 数据冻结与“公开、可取得、PBMC/血液、donor-level SLE-vs-HC、样本量足以评估”的纳入标准下，公开数据已基本见底。
2. 三个主 benchmark 队列固定为 `GSE135779`、`GSE285773`、`GSE174188`。
3. `GSE137029` 只进入 donor-overlap 审计；在与 `GSE174188` 完成去重前，不作为独立外部验证队列。
4. `GSE250024` 已解析：仅 3 名 SLE donor、每人 3 个疫苗时间点，共 9 个 scRNA-seq 样本；不适合作为主病例-健康对照 benchmark，可作为纵向/疫苗扰动辅助任务。
5. 小样本、皮肤组织或非 SLE 对照来源只保留 accession 级 provenance，不投入主计算。
6. “见底”是带日期和纳入标准的项目结论，不是永久断言；SLECA 数据冻结后的新公开数据须通过 amendment 登记。

机器可读版本见：`configs/dataset_scope.yaml`。项目总纲见 `docs/00_project_charter.md` 和 `configs/project_charter.yaml`。

---

## 1. 存储布局

系统工作区 `/root/autodl-tmp` 只有 50GB，仅放代码、短期日志和软链接。所有环境、模型、数据、缓存与结果放 `/autodl-fs/data`。

建议根目录：

```bash
export RL_ROOT=/autodl-fs/data/rheumlens
mkdir -p "$RL_ROOT"/{
  code,envs,tools,
  data/raw,data/processed,data/external,data/manifests,
  metadata,immutable,splits,
  models/scgpt,models/geneformer,models/references,
  tokenized,embeddings,representations,
  results,logs,tmp,
  cache/huggingface,cache/torch,cache/pip
}
```

建议将代码仓库放在：

```text
/autodl-fs/data/rheumlens/code/RheumLens_v3_code
```

创建环境变量文件：

```bash
cat > "$RL_ROOT/env.sh" <<'EOF'
export RL_ROOT=/autodl-fs/data/rheumlens
export HF_HOME=$RL_ROOT/cache/huggingface
export HUGGINGFACE_HUB_CACHE=$RL_ROOT/cache/huggingface/hub
export TRANSFORMERS_CACHE=$RL_ROOT/cache/huggingface/transformers
export TORCH_HOME=$RL_ROOT/cache/torch
export PIP_CACHE_DIR=$RL_ROOT/cache/pip
export TMPDIR=$RL_ROOT/tmp
export PYTHONNOUSERSITE=1
export TOKENIZERS_PARALLELISM=false
EOF
source "$RL_ROOT/env.sh"
```

### 建议空间预算

| 类别 | 初始预算 |
|---|---:|
| 现有/下载的 raw 与 processed 数据 | 150–250GB |
| Geneformer tokenized 数据 | 80–180GB |
| scGPT/Geneformer embeddings | 10–40GB |
| 多折 representation cache | 80–180GB |
| 结果、bootstrap、permutation、日志 | 50–100GB |
| 模型、环境和缓存 | 30–80GB |

不下载 FASTQ/SRA 原始 reads，除非后续明确需要重新比对。

---

## 2. 环境分层

不要把 benchmark、scGPT、Geneformer 全塞进一个 Python 环境。使用四个隔离环境。

### 2.1 `rheumlens-core`

用途：数据审计、fold-contained PCA/HVG、所有固定聚合、统计、bootstrap/permutation、绘图、CPU/GPU 小模型。

- Python 3.11
- PyTorch 2.5.1 + CUDA 12.4
- NumPy 1.26.4
- pandas 2.2.3
- SciPy 1.13.1
- scikit-learn 1.5.2
- anndata 0.10.9
- scanpy 1.11.5
- pyarrow、h5py、zarr、statsmodels
- POT、geomloss（OT/距离方法）
- pytest、ruff、mypy（开发检查）

环境文件：`envs/environment-core.yml`。

### 2.2 `rheumlens-geneformer`

用途：官方 Geneformer tokenizer 和 V2 embedding 提取。

- Python 3.10
- PyTorch 2.5.1 + CUDA 12.4
- Geneformer 固定 revision：`04c2b2e...`
- V2-104M 作为已验证 checkpoint
- V2-316M 作为后续扩展，不取代已验证主 sensitivity
- 官方 token dictionary、gene median 文件、rank-value tokenizer

所有 clone、checkpoint 和 tokenizer 文件必须记录 Git commit 与 SHA256。

### 2.3 `rheumlens-scgpt`

用途：复现已有 scGPT embedding 或补提取新的 cell subset。

- Python 3.10
- scGPT 0.2.4
- 使用已验证 whole-human checkpoint
- 建议优先复用已有有效 embedding；只有 cell inclusion 或 extraction spec 改变时才重提取

### 2.4 `rheumlens-r`

用途：标准 coloc 与 count-based pseudobulk 复核。

- R 4.3/4.4
- coloc 5.2.3
- DESeq2、edgeR、limma
- data.table、arrow、jsonlite

---

## 3. 系统工具

建议一次安装：

```bash
apt-get update
apt-get install -y \
  tmux git git-lfs rsync curl wget aria2 pigz parallel jq tree htop \
  build-essential cmake pkg-config \
  libhdf5-dev libopenblas-dev liblapack-dev \
  libcurl4-openssl-dev libssl-dev libxml2-dev

git lfs install
```

`tmux` 用于长任务；`aria2c` 用于断点下载；`pigz` 用于多核压缩；`parallel` 用于 CPU 任务调度。

---

## 4. 启动前必须找到或下载的数据材料

### 4.1 项目不可变资产（必须先从旧项目复制）

以下文件比重新下载数据更优先：

- 三队列权威 donor inclusion 表；
- 三队列 disease labels；
- 权威 primary folds；
- 30-seed repeated-CV split definitions；
- 原始 counts 或明确的 `layers['counts']`；
- 已验证 scGPT cell embeddings；
- GSE135779 matched-500 cell ID 清单；
- 有效 Geneformer V2-104M 的 tokenized/embedding manifest；
- 现有 OOF predictions、bootstrap arrays、permutation nulls；
- 公共元数据结构化表；
- `FINAL_MANIFEST_SHA256.csv`。

复制后放到：

```text
$RL_ROOT/immutable/
$RL_ROOT/metadata/
$RL_ROOT/splits/
```

任何脚本不得修改 `immutable/`。

### 4.2 三个主队列

#### GSE135779

必须准备：

- QC 前或 QC 后可追溯的 raw integer count matrix；
- cell barcode、gene symbol/Ensembl ID；
- donor ID、SLE/HC label；
- pediatric primary donor list；
- adult transfer donor list；
- Batch、Age、Gender、Race、Ethnicity、collection year；
- UMI、genes/cell、mitochondrial fraction；
- 原作者 donor-level 20-cluster counts；
- 项目 marker-defined T-lineage/monocyte flags；
- matched-500 cell IDs。

#### GSE285773

必须准备：

- 26 donor raw counts/processed H5AD；
- donor ID 与 SLE/HC label；
- QC summaries；
- GEO SOFT/MINiML 与 BioSample metadata；
- 公开记录中缺失的 age/sex/ancestry/treatment/site/batch 保持 NA，不推断。

#### GSE174188

必须准备：

- 项目已使用的 H5AD/CELLxGENE-derived CD4 subset；
- donor/sample/library/suspension mapping；
- disease state、sex、ethnicity、processing cohort；
- managed/control-only 与 single-sample-donor flags；
- 权威 donor inclusion 和 folds。

GEO 本身不公开 raw/processed 文件，项目继续以已经归档的合法可用对象为准。

### 4.3 只下载 metadata 或辅助材料的 accession

#### GSE137029

第一阶段只下载：

- GEO SOFT/MINiML；
- sample names、pool IDs、donor aliases；
-论文补充 donor table；
- 与 GSE174188 的 overlap audit 所需字段。

只有完成去重并明确任务后，才下载约 2.2GB 的公开 SLE MTX。

#### GSE250024

记录为辅助纵向数据：

- 3 名 SLE donor；
- day0、day1–2、day23；
- 9 个样本；
- raw TAR 约 472MB。

主 benchmark 不需要下载。若后续验证 FOCUS 对疫苗诱导状态变化的灵敏度，再下载。

#### 小样本/非 PBMC/非 SLE来源

只保存 GEO accession 页面、SOFT、论文和排除原因，不下载大矩阵。

---

## 5. 模型和参考资源

### 5.1 scGPT

需要：

- whole-human checkpoint；
- model config；
- gene vocabulary；
- scGPT 0.2.4 source/version record；
- checkpoint SHA256；
- extraction parameters：`max_length=1200`、CLS、frozen weights。

### 5.2 Geneformer

需要：

- Geneformer source repository at fixed revision `04c2b2e...`；
- V2-104M checkpoint；
- V2 token dictionary；
- gene median file；
- tokenizer config；
- model/checkpoint/tokenizer SHA256；
- successful matched-500 extraction manifest。

### 5.3 Gene annotation

需要固定版本：

- GRCh38 Ensembl GTF；
- gene symbol ↔ Ensembl ID mapping；
- 去版本号规则；
- duplicate symbol/Ensembl collapse log；
- 每个队列 vocabulary match summary。

优先复用有效 Geneformer 重跑时已经生成的映射表，避免重新在线映射造成漂移。

### 5.4 机制和细胞状态资源

第一阶段准备：

- 固定 15-gene ISG panel；
- Reactome/GO 开放基因集；
- FOCUS query bank YAML；
- 每个 query 的名称、描述、正向基因、负向基因、适用 cell type、来源和版本。

可选：CellTypist/Azimuth 等冻结 reference model，用于软 cell-type assignment；必须记录模型版本。

---

## 6. 下载优先级

### P0：不开 GPU 前必须完成

1. 旧项目不可变资产复制与 SHA256；
2. 三个主队列 raw counts/processed objects 定位；
3. 权威 donor/fold 文件定位；
4. scGPT 和有效 Geneformer 模型资产定位；
5. Gene mapping 与 metadata 表定位；
6. 磁盘容量和文件权限检查。

### P1：环境安装时同步完成

1. SLECA preprint 与 supplementary provenance；
2. 五个主要 accession 的 GEO SOFT/MINiML；
3. GSE137029 overlap-audit metadata；
4. GSE250024 metadata（raw TAR 暂不下）；
5. Reactome/GO gene sets；
6. 论文 PDF 和补充表。

### P2：仅在对应分析启动时下载

- GSE137029 大 MTX；
- GSE250024 raw TAR；
- V2-316M checkpoint；
- 外部 cell reference model；
- 额外疾病 benchmark 数据。

---

## 7. 数据验收标准

每个主队列必须生成 `dataset_audit.json`，至少包含：

- cells × genes；
- donor 数、case/control 数；
- 每 donor cell count 分布；
- raw counts 是否非负且近整数；
- obs/var 必需字段；
- donor-label 一致性；
- duplicate cell IDs；
- duplicate donor aliases；
- gene symbol/Ensembl mapping rate；
- scGPT/Geneformer vocabulary overlap；
- cell type/lineage 可用性；
- metadata missingness；
- 文件 SHA256。

任何主数据不通过审计，不进入 embedding 或 benchmark。

---

## 8. 环境验收标准

```bash
nvidia-smi
python - <<'PY'
import torch
print(torch.__version__, torch.version.cuda)
print(torch.cuda.is_available())
print(torch.cuda.get_device_name(0))
print(torch.cuda.get_device_properties(0).total_memory / 1024**3)
PY
```

要求：

- PyTorch 能识别 A800；
- 单卡总显存约 80GB；
- 1000-cell scGPT/Geneformer smoke test 可完成；
- benchmark 单元测试通过；
- `/root/autodl-tmp` 不被模型缓存占满；
- 所有缓存实际写入 `/autodl-fs/data/rheumlens/cache`。

---

## 9. 前期完成标志

只有以下文件齐全，才开始 240 小时正式计算：

```text
immutable/authoritative_donors.csv
immutable/authoritative_folds.csv
immutable/repeated_cv_splits.csv
immutable/FINAL_MANIFEST_SHA256.csv
metadata/cohort_metadata.parquet
data/manifests/datasets.tsv
models/manifests/models.tsv
configs/dataset_scope.yaml
configs/prereg.yaml
results/preflight/environment_report.json
results/preflight/data_audit_summary.csv
```

前期目标不是立即获得 AUC，而是确保之后的每个结果都能追溯到固定输入、固定模型和固定 estimand。
