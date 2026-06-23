# RheumLens A800 前期环境、材料与项目章程包

版本：2026-06-20

这个压缩包不是单纯的环境安装包，而是 RheumLens 项目在 A800 服务器上的**启动上下文**。它包含：

1. 项目总目标、科学问题与核心假设；
2. 当前论文资产、最新权威结论与项目边界；
3. SLE 公开数据版图和纳入/排除策略；
4. A800 环境、存储和软件准备；
5. 需要复制、下载和核查的材料；
6. 分阶段执行路线、成功标准、临床意义和最终交付物；
7. 机器可读配置、环境 YAML、下载和验收脚本。

## 建议阅读顺序

- `docs/00_project_charter.md`：项目总纲与研究目标
- `docs/01_scientific_questions_and_hypotheses.md`：核心科学问题与假设
- `docs/02_scope_data_and_authoritative_state.md`：数据范围与当前权威状态
- `docs/03_environment_and_materials.md`：服务器环境、目录和材料清单
- `docs/04_execution_phases_and_deliverables.md`：执行阶段、任务和交付物
- `docs/05_success_criteria_claims_and_non_goals.md`：成功标准、可主张结论与非目标
- `docs/06_clinical_relevance.md`：临床转化意义和证据边界

## 机器可读配置

- `configs/project_charter.yaml`
- `configs/dataset_scope.yaml`

## 一键准备

```bash
bash scripts/bootstrap_a800_env.sh
bash scripts/download_public_metadata.sh
python scripts/check_materials.py \
  --manifest templates/materials_manifest.tsv \
  --root /autodl-fs/data/rheumlens
```

任何正式计算前，先冻结：donor inclusion、labels、authoritative folds、input hashes、model revisions 和 environment manifests。
