# Supervisor 控制层

从 `00_SUPERVISOR_PROMPT.md` 开始。`01_CURRENT_STATE.md` 是当前进度入口；原始 `docs/`、`configs/`、`envs/`、`scripts/` 和 `templates/` 来自最初的项目总纲/环境准备包，用于背景与基础约束。

本包不包含服务器密码、大型数据、embedding、失效实例环境或未经审查的临时并行 wrapper。当前凭据由用户通过独立安全渠道提供。

建议 Supervisor 首次创建：

```text
control/STATE_LEDGER.csv
control/jobs/<job_id>/JOB_SPEC.yaml
control/incidents/
control/decisions/
```

正式结果仍应落在项目规定的持久化结果目录；控制层只保存小型状态、规格、日志索引和manifest。
