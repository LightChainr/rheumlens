# RheumLens Supervisor 交接包

请从以下文件开始：

```text
SUPERVISOR/00_SUPERVISOR_PROMPT.md
SUPERVISOR/01_CURRENT_STATE.md
SUPERVISOR/03_ARCHITECTURE_AND_ROLES.md
SUPERVISOR/04_JOB_AND_RESOURCE_PROTOCOL.md
SUPERVISOR/05_P8_P10_EXECUTION_PLAN.md
```

本包以最初的“项目总纲、环境与材料准备包”为基础，新增了截至2026-06-21的项目进度、三次实例失联复盘、P8.2恢复方案和Supervisor–Worker控制协议。

当前最重要任务是：从持久盘复核状态，隔离失败的48-`spawn` wrapper，保留有效的 `scgpt_mean` 正式 permutation，并用通过审查、可checkpoint恢复的runner完成另外两个P8.2 primary methods。

不要从头重跑P4–P7，不要使用旧服务器凭据或PID，不要把临时Worker汇报直接视为正式验收。
