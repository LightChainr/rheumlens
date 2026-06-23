# Supervisor–Worker 架构与职责

## Supervisor

- 维护唯一 `STATE_LEDGER.csv`；
- 冻结 job spec、输入/配置 hash 和结果命名空间；
- 指派 Worker，不直接接受口头“已完成”；
- 审查代码与资源模型；
- 控制作业并发和一写者原则；
- 验证 smoke、checkpoint恢复和正式输出；
- 决定 GO/NO-GO；
- 汇总失败、跳过和 provenance。

## Worker 类型

### Audit Worker

只读核验数据、folds、hash、OOF和阶段状态；不得修改正式结果。

### Implementation Worker

根据冻结 job spec 实现/修复 runner，提交代码、测试和资源估计；默认不得直接启动正式运行。

### Execution Worker

运行已批准 runner，维护 heartbeat/checkpoint，报告 PID、进度和资源；不得修改科学定义。

### Validation Worker

独立验证结果完整性、serial/parallel一致性、统计公式和 manifest；不得使用实现者的内存结论替代文件证据。

### Monitoring Worker

低开销记录系统与job状态；不得高频SSH轮询或执行递归磁盘扫描。

一个 Agent 可以在不同时段承担不同角色，但同一个 job 的 implementation approval 和 final validation 必须形成明确的两个步骤，不能“写完即宣布通过”。

## Job 状态机

```text
DRAFT
→ SPEC_FROZEN
→ IMPLEMENTED
→ CODE_REVIEWED
→ SMOKE_PASSED
→ RECOVERY_TEST_PASSED
→ RUNNING
→ COMPLETED_TECHNICAL
→ VALIDATED
→ ACCEPTED
```

异常状态：

```text
BLOCKED_INPUT
BLOCKED_SCIENTIFIC
FAILED_CODE
FAILED_RESOURCE
INTERRUPTED_RECOVERABLE
QUARANTINED
SKIPPED_MISSING_ASSET
SKIPPED_NOT_IMPLEMENTED
```

## 一写者原则

每个 `cohort/method/analysis_version` 命名空间同时只有一个 Execution Worker。其他 Worker只能读。发布采用：

```text
local workdir → validation → shared .partial → fsync/hash → atomic rename
```
