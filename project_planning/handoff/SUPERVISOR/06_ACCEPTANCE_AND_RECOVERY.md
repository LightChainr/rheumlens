# 验收、恢复与停止条件

## 科学验收

- donor-level独立性；
- train/test donor无重叠；
- 每donor OOF恰好一次；
- preprocessing/HVG/PCA/scaling/监督选择只在outer-train拟合；
- permutation重跑完整fold-contained监督管线，不能只打乱最终prediction；
- bootstrap使用同一cohort内配对且outcome-stratified donor indices；
- 缺失协变量保持NA，不推断；
- target labels不参与transfer训练或选择。

## 工程验收

- PID、命令、时间、退出码和资源峰值齐全；
- config/input/code/folds/output hashes齐全；
- checkpoint恢复测试通过；
- 没有引用QUARANTINED；
- 无重复job writer；
- 输出可由manifest定位和重算；
- failed/skipped有明确状态。

## 异常恢复

### 仅SSH断开

tmux/background job继续；外部控制面停止高频重试，指数退避。恢复后从状态文件核验，不凭聊天记忆。

### 实例崩溃/回收

本地PID、tmux和local checkpoint作废。从持久盘 `state.json + rep checkpoints + hashes` 恢复。不能确定完整性的本地残留一律隔离。

### Worker失忆

读取 `STATE_LEDGER.csv`、job spec、checkpoint、latest report和manifest；不得重跑已验收阶段。

### 作业异常

先停止领取新任务，flush当前rep，记录traceback和资源；禁止静默NaN、无限自动重试或直接kill -9。

## 强制停止条件

- OOM/cgroup OOM、GPU Xid；
- 多个相关PID持续D-state超过5分钟；
- 本地/共享盘空间或inode接近耗尽；
- io PSI持续full stall；
- SSH连接风暴或控制面持续失联；
- input/folds/hash冲突；
- OOF donor重复、缺失或泄漏；
- 需要修改冻结科学定义；
- 结果命名空间已有其他writer。

资源异常停止后，Supervisor可以降低worker数并从有效checkpoint恢复；不必为普通资源调优等待用户。
