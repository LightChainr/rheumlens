# 基础设施事件与并行 wrapper 复盘

## 结论

第三次失联与错误并行实现具有较强因果关联，但没有内核日志时不能宣称已100%证明。前两次事件可能涉及共享盘 I/O、实例/入口平台或其他重任务，证据不足。

## 第三次事件链

1. 正式 permutation 原实现单核运行，PCA 约33.6秒/rep。
2. 为提高利用率，临时生成48-worker wrapper。
3. wrapper 使用 `multiprocessing.get_context("spawn")`。
4. 48个 worker 各自加载完整 NPZ，而非父进程一次加载或 memmap。
5. 出现本地盘读取、NPZ解压、页错误和内存分配峰值。
6. SSH 进入 `Connection timed out during banner exchange`。
7. 代码复查发现 checkpoint、隔离和完整性断言存在确定性缺陷。

## 不能接受的归因方式

- “端口开着，所以 sshd 正常”；
- “48 workers 加载导致 banner 慢是预期行为”；
- “三次都是 sshd 自己僵死”；
- “换GPU型号必然解决”；
- “cron定时重启sshd可以预防宿主机/I/O问题”。

## 下一实例的必要证据

Supervisor 应让独立监控 Worker 每60秒写本地采样，并每5分钟写小型持久快照：

- PID/PPID/state/wchan/RSS/CPU；
- D-state 与 zombie 分开统计；
- CPU/memory/io PSI；
- MemAvailable；
- GPU温度、显存、Xid；
- 本地/共享盘空间和 inode；
- SSH连接数；
- 当前 job_id、rep进度和最近 checkpoint。

发生异常时按证据分类为 shared I/O、local I/O/memory、GPU、SSH并发、项目代码或平台/宿主机，不允许只凭现象命名根因。

## 成本策略

P8 permutation、bootstrap、repeated CV 和大部分敏感性分析是 CPU 任务。除非需要重新提取 foundation-model embedding，不应持续租用昂贵 A800。优先使用高核/大内存 CPU 实例；GPU任务与CPU统计任务分离。
