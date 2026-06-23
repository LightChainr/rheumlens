# Supervisor overnight autonomy directive

Authorization token:

```text
OVERNIGHT_AUTONOMY_P82_V2
```

用户将离线休息。你可以在下列冻结范围内自主授权Worker连续推进，不需要逐gate等待用户回复。目标是安全完成迁移收尾、P8.2两项剩余formal及联合验收；不要用无意义负载消耗服务器。

## 已知状态

- B→A单流rsync在后台运行；不要启动第二条传输。
- B侧finalizer将在rsync结束后生成源关键manifest并传到A侧，然后写：

```text
/autodl-fs/data/rheumlens.incoming_20260621/.SOURCE_KEY_MANIFEST.tsv
/autodl-fs/data/rheumlens.incoming_20260621/.SOURCE_TRANSFER_STATUS
```

- 只有状态同时满足：

```text
source_exit_code=0
manifest_ready=1
```

才进入目标manifest验证。

- v2代码包已在A侧：

```text
/root/autodl-tmp/rheumlens_worker_codepack_v2_20260622.zip
SHA256 29e12a7bccef696a958ea7c6432abdf3f8a22e725c24cc36b99f6db2ef9bfcab
/root/autodl-tmp/supervisor_codepack_v2/rheumlens_worker_codepack_v2
```

## 自主推进流水线

按顺序执行；每个gate通过后自动进入下一个，不询问用户。

1. 每10分钟最多一次监控迁移；普通RUNNING不汇报。
2. marker完成后，按 `.SOURCE_KEY_MANIFEST.tsv` 的relpath在incoming计算目标size/SHA256并严格diff。
3. diff为空且 `/autodl-fs/data/rheumlens` 不存在时，发布incoming；否则停止并写BLOCKED报告。
4. 删除一次性transfer public-key授权；不要删除项目资产。
5. 将v2包部署到canonical root的 `supervisor_control/`，再次校验包内manifest。
6. 在本地盘重建包内冻结 `rheumlens-core`；要求compileall、pytest 11/11、smoke 9/9。
7. 生成P8.2 input manifest；严格审计scgpt_mean 1000 reps。
8. 隔离明确属于失败48-spawn wrapper的shared产物；不删除、不扩大匹配范围。
9. 对PCA执行serial/parallel、真实TERM/resume、bad-hash拒绝和16-worker resource gate。
10. Gate全通过后，自动授权：

```text
AUTHORIZATION: PCA_FORMAL_1000
```

默认32 workers；若16-worker smoke出现资源压力则24或16。不要单核运行。
11. PCA完成后运行strict validator。只有1000/1000、exit=0、errors空、identity/seed完整才记录ACCEPTED。
12. PCA通过后对KME执行同样gates，并自动授权：

```text
AUTHORIZATION: KME_FORMAL_1000
```

KME默认24 workers；资源证据良好可32。
13. KME完成并验收后，生成P8.2 combined report、结果manifest和state ledger。
14. P8.2全通过后，不临时发明P8.3–P8.8实现。只从迁移项目中查找已有冻结脚本/config：
    - 若存在明确、已实现且无需科学选择的P8.3入口，可完成小型smoke和执行计划，但不未经代码审查启动完整P8.3。
    - 若不存在，生成 `P8_3_NEXT_EXECUTION_SPEC.md` 并停止重任务。
15. 写 `SAFE_TO_STOP_EXPENSIVE_INSTANCES.md`，明确source A800是否可关停、新节点是否仍有任务。

## 自动停止条件

出现以下任一情况停止后续重任务，但保留证据：

- 迁移exit非0或manifest diff；
- canonical root冲突；
- 关键P4–P8 hash不符；
- 环境测试未达11/11或9/9；
- scgpt seed/结果审计失败；
- runner gate失败；
- OOM、持续D-state、严重I/O stall、SSH控制面明显迟滞；
- 需要改变folds、seed表、observed AUC、primary methods；
- 需要删除/覆盖已验收结果。

普通依赖路径修复、worker数从32降到24/16、从有效checkpoint恢复可自主处理。

## 资源纪律

- 一次一个重任务；
- P8强制CPU，GPU空闲正常；
- input只从共享盘串行stage一次；
- worker不并发读取共享盘；
- 不并发下载、hash、解压和formal；
- 不启动keepalive负载；
- 完成后主动建议停止昂贵实例。

最终向用户提交一份晨间报告：迁移、环境、scgpt、PCA、KME、P8.2、资源峰值、失败/恢复、关键SHA256、是否可关机。
