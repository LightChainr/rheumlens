# Worker overnight execution directive

Supervisor delegation token:

```text
OVERNIGHT_AUTONOMY_P82_V2
```

用户将离线。你可以在本提示词定义的P8.2范围内连续自主推进；gate通过后无需等待逐条确认。你不能改变科学定义，也不能超出Supervisor流水线。

## 当前任务顺序

```text
迁移完成与manifest验证
→ canonical root发布
→ v2代码包部署
→ rheumlens-core 11/11 + 9/9
→ scgpt_mean严格审计
→ PCA完整gates
→ PCA formal 1000 + strict validation
→ KME完整gates
→ KME formal 1000 + strict validation
→ P8.2 combined report/manifest
→ P8.3已有实现的只读调查与spec
```

## 迁移

每10分钟最多检查一次：

```text
/autodl-fs/data/rheumlens.incoming_20260621/.SOURCE_TRANSFER_STATUS
```

只有 `source_exit_code=0` 且 `manifest_ready=1` 才继续。使用同目录 `.SOURCE_KEY_MANIFEST.tsv` 的relpath逐项计算目标size/hash并diff。不能访问B侧不是阻塞，source manifest会自动传到A侧。

diff为空且canonical root不存在才发布。存在冲突立即停止，不得移动/覆盖已有root。

## 代码包

只使用：

```text
/root/autodl-tmp/rheumlens_worker_codepack_v2_20260622.zip
SHA256 29e12a7bccef696a958ea7c6432abdf3f8a22e725c24cc36b99f6db2ef9bfcab
```

旧Markdown包、旧spawn wrapper和旧partial全部禁止。

发布后将v2包部署到：

```text
/autodl-fs/data/rheumlens/supervisor_control/rheumlens_worker_codepack_v2
```

先读 `SUPERVISOR_DEPLOYMENT.md`、`README.md`、`TEST_REPORT.md`，校验内部manifest。

## 环境与审计

- 项目环境安装在本地盘；不修改系统Python。
- 使用包内 `environment-core-p82.yml`。
- compileall、pytest必须11/11、project smoke必须9/9；禁止掩盖失败。
- P8设置BLAS线程=1、`CUDA_VISIBLE_DEVICES=""`。
- 生成input manifest。
- 使用包内strict audit审计scgpt_mean；不能只grep字符串。

## PCA

先运行包内gate controller。要求：

```text
serial_parallel PASS
real TERM/resume PASS
bad-hash resume rejected
16-worker resource smoke PASS
```

通过即视为获得：

```text
AUTHORIZATION: PCA_FORMAL_1000
```

默认32 workers，资源异常降至24/16；禁止单核。使用tmux/nohup和v3 checkpoint。SSH断开时先检查PID/checkpoint，不重复启动。

PCA结束后必须运行strict validator。只有完整1000行、准确seed、finite、exit0、errors空才能报告 `ACCEPT_RECOMMENDED`。

## KME

只有PCA strict validation通过才进入。重复全部gates。通过即获得：

```text
AUTHORIZATION: KME_FORMAL_1000
```

默认24 workers，资源稳定可32。不得与PCA并发。

## 完成

生成：

- P8.2 combined report；
- PCA/KME validation JSON；
- result和environment manifests；
- failed/skipped registry；
- state ledger更新；
- `SAFE_TO_STOP_EXPENSIVE_INSTANCES.md`。

P8.2完成后只调查已有P8.3冻结入口并写spec/smoke计划；没有经过审查的现成实现时，不临时写新formal runner。

## 停止条件

manifest/hash冲突、root冲突、测试失败、scgpt审计失败、gate失败、OOM、持续D-state、严重I/O或科学定义歧义时停止重任务并写BLOCKED报告。不要无限重试，不要删数据，不要用keepalive浪费费用。

普通路径/依赖修复、减少worker数、从有效checkpoint恢复可以自主处理并记录。

晨间汇报必须说明：完成到哪一步、当前PID、P8.2三方法状态、资源峰值、恢复次数、关键SHA256、剩余阻塞、两台昂贵实例是否可以关闭。
