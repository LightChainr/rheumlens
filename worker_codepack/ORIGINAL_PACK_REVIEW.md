# RheumLens Worker Code Pack 20260621 审查报告

## 结论

```text
Verdict: REFERENCE_ONLY / FORMAL_EXECUTION_NOT_AUTHORIZED
```

压缩包完整，`MANIFEST.json` 中14个文件的大小和SHA256均匹配；压缩包SHA256为：

```text
63a05bcfd7ea53f0ccb1432282b61d1c2626c4345c2d583194ae9377f9e35e28
```

但该包只包含Markdown操作手册与代码片段，没有实际 `perm_parallel_v2.py`。它可以用于讨论和修订，不能单独证明runner可用，也不能授权正式1000-rep permutation。

## P0：必须修复后才能运行formal

### 1. 环境恢复定位错误且会掩盖测试失败

`02_blackwell_host_env_core.md`：

- 项目真实候选目录是 `code/RheumLens_A800_full_code`，脚本只检查 `code/rheumlens`，可能在错误根目录执行compile/pytest。
- 创建的是未锁版本的“最小环境”，不符合冻结复现要求。
- `pytest -q || true` 会把测试失败吞掉。
- 没有执行要求中的9/9 synthetic smoke。
- 文档允许“非P8失败时继续”，与正式gate冲突。

要求：从打包的environment/explicit manifest重建本地环境；明确项目目录；compileall、11/11 tests、9/9 smoke任一失败都不得进入formal。

### 2. scgpt_mean审计可产生假阳性

`03_scgpt_mean_audit.md` 中：

```python
has_method = "scgpt_mean" in text_blob or True
```

该表达式永远为True。审计还存在：

- 仅搜索字符串 `0.000999`，没有从null AUC和observed AUC重算P值；
- 没有验证rep集合恰为0..999；
- 没有逐项验证冻结seed表；
- 没有验证完整fold-contained pipeline、fold/config/input/code identity；
- 没有排除重复行和错误方法结果。

要求：程序化读取唯一指定结果文件；严格验证1000行、rep 0..999、seed表、finite、method/config/folds/input identity，并重算P值。

### 3. 新seed算法可能改变正式统计定义

`06_runner_hardening_reference.md` 使用 `SeedSequence.spawn` 生成rep seeds；历史wrapper和现有状态使用过 `default_rng(base_seed).integers(...)`。两者生成完全不同的置换序列。

要求：不得临时更换seed算法。优先从已完成 `scgpt_mean` 正式结果恢复并冻结实际1000行seed表；PCA和KME复用同一表。若无法恢复，Supervisor必须明确决定，而不是由Worker自行选择。

### 4. resume identity没有校验fold/config文件内容

`make_run_identity` 只保存folds和config的路径字符串，没有SHA256。文件内容变化时resume仍可被接受。

要求：identity必须包括method、n_reps、完整seed-table hash、runner code hash、config hash、folds hash、每个input hash、method parameters hash和scientific pipeline revision。

### 5. kill/resume gate不是真实中断测试

`07_pca_gate_harness.md` 使用 `--stop-after-reps 7`，这是runner主动正常退出，不是SIGTERM/进程丢失恢复。并且：

- `|| true` 会掩盖初始运行的任何错误；
- 没有验证错误config/input hash的resume必须被拒绝；
- 没有验证checkpoint已写到持久盘；
- comparator只比较8-worker interrupted/uninterrupted，不包含真正的1-worker serial vs parallel测试。

要求：分别完成：1-worker vs 4-worker相同rep逐项比较；真实SIGTERM恢复；从持久checkpoint恢复；篡改config/hash后resume拒绝。

### 6. formal launcher可能泄露凭据

`08_formal_execution_controller.md` 将：

```bash
env | sort
```

完整写入共享盘 `environment.txt`，可能泄露token、密码、代理凭据和云服务环境变量。

要求：使用环境变量allowlist，仅记录Python/package/thread/CUDA/path等非敏感字段。

### 7. checkpoint位置违反恢复目标

formal launcher将 `RUN` 和 `CKPT` 都直接设在共享盘。该设计既可能造成高频共享盘I/O，也没有实现“本地高速checkpoint + 周期性紧凑持久快照”的分层策略。

要求：rep shard首先原子写本地；协调器每固定rep数或5分钟生成紧凑持久快照；resume从持久快照恢复。共享盘同时只能有一个发布者。

### 8. formal validator会选错结果文件并错误ACCEPT

`09_formal_validation_pvalue.md`：

- `find_result_csv` 选择目录中最大的CSV，可能选到seed表、诊断表或其他CSV；
- 找不到metric列时选择“最后一个numeric column”，可能把seed当作AUC；
- observed可以缺失；即使P值未计算，recommendation仍固定为 `ACCEPT_RECOMMENDED`；
- 没有严格检查rep集合、行数、重复、seed、method identity、exit code、errors.jsonl或run identity；
- observed可能从任意递归JSON中读取到陈旧值。

要求：runner在state.json中声明唯一正式结果文件；validator只读取该文件和冻结identity；observed缺失、错误日志非空、退出码非0或任何完整性条件失败都必须FAIL。

### 9. combined report可能把失败结果汇总为通过

`11_p8_2_combined_report.md` 只要找到PCA/KME validation JSON就可能给出P8.2接受建议，没有要求其recommendation为PASS。它还通过文本子串判断scGPT状态，容易误判。

要求：combined report读取结构化Supervisor decision ledger，只接受明确 `ACCEPTED` 的run ID/hash；不允许按最新mtime自动挑选结果。

### 10. manifest自包含导致不稳定

`08` 和 `11` 使用重定向创建manifest，同时 `find . -type f` 会把正在写入的manifest自身纳入hash，生成不可复核条目。

要求：显式排除manifest自身，先写临时文件，完成后原子rename。

## P1：重要工程问题

### 11. 原子JSONL并不支持多进程并发

`atomic_append_jsonl` 只是普通append，没有锁，也不是“atomic”文件更新。多个worker写同一个errors.jsonl可能交错。

要求：每worker独立error shard，最后由单一协调器合并；或使用可靠文件锁/队列单写者。

### 12. 完整性断言不够严格

`assert_complete_results` 只检查unique和finite数量，没有要求：

```text
len(df) == n_reps
rep_id集合 == {0,...,n_reps-1}
rep_id无重复
seed与冻结seed表逐项一致
```

必须补齐。

### 13. gate脚本失败时可能没有exit_code

`run_gate_method.sh` 使用 `set -euo pipefail` 后直接运行 `python | tee`。命令失败时shell会在写exit_code之前退出。

要求：围绕runner使用 `set +e`、保存 `${PIPESTATUS[0]}`、恢复 `set -e`，无论成功失败都写exit_code。

### 14. local staging标为optional且无校验

正式多进程运行时本地staging应为强制；复制后必须核验size/hash。还应只stage当前方法需要的输入，而不是每个方法都复制raw/lognorm/scGPT全集。

### 15. promotion过于激进

`01_migration_verify_promote.md` 如果发现已有 `/autodl-fs/data/rheumlens` 会自动移动它。Worker不应自行移动未知现有权威根。

要求：已有正式根时停止并请求Supervisor决策；只有明确确认不存在/是空壳时才能promote。

### 16. A侧Worker无法执行源侧步骤

迁移文档要求在源服务器读取 `/root/autodl-tmp/rheumlens_migrate.exit` 和生成source manifest，但当前Worker只能连接A侧。

要求：使用目标已存在的 `.SOURCE_TRANSFER_STATUS`，并由迁移控制方提供source key manifest；A侧Worker不得把无法执行的步骤标为通过。

### 17. runner候选选择不安全

`head -1 RUNNER_CANDIDATES.txt` 会任意选择第一个同名文件，可能是旧版或失败版。

要求：Supervisor显式指定runner绝对路径和批准SHA256；Worker不得自动挑选。

### 18. 通用label shuffle不能替代项目管线

reference中的 `permute_labels(y)` 只是普通数组洗牌。正式P8.2必须调用已审计的donor-level permutation和完整fold-contained方法执行函数，并证明expression_data/embedding labels一致处理。

## 可保留的部分

- 包级MANIFEST完整；
- 总体阶段顺序正确；
- 禁止重跑P4–P7、禁用旧spawn/partial的约束正确；
- CPU线程限制与渐进式4→8→16→32策略合理；
- Blackwell GPU仅做基础验收、P8强制CPU合理；
- 一次只跑一个formal、明确run目录、保存日志/资源信息的方向正确；
- `atomic_write_text` / `atomic_save_npz` 的基本思路可复用，但需补目录fsync和并发策略。

## Supervisor下一步要求

1. 将本包标记为 `REFERENCE_ONLY`。
2. 不授权任何formal命令。
3. 提供实际 `perm_parallel_v2.py` 及其依赖的项目revision进行独立代码审查。
4. 先修复P0问题并提交第二版code pack。
5. 第二版必须包含可执行runner、单元测试、固定seed表格式、job spec和明确结果schema，而不只是Markdown代码片段。
6. 通过：syntax/unit tests、serial-parallel、真实kill/resume、bad-hash resume rejection和8/16-worker resource smoke后，才可授权1000-rep运行。
