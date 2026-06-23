# Job、资源与持久化协议

## Job spec 必填字段

- job_id、阶段、cohort、method_id、analysis_version；
- 代码 revision/hash；
- input paths + SHA256；
- folds path + SHA256；
- config + SHA256；
- n_reps、base_seed、rep seed算法；
- 统计单位与防泄漏说明；
- 本地和共享输出目录；
- worker数、每worker线程数、RAM预算、I/O策略；
- checkpoint频率与恢复流程；
- smoke和验收断言；
- stop conditions。

## 资源策略

对144核、1TB RAM主机，允许高利用率，但必须渐进验证：

1. 4 workers做serial/parallel数值一致性；
2. 8 workers做恢复测试；
3. 16 workers做2–5分钟资源 smoke；
4. 无D-state、PSI、RSS或SSH异常时升至32；
5. 48 workers仅在输入单次加载/memmap且32-worker稳定后使用。

始终保留至少8个逻辑核给控制面。每个进程：

```text
OMP_NUM_THREADS=1
MKL_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
```

这些变量必须在数值库 import 前设置。禁止 `workers × BLAS threads` 嵌套扩张。

## I/O策略

- 共享盘串行复制到本地 staging，流式计算 SHA256；禁止 `read_bytes()` 读取大文件做hash。
- 只复制方法真正需要的输入。
- 并行阶段不得让每个 worker 各自读取共享盘完整数据。
- 优先父进程加载后只读 `fork`，或转换为本地只读 memmap。
- 不并发下载、解压、hash和正式分析。

## Checkpoint协议

- method/config/run完全隔离目录；
- 每个rep由唯一 `rep_id` 和确定性seed标识；
- 每5–10 reps更新本地checkpoint；
- 每25 reps或5分钟同步紧凑checkpoint到持久盘；
- 每个worker写独立shard，不能并发写同一文件；
- 最终只从rep-level checkpoint按rep_id合并；
- resume必须拒绝错误config/input hash和重复rep；
- SIGTERM/INT处理器应停止领取新rep、完成当前rep、flush并退出。

原子Numpy保存示例：

```python
def atomic_save_npy(path, array):
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "wb") as handle:
        np.save(handle, array)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, path)
```

不要把字符串 `*.tmp` 直接传给 `np.save` 后再 rename 原名，因为 NumPy会自动追加 `.npy`。

## 正式完成断言

```text
n_requested == n_finite == 1000
rep_id == 0..999，唯一且无缺失
seed与冻结seed表逐项一致
无静默异常、NaN或Inf
observed metric与权威OOF重算一致
serial/parallel smoke逐rep一致
empirical P公式正确
结果、日志、config、input和code均有hash
```
