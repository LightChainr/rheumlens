# RheumLens Supervisor 增量更新：RTX PRO 6000 Blackwell 新节点

请将本文件视为对既有 Supervisor 交接包和旧 Worker 提示词的**增量覆盖指令**。发生冲突时，以本文件为准。不要因为上下文中仍保留旧服务器信息而恢复旧节点、旧PID、旧tmux或旧本地环境。

## 1. Supervisor 可能不知道的最新基础设施状态

前三台 A800 实例均已停止使用。新的有效节点位于不同集群：

```text
cluster: singapore-a
GPU: NVIDIA RTX PRO 6000 Blackwell Server
VRAM: approximately 96 GB
OS: Ubuntu 22.04
base image: PYTORCH_27_CUDA128_UBUNTU2204
system Python: 3.12.3
system PyTorch: 2.7.0+cu128
CUDA runtime: 12.8
cuDNN: 9.7.1
idle shutdown policy: 24h
expected shared root: /autodl-fs/data/rheumlens
```

SSH地址与密码由用户通过安全渠道另行提供。不得把密码写入日志、报告、脚本、shell历史或共享盘。

旧信息全部失效：

- A800、singapore-b；
- PyTorch 2.5.1+cu124 作为“系统基础镜像”的假设；
- 旧实例的 PID、tmux、`/root/autodl-tmp`、micromamba环境和SSH配置；
- 对第三台失败48-worker任务仍在运行的假设。

共享盘上的已验收资产、结果和manifest仍可能有效，必须按实际文件和hash恢复。

## 2. Blackwell兼容性的正确处理

新镜像已是 PyTorch 2.7.0+cu128。PyTorch 2.7是正式支持Blackwell/CUDA 12.8的版本，因此不应把新节点误判为“PyTorch 2.5.1不支持Blackwell”。

但“版本号正确”仍不等于GPU已验收。首次连接后进行一次实际CUDA运算：

```bash
python3 - <<'PY'
import sys, torch
print("python", sys.version)
print("torch", torch.__version__)
print("cuda build", torch.version.cuda)
print("available", torch.cuda.is_available())
print("device", torch.cuda.get_device_name(0))
print("capability", torch.cuda.get_device_capability(0))
print("arch list", torch.cuda.get_arch_list())
x = torch.randn(2048, 2048, device="cuda")
y = x @ x
torch.cuda.synchronize()
print("cuda_matmul_ok", float(y.mean()))
PY
```

成功后记录 `BLACKWELL_BASE_ACCEPTED`。失败时只阻塞未来GPU任务，不阻塞P8 CPU统计；保留完整错误，禁止在系统Python中临时降级/覆盖torch。

## 3. 系统环境与项目环境必须隔离

系统Python 3.12.3和PyTorch 2.7.0+cu128只作为Blackwell基础环境。不得直接在系统Python中执行项目级 `pip install`、降级torch或安装Geneformer/scGPT依赖。

项目环境安装在新节点本地盘，而不是共享盘：

```text
rheumlens-core        冻结的CPU统计与项目测试环境
rheumlens-scgpt       Python 3.10/3.11，未来确需scGPT GPU运行时使用
rheumlens-geneformer  Python 3.10/3.11，未来确需Geneformer GPU运行时使用
rheumlens-r           原冻结R环境
```

不要假定共享盘 `env.sh` 中的旧本地路径仍有效。先读取环境manifest/yaml，再在新节点本地重建必要环境。当前只需优先建立 `rheumlens-core`；scGPT/Geneformer embedding 已存在且此前通过provenance审计，不应为了环境验证重新提取。

环境重建不等于重跑P4–P7。`rheumlens-core` 完成后必须达到：

```text
compileall PASS
pytest 11/11
synthetic smoke 9/9
```

## 4. 当前分析不需要GPU

当前恢复目标是 P8.2 formal permutation：

```text
scgpt_mean
donor_expression_pca
kme_multiscale@scgpt
```

它们读取已经生成的embedding/表达数据，主要使用NumPy/scikit-learn/PCA/逻辑回归/KME，是CPU统计任务。P8运行时应显式：

```bash
export CUDA_VISIBLE_DEVICES=""
```

GPU空闲是预期现象，不得因为租用了96GB Blackwell就擅自把CPU统计改写成GPU实现。只有未来确需重新提取embedding或重新训练神经网络时，才启用隔离的Blackwell GPU环境。

## 5. 最新项目状态，不得退回旧阶段

已验收或应只读复核：

```text
P4 Geneformer: ACCEPTED
P5 GSE135779 matched-500: ACCEPTED, 27/27
P6 GSE285773: ACCEPTED, 11/11
P6 GSE174188_CD4: Agent已完成Gate A/B和17/17；从共享盘复核，不重跑固定benchmark
P6三队列总结: 已报告生成
P7A: 已报告5个learnable method_id，需从实际状态表恢复名称
P7B: 6个originals，产物位于P6目录，不复制第二套OOF
P8.1: 10,000 paired bootstrap已报告生成
```

P8.2准确状态：

- `scgpt_mean`：报告1000/1000完成，p=0.000999；必须从共享盘程序化复核后保留，不能因为换服务器而重跑。
- `donor_expression_pca`：旧串行partial无checkpoint，后续错误48-worker运行无效；正式结果需要从冻结seed序列重新运行。
- `kme_multiscale@scgpt`：正式1000次尚未可靠完成，需要重新运行。
- `raw_pseudobulk` 和 `donorclr` 的100次结果只是Gate B sanity，不是P8.2正式primary permutation。
- P8.3–P8.8、P9、P10尚未完成，P8.2结束后不能直接跳P9。

## 6. 第三台服务器事故的确定性代码事实

旧48-worker `spawn` wrapper不得复制到新节点、不得恢复、不得合并其partial。已确认的问题包括：

1. 48个spawn worker分别加载完整数据，造成并发读取、解压、内存分配和page-fault高峰；
2. `np.save("*.tmp", ...)` 实际生成 `*.tmp.npy`，随后rename原路径必然失败；
3. PCA与KME共用checkpoint/shard目录，会跨方法污染；
4. checkpoint只在旧实例本地盘，不能抵抗实例失效；
5. 异常被静默转换为NaN；
6. 未强制 `n_finite == n_requested`；
7. 未验证serial/parallel相同rep结果一致；
8. 旧KME smoke存在可疑重复/错误调用。

新节点本地盘不会包含旧实例 `/root/autodl-tmp/rheumlens/perm_par`。不要浪费时间寻找旧本地目录。只需清点共享盘中可能留下的失败结果，并移动到明确的quarantine命名空间；不得删除。

## 7. Candidate runner不能只做grep审计

任何名为 `perm_parallel_v2.py` 的候选文件都不是自动可信资产。Supervisor必须让Implementation/Audit角色完整审查源码，至少证明：

- 父进程一次加载只读数据后fork，或使用新节点本地memmap；
- worker不从共享盘独立重复加载完整数据；
- PCA、KME和不同config使用独立run/checkpoint目录；
- 原子保存实现不会触发NumPy自动扩展名问题；
- 每个rep具有唯一rep_id和确定性seed；
- 异常写入errors.jsonl和traceback，不静默接受NaN；
- checkpoint周期性同步到持久共享盘；
- resume拒绝错误code/config/input/folds hash；
- 最终严格断言1000个finite、唯一且完整rep；
- observed AUC与权威OOF重算一致；
- serial和parallel对相同rep seeds逐项一致；
- 输出目录由显式 `--run-dir`/`--checkpoint-dir` 控制，不使用隐式默认目录。

仅进行静态grep、计算脚本hash或看到“使用fork”字样，均不足以通过代码审查。

## 8. Resume是正式运行的硬门槛

旧Worker提示词中如存在：

```text
RESUME_NOT_SUPPORTED; FORMAL RUN MUST COMPLETE IN ONE SESSION
```

该规则现已作废。新的硬规则是：

```text
kill/resume test不通过，禁止启动1000-rep formal。
```

至少完成：

1. 20-rep测试；
2. 中途优雅终止；
3. 从持久checkpoint恢复；
4. 最终20/20，无重复、无缺失；
5. 与不中断的同seed参考结果一致；
6. 修改config/hash后必须拒绝错误resume。

## 9. 新节点资源策略

不要再次从4-worker smoke直接跳到48个spawn。也不要为了“保守”长期单核运行。

正式放大顺序：

```text
4 workers: serial/parallel数值一致性
8 workers: kill/resume测试
16 workers: 2–5分钟资源smoke
32 workers: 默认formal规模
48 workers: 仅在32稳定、输入单次加载、总RSS和I/O有余量时允许
```

所有worker在导入NumPy/SciPy/sklearn前设置：

```text
OMP_NUM_THREADS=1
MKL_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
```

至少保留8个逻辑核给控制面。大型输入串行stage到本地并流式hash；正式计算期间不并发下载、解压、全量hash或第二个重任务。

## 10. 审计输出的落盘顺序

新节点刚接入时，先确认共享盘和本地盘类型。初始health日志先写本地控制目录，确认共享盘正常后再原子同步。不要在尚未验证共享盘时把每条监控记录直接写共享盘。

正式runner必须显式接收：

```text
--run-dir
--checkpoint-dir
--method
--n-reps
--n-workers
--base-seed/seed-table
--input-manifest
--folds
--config
```

每个tmux job写：

```text
state.json
heartbeat.json
command.txt
environment.txt
exit_code.txt
errors.jsonl
MANIFEST_SHA256.json
```

## 11. Supervisor现在应分派的第一批任务

### Worker A：Host/Environment Audit

- 验证Blackwell实际CUDA矩阵运算；
- 识别本地盘、共享盘、CPU核数、RAM、inode；
- 检查系统Python但不修改；
- 在本地重建/验证 `rheumlens-core`；
- compileall、11/11、9/9。

### Worker B：Shared-State Recovery Audit

- 只读恢复P4–P8.1状态；
- 程序化验收 `scgpt_mean` 1000 reps；
- 清点并隔离共享盘失败parallel产物；
- 不读取/恢复旧实例本地partial。

### Worker C：Permutation Runner Review

- 完整审查候选runner；
- 修复后提交code hash、单元测试和资源模型；
- 不直接启动formal。

### Worker D：Independent Validation

- 对修复runner执行serial/parallel、kill/resume、错误hash resume测试；
- 给Supervisor提交PASS/FAIL证据。

Supervisor只有在A/B/C/D相应gate通过后，才授权Execution Worker启动 `donor_expression_pca` 的正式1000次。随后独立验收，再执行KME。

## 12. 首次更新汇报

首次连接并完成最小只读检查后，Supervisor应报告：

```markdown
# Blackwell Supervisor Update
- host / cluster:
- system Python / torch / CUDA:
- CUDA matmul: PASS / FAIL
- local/shared storage:
- rheumlens-core status:
- shared authoritative state recovered:
- scgpt_mean formal audit:
- failed wrapper shared artifacts:
- candidate runner review status:
- next authorized gate:
```

不要在首次汇报前启动正式permutation。完成必要gate后可自主推进，无需逐条等待用户；涉及科学定义冲突、已验收资产覆盖、无法解释的hash不一致或新增大规模付费资源时才停止等待用户。
