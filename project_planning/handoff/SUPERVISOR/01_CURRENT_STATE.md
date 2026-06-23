# 当前权威状态（截至 2026-06-21）

本文件区分“已通过关键哈希/结构核验”和“由执行 Agent 最新汇报、需从共享盘复核”的状态。

## 已验收

| 阶段 | 状态 | 关键结论 |
|---|---|---|
| P0 | ACCEPTED | compileall；tests 11/11；smoke 9/9；CUDA 可用 |
| P4 Geneformer | ACCEPTED | V2-104M，geneformer_mean OOF AUC 0.9256，有效 embedding hash 已匹配 |
| P5 GSE135779 matched-500 | ACCEPTED | 27/27 方法；关键结果文件 hash 已匹配 |
| P6 GSE285773 | ACCEPTED | 11/11 方法；26 donors；amendment folds 明确标记且 hash 匹配 |

P4 有效 embedding：

```text
embeddings/geneformer/v2_20260620T0915/GSE135779_matched500.npz
SHA256 87276833b0e7aecded2692915c66c2582ffabd22ba044200c23923845a482163
```

P5 关键 hash：

```text
all_oof.csv 742b321ec372ad459fd9144a5bb3de02f50942bc7c7af71951d8c3d7a053fe24
method_summary.csv 26190064666463c8c3e57f1d4938ad10b76e665c279ff51395278390bf3678d1
```

GSE285773 amendment folds hash：

```text
90fe8a0cb471725f2e3c17193277ece83f9d51ca58a1da6953bfc64f45da775b
```

## 第三台服务器上已报告完成，Supervisor 应从共享盘复核

| 项目 | 最新状态 |
|---|---|
| Gate A GSE174188 scGPT provenance | PROVENANCE_RECOVERED；audit + manifest 已生成 |
| Gate B GSE174188 P6 | ACCEPTED；100-rep sanity 对3方法均为 p=0.0099 |
| P6 GSE174188_CD4 | 17/17；261 OOF rows/method；0 NaN；AUC 0.9355–0.9879 |
| P6 三队列总结 | 已生成；汇报为58个 cohort-method entries，含4 CSV和summary |
| P7A learnable models | 汇报为5个 method_id；必须从实际状态表恢复第五个准确名称 |
| P7B originals | 6个，产物位于 P6 目录；包含 `cckme_u_weighted`；不得复制第二份 OOF |
| P8.1 | 10,000 paired stratified bootstrap 已生成 |

GSE174188 authoritative folds 是5 folds、261 donors，不是5 donors。原 JSON hash：

```text
8fa9755e702ab79d2c8ff3639282b2589379c93ebfd1e719bf07aa43a435e13e
```

## P8.2 当前准确状态

正式 primary methods：

```text
scgpt_mean
donor_expression_pca
kme_multiscale@scgpt
```

| 方法 | 状态 | 处理 |
|---|---|---|
| scgpt_mean | Agent 报告正式1000次完成，p=0.000999 | 核验 n=1000、seeds、完整管线、hash 后保留 |
| donor_expression_pca | 单核运行到约225/1000后被终止，无 checkpoint | 旧部分不使用；由合格并行 runner 从冻结 seed 序列重跑 |
| kme_multiscale@scgpt | 正式作业未可靠完成 | 由合格并行 runner执行 |

## 明确无效的运行

第三台服务器曾启动一个 48-worker `spawn` wrapper。该 wrapper 代码审计失败，其产物不得进入正式结果：

- 每个 worker 独立加载完整数据，引发并发 I/O/解压/内存分配高峰；
- `np.save("*.tmp", ...)` 自动生成 `*.tmp.npy`，随后 rename 原路径必然失败；
- PCA/KME 共用 checkpoint 和 shard 目录，会跨方法污染；
- checkpoint 仅在实例本地盘，实例失效后不可恢复；
- worker 异常被静默转换为 NaN；
- 未强制 `n_finite == n_requested`；
- 未验证 parallel 与 serial smoke 逐 rep 一致；
- KME smoke 中存在可疑/无意义的重复调用。

所有相关目录、日志和产物应标为：

```text
QUARANTINED_FAILED_PARALLEL_WRAPPER_20260621
```

## 未开始

```text
P8.3 repeated CV
P8.4 cell-count sensitivity
P8.5 surrogate nulls
P8.6 covariate sensitivity
P8.7 kernel diagnostics
P8.8 query/removal controls
P9 cross-cohort
P10 final integration
```

## 基础设施事件

三台 A800 实例先后失联。第三次与错误的48-`spawn`加载风暴时间高度一致；前两次缺少充分日志，不能统一归因为 sshd。TCP端口开放但 SSH banner timeout 只说明入口接受连接，不能证明后端 sshd 健康。

旧实例地址、密码、PID、tmux 和本地环境均不属于权威状态。本包不保存凭据；用户通过独立安全渠道提供当前访问方式。
