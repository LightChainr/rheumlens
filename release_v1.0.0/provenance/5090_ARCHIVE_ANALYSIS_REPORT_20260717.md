# RTX 5090 cap500 归档本地分析报告

**分析日期：** 2026-07-17  
**归档：** `patient_generalization_5090_cap500_complete_20260717.tar.gz`  
**归档 SHA256：** `5f08f97007f4358039afb85185816d298ce222420c0e221f2732d6d3c5f32c5d`

## 1. 归档识别

服务器结果以 `patient_generalization_5090_cap500_complete_20260717.tar.gz` 作为不可变原始输入，并在独立目录中解包后完成校验和二次分析。本公开记录仅保留可复现的归档标识符、哈希和分析结果，不包含本机目录或服务器连接信息。

## 2. 完整性与复现

- 归档共 342 个条目，解包后约 629 MB。
- 两个数据集共核验 89 个 embedding shard；metadata 和 NumPy 文件 SHA256 全部通过。
- `GSE174188`：261 donors，127,129 cells，1152 dimensions，cap500，seed 1。
- `GSE285773`：26 donors，13,000 cells，1152 dimensions，cap500，seed 1。
- Geneformer checkpoint：`ctheodoris/Geneformer`，revision `04c2b2e84da7c0f385c3f9ad8f3ec24bab6650e5`。
- Model weight SHA256：`965ceccea81953d362081ef3843560a0e4fef88d396c28017881f1e94b1246f3`。
- 服务器环境：Python 3.12.3、PyTorch 2.8.0+cu128、Transformers 4.53.2、scikit-learn 1.7.0。
- 本地使用 Python 3.13.7、scikit-learn 1.9.0 重新运行完整评估；gate decisions 完全一致。
- 本地与服务器 AUC/Brier/ECE 的最大数值差约 `1.6e-5`，属于软件版本和浮点实现差异，不改变任何判断。

服务器日志中的 “pooler weights newly initialized” 警告不影响当前结果：提取代码使用 `last_hidden_state[:, 0, :]`，并未使用随机初始化的 BERT `pooler_output`。

## 3. 必须排除的旧输出

归档同时保留了：

`outputs/distributional_pooling_cap500_sampling_order_mismatch/`

该目录的 historical reproduction 标志为 `false`，属于已经发现并纠正的采样顺序错误。论文、图件和 Zenodo 均只能使用：

`outputs/distributional_pooling_cap500/`

后者在两个队列均通过历史 cap500 donor-mean reproduction。

## 4. 固定 pooling 压力测试

原服务器协议比较 coordinate mean 与七种固定聚合：median、10% trimmed mean、mean+SD、source-PCA quantiles、source-PCA moments 和两个 donor-PCA mean。

主要结果保持不变：

- 小目标 `GSE285773 (n=26)` 中，多种方法表面优于 mean；最高为 source-PCA8 donor mean，AUC 0.944。
- 大目标 `GSE174188 (n=261)` 中，这些优势普遍反转；mean 为 0.884，median 为 0.858，source-PCA8 donor mean 为 0.834。
- 七种候选均未通过预设双向 continuation gate，因而无需继续 GPU cap1000 扩展。
- 原协议 14 项 DeLong 检验经 BH 校正后，大目标中 median、PCA32 quantiles、PCA8 donor mean 和 PCA16 donor mean 显著低于 mean。

这支持当前正文的核心表述：**单一方向或小目标队列中的提升不能证明 pooling 改进具有可迁移性。**

## 5. 新发现：AUC 保留但概率尺度崩塌

大源队列到 26 人小目标的 coordinate mean 结果为：

- ROC-AUC：0.8625
- Brier：0.6154
- 10-bin ECE：0.6154
- 26/26 供体预测概率均低于 0.01

也就是说，模型仍可按病例风险排序，但把所有目标供体都置于几乎为零的概率尺度。AUC 在这里不能被解释为可用的疾病概率或临床风险。

source-PCA8 donor mean 在同一目标中达到 AUC 0.950、Brier 0.148，但在 261 人反向目标中降至 AUC 0.834。该结果进一步区分了两个问题：

1. 某种 source-fitted transform 可偶然修复一个方向的概率尺度；
2. 它并不因此成为跨方向稳定的患者表征。

建议在 Results 增加一句，在 Discussion 增加一个集中段落，明确“transported ranking”与“transported probability”不同。DCA 不应基于这些未经外部校准的概率作为主证据。

## 6. 新发现：cohort shift 主导供体嵌入几何

在不使用疾病标签拟合方向的联合 descriptive PCA 中：

- Pooled donor PC1 解释 45.4% 方差；PC2 解释 22.2%。
- 两个队列在 donor-mean Geneformer 空间中几乎完全分离。
- 队列质心距离约为 pooled within-cohort RMS dispersion 的 2.95 倍。
- 同一指标下，病例/对照质心分离仅为 0.59（GSE174188）和 0.73（GSE285773）。

这不是新的监督预测结果，也不能解释为批次的单一来源；两个队列同时存在年龄、招募、实验和处理差异。但它直接可视化了为什么 within-cohort 高 AUC 不足以选择可迁移表示：**队列位移大于疾病轴。**

建议将该图作为主文或高位补充图，并在图注中明确 PCA 仅用于描述、未进入迁移模型。

## 7. 新发现：donor mean 后的有效维度收缩

供体均值的完整 centered spectrum 显示：

| Dataset | Embedding dimensions | Entropy effective rank | Participation rank | PC1 fraction |
|---|---:|---:|---:|---:|
| GSE174188 | 1152 | 6.83 | 3.64 | 48.1% |
| GSE285773 | 1152 | 6.69 | 4.36 | 38.9% |

固定种子抽取 10,000 cells、使用 256-component randomized spectrum 时：

- GSE174188 cell-level entropy rank：51.4；top-256 捕获 94.9% 方差。
- GSE285773 cell-level entropy rank：31.8；top-256 捕获 96.6% 方差。
- cell-level 到 donor mean 的 rank ratio 约为 7.53 和 4.75。

这为“均值池化为何未兑现细胞基础模型的高维信息”提供了直接几何线索：细胞层存在的多维变化在供体均值后集中到约 7 个有效方向。

当前结果仍属于探索性几何证据。进入正文前建议用 5–10 个 cell-sampling seeds 重复 top-256 spectrum，并预先固定 effective-rank 定义。不能直接写成信息论证明或因果机制。

## 8. source-CV 公平性敏感性

服务器 pooling 脚本和原主迁移脚本均为不同方法使用不同的 source-CV 随机折来选择 logistic-regression `C`。这不涉及目标标签泄漏，但削弱了方法间调参条件的严格配对。

已完成三种敏感性：

1. 共享 source-CV 折、原 full-source scaler；
2. 原方法特异折、fold-contained classifier scaler；
3. 共享 source-CV 折、fold-contained classifier scaler。

结果：

- 三种协议下，七种 pooling 候选都没有通过双向 AUC gate。
- mean pooling 在原折下采用 fold-contained scaler 时精确复现原结果，说明 mean 的变化不是 scaler 泄漏导致。
- 26 人 source 中，单次 CV 折可使选择的 `C` 从 `1e-4/1e-3` 跳到 `1e-2`，并使大目标 Geneformer AUC 变化约 0.025。
- 个别本已失败的 PCA 聚合算子对调参协议更敏感，最大 AUC 变化约 0.051。

### 主三表征共享折候选

在 Geneformer、HVG pseudobulk 和 PCA pseudobulk 间共享同一 source-CV 折后：

| Direction | Geneformer | HVG pseudobulk | PCA pseudobulk |
|---|---:|---:|---:|
| GSE174188 -> GSE285773 | 0.900 | 0.944 | 0.981 |
| GSE285773 -> GSE174188 | 0.890 | 0.916 | 0.916 |

在 261 人大目标中，pseudobulk 仍优于 Geneformer：

- HVG minus Geneformer：0.0264；BH-adjusted paired DeLong `q=0.029`。
- PCA minus Geneformer：0.0267；BH-adjusted paired DeLong `q=0.029`。

因此共享折使效果量更保守，但核心结论依然成立。该版本在方法学上优于当前方法特异随机折，建议作为下一版候选主结果。

## 9. cap500 与 cap1000：真正敏感的是正则化

两个 cap 的 donor embeddings 本身非常接近：

- raw cosine median：0.999989（GSE174188）和 0.999964（GSE285773）；
- centered cosine median：0.9982 和 0.9935；
- donor-distance matrix correlation：0.9990 和 0.9967。

固定同一个 `C` 后，cap500 与 cap1000 的外部 AUC 曲线几乎重合。表面出现的 0.859 对 0.890 差异主要来自单次 source-CV 选择了不同 `C`，不能归因于多使用了 500 个细胞。

这要求修改 cell-budget 解释：

> Donor embedding geometry was stable between 500 and 1000 cells, whereas externally transferred probabilities remained sensitive to regularization selected from the small source cohort.

比“500 cells are sufficient”更准确，也形成新的、可引用的方法学结论：**细胞采样稳定不等于下游模型选择稳定。**

## 10. 归档缺口

### 必须修复

1. `archive_manifest.json` 列出 `models/ctheodoris_Geneformer/download_manifest.json`，但压缩包中缺失该文件。
2. GSE285773 输入是目录，embedding manifest 中 `input_sha256=null`；Zenodo 前应保存 GEO raw archive 或逐文件校验清单。
3. 作废的 `sampling_order_mismatch` 目录不应进入公开主结果路径，应移入明确的 `invalidated/` 或仅保存在内部归档。
4. shard JSON 保存的是服务器绝对路径；统计重放不受影响，但本地 resumability 会失败。公开版应改为归档相对路径。

### 已补救

- 已根据锁定的 `download_geneformer.py` 生成 `reconstructed_model_provenance.json`，明确标为 reconstructed，不伪造原始文件大小或时间戳。
- 论文当前已经记录 checkpoint revision、model SHA256 和服务器软件版本。

## 11. 对论文的建议

### 可直接纳入

1. 共享 source-CV 折的主三表征迁移结果，作为更保守、更公平的候选主分析。
2. AUC 与概率校准分离的结果，强化“不应把迁移 AUC 当作临床概率”的论点。
3. 双向 pooling rank reversal 图，支持“单向改进不具可迁移性”。
4. cap500/cap1000 fixed-C path，支持“embedding stability does not imply tuning stability”。

### 完成小型本地敏感性后纳入

1. cohort-shift donor PCA。
2. cell-to-donor effective-rank contraction。
3. 重复 source-CV 正则化选择，而不是依赖一次五折选择 `C`。

### 暂不纳入主结论

- 把低秩收缩写成均值池化失败的理论证明；目前只有两个队列和一个 embedding model。
- 基于未校准外部概率的 DCA。
- 将队列 PCA 分离归结为单一技术批次或单一人口学因素。

## 12. 当前判断

这批文件的价值超过“补齐服务器归档”。它提供了论文此前缺少的解释层：

> Frozen Geneformer preserves disease ranking across cohorts, but donor averaging concentrates cell-level variation into a low-dimensional space dominated by cohort shift; under small-source regularization uncertainty, AUC can remain high even when the transported probability scale fails.

这使论文能够从简单的 pseudobulk-vs-Geneformer benchmark，推进到更有增量的主题：**patient-level representation under cohort shift requires simultaneous control of geometry, calibration, and source-sample tuning stability.**

所有 GPU 密集步骤已经完成。下一步的 repeated source-CV、effective-rank 多种子、图件重排和正文更新均可在本机完成，当前不需要继续租用服务器。
