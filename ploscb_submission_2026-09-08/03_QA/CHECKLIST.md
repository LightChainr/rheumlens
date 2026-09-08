# 投稿前自查 · v5（落实 issue #9、补充审稿意见与「最小改动版」精修意见后）

日期：2026-09-08（第二轮）

## 统计定义

| 项 | 状态 | 依据 |
|---|---|---|
| 设计矩阵无标签泄漏 | ✅ | `audit/build_design_manifest.py` 由真实 `design_matrix()` 生成清单；标签列进入任一块的次数 **0** |
| 清单与实际模型矩阵一致 | ✅ | 逐块列宽比对 **34/34** 一致，机器校验表随包发布 |
| V_D 与其 p 值是同一统计量 | ✅ | 置换标签推过同一组折；样本内那一对独立命名为 `p_V_D_insample`，不混用 |
| 置换检验命名与条件集一致 | ✅ | 全文 `collection-preserving`；`collection_strata()` 返回并写出自己的分层定义 |
| 置换检验的失效模式已量化 | ✅ | §3.4 / Figure S1 / Table S1：γ=0 时校准（6.7%，95% CI 4.4–10.1，名义 5%），层内 \|corr\|=0.61 时假阳性率升到 89.3%，与 free 检验无差别 |
| 推断流程唯一 | ✅ | 冻结流程为全文唯一推断统计量；调参 AUC 为描述量，差值作为列发布 |
| 0.20 任意阈值已移除 | ✅ | `not estimable` 只保留无混合层与完全共线两种真正不可识别情形 |
| 残差化损失的机制解释已检验 | ✅ | §3.5 / Figure S2 / Table S1：同宽度（92 列）0.045、CMV 形状 0.091、峰值 0.229（48–64 列，非单调），均小于观察到的 0.282；差额写入 Limitations |
| 模拟连续臂真的预测连续结局 | ✅ | ridge + concordance index（二元时等于 AUC） |
| 模拟不平衡臂已分离患病率与阈值错配 | ✅ | `imbalanced` 与 `imbalanced_offset` 两臂 |

## 事实一致性

| 项 | 状态 | 依据 |
|---|---|---|
| CMV 结构描述属实 | ✅ | 36 batch / 46 pool / 46 层，无供者内配对；§4.2 逐项写明 |
| Table 1 verdict 由代码判定生成 | ✅ | `tools/build_table1.py`，不再手工转录 |
| Figure 7A 队列标注 | ✅ | 已由 GSE174188 更正为 GSE135779 |
| §5.5 迁移方向与图对应 | ✅ | 261→26 对应 Fig 8；26→261 只在正文给数值 |
| §9.1 数据集计数 | ✅ | 五个进主筛查、两个只用于 §5 下游分析 |
| 数据集原始文献已引用 | ✅ | 新增 [39-42]；[5]（scVI）移到 §9.12 |
| CMV 设计列数全文统一 | ✅ | 89；walkthrough 改为 import screen 的 `design_matrix()`，不再自行实现 |
| 参考文献编号闭合 | ✅ | 脚本双向校验 |
| 所有 panel 在正文被引用 | ✅ | 脚本校验 |
| Table 1 与筛查输出逐值一致 | ✅ | `tools/build_table1.py` 直接生成；`tools/check_manuscript_numbers.py` 逐条比对 **21/21** |

## 表述

| 项 | 状态 | 依据 |
|---|---|---|
| 正文无 rebuttal / 审稿史语言 | ✅ | `grep -c "referee\|reviewer\|earlier version\|previous version"` = 0 |
| 绝对化措辞已收紧 | ✅ | "Nothing used it"、"face value"、"shows nothing at all" 三处改写 |
| 三组变量不再等同于实验设计 | ✅ | §4.1 分别定义；collection 为主，另两组为预设敏感性分析 |
| 标题计数单位 | ✅ | comparisons，不是 cohorts |
| learned pooling 移入补充 | ✅ | Figure S6 / Table S5 |

## 可复现

| 项 | 状态 | 依据 |
|---|---|---|
| 随机种子跨进程稳定 | ✅ | `hashlib.blake2b` 派生；每格 seed 写入结果表 |
| 无本机绝对路径 | ✅ | `grep -r "/Volumes/Mac Data"` 在 sim/walkthrough/figures/tools/audit 下为 0 |
| 重复脚本已删 | ✅ | `run_design_screen_seed20260908.py`（`--seed` 已覆盖） |
| 全部交付物可由一条命令重建 | ✅ | `bash tools/build_all.sh`：表 → Table 1 → 图 → HTML → 两个一致性检查，任一失败即非零退出 |
| 投稿包由脚本组装 | ✅ | `tools/build_submission.py`；图号取自 HTML builder 的映射，该映射被结构检查绑定到正文图注 |
| 每张补充表都有生成脚本 | ✅ | 手工维护的两张（队列登记表、超参表）移入 `supplementary_src/`，由 `stage_static_tables.py` 复制，不再手工放进 `supplementary/` |

## PLOS 格式

| 项 | 状态 | 依据 |
|---|---|---|
| Abstract ≤ 300 词、无小标题 | ✅ | **292 词**（留 8 词余量），脚本校验 |
| Short title ≤ 70 字符 | ✅ | `Study metadata and single-cell patient classification`，53 字符，脚本校验 |
| Cover letter 约一页 | ✅ | 568 词（原 728），脚本校验 ≤650 |
| Author Summary 150–200 词、位于 Abstract 之后 | ✅ | 199 词 |
| Results / Discussion / Methods 分区可识别 | ✅ | 新增分区标题 |
| 主图 8（决策树）与正文一致 | ✅ | 重画：删除已撤回的 `tuned ≈ frozen` 门槛与 `model-selection artefact` 停止条件；Route A 由「残差化是主结论」改为「几乎无可移除，调整须报告代价」，与 §3.5 / Fig 6C 一致 |
| 置换次数正文=代码 | ✅ | `N_PERM=1000`（floor 1/1001）、V_D 与 design-only AUC null 200（floor 1/201）；`2,000` / `0.0005` 全部清除，脚本比对代码常量 |
| 模拟 regime 数正文=输出 | ✅ | 7 regimes / 19,600，脚本从 `extended_simulation_summary.tsv` 算出核对 |
| 补充表列名=正文词汇 | ✅ | `tools/si_names.py` 单一映射，导出时翻译；表头带旧词即导出失败 |
| 队列登记表（Table S6）完整 | ✅ | 由 `build_registry_table.py` 从 registry.yaml + 真实 `donor_labels.tsv` 生成；7 个数据集、字段无空值、发布供者数与实际建模数分列 |
| 图 / 表编号与首次引用顺序一致 | ✅ | 正文 Fig 1–8、Fig S1–S7、Table S1–S9 全部按引用顺序重排；脚本校验引用顺序 = 声明顺序 |
| 每个声明的图表在正文至少被引用一次 | ✅ | 脚本校验；本轮补上了此前从未引用的 Table S1/S6/S7/S9（旧编号）|
| AI 声明具体到模型、环节、核验方式 | ✅ | §9.18 |
| **行号** | ❌ | 需在最终 PDF 生成环节加 |
| 通讯作者确认 Cover Letter 与署名分工 | ❌ | 待人工 |
