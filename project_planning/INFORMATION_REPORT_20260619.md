# RheumLens 论文完整信息报告

**日期**: 2026-06-19  
**审计 Worker Agent**  
**来源**: 项目文件审计 + GEO 联网核验 + 历史对话记录 + Notion 文档 + A800 服务器执行记录  

---

## 一、三队列正式身份

### GSE135779 — 儿科发现队列

| 字段 | 值 | 来源 |
|------|-----|------|
| GEO accession | [GSE135779](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135779) | GEO |
| 论文 | **Nehar-Belaid et al., Nat Immunol 2020** | PMID 32747814 |
| 标题 | "Mapping systemic lupus erythematosus heterogeneity at the single-cell level" | PubMed |
| 第一作者 | Djamel Nehar-Belaid & Seunghee Hong | PubMed |
| 通讯作者 | Jacques F Banchereau, Virginia Pascual, Paul Robson | PubMed |
| 提交者 | William F Flynn, The Jackson Laboratory | GEO |
| 公开日期 | Sep 02, 2020 | GEO |
| 物种 | Homo sapiens | GEO |
| 平台 | **Illumina HiSeq 4000** (GPL20301) | GEO |
| 文库 | 10x Genomics single-cell RNA-seq | GEO |
| 细胞类型 | **全 PBMC**（非分选） | GEO |
| 样本数 | 56（33 cSLE + 11 cHD + 8 aSLE + 4 aHD，部分缺失） | GEO |
| 儿童队列分析用 | **33 cSLE + 11 cHD** = 44 donors | 项目约定 |
| 成人子集 | **7 aSLE + 5 aHD** = 12 donors（age 24-63） | 项目约定 |
| 原始数据 | dbGaP 受控访问 (phs002048.v1.p1) | GEO |
| ISG 方法 | **IFN 模块 M1.2/M3.4/M5.12 评分**, ~100 基因 | 原文 |
| 细胞数 (项目用) | 281,402 cells, 24,069 genes | H5AD |
| Cell-type annotation | **无** — 全 PBMC，无 per-cell 标签列 | 项目文件 |

### GSE285773 — 儿科复现队列

| 字段 | 值 | 来源 |
|------|-----|------|
| GEO accession | [GSE285773](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285773) | GEO |
| 论文 | **未发表**（"Citation missing"） | GEO |
| 标题 | "Single cell RNA profiling of blood CD4+ T cells identifies distinct helper and dysfunctional regulatory clusters in children with SLE" | GEO |
| 提交者 | Simone Caielli, Weill Cornell Medicine | GEO |
| 公开日期 | Sep 05, 2025 | GEO |
| 最后更新 | Sep 06, 2025 | GEO |
| 物种 | Homo sapiens | GEO |
| 平台 | **Illumina NovaSeq 6000** (GPL24676) | GEO |
| 细胞类型 | **FACS-sorted CD4+ T cells**（物理分选，非 computational subset） | GEO 标题 |
| 样本数 | **26**（10 HD + 16 SLE） | GEO |
| BioProject | PRJNA1205970 | GEO |
| 补充文件 | GSE285773_RAW.tar (1.7 GB) | GEO |
| 细胞数 (项目用) | 262,908 cells, 19,606 genes | H5AD |
| 年龄/性别/ancestry | **全部缺失** — archived H5AD 无这些字段 | donor_covariates.csv |
| batch/site/treatment | **100% 缺失** | donor_covariates.csv |
| Cell-type annotation | **无** — 全 CD4+，无 per-cell 亚型标签 | 项目文件 |
| ISG 基因可用性 | **不可用**（run_metadata: `ifn_genes_available: []`） | run_metadata.json |

### GSE174188 — 成人复现队列

| 字段 | 值 | 来源 |
|------|-----|------|
| GEO accession | [GSE174188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188) | GEO |
| 论文 | **Perez et al., Science 2022** | PMID 35389781 |
| 标题 | "Multiplexed scRNA-seq reveals the cellular and genetic correlates of systemic lupus erythematosus" | GEO |
| 第一作者 | Richard K Perez | GEO |
| 提交者 | Yang Sun, UCSF | GEO |
| 公开日期 | Dec 01, 2021 | GEO |
| 物种 | Homo sapiens | GEO |
| 平台 | **Illumina NovaSeq 6000** (GPL24676) | GEO |
| 文库 | multiplexed scRNA-seq (mux-seq) | GEO |
| 细胞类型 | **全 PBMC → 项目筛选为 CD4+ T cells** | CELLxGENE |
| 样本数 | 88 GSM (162 SLE + 99 HC = 261 donors) | GEO |
| 原始数据 | dbGaP 受控访问 | GEO |
| 项目数据源 | **CELLxGENE Census** (snapshot 2025-11-08), dataset_id 218acb0f | 项目文件 |
| CELLxGENE 文件 | SHA256 `7ba85edbdc033a9aeeecabe7a1cfcb0660282cc24b04f49fb52dd3825052c620` (12.2 GB) | AUDIT_AND_RESULTS.md |
| 细胞数 (项目用) | 380,275 CD4+ cells after QC | AUDIT_AND_RESULTS.md |
| Sex/Ethnicity/Processing cohort | **可用** | donor_covariates.csv |
| 年龄/sequencing lane/date/site | **缺失** | donor_covariates.csv |
| Cell-type annotation | **可用**（CELLxGENE `obs["cell_type"]`） | 项目文件 |
| ISG 基因可用性 | **不可用**（run_metadata: `ifn_genes_available: []`） | run_metadata.json |
| Disease state | **managed**(146) / **flare**(11) / **treated**(5) / **na**(99 HC) | donor_covariates.csv |
| 多样本 donor | 11 个（n_samples=2-3），FLARE 系列可能有重复采样 | donor_covariates.csv |
| 多 library | 大量 donor n_libraries=4-16，n_suspensions=1-3 | donor_covariates.csv |

---

## 二、平台交叉验证

| 队列 | 平台 | 文库 |
|------|------|------|
| GSE135779 | **HiSeq 4000** | 10x Genomics |
| GSE285773 | **NovaSeq 6000** | 10x Genomics |
| GSE174188 | **NovaSeq 6000** | mux-seq (10x) |

**结论**: GSE135779 (HiSeq) ≠ GSE285773/GSE174188 (NovaSeq)。**可以声称 cross-platform**（HiSeq vs NovaSeq），手稿中可注明为 "two sequencing platforms"。

---

## 三、权威主结果（fold-contained statistical audit, 2026-06-19）

**目录**: `/Users/lc/Documents/RheumLens/results/fold_contained/`  
**协议**: donor 为统计单位, training-fold-only HVG/PCA/scaling, 10,000 paired outcome-stratified donor bootstrap

| 队列 | scGPT AUC (95% CI) | Expr-PCA AUC (95% CI) | HVG-pseudobulk AUC (95% CI) | PCA − scGPT (95% CI) |
|------|-------------------|----------------------|---------------------------|---------------------|
| GSE135779 | 0.854 (0.733–0.950) | 0.948 (0.871–1.000) | 0.898 (0.755–0.995) | +0.094 (+0.028–+0.182) |
| GSE285773 | 0.950 (0.850–1.000) | 0.988 (0.944–1.000) | 0.975 (0.913–1.000) | +0.038 (−0.038–+0.144) |
| GSE174188 | 0.978 (0.963–0.990) | 0.981 (0.967–0.993) | 0.984 (0.970–0.995) | +0.003 (−0.008–+0.015) |

### IFN 残差（仅 GSE135779，training-fold-only）

- 15-gene donor ISG score 与 scGPT OOF 预测相关 (r=0.521, P=0.000289)
- SLE vs HC ISG score Cohen's d=1.258
- 残差 scGPT AUC: **0.760** (95% CI 0.601–0.895)
- 置换检验: **P=0.020**（1,000 次 donor-label permutation）→ **残差仍显著** ← 逆转旧结论

### 旧 GPU 服务器结果（已 DEPRECATED）

| 指标 | 旧值 | 新值 | 偏差 |
|------|------|------|------|
| GSE135779 scGPT | 0.890 | 0.854 | −0.036 |
| GSE285773 scGPT | 0.894 | 0.950 | +0.056 |
| GSE174188 scGPT | 0.852 | 0.978 | +0.126 |
| IFN 残差 | 0.609 | 0.760 | +0.151 |
| 残差 P | 0.184 (NS) | 0.020 (显著) | 结论反转 |

**旧值来源**: A800 GPU 服务器, pre-fold HVG/PCA（信息泄漏），单次 5-fold 无 bootstrap。**不可用于论文**。

---

## 四、15-gene ISG Panel 完整溯源

### 权威定义（来自 fold-contained audit）

```
ISG15, IFI6, MX1, OAS1, OAS2, OAS3, IFIT1, IFIT3,
IFI44, IFI44L, STAT1, RSAD2, IFITM1, IFITM3, HERC5
```

**首次记录**: `GSE135779_run_metadata.json:14-29` (2026-06-19)

### 历史演变

| 时间 | 来源 | 基因数 | 说明 |
|------|------|--------|------|
| 06-16 | Notion 早期交接 §9 | 7 | ISG15/IFI6/MX1/OAS1/IFIT1/IFIT3/STAT1 |
| 06-18 01:06 | GF 诊断(b) | 19 (混合) | ISG + GWAS 基因 + DE 基因混列表 |
| 06-19 | fold-contained audit | **15** | 首次作为纯 ISG panel 固定 |
| ~~早期~~ | "24-gene" | 24 | 不在任何当前文件中，已被手稿 v0.3 删除 |

### GSE135779 原文的 ISG 方法

Nehar-Belaid et al. (Nat Immunol 2020, PMID 32747814) **不使用基因列表评分**。原文使用：
- **IFN 模块 M1.2, M3.4, M5.12**（来自 Chaussabel/Banchereau 模块框架）计算模块评分
- 约 **100 个 IFN 相关基因**（来自 Gene Ontology + 上述模块）用于热图展示
- 蛋白水平仅验证了 **ISG15** 一个基因

**论文中的 15-gene panel 不是 GSE135779 原文方法的复现**，而是 RheumLens 项目的内部预设。方法部分需标注为 "curated ISG panel"，并引用原文模块来源。

---

## 五、二次结果溯源

### 已确认（来自 fold-contained 或独立核验）

| 结果 | 值 | 来源 | 状态 |
|------|-----|------|------|
| scGPT ARI vs cell type | 0.762 | Desktop/step0_audit.py | CONFIRMED |
| Geneformer ARI vs cell type | −0.007 | Desktop/gf_diag_c_fast.py | CONFIRMED |
| Coloc PP4 (STAT4/IFIH1/PRDM1/TET2) | 0.0 | coloc_final/a1_coloc_results.json | CONFIRMED |
| PRDM1 PP3 | 0.839 | coloc_final | CONFIRMED |
| Coloc nsnps (STAT4) | 3558 | coloc_final | CONFIRMED (GWAS p≈1e-65 正常) |
| eQTL 源 | QTD000479 = DICE naive CD4+ T (n=176) | coloc_run.py | CONFIRMED |

### 来自旧 GPU 运行，需 fold-contained 重跑

| 结果 | 值 | 手稿位置 | 状态 |
|------|-----|---------|------|
| T-cell-only AUC | 0.904 | v0.3 line 100 | NEEDS_RERUN |
| Monocyte-only AUC | 0.871 | v0.3 line 100 | NEEDS_RERUN |
| Cell-composition AUC | 0.796 | v0.3 line 100 | NEEDS_RERUN |
| 322 DE genes | 322 | v0.3 line 106 | NEEDS_RERUN |
| GSE135779→GSE285773 transfer (legacy) | 0.950 (硬编码) | p3_transfer.py (已删除) | UNSUPPORTED |

### 跨队列 transfer（从 evidence_package 重跑）

| 方向 | 缩放 | AUC | 说明 |
|------|------|-----|------|
| GSE135779 PBMC → GSE285773 CD4 | source_scaled | 0.631 | exploratory, unmatched cell type |
| GSE135779 PBMC → GSE285773 CD4 | unscaled | 0.938 | exploratory, preprocessing-sensitive |

**结论**: 无 matched CD4→CD4 transfer（GSE135779 无 cell_type 列，无法切 CD4）。**0.950 是 GSE285773 within-cohort OOF AUC，不是 transfer**。

---

## 六、缺失信息清单

### 元数据缺失

| 项目 | GSE135779 | GSE285773 | GSE174188 |
|------|-----------|-----------|-----------|
| Batch / site | ❌ | ❌ | ❌ (processing_cohort 作代理) |
| Age | ✅ (child/adult) | ❌ | ⚠️ (development_stage 有区间) |
| Sex | ❌ | ❌ | ✅ |
| Ancestry | ❌ | ❌ | ✅ |
| Treatment / drug | ❌ | ❌ | ⚠️ (disease_state: managed/flare/treated) |
| Sample / visit | ❌ | ❌ | ⚠️ (n_samples, 无 timepoint) |
| Library / lane | ❌ | ❌ | ⚠️ (n_libraries 4-16, 无 lane) |
| Cell-type annotation | ❌ | ❌ | ✅ |

### 统计缺失

| 项目 | 状态 |
|------|------|
| GSE174188 三个主模型置换检验 | **MISSING** |
| GSE285773 三个主模型置换检验 | **MISSING** |
| Repeated-CV sensitivity | **MISSING**（仅单一 5-fold） |
| Leave-processing-cohort-out | **MISSING** |
| GSE285773 IFN residual | **MISSING**（ISG 基因不可用） |
| GSE174188 IFN residual | **MISSING**（ISG 基因不可用） |
| T-cell/monocyte/composition 重跑 | **NEEDS_RERUN** |
| 322 DE genes 重跑 | **NEEDS_RERUN** |

### 可复现性缺失

| 项目 | 状态 |
|------|------|
| 环境锁文件 (conda-lock / pip freeze) | **MISSING** |
| 公共代码仓库 URL / DOI (GitHub/Zenodo) | **MISSING** |
| 引用审计 | **MISSING** |
| 作者信息 / 通讯作者 / 资金来源 | **MISSING** |

---

## 七、Fold 与模型开发历史

| 队列 | Fold 文件 | 创建方式 | 试过多套 Folds？ |
|------|----------|---------|----------------|
| GSE135779 | `gse135779_child_donor_folds.json` | 固定 5-fold 划分（早期创建，来源未记录） | **无证据** |
| GSE285773 | `gse285773_donor_folds.json` | `StratifiedKFold(5, shuffle=True, random_state=42)` | **单次生成** |
| GSE174188 | `gse174188_cd4_donor_folds.json` | `StratifiedKFold(5, shuffle=True, random_state=42)` | **单次生成** |

**证据**: 无 multiple fold set 记录，无 fold selection 依据 OOF AUC，**无 repeated-CV sensitivity**。fold-contained benchmark 脚本仅调用单次 StratifiedKFold。→ 建议在 Discussion/Limitations 中注明"folds 创建后固定未替换"，若审稿人要求可补 repeated-CV。

---

## 八、A1 Coloc 完整配置

| 参数 | 值 |
|------|-----|
| eQTL 数据 | QTD000479 (eQTL Catalogue), DICE naive CD4+ T cells |
| eQTL n | 176 |
| GWAS 数据 | Bentham 2015 SLE GWAS (GCST003156) |
| 基因数 | 6 (STAT4/IFIH1/PRDM1/TET2/BCL2/FCGR2A) |
| 窗口 | ±200kb per gene |
| 合并方式 | rsid matching |
| 方法 | coloc.abf (custom Python implementation) |
| R coloc 版本 | v5.2.3（定稿版） |
| 脚本 | `coloc/a1_coloc_run.py` + `coloc/run_coloc.R` + `coloc/merge_tables.py` |
| 输出 | `a1_coloc_results.json` + 6 × `*_merged.tsv` |
| 持久化 | `/autodl-fs/data/rheumlens/coloc_final/` |
| 结论 | PP4=0.0 所有基因，有信息量阴性 |

---

## 九、手稿状态

| 版本 | 文件 | 状态 |
|------|------|------|
| v0.3 (当前) | `RheumLens_manuscript_working_v0.3.md` | 已使用所有 fold-contained 修正值 |
| v0.3 清除的旧值 | 0.890/0.894/0.852/0.609/P=0.184/24-gene | **全部已删除** ✅ |
| v0.2 | docx | 旧值版本 |

---

## 十、Mac 本地备份清单

| 文件 | 大小 | SHA256 |
|------|------|--------|
| `rheumlens_sub_20260618_0152.tar.gz` | 2.1 GB | e08ad073... |
| `rheumlens_scgpt_emb_20260618_1107.tar.gz` | 5.9 GB | af63a398... |
| `rheumlens_figdata_full_20260618_1610.tar.gz` | 30 MB | 90751093... |
| `rheumlens_figdata_20260618_1547.tar.gz` | 18 MB | 6a2fd3af... (被 full 取代) |

---

## 十一、关键待修复声明

| # | 当前声明 | 应修正为 | 原因 |
|---|---------|---------|------|
| 1 | GSE285773 "computational CD4 subset" | **"FACS-sorted CD4+ T cells"** | GEO 标题明确注明 |
| 2 | 15-gene ISG 来自原文方法 | **"curated panel; original study used IFN modules M1.2/M3.4/M5.12"** | Nehar-Belaid 2020 使用模块评分 |
| 3 | GSE285773 transfer AUC 0.950 (旧版) | **删除** | 为 GSE285773 within-cohort OOF |
| 4 | "IFN removal eliminates discrimination" (旧版) | **"IFN explains substantial component but residual remains significant (P=0.020)"** | fold-contained 结果 |
| 5 | "no cross-platform claim" | **可声称 "two sequencing platforms (HiSeq 4000 / NovaSeq 6000)"** | 已验证 |
| 6 | GSE285773 unpublished | **引用 GEO accession GSE285773 + BioProject PRJNA1205970** | 无 PMID |
