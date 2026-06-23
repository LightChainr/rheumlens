# RheumLens local data index
Generated: 2026-06-23 Asia/Shanghai
Root: `/Volumes/Mac Data/Research/RheumLens_20260622/server_snapshot/rheumlens`
Critical tar: `/Volumes/Mac Data/Research/RheumLens_20260622/relay_download/rheumlens_critical_20260622.tar`
Inventory CSV: `/Volumes/Mac Data/Research/RheumLens_20260622/data_index/RHEUMLENS_LOCAL_DATA_INVENTORY.csv`
Key assets CSV: `/Volumes/Mac Data/Research/RheumLens_20260622/data_index/RHEUMLENS_KEY_ASSETS.csv`

## Snapshot summary
- Total files: 1506
- Indexed data/report files: 1381
- Total snapshot size: 2.3GB

## Top-level directories
- `data/`: 1.5GB; indexed files=4
- `embeddings/`: 726.0MB; indexed files=4
- `results/`: 103.6MB; indexed files=1321
- `splits/`: 197.2KB; indexed files=40
- `configs/`: 34.0KB; indexed files=11
- `metadata/`: 5.1KB; indexed files=1

## Key assets
| Asset | Exists | Size | SHA256 |
|---|---:|---:|---|
| `embeddings/scgpt/GSE135779_matched500.npz` | True | 20.0MB | `00e4427ed0f3a0c361c4a8722934a4fadb6c5435b66e5952db71f5039162d31f` |
| `embeddings/scgpt/GSE174188_CD4_v1/GSE174188_CD4_scgpt.npz` | True | 442.0MB | `68eea032a2d5da4379b260b293da9315ffc9878ac91f30dc5dd487c0c7c64a9a` |
| `embeddings/scgpt/GSE285773_v1/GSE285773_scgpt.npz` | True | 238.2MB | `6d9203333920235449165b8408e7a4ad3a669f494ed08c4517c99d85bc978471` |
| `embeddings/geneformer/v2_20260620T0915/GSE135779_matched500.npz` | True | 25.8MB | `87276833b0e7aecded2692915c66c2582ffabd22ba044200c23923845a482163` |
| `data/processed/GSE135779_matched500/lognorm.npz` | True | 51.2MB | `ed4f5423f4bef6afec62818db03ad8be85526301b9b14851af29a0876eca196c` |
| `data/processed/GSE135779_matched500/raw_counts.npz` | True | 49.4MB | `19c23e18040bc3cc72b4de81f0b58be5184e84f3ddc4e816017863ab5d08fb69` |
| `data/processed/GSE174188_CD4/lognorm.npz` | True | 538.3MB | `fb236c1f4169db49645c33010dc815cc2b3e171bbaf08595c1ff07f4079a9b46` |
| `data/processed/GSE285773/lognorm.npz` | True | 903.7MB | `aac6edc2e206c916435360bde8314993038a6d70d77255e29d19addad7033c90` |
| `splits/authoritative_primary/GSE135779.csv` | True | 4.1KB | `bf30ae4eb8807ae3e2b8eb4fc98958dbae548ac669b4cef8fce7efe730e9ed35` |
| `splits/authoritative_primary/GSE174188_CD4.csv` | True | 21.9KB | `d06b67ccf47d43e42b540aa1a3a134046940996fb9679b2474dd7828d996c869` |
| `splits/authoritative_primary/GSE285773.amendment_v1.csv` | True | 2.0KB | `90fe8a0cb471725f2e3c17193277ece83f9d51ca58a1da6953bfc64f45da775b` |
| `results/P5_matched500_v2_20260620T113549Z/final_tables/method_summary.csv` | True | 3.5KB | `26190064666463c8c3e57f1d4938ad10b76e665c279ff51395278390bf3678d1` |
| `results/P6_GSE174188_v1/final_tables/method_summary.csv` | True | 3.0KB | `4040bcf8791e9c14c4e6b7ebd4ea68536b670c1d09f4a433a4c29bb7f650fe89` |
| `results/P8_6_covariate_sensitivity_v3/COMPLETION_REPORT_V3_1.md` | True | 2.8KB | `f7032bca7d27d8fec813d6d4a9bcebec79106024e8fe68eed14719294f87e66b` |
| `results/P8_7_kernel_diagnostics/P8_7_COMPLETION_REPORT.md` | True | 3.1KB | `a75327cc5b1cf0e1ae820cac6bcf01f95968dc50b77f84864b3745ae300a4d75` |

## Practical interpretation
- 本地已有 P9/P10 收尾所需的核心数据：三队列 splits、scGPT/Geneformer embedding、lognorm/raw_counts 子集、P5/P6/P8 结果表与审计报告。
- 大体积原始全量数据没有全部下到 Mac；这份 snapshot 是 critical package，不是完整服务器镜像。
- 后续正式分析应优先使用 `server_snapshot/rheumlens` 作为只读输入，把新结果写入 `mac_workspace` 或新的 `results_mac_final`。
