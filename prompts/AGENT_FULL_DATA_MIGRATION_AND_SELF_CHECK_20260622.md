# RheumLens 全量数据迁移与自检任务

你负责把服务器上的完整 RheumLens 数据复制到你自己的持久工作区；源服务器只读，禁止删除、移动或覆盖。

## 源端

- SSH：使用用户提供的当前连接信息，不把密码写入脚本、日志或报告。
- 主项目：`/autodl-fs/data/rheumlens`（约 54.1 GB）
- Geneformer 历史环境/材料：`/autodl-fs/data/rheumlens_20260619`（约 10.8 GB）
- 两棵目录都必须复制，保留符号链接、mtime、权限与目录结构。

## 执行规范

1. 先检查目标持久盘剩余空间，要求至少 100 GB；确认目标不是临时目录。
2. 使用单个 rsync 进程顺序复制，不并发扫描共享盘。推荐：
   `rsync -a --partial --append --inplace --timeout=180 -e 'ssh ...' SOURCE DEST`。
3. SSH 中断后自动等待 30 秒并续传；不得因一次断线重新从零下载。
4. 复制期间每 10 分钟最多检查一次：rsync PID、目标字节增长、磁盘余量、D-state 数；禁止高频轮询。
5. 首轮完成后执行一次 `rsync -aRcni` checksum dry-run。输出必须为空；如有差异，用
   `rsync -aRc --inplace` 修复后再次 dry-run。
6. 不执行 GPU 训练、embedding 提取或大规模解压，直到迁移验收通过。

## 强制自检

分别对两棵目录记录：

- 源端与目标端普通文件数；
- 源端与目标端普通文件总字节数；
- 符号链接数及断链列表；
- checksum dry-run 差异数；
- 可用磁盘空间；
- 关键资产存在性、大小和 SHA256。

关键资产至少包括：

- GSE174188_CD4、GSE285773 的 lognorm 与 scGPT NPZ；
- GSE135779 matched-500 scGPT 与已验收 Geneformer NPZ；
- `splits/authoritative_primary/`；
- `configs/`；
- P4/P5/P6/P8 系列 `results/`；
- Geneformer checkpoint、token dictionary、gene median dictionary及 gene mapping。

再用只读 Python 检查每个关键 NPZ：可打开、矩阵 shape、donor 数、标签分布、NaN/Inf、cell-ID 唯一性；
检查 GSE174188/GSE285773 的 embedding、lognorm、folds donor 集合一致。不得把 folds 的 `role` 当临床标签。

## 验收输出

生成：

- `MIGRATION_REPORT.md`
- `SOURCE_TARGET_COUNTS.tsv`
- `CRITICAL_ASSET_SHA256.tsv`
- `NPZ_STRUCTURAL_AUDIT.tsv`
- `BROKEN_SYMLINKS.txt`
- `CHECKSUM_DRYRUN.txt`
- `MIGRATION_COMPLETE.json`

只有文件数、总字节数匹配且 checksum dry-run 为零差异时，状态才能写 `ACCEPTED`。否则写
`INCOMPLETE_RESUMABLE` 并继续断点续传。完成后不要自行启动科学分析，只汇报目标路径、总大小、验收状态和上述报告路径。
