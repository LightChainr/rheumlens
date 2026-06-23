# RheumLens 一小时关键包下载任务

从 Supervisor 获取 `URL`、`TOKEN` 和 `EXPIRES`。立即执行：

```bash
mkdir -p "$HOME/rheumlens_incoming"
cd "$HOME/rheumlens_incoming"
curl --fail --location --retry 20 --retry-delay 5 --continue-at - \
  -H "Authorization: Bearer $TOKEN" \
  "$URL/rheumlens_critical_20260622.tar" \
  -o rheumlens_critical_20260622.tar
curl --fail --location --retry 10 \
  -H "Authorization: Bearer $TOKEN" \
  "$URL/SHA256SUMS.txt" -o SHA256SUMS.txt
sha256sum -c SHA256SUMS.txt
```

如果一小时内未完成，不删除部分 tar；Supervisor 再开一个新的一小时链接后，使用同一条
`curl --continue-at -` 命令续传。SHA256 通过前不得解包。

校验通过后：

```bash
mkdir -p rheumlens_critical
tar -xf rheumlens_critical_20260622.tar -C rheumlens_critical
```

自检：关键 NPZ 存在且可读；GSE174188/GSE285773 的 embedding 与 lognorm donor 集合一致；folds 中
`role` 只表示 train/test，不是临床标签。汇报 tar 字节数、SHA256、解包目录及自检结果。不要启动分析。
