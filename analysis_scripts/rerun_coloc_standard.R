#!/usr/bin/env Rscript

.libPaths(c(".Rlib", .libPaths()))
suppressPackageStartupMessages(library(coloc))

root <- normalizePath(".")
input_dir <- file.path(root, "results", "evidence_package", "historical_recovery",
                       "rheumlens_sub_20260618_0152", "coloc")
out_dir <- file.path(root, "results", "evidence_package", "coloc_standard_rerun")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genes <- c("STAT4", "IFIH1", "PRDM1", "TET2", "BCL2", "FCGR2A")
prior_grid <- c(1e-6, 1e-5, 1e-4)
n_gwas <- 7219 + 15991
case_fraction <- 7219 / n_gwas
rows <- list()

for (gene in genes) {
  path <- file.path(input_dir, paste0(gene, "_merged.tsv"))
  d <- read.delim(path, check.names = FALSE)
  d <- d[!duplicated(d$rsid), , drop = FALSE]
  required <- c("rsid", "eqtl_beta", "eqtl_se", "eqtl_maf", "eqtl_an",
                "gwas_beta", "gwas_se", "gwas_p")
  if (nrow(d) == 0) {
    for (p12 in prior_grid) {
      rows[[length(rows) + 1]] <- data.frame(
        gene = gene, n_snps = 0, p1 = 1e-4, p2 = 1e-4, p12 = p12,
        PP0 = NA, PP1 = NA, PP2 = NA, PP3 = NA, PP4 = NA,
        status = "no_shared_variants"
      )
    }
    next
  }
  d <- d[complete.cases(d[, required]), , drop = FALSE]
  d <- d[d$eqtl_se > 0 & d$gwas_se > 0 & d$eqtl_maf > 0 & d$eqtl_maf < 1, , drop = FALSE]
  ne <- median(as.numeric(d$eqtl_an)) / 2
  gwas <- list(
    beta = as.numeric(d$gwas_beta),
    varbeta = as.numeric(d$gwas_se)^2,
    pvalues = as.numeric(d$gwas_p),
    snp = as.character(d$rsid),
    type = "cc", s = case_fraction, N = n_gwas
  )
  eqtl <- list(
    beta = as.numeric(d$eqtl_beta),
    varbeta = as.numeric(d$eqtl_se)^2,
    MAF = as.numeric(d$eqtl_maf),
    snp = as.character(d$rsid),
    type = "quant", N = ne
  )
  for (p12 in prior_grid) {
    fit <- coloc.abf(gwas, eqtl, p1 = 1e-4, p2 = 1e-4, p12 = p12)
    sm <- fit$summary
    rows[[length(rows) + 1]] <- data.frame(
      gene = gene, n_snps = nrow(d), p1 = 1e-4, p2 = 1e-4, p12 = p12,
      PP0 = unname(sm[["PP.H0.abf"]]), PP1 = unname(sm[["PP.H1.abf"]]),
      PP2 = unname(sm[["PP.H2.abf"]]), PP3 = unname(sm[["PP.H3.abf"]]),
      PP4 = unname(sm[["PP.H4.abf"]]), status = "standard_coloc_abf"
    )
  }
}

result <- do.call(rbind, rows)
write.csv(result, file.path(out_dir, "coloc_prior_sensitivity.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "R_sessionInfo.txt"))
writeLines(c(
  "Input limitation: merged tables contain GWAS alleles but no eQTL effect/other alleles.",
  "coloc.abf uses squared association evidence, so sign alignment does not alter ABFs,",
  "but allele harmonization and palindromic-SNP auditing cannot be reconstructed from these merged files.",
  "BCL2 and FCGR2A have zero shared rows and are not PP4=0 results."
), file.path(out_dir, "INPUT_LIMITATIONS.txt"))
print(result)
