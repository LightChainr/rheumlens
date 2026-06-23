#!/usr/bin/env Rscript

.libPaths(c(".Rlib", .libPaths()))
suppressPackageStartupMessages(library(edgeR))

root <- normalizePath(".")
count_path <- file.path(root, "data", "secondary_features", "GSE135779_project_marker_T_raw_count_sums.csv.gz")
public <- file.path(root, "results", "evidence_package", "public_metadata", "structured")
out <- file.path(root, "results", "evidence_package", "secondary_foldcontained")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

counts <- read.csv(gzfile(count_path), row.names = 1, check.names = FALSE)
clinical <- read.csv(file.path(public, "GSE135779_ST1b_donor_clinical.csv"), check.names = FALSE)
mapping <- read.csv(file.path(public, "GSE135779_study_name_to_donor_id.csv"), check.names = FALSE)
meta <- merge(clinical, mapping[, c("study_name", "donor_id")], by.x = "Names", by.y = "study_name")
meta <- meta[match(colnames(counts), meta$donor_id), ]
stopifnot(all(meta$donor_id == colnames(counts)))
group <- factor(ifelse(meta$Groups == "cSLE", "SLE", "HD"), levels = c("HD", "SLE"))

fit_model <- function(design, coef_name, label) {
  y <- DGEList(counts = counts, group = group)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = which(colnames(design) == coef_name))
  tab <- topTags(test, n = Inf, sort.by = "PValue")$table
  tab$gene <- rownames(tab)
  tab$model <- label
  write.csv(tab, file.path(out, paste0("GSE135779_project_marker_T_DE_", label, ".csv")), row.names = FALSE)
  data.frame(model = label, tested_genes = nrow(tab), fdr_lt_005 = sum(tab$FDR < .05),
             up_fdr_lt_005 = sum(tab$FDR < .05 & tab$logFC > 0),
             down_fdr_lt_005 = sum(tab$FDR < .05 & tab$logFC < 0))
}

design_unadjusted <- model.matrix(~ group)
results <- list(fit_model(design_unadjusted, "groupSLE", "disease_only"))

meta$Batch <- factor(meta$Batch)
meta$Gender <- factor(meta$Gender)
meta$Age <- as.numeric(meta$Age)
design_adjusted <- model.matrix(~ Batch + Age + Gender + group, data = meta)
if (qr(design_adjusted)$rank == ncol(design_adjusted)) {
  results[[2]] <- fit_model(design_adjusted, "groupSLE", "batch_age_sex_adjusted")
} else {
  results[[2]] <- data.frame(model = "batch_age_sex_adjusted", tested_genes = NA,
                             fdr_lt_005 = NA, up_fdr_lt_005 = NA, down_fdr_lt_005 = NA)
  writeLines("Adjusted design is rank deficient; result not estimated.",
             file.path(out, "GSE135779_project_marker_T_DE_adjusted_design_error.txt"))
}

summary <- do.call(rbind, results)
summary$annotation <- "project_fixed_marker_score_not_author_cell_labels"
write.csv(summary, file.path(out, "GSE135779_project_marker_T_DE_summary.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "GSE135779_project_marker_T_DE_R_sessionInfo.txt"))
print(summary)
