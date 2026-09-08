# Figures 5-8 and S4: the two lupus cohorts.
# Figure numbers live in the manuscript and in tools/build_html.py, never here.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(dplyr);library(tidyr)})
# Project root: the directory holding figures/, sim/, data/. Resolved from this
# script's own location, or from RHEUMLENS_ROOT, so the script runs anywhere.
W <- Sys.getenv("RHEUMLENS_ROOT", unset = NA)
if (is.na(W)) {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  W <- if (length(f)) normalizePath(file.path(dirname(f), "..", ".."))
       else normalizePath(getwd())
}
setwd(W); L <- "data/legacy/"

base <- theme_classic(base_size=8) +
  theme(axis.text=element_text(colour="grey20",size=6.6), axis.title=element_text(size=7.5),
        plot.title=element_text(size=8,face="bold",margin=margin(b=1)),
        plot.subtitle=element_text(size=6.3,colour="grey35",margin=margin(b=4)),
        plot.title.position="plot", legend.key.size=unit(6,"pt"),
        legend.text=element_text(size=6.4), legend.title=element_text(size=6.8),
        legend.margin=margin(0,0,0,0), strip.background=element_blank(),
        strip.text=element_text(size=6.8, face="bold"))
CASE <- "#C4633E"; CTRL <- "#3C6E9F"
repl <- c(geneformer="Geneformer", hvg_pseudobulk="HVG pseudobulk",
          pca_pseudobulk="PCA pseudobulk")
rpal <- c(Geneformer="#7B6BA8", `HVG pseudobulk`="#3C6E9F", `PCA pseudobulk`="#4E8C6E")

## ======================== Fig_lupus_design =================================
wv <- read.delim(paste0(L,"raw_h5ad_wave_by_disease_cd4_only.tsv"), check.names=FALSE)
names(wv) <- c("level","control","case"); wv$panel <- "GSE174188 · processing wave"
wv$level <- paste0("W", as.integer(wv$level))
bb <- read.delim(paste0(L,"gse135779_batch_by_label.tsv"))
bb <- data.frame(level=bb$batch, control=bb$control, case=bb$case,
                 panel="GSE135779 · sequencing batch")
yy <- read.delim(paste0(L,"gse135779_year_by_label.tsv"))
yy <- data.frame(level=as.character(yy$collection_year), control=yy$control,
                 case=yy$case, panel="GSE135779 · collection year")
comp <- bind_rows(wv, bb, yy) |>
  filter(!grepl("^(Total|All|total)$", level)) |>
  mutate(panel=factor(panel, levels=c("GSE174188 · processing wave",
                                      "GSE135779 · sequencing batch",
                                      "GSE135779 · collection year")),
         pure = (control==0 | case==0)) |>
  pivot_longer(c(case, control), names_to="class", values_to="n")

p5a <- ggplot(comp, aes(level, n, fill=class)) +
  geom_col(width=.72) +
  geom_point(data=comp |> group_by(panel) |>
               mutate(off=-0.045*max(tapply(n, level, sum))) |> ungroup() |>
               distinct(panel, level, pure, off) |> filter(pure),
             aes(level, y=off), inherit.aes=FALSE, shape=17, size=1.2, colour="#8C3B26") +
  facet_wrap(~panel, scales="free", nrow=1) +
  scale_fill_manual(values=c(case=CASE, control=CTRL), name=NULL) +
  labs(x=NULL, y="donors", title="A   Cases and controls are not spread evenly over the design",
       subtitle="▲ marks a stratum containing only one kind of donor") +
  base + theme(legend.position="bottom")

s <- read.delim(paste0(L,"summary.tsv"))
cov <- s |> filter(adjustment=="covariates_only") |>
  mutate(block=recode(method, covariates_all="all", covariates_batch="batch",
                      covariates_collection_year="collection year",
                      covariates_demographic="demographic", covariates_qc="sample quality"),
         cohort=recode(cohort, GSE174188_CD4="GSE174188 CD4"))
cov$block <- factor(cov$block, levels=c("sample quality","demographic","collection year",
                                        "batch","all"))
# GSE174188 has no recorded collection year; insert an explicit gap so the line breaks
# there rather than implying a value.
cov <- bind_rows(cov, tidyr::tibble(cohort="GSE174188 CD4", block=factor("collection year",
                 levels=levels(cov$block)), auc_mean=NA_real_, auc_p025=NA_real_,
                 auc_p975=NA_real_))
p5b <- ggplot(cov, aes(block, auc_mean, colour=cohort, group=cohort)) +
  geom_hline(yintercept=.5, linetype="22", colour="grey60", linewidth=.3) +
  geom_errorbar(aes(ymin=auc_p025, ymax=auc_p975), width=.12, linewidth=.4) +
  geom_line(linewidth=.5) + geom_point(size=2) +
  scale_colour_manual(values=c(GSE135779="#B4913C", `GSE174188 CD4`="#3C6E9F"), name=NULL) +
  scale_y_continuous(limits=c(.4,1)) +
  labs(x=NULL, y="design-only AUC",
       title="B   Recorded metadata alone predicts the diagnosis in both cohorts",
       subtitle="bars are 2.5th-97.5th percentiles over 20 repeated donor splits. GSE174188 has no recorded collection year (gap). In GSE135779 batch alone is at chance:\nthe label is spread across year, sample quality and demographics instead") +
  base + theme(axis.text.x=element_text(angle=20, hjust=1), legend.position=c(.16,.87),
               legend.background=element_blank())

F5 <- p5a / p5b + plot_layout(heights=c(1,1.1))
ggsave("figures/out/Fig_lupus_design.svg", F5, width=180, height=132, units="mm",
       device=svglite::svglite)
ragg::agg_png("figures/out/Fig_lupus_design.png", width=180, height=132, units="mm", res=400)
print(F5); invisible(dev.off())

## ======================== Fig_exposure_vs_use ==============================
bp <- read.delim(paste0(L,"representation_batch_predictability_summary.tsv"))
best <- bp |> group_by(cohort, method) |> slice_max(auc_mean, n=1) |> ungroup() |>
  transmute(cohort, representation=repl[method], batch_auc=auc_mean, batch=batch)
dis <- s |> filter(adjustment=="unadjusted") |>
  transmute(cohort, representation=repl[method], disease_auc=auc_mean)
m6 <- inner_join(best, dis, by=c("cohort","representation")) |>
  mutate(cohort=recode(cohort, GSE174188_CD4="GSE174188 CD4"))
m6l <- m6 |> pivot_longer(c(disease_auc, batch_auc))
m6l$name <- factor(m6l$name, levels=c("disease_auc","batch_auc"),
                   labels=c("diagnosis","best single batch"))
p6a <- ggplot(m6l, aes(value, representation, colour=representation, shape=name)) +
  geom_line(aes(group=representation), colour="grey82", linewidth=.5) +
  geom_point(size=2.2, stroke=.6, fill="white") +
  facet_wrap(~cohort, nrow=1) +
  scale_shape_manual(values=c(diagnosis=16, `best single batch`=21), name=NULL) +
  scale_colour_manual(values=rpal, guide="none") +
  scale_x_continuous(limits=c(.8,1.02)) +
  labs(x="cross-validated AUC", y=NULL,
       title="A   The representations recognise the batch as well as the diagnosis",
       subtitle="one-versus-rest AUC for the most predictable recorded batch") +
  base + theme(legend.position="bottom")

p6b <- ggplot(data.frame(x=c(0.9930,0.4992),
                         y=factor(c("representation → batch","batch → diagnosis"),
                                  levels=c("batch → diagnosis","representation → batch"))),
              aes(x,y)) +
  geom_vline(xintercept=.5, linetype="22", colour="grey60", linewidth=.3) +
  geom_point(size=2.6, colour="#8C3B26") +
  geom_text(aes(label=sprintf("%.3f", x)), nudge_y=.22, size=2.3, colour="grey25") +
  scale_x_continuous(limits=c(.42,1.04)) +
  labs(x="AUC", y=NULL,
       title="B   The batch is legible but uninformative",
       subtitle="GSE135779: the batch is recovered almost perfectly,\nbut says nothing about the diagnosis") +
  base

delta <- s |> filter(cohort=="GSE135779", adjustment %in% c("unadjusted","residual_batch")) |>
  select(method, adjustment, auc_mean) |>
  pivot_wider(names_from=adjustment, values_from=auc_mean) |>
  mutate(representation=repl[method], d=residual_batch-unadjusted)
p6c <- ggplot(delta, aes(d, representation, colour=representation)) +
  geom_vline(xintercept=0, colour="grey45", linewidth=.35) +
  geom_segment(aes(x=0, xend=d, yend=representation), linewidth=.5) +
  geom_point(size=2.2) +
  geom_text(aes(label=sprintf("%+.3f", d)), nudge_x=ifelse(delta$d<0, -0.006, 0.006),
            hjust=ifelse(delta$d<0, 1, 0), size=2.3, colour="grey25") +
  scale_colour_manual(values=rpal, guide="none") +
  scale_x_continuous(limits=c(-.03,.03)) +
  scale_y_discrete(limits=rev) +
  labs(x="change in AUC for diagnosis after removing batch", y=NULL,
       title="C   Removing it changes nothing",
       subtitle="all three shifts are smaller than 0.01\n") +
  base

F6 <- p6a / (p6b | p6c) + plot_layout(heights=c(1,.8))
ggsave("figures/out/Fig_exposure_vs_use.svg", F6, width=180, height=112, units="mm",
       device=svglite::svglite)
ragg::agg_png("figures/out/Fig_exposure_vs_use.png", width=180, height=112, units="mm", res=400)
print(F6); invisible(dev.off())
cat("Fig_lupus_design, Fig_exposure_vs_use ok\n")
