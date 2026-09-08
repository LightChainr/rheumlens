# Figures 7, 8 and S4.
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
        strip.text=element_text(size=6.8,face="bold"))
repl <- c(geneformer="Geneformer", hvg_pseudobulk="HVG pseudobulk", pca_pseudobulk="PCA pseudobulk")
rpal <- c(Geneformer="#7B6BA8", `HVG pseudobulk`="#3C6E9F", `PCA pseudobulk`="#4E8C6E")

## ======================== Fig_adjustment ===================================
rs <- read.delim(paste0(L,"residualisation_sensitivity_summary.tsv")) |>
  mutate(representation=repl[representation],
         adj=recode(adjustment, unadjusted="no adjustment",
                    batch_location_adjustment="batch means removed",
                    overlap_weighted_full_design="overlap weighting",
                    random_forest_full_design="random-forest residuals",
                    ridge_full_design="ridge residuals"))
rs$adj <- factor(rs$adj, levels=c("no adjustment","batch means removed","overlap weighting",
                                  "random-forest residuals","ridge residuals"))
p7a <- ggplot(rs, aes(mean, adj, colour=representation)) +
  geom_vline(xintercept=.5, linetype="22", colour="grey60", linewidth=.3) +
  geom_errorbarh(aes(xmin=p025, xmax=p975), height=.16, linewidth=.4,
                 position=position_dodge(width=.55)) +
  geom_point(size=1.9, position=position_dodge(width=.55)) +
  scale_colour_manual(values=rpal, name=NULL) +
  scale_x_continuous(limits=c(.38,1)) +
  scale_y_discrete(limits=rev) +
  labs(x="AUC for diagnosis", y=NULL,
       title="A   How much the answer changes depends on how you remove the metadata",
       subtitle="GSE135779. The same cohort scores 0.94 or 0.52 depending only on which adjustment is chosen") +
  base + theme(legend.position="bottom")

pw <- read.delim(paste0(L,"pure_wave_stratum.tsv")) |>
  mutate(representation=repl[representation],
         grp=ifelse(block=="C2_pure_wave","wave 4 only\n(n = 66)",
                    "random subset\nmatched for size (n = 66)"))
p7b <- ggplot(pw, aes(auc_mean, grp, colour=representation)) +
  geom_errorbarh(aes(xmin=auc_p025, xmax=auc_p975), height=.14, linewidth=.4,
                 position=position_dodge(width=.5)) +
  geom_point(size=1.9, position=position_dodge(width=.5)) +
  scale_colour_manual(values=rpal, guide="none") +
  scale_x_continuous(limits=c(.82,1)) +
  labs(x="AUC for diagnosis", y=NULL,
       title="B   Restricting to one processing wave costs much less",
       subtitle="GSE174188 CD4. The signal survives inside a single wave,\nso it is not only the wave being read") +
  base

# Read the CMV numbers from the walkthrough output rather than hardcoding them:
# the design matrix width changed once already and a pasted number went stale.
wt  <- read.delim("walkthrough/results/decision_tree_walkthrough.tsv")
cmv <- wt[wt$cohort == "CMV_HIHA", ]
stopifnot(nrow(cmv) == 1)
nc <- data.frame(step=factor(c("before\nadjustment","after ridge\nresiduals"),
                             levels=c("before\nadjustment","after ridge\nresiduals")),
                 auc=c(cmv$disease_auc, cmv$residualised_auc))
p7c <- ggplot(nc, aes(step, auc)) +
  geom_hline(yintercept=.5, linetype="22", colour="grey60", linewidth=.3) +
  geom_line(aes(group=1), linewidth=.5, colour="#8C3B26") +
  geom_point(size=2.4, colour="#8C3B26") +
  geom_text(aes(label=sprintf("%.3f", auc), vjust=ifelse(auc > .7, -1.1, 1.9)),
            size=2.3, colour="grey25") +
  annotate("segment", x=1.5, xend=1.5, y=cmv$disease_auc, yend=cmv$residualised_auc,
           linewidth=.3, colour="grey45", arrow=arrow(length=unit(2,"pt"), ends="both")) +
  annotate("text", x=1.55, y=(cmv$disease_auc+cmv$residualised_auc)/2,
           label=sprintf("\u2212%.3f", cmv$residualisation_loss),
           hjust=0, size=2.4, colour="grey25") +
  scale_y_continuous(limits=c(.45,.95)) +
  scale_x_discrete(expand=expansion(add=.55)) +
  labs(x=NULL, y="AUC for diagnosis",
       title="C   The same adjustment breaks an unconfounded cohort",
       subtitle=sprintf(paste0("CMV: design-only AUC is %.3f, so there is nothing to remove.\n",
                               "%d design columns for %d donors take the diagnosis with them."),
                        cmv$design_only_auc, cmv$n_design_col, cmv$n_donor)) +
  base
F7 <- (p7a / (p7b | p7c + plot_layout(widths=c(1,1)))) + plot_layout(heights=c(1,.82))
ggsave("figures/out/Fig_adjustment.svg", F7, width=180, height=126, units="mm", device=svglite::svglite)
ragg::agg_png("figures/out/Fig_adjustment.png", width=180, height=126, units="mm", res=400)
print(F7); invisible(dev.off())

## ======================== Fig_transfer =====================================
dt <- read.delim(paste0(L,"deconfounded_transfer.tsv")) |>
  mutate(representation=repl[representation],
         src=recode(source_set, full="all 261 donors", wave4="wave 4 only (n = 89)",
                    random_n89="random n = 89"))
dt$src <- factor(dt$src, levels=c("all 261 donors","random n = 89","wave 4 only (n = 89)"))
p8a <- ggplot(dt, aes(target_auc, src, colour=representation)) +
  geom_vline(xintercept=.5, linetype="22", colour="grey60", linewidth=.3) +
  geom_errorbarh(aes(xmin=ci_low, xmax=ci_high), height=.14, linewidth=.4,
                 position=position_dodge(width=.55)) +
  geom_point(size=1.9, position=position_dodge(width=.55)) +
  scale_colour_manual(values=rpal, name=NULL) +
  scale_x_continuous(limits=c(.55,1.02)) + scale_y_discrete(limits=rev) +
  labs(x="AUC in GSE285773 (26 donors), trained on GSE174188", y=NULL,
       title="Training inside one wave does not improve transfer to a new cohort",
       subtitle="the wave-4 model is no better than a random subset of the same size, and worse than using all donors. Bars are DeLong intervals, except for the random subset, where they are the 2.5th-97.5th percentiles over 20 draws") +
  base + theme(legend.position="bottom")

F8 <- p8a + labs(title="Training inside one wave does not improve transfer to a new cohort")
ggsave("figures/out/Fig_transfer.svg", F8, width=180, height=76, units="mm", device=svglite::svglite)
ragg::agg_png("figures/out/Fig_transfer.png", width=180, height=76, units="mm", res=400)
print(F8); invisible(dev.off())

## ================ Fig_S5_celltype_composition =============================
cs <- read.delim(paste0(L,"composition_summary.tsv")) |>
  mutate(v=recode(variant, raw_proportions="cell-type proportions",
                  raw_plus_log_cells="proportions + cell count",
                  clr_composition="CLR-transformed", clr_plus_log_cells="CLR + cell count"),
         a=recode(adjustment, unadjusted="no adjustment",
                  residual_processing_wave="wave removed"))
s4a <- ggplot(cs, aes(auc_mean, v, colour=a)) +
  geom_vline(xintercept=.5, linetype="22", colour="grey60", linewidth=.3) +
  geom_errorbarh(aes(xmin=auc_p025, xmax=auc_p975), height=.14, linewidth=.4,
                 position=position_dodge(width=.5)) +
  geom_point(size=1.9, position=position_dodge(width=.5)) +
  scale_colour_manual(values=c(`no adjustment`="#3C6E9F", `wave removed`="#C4633E"), name=NULL) +
  scale_x_continuous(limits=c(.45,1)) +
  labs(x="AUC for diagnosis", y=NULL,
       title="A   Cell-type proportions alone separate cases from controls",
       subtitle="GSE174188. Removing the processing wave costs about 0.20 AUC in every encoding") +
  base + theme(legend.position="bottom")

cr <- read.delim(paste0(L,"composition_restriction_metrics.tsv")) |>
  mutate(v=recode(variant, raw_proportions="proportions", raw_plus_log_cells="prop + count",
                  clr_composition="CLR", clr_plus_log_cells="CLR + count"),
         st=recode(stratum, pure_wave_4="wave 4 only", dominant_wave_4="wave 4 dominant"),
         ct=recode(control_type, observed="observed stratum",
                   size_label_matched="size- and label-matched random"))
s4b <- ggplot(cr, aes(roc_auc, v, colour=ct)) +
  geom_boxplot(outlier.size=.4, linewidth=.32, width=.55,
               position=position_dodge(width=.7)) +
  facet_wrap(~st, nrow=1) +
  scale_colour_manual(values=c(`observed stratum`="#C4633E",
                               `size- and label-matched random`="#3C6E9F"), name=NULL) +
  labs(x="AUC for diagnosis", y=NULL,
       title="B   Inside one wave the composition signal is still there",
       subtitle="20 repeated splits per box. Restricting to a wave costs 0.03-0.04 AUC on average — far less than the 0.20 that removing the wave costs in panel A") +
  base + theme(legend.position="bottom")
S4 <- s4a / s4b + plot_layout(heights=c(.85,1))
ggsave("figures/out/Fig_S5_celltype_composition.svg", S4, width=180, height=115, units="mm", device=svglite::svglite)
ragg::agg_png("figures/out/Fig_S5_celltype_composition.png", width=180, height=115, units="mm", res=400)
print(S4); invisible(dev.off())
cat("Fig_adjustment, Fig_transfer, Fig_S5_celltype_composition ok\n")
