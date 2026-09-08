# Figures S2, S6 and S7.
#   S5  residualisation loss when the design is independent of the diagnosis
#   S6  type-I error of the collection-preserving permutation with a within-stratum
#       confounder and no biological signal
#   S7  learned pooling, moved out of the main text
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(dplyr);library(tidyr)})
W <- Sys.getenv("RHEUMLENS_ROOT", unset = NA)
if (is.na(W)) {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  W <- if (length(f)) normalizePath(file.path(dirname(f), "..", ".."))
       else normalizePath(getwd())
}
setwd(W)
base <- theme_classic(base_size=8) +
  theme(axis.text=element_text(colour="grey20",size=6.6), axis.title=element_text(size=7.5),
        plot.title=element_text(size=8,face="bold",margin=margin(b=1)),
        plot.subtitle=element_text(size=6.3,colour="grey35",margin=margin(b=4)),
        plot.title.position="plot", legend.key.size=unit(6,"pt"),
        legend.text=element_text(size=6.4), legend.title=element_text(size=6.8),
        legend.margin=margin(0,0,0,0), strip.background=element_blank(),
        strip.text=element_text(size=6.8,face="bold"))

## ---------------------------------------------------------------- Figure S2
b <- read.delim("sim/results/calibration_arm_b_raw.tsv")
bal <- b |> filter(layout=="balanced")
rag <- b |> filter(layout=="cmv_ragged")
sumb <- bal |> group_by(n_design_col) |>
  summarise(loss=mean(residualisation_loss),
            lo=quantile(residualisation_loss,.025),
            hi=quantile(residualisation_loss,.975), .groups="drop")
ragl <- mean(rag$residualisation_loss)
ragci <- quantile(rag$residualisation_loss, c(.025,.975))

cmv_loss <- read.delim("walkthrough/results/decision_tree_walkthrough.tsv") |>
  subset(cohort == "CMV_HIHA") |> getElement("residualisation_loss")
stopifnot(length(cmv_loss) == 1)

s5 <- ggplot(sumb, aes(n_design_col, loss)) +
  geom_hline(yintercept=0, linetype="22", colour="grey60", linewidth=.3) +
  geom_ribbon(aes(ymin=lo, ymax=hi), fill="#3C6E9F", alpha=.15) +
  geom_line(colour="#3C6E9F", linewidth=.6) +
  geom_point(colour="#3C6E9F", size=1.6) +
  annotate("errorbar", x=mean(rag$n_design_col), ymin=ragci[1], ymax=ragci[2],
           width=2, colour="#B4913C", linewidth=.4) +
  annotate("point", x=mean(rag$n_design_col), y=ragl, colour="#B4913C", size=2.4,
           shape=17) +
  annotate("text", x=mean(rag$n_design_col)-3, y=ragl, hjust=1, size=2.3,
           colour="#8A6E2A", label=sprintf("CMV-shaped design\n%.3f", ragl)) +
  # read the observed loss rather than pasting it: it moved once already
  geom_hline(yintercept=cmv_loss, linetype="21", colour="#8C3B26", linewidth=.45) +
  annotate("text", x=2, y=cmv_loss, vjust=-0.6, hjust=0, size=2.4, colour="#8C3B26",
           label=sprintf("observed CMV loss %.3f", cmv_loss)) +
  labs(x="one-hot design columns (108 donors)", y="AUC lost to residualisation",
       title="Figure S2  Residualisation loss with no confounding present",
       subtitle=paste0("Design generated independently of the diagnosis; biological effect ",
                       "calibrated so unadjusted AUC matches the CMV cohort.\n",
                       "Ribbon is the 2.5th-97.5th percentile over 60 draws. ",
                       "The triangle reuses the real batch-by-pool level sizes.")) +
  base
ggsave("figures/out/Fig_S2_residualisation_width.svg", s5, width=180, height=88,
       units="mm", device=svglite::svglite)
ragg::agg_png("figures/out/Fig_S2_residualisation_width.png", width=180, height=88,
              units="mm", res=400); print(s5); invisible(dev.off())

## ---------------------------------------------------------------- Figure S1
a <- read.delim("sim/results/calibration_arm_a_raw.tsv")
rej <- a |> group_by(gamma) |>
  summarise(within=mean(within_stratum_abs_corr),
            free=mean(p_unstratified<=.05),
            coll=mean(p_collection_preserving<=.05, na.rm=TRUE),
            n=n(), .groups="drop") |>
  pivot_longer(c(free, coll), names_to="test", values_to="rate") |>
  mutate(test=recode(test, free="free permutation",
                     coll="collection-preserving permutation"))
# Wilson intervals: the claim at gamma = 0 is that the rate is *consistent with*
# 0.05, which needs an interval rather than a point.
wilson <- function(k, n, z=1.959964) {
  ph <- k/n; d <- 1 + z^2/n
  c((ph + z^2/(2*n) - z*sqrt(ph*(1-ph)/n + z^2/(4*n^2)))/d,
    (ph + z^2/(2*n) + z*sqrt(ph*(1-ph)/n + z^2/(4*n^2)))/d)
}
ci <- t(mapply(function(r, n) wilson(round(r*n), n), rej$rate, rej$n))
rej$lo <- ci[,1]; rej$hi <- ci[,2]
s6a <- ggplot(rej, aes(gamma, rate, colour=test)) +
  geom_hline(yintercept=.05, linetype="22", colour="grey55", linewidth=.35) +
  geom_line(linewidth=.6) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=.05, linewidth=.35) +
  geom_point(size=1.9) +
  annotate("text", x=max(rej$gamma), y=.05, vjust=-0.7, hjust=1, size=2.3,
           colour="grey40", label="nominal 0.05") +
  scale_colour_manual(values=c(`free permutation`="#3C6E9F",
                               `collection-preserving permutation`="#C4633E"), name=NULL) +
  scale_y_continuous(limits=c(0,1)) +
  labs(x="strength of the within-stratum sample-quality → diagnosis association (γ)",
       y="fraction of runs with p ≤ 0.05",
       title="A   Neither test protects against a confounder that survives inside the strata",
       subtitle="No disease effect exists in the generating model, so every rejection here is a false positive") +
  base + theme(legend.position="bottom")
s6b <- ggplot(rej |> distinct(gamma, within), aes(gamma, within)) +
  geom_line(colour="grey35", linewidth=.5) + geom_point(size=1.7, colour="grey25") +
  labs(x="γ", y="within-stratum |corr|",
       title="B   What survives the restriction: the association the permutation never conditions on") +
  base + theme(plot.title=element_text(size=7.6))
S6 <- s6a / s6b + plot_layout(heights=c(1,.55)) +
  plot_annotation(title="Figure S1  Type-I error of the collection-preserving permutation",
                  theme=theme(plot.title=element_text(size=8.5, face="bold")))
ggsave("figures/out/Fig_S1_permutation_calibration.svg", S6, width=180, height=125,
       units="mm", device=svglite::svglite)
ragg::agg_png("figures/out/Fig_S1_permutation_calibration.png", width=180, height=125,
              units="mm", res=400); print(S6); invisible(dev.off())

## ---------------------------------------------------------------- Figure S6
pt <- read.delim("data/legacy/learned_pooling_paired_tests.tsv") |>
  mutate(dir=recode(direction, GSE174188_to_GSE285773="GSE174188 → GSE285773 (26 targets)",
                    GSE285773_to_GSE174188="GSE285773 → GSE174188 (261 targets)"),
         m=recode(method, deepsets="DeepSets", gatedmil="gated MIL", pma="PMA"),
         b=recode(baseline, mean_pool_1152d="vs mean pooling",
                  pca32_donor_mean="vs 32-PC donor mean", hvg_pseudobulk="vs HVG pseudobulk",
                  pca_pseudobulk="vs PCA pseudobulk"),
         sig=delong_bh_q < 0.05)
pt$b <- factor(pt$b, levels=c("vs mean pooling","vs 32-PC donor mean",
                              "vs HVG pseudobulk","vs PCA pseudobulk"))
s7 <- ggplot(pt, aes(auc_diff, b, colour=m, shape=sig)) +
  geom_vline(xintercept=0, colour="grey45", linewidth=.35) +
  geom_errorbarh(aes(xmin=diff_ci_low, xmax=diff_ci_high), height=.15, linewidth=.35,
                 position=position_dodge(width=.62)) +
  geom_point(size=1.8, fill="white", stroke=.6, position=position_dodge(width=.62)) +
  facet_wrap(~dir, nrow=1) +
  scale_shape_manual(values=c(`TRUE`=16, `FALSE`=21), guide="none") +
  scale_colour_manual(values=c(DeepSets="#7B6BA8", `gated MIL`="#B4913C", PMA="#4E8C6E"),
                      name=NULL) +
  scale_y_discrete(limits=rev) +
  labs(x="difference in AUC (attention pooling − baseline)", y=NULL,
       title="Figure S6  Learned pooling never beats averaging the cells",
       subtitle="Paired DeLong intervals. Filled points are significant after Benjamini-Hochberg adjustment; every one of them is negative.") +
  base + theme(legend.position="bottom")
ggsave("figures/out/Fig_S6_learned_pooling.svg", s7, width=180, height=88, units="mm",
       device=svglite::svglite)
ragg::agg_png("figures/out/Fig_S6_learned_pooling.png", width=180, height=88, units="mm",
              res=400); print(s7); invisible(dev.off())
cat("Fig_S1, Fig_S2, Fig_S6 ok\n")
