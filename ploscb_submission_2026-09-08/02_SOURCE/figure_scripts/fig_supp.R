# Supplementary Figures S3-S3
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(dplyr);library(tidyr);library(grid)})
# Project root: the directory holding figures/, sim/, data/. Resolved from this
# script's own location, or from RHEUMLENS_ROOT, so the script runs anywhere.
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
        legend.text=element_text(size=6.4), legend.title=element_text(size=6.8))

# ---------------- S1: admitting a disease-caused covariate --------------------
m <- read.delim("sim/results/mediator_arm_summary.tsv")
m$admitted <- factor(m$admitted, levels=c("False","True"),
                     labels=c("design variables only",
                              "plus one disease-caused covariate"))
if (all(is.na(m$admitted))) {                     # R may read TRUE/FALSE as logical
  m <- read.delim("sim/results/mediator_arm_summary.tsv")
  m$admitted <- factor(ifelse(as.logical(m$admitted),
                              "plus one disease-caused covariate",
                              "design variables only"),
                       levels=c("design variables only",
                                "plus one disease-caused covariate"))
}
long <- m |>
  select(rho, admitted, `design-only AUC`=design_auc, `V_D (cross-fitted)`=I_D_cv,
         `runs flagged significant`=sig, `AUC after residualisation`=residualised_auc) |>
  pivot_longer(-c(rho, admitted))
long$name <- factor(long$name, levels=c("design-only AUC","V_D (cross-fitted)",
                                        "runs flagged significant","AUC after residualisation"))
pal2 <- c("design variables only"="#3C6E9F",
          "plus one disease-caused covariate"="#C4633E")
S1 <- ggplot(long, aes(rho, value, colour=admitted)) +
  geom_line(linewidth=.6) + geom_point(size=1.4) +
  facet_wrap(~name, nrow=1, scales="free_y") +
  scale_colour_manual(values=pal2, name=NULL) +
  labs(x="true design-diagnosis association ρ", y=NULL,
       title="Figure S3  One disease-caused covariate is enough to manufacture the finding",
       subtitle="At ρ = 0 design and diagnosis are unrelated by construction. Admitting a severity variable flags the cohort in 100% of runs.") +
  base + theme(legend.position="bottom",
               strip.background=element_blank(),
               strip.text=element_text(size=6.8, face="bold"))
ggsave("figures/out/Fig_S3_mediator.svg", S1, width=180, height=68, units="mm",
       device=svglite::svglite)
ragg::agg_png("figures/out/Fig_S3_mediator.png", width=180, height=68, units="mm", res=400)
print(S1); invisible(dev.off())

# ---------------- S2: design composition of every cohort ----------------------
s <- read.delim("data/strata_composition.tsv")
lab <- c(CMV_HIHA="CMV · HIHA", COVID_REN_ASSAY_10x_5_v2="COVID-19 · Ren (one assay)",
         COVID_STEPHENSON="COVID-19 · Stephenson", COVID_REN="COVID-19 · Ren",
         COMBAT_COVID="COVID-19 · COMBAT", SLE_GSE174188_ALLCELL="SLE · GSE174188 all",
         SLE_GSE174188_CD4="SLE · GSE174188 CD4", COMBAT_CROSS="Sepsis vs COVID-19",
         COMBAT_INFLUENZA="Influenza · COMBAT")
ordc <- c("CMV_HIHA","COVID_REN_ASSAY_10x_5_v2","COVID_STEPHENSON","COVID_REN",
          "COMBAT_COVID","SLE_GSE174188_ALLCELL","SLE_GSE174188_CD4",
          "COMBAT_CROSS","COMBAT_INFLUENZA")
s$cohort <- factor(s$cohort, levels=ordc, labels=lab[ordc])
sl <- s |> pivot_longer(c(cases, controls), names_to="class", values_to="n")
S2 <- ggplot(sl, aes(factor(stratum), n, fill=class)) +
  geom_col(width=.78) +
  geom_point(data=subset(s, mixed==1), aes(factor(stratum), y=-1.5),
             inherit.aes=FALSE, shape=17, size=.9, colour="grey25") +
  facet_wrap(~cohort, scales="free", ncol=3) +
  scale_fill_manual(values=c(cases="#C4633E", controls="#3C6E9F"), name=NULL) +
  scale_x_discrete(labels=NULL, breaks=NULL) +
  labs(x="collection stratum (ordered by size; identity not meaningful)", y="donors",
       title="Figure S4  Design composition of every cohort",
       subtitle="Bars are donors per collection stratum. A triangle marks a stratum holding both cases and controls — only those can be permuted within, or restricted to.") +
  base + theme(legend.position="bottom", strip.background=element_blank(),
               strip.text=element_text(size=6.6, face="bold"),
               axis.ticks.x=element_blank())
ggsave("figures/out/Fig_S4_composition.svg", S2, width=180, height=150, units="mm",
       device=svglite::svglite)
ragg::agg_png("figures/out/Fig_S4_composition.png", width=180, height=150, units="mm", res=400)
print(S2); invisible(dev.off())

# ---------------- S3: seed stability ------------------------------------------
# Select columns BY NAME. An earlier version renamed by position against a table
# that had since gained columns, which silently mislabelled every panel.
sd <- read.delim("data/seed_stability_all.tsv", check.names=FALSE)
lab <- c(CMV_HIHA                 = "CMV · HIHA",
         COVID_REN_ASSAY_10x_5_v2 = "COVID-19 · Ren (one assay)",
         COVID_STEPHENSON         = "COVID-19 · Stephenson",
         COVID_REN                = "COVID-19 · Ren",
         COMBAT_COVID             = "COVID-19 · COMBAT",
         SLE_GSE174188_ALLCELL    = "SLE · GSE174188 all cells",
         SLE_GSE174188_CD4        = "SLE · GSE174188 CD4",
         COMBAT_CROSS             = "Sepsis vs COVID-19 · COMBAT",
         COMBAT_INFLUENZA         = "Influenza · COMBAT")
stopifnot(all(sd$cohort %in% names(lab)))
sd$cohort <- factor(unname(lab[sd$cohort]), levels = unname(lab))

auc_lv <- c("design-only AUC","V_D","diagnosis AUC")
p_lv   <- c("p, design-only AUC","p, V_D","p, free permutation",
            "p, collection-preserving")
long <- sd |>
  select(cohort, seed,
         `design-only AUC`          = design_auc_linear,
         `V_D`                      = I_D_cv,
         `diagnosis AUC`            = disease_auc,
         `p, design-only AUC`       = p_design_auc,
         `p, V_D`                   = p_I_D,
         `p, free permutation`      = p_standard,
         `p, collection-preserving` = p_collection_preserving) |>
  pivot_longer(-c(cohort, seed))

panel <- function(d, lv, logy) {
  d$name <- factor(d$name, levels = lv)
  g <- ggplot(d, aes(factor(seed), value, colour=cohort, group=cohort)) +
    geom_line(linewidth=.4, alpha=.85) + geom_point(size=1.1) +
    facet_wrap(~name, nrow=1, scales="free_y") +
    scale_colour_manual(values=c("#4E8C6E","#E0A07A","#D4794E","#7A3B22","#A8452A",
                                 "#6E9BC4","#3C6E9F","#B4913C","#7B6BA8"), name=NULL,
                        drop=FALSE) +
    base + theme(strip.background=element_blank(),
                 strip.text=element_text(size=6.8,face="bold"),
                 axis.text.x=element_text(angle=90, hjust=1, size=4.8))
  if (logy) g <- g +
    geom_hline(yintercept=0.05, linetype="22", colour="grey55", linewidth=.3) +
    scale_y_log10(breaks=c(0.005,0.05,0.5,1),
                  labels=c("0.005","0.05","0.5","1")) +
    coord_cartesian(ylim=c(0.004,1))
  g
}
# The p-value floor is 1/(n_perm+1) = 0.005; a log axis is the only way to see
# that most comparisons sit on it rather than merely "near zero".
S3 <- (panel(filter(long, name %in% auc_lv), auc_lv, FALSE) +
         labs(x=NULL, y=NULL)) /
      (panel(filter(long, name %in% p_lv), p_lv, TRUE) +
         labs(x="seed", y=NULL)) +
  plot_layout(guides="collect", heights=c(1,1)) &
  theme(legend.position="bottom", legend.text=element_text(size=5.8),
        legend.key.height=grid::unit(3,"mm"))
S3 <- S3 + plot_annotation(
  title="Figure S7  Seed stability, nine comparisons x five seeds",
  subtitle=paste("Every per-seed value is released rather than summarised.",
                 "Dashed line marks 0.05; the axis floor 0.005 is the smallest",
                 "value 200 permutations can return."),
  theme=theme(plot.title=element_text(size=8.5,face="bold"),
              plot.subtitle=element_text(size=6.4,colour="grey35")))
ggsave("figures/out/Fig_S7_seeds.svg", S3, width=180, height=118, units="mm",
       device=svglite::svglite)
ragg::agg_png("figures/out/Fig_S7_seeds.png", width=180, height=118, units="mm", res=400)
print(S3); invisible(dev.off())
cat("Fig_S3, Fig_S4, Fig_S7 ok\n")
