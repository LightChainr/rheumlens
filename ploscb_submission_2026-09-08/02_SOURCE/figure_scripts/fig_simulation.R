# Fig_simulation: the extended simulation, seven regimes x seven rho x two
# label types.
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
setwd(W)
d <- read.delim("sim/results/extended_simulation_raw.tsv")

armlab <- c(base="additive linear", nonlinear_design="nonlinear design effect",
            heteroscedastic="heteroscedastic noise",
            design_by_celltype="design × cell-type interaction",
            imbalanced="20% cases, matched thresholds",
            imbalanced_offset="20% cases, offset collection boundary",
            mechanism_shift="target acquired differently")
apal <- c(base="#2F3A45", nonlinear_design="#3C6E9F", heteroscedastic="#4E8C6E",
          design_by_celltype="#7B6BA8", imbalanced="#C4633E",
          imbalanced_offset="#8C3B26", mechanism_shift="#B4913C")
base_t <- theme_classic(base_size=8) +
  theme(axis.text=element_text(colour="grey20",size=6.6), axis.title=element_text(size=7.5),
        plot.title=element_text(size=8,face="bold",margin=margin(b=1)),
        plot.subtitle=element_text(size=6.3,colour="grey35",margin=margin(b=4)),
        plot.title.position="plot", legend.key.size=unit(6,"pt"),
        legend.text=element_text(size=6.3), legend.title=element_text(size=6.6),
        legend.margin=margin(0,0,0,0))

agg <- d %>% group_by(arm,label_type,rho) %>%
  summarise(across(c(design_c,I_D_cv,disease_c,residualised_c,
                     restricted_c,external_c), mean, na.rm=TRUE),
            p_sig=mean(p_I_D<=0.05, na.rm=TRUE), .groups="drop")

# A: does p(V_D) calibrate?
A <- filter(agg, label_type=="binary")
pA <- ggplot(A, aes(rho, p_sig, colour=arm)) +
  geom_hline(yintercept=c(0.05,1), linetype="22", colour="grey65", linewidth=.3) +
  geom_line(linewidth=.55) + geom_point(size=1.2) +
  scale_colour_manual(values=apal, labels=armlab, name=NULL) +
  scale_y_continuous(limits=c(0,1.02), breaks=c(0,.25,.5,.75,1)) +
  labs(x="design-label association ρ", y="fraction of runs with p(V_D) ≤ 0.05",
       title="A   The screen's null calibrates",
       subtitle="at ρ = 0 the false-positive rate sits at the nominal level in every regime") +
  base_t + theme(legend.position="none")

# B: residualisation and restriction diverge
B <- filter(agg, label_type=="binary", arm=="base") %>%
  select(rho, unadjusted=disease_c, residualised=residualised_c,
         restricted=restricted_c, external=external_c) %>%
  pivot_longer(-rho)
B$name <- factor(B$name, levels=c("unadjusted","external","restricted","residualised"))
pB <- ggplot(B, aes(rho, value, colour=name)) +
  geom_hline(yintercept=.5, linetype="22", colour="grey65", linewidth=.3) +
  geom_line(linewidth=.6) + geom_point(size=1.2, na.rm=TRUE) +
  scale_colour_manual(values=c(unadjusted="#2F3A45", external="#3C6E9F",
                               restricted="#4E8C6E", residualised="#C4633E"), name=NULL) +
  scale_y_continuous(limits=c(.3,.95)) +
  labs(x="design-label association ρ", y="AUC",
       title="B   Residualisation and restriction separate",
       subtitle="additive linear regime; restriction is undefined at ρ = 1 (no mixed stratum)") +
  base_t + theme(legend.position=c(.16,.28), legend.background=element_blank())

# C: the gap is regime-independent, except under label imbalance
C <- filter(agg, label_type=="binary") %>%
  mutate(gap = restricted_c - residualised_c)
pC <- ggplot(C, aes(rho, gap, colour=arm)) +
  geom_hline(yintercept=0, colour="grey65", linewidth=.3) +
  geom_line(linewidth=.55) + geom_point(size=1.2, na.rm=TRUE) +
  scale_colour_manual(values=apal, labels=armlab, name=NULL) +
  guides(colour=guide_legend(nrow=3)) +
  labs(x="design-label association ρ", y="restricted AUC − residualised AUC",
       title="C   The gap is not an artefact of the linear regime",
       subtitle="it opens in every generative model tested.\nRestriction is undefined at ρ = 1 except under imbalance") +
  base_t + theme(legend.position="bottom")

# D: an external target sharing the design mechanism confirms a design-driven model
D <- filter(agg, label_type=="binary", arm %in% c("base","mechanism_shift")) %>%
  select(arm, rho, internal=disease_c, external=external_c) %>%
  pivot_longer(c(internal,external))
pD <- ggplot(D, aes(rho, value, colour=arm, linetype=name)) +
  geom_line(linewidth=.6) + geom_point(size=1.2) +
  scale_colour_manual(values=apal[c("base","mechanism_shift")],
                      labels=c(base="target shares the design mechanism",
                               mechanism_shift="target acquired differently"), name=NULL) +
  scale_linetype_manual(values=c(internal="solid", external="22"), name=NULL,
                        labels=c(internal="internal AUC", external="external AUC")) +
  scale_y_continuous(limits=c(.75,.95)) +
  guides(colour=guide_legend(nrow=2, order=1),
         linetype=guide_legend(nrow=2, order=2,
                               override.aes=list(colour="grey30"))) +
  labs(x="design-label association ρ", y="AUC",
       title="D   External validation can confirm a design-driven model",
       subtitle="external AUC rises with confounding when the target\nshares the acquisition process") +
  base_t + theme(legend.position="bottom", legend.box="vertical",
                 legend.spacing.y=unit(0,"pt"))

fig <- (pA|pB)/(pC|pD) + plot_layout(heights=c(1,1.18))
ggsave("figures/out/Fig_simulation.svg", fig, width=180, height=150, units="mm",
       device=svglite::svglite)
ragg::agg_png("figures/out/Fig_simulation.png", width=180, height=150, units="mm", res=400)
print(fig); invisible(dev.off()); cat("ok\n")
