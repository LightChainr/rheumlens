# Figure: the design-label confounding spectrum across 9 contrasts in 4 diseases
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
})
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

# All four inputs are produced by tools/stage_screen_results.py from the screen's
# own output, so no number in this figure is transcribed.
blocks <- read.delim("data/design_screen_all_blocks.tsv", check.names = FALSE)
spec   <- read.delim("data/FINAL_confounding_spectrum.tsv", check.names = FALSE)
coll   <- read.delim("data/collection_block.tsv", check.names = FALSE)
seed   <- read.delim("data/seed_stability_all.tsv", check.names = FALSE)
# Frozen pipeline throughout: it is the statistic every p-value here was
# computed against, and the one Table 1 prints. Plotting the tuned nested-CV
# AUC instead put two different numbers under one name - for COMBAT influenza
# the diagnosis classifier read 0.996 here and 0.475 in Table 1.
names(seed)[names(seed) == "design_auc_frozen"] <- "design_auc"
spec$design_auc     <- spec$design_auc_frozen
spec$disease_auc    <- spec$observed_frozen_auc

lab <- c(CMV_HIHA                 = "CMV · HIHA",
         COVID_REN_ASSAY_10x_5_v2 = "COVID-19 · Ren (one assay)",
         COVID_STEPHENSON         = "COVID-19 · Stephenson",
         COVID_REN                = "COVID-19 · Ren",
         COMBAT_COVID             = "COVID-19 · COMBAT",
         SLE_GSE174188_ALLCELL    = "SLE · GSE174188 all cells",
         SLE_GSE174188_CD4        = "SLE · GSE174188 CD4",
         COMBAT_CROSS             = "Sepsis vs COVID-19 · COMBAT",
         COMBAT_INFLUENZA         = "Influenza · COMBAT")
dis <- c(CMV_HIHA="CMV", COVID_REN_ASSAY_10x_5_v2="COVID-19",
         COVID_STEPHENSON="COVID-19", COVID_REN="COVID-19", COMBAT_COVID="COVID-19",
         SLE_GSE174188_ALLCELL="SLE", SLE_GSE174188_CD4="SLE",
         COMBAT_CROSS="Sepsis vs COVID-19", COMBAT_INFLUENZA="Influenza")
pal <- c(SLE="#3C6E9F", "COVID-19"="#C4633E", Influenza="#7B6BA8",
         CMV="#4E8C6E", "Sepsis vs COVID-19"="#B4913C")

base <- theme_classic(base_size = 8) +
  theme(axis.text     = element_text(colour = "grey20", size = 6.6),
        axis.title    = element_text(size = 7.5),
        plot.title    = element_text(size = 8, face = "bold", hjust = 0,
                                     margin = margin(b = 1)),
        plot.subtitle = element_text(size = 6.3, colour = "grey35",
                                     margin = margin(b = 4)),
        plot.title.position = "plot",
        legend.key.size = unit(6, "pt"),
        legend.text   = element_text(size = 6.4),
        legend.title  = element_text(size = 6.8),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(3, "pt"))

# strongest confounding at the top
ord <- spec$cohort[order(spec$design_auc)]        # ascending
ylev <- ord                                       # ggplot draws level 1 at the bottom

# ---------------- Panel A ----------------------------------------------------
A <- spec; A$disease <- dis[A$cohort]
rng <- seed |> group_by(cohort) |>
  summarise(lo = min(design_auc), hi = max(design_auc), .groups = "drop")
# EVERY join must happen before the factor call: a left_join on a factor key
# coerces it back to character and the row order silently reverts.
A <- left_join(A, rng, by = "cohort")
A <- left_join(A, coll[c("cohort","design_auc_frozen")] |>
                 rename(collection_auc = design_auc_frozen), by = "cohort")
A$cohort <- factor(A$cohort, levels = ylev)
stopifnot(!any(is.na(A$cohort)), is.factor(A$cohort))
# Verdict comes from the screen's own gate columns, not from a hand-kept list.
A$verdict <- ifelse(A$degenerate, "cannot be answered",
             ifelse(A$underpowered, "underpowered", NA))

pA <- ggplot(A, aes(y = cohort)) +
  geom_vline(xintercept = .5, linetype = "22", colour = "grey65", linewidth = .3) +
  geom_segment(aes(x = design_auc, xend = disease_auc, yend = cohort),
               colour = "grey82", linewidth = .5) +
  geom_point(aes(x = disease_auc), shape = 21, fill = "white",
             colour = "grey45", size = 1.8, stroke = .45) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                linewidth = .9, colour = "grey25", na.rm = TRUE) +
  geom_point(aes(x = design_auc, fill = disease), shape = 21, size = 2.5,
             colour = "white", stroke = .5) +
  geom_point(aes(x = collection_auc), shape = 124, size = 2.4, colour = "grey20",
             na.rm = TRUE) +
  geom_text(aes(x = 0.335, label = verdict), hjust = 0, size = 1.95,
            colour = "grey35", fontface = "italic", na.rm = TRUE) +
  scale_y_discrete(labels = lab) +
  scale_x_continuous(breaks = seq(.4, 1, .1), expand = c(0, 0)) +
  coord_cartesian(xlim = c(.325, 1.0), clip = "off") +
  scale_fill_manual(values = pal, name = NULL) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "cross-fitted AUC", y = NULL, title = "A   How strongly recorded metadata predicts the diagnosis",
       subtitle = "filled = all recorded metadata   | = collection only\nopen = full diagnosis classifier   bar = 5-seed range") +
  base + theme(legend.position = "bottom",
               plot.margin = margin(4, 6, 2, 2))

# ---------------- Panel B ----------------------------------------------------
grid <- expand.grid(cohort = ylev,
                    block  = c("quality","demographic","collection","all"),
                    stringsAsFactors = FALSE)
B <- left_join(grid, blocks[c("cohort","block","design_auc_frozen","p_design_auc")],
               by = c("cohort","block"))
B$cohort <- factor(B$cohort, levels = ylev)
B$block  <- factor(B$block, levels = c("quality","demographic","collection","all"))
B$sig  <- !is.na(B$p_design_auc) & B$p_design_auc <= .05
B$txt  <- ifelse(is.na(B$design_auc_frozen), "n/a",
                 sprintf("%.2f", B$design_auc_frozen))

pB <- ggplot(B, aes(x = block, y = cohort)) +
  geom_tile(aes(fill = design_auc_frozen), colour = "white", linewidth = 1.1) +
  geom_tile(data = subset(B, sig), fill = NA, colour = "grey15", linewidth = .45) +
  geom_text(aes(label = txt,
                colour = !is.na(design_auc_frozen) & design_auc_frozen > .78),
            size = 2.05, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20")) +
  scale_fill_gradient(low = "#F4F1EA", high = "#8C3B26", limits = c(.3, 1),
                      na.value = "grey94", name = "design AUC",
                      breaks = c(.4, .7, 1.0)) +
  scale_y_discrete(labels = NULL, breaks = NULL) +
  guides(fill = guide_colourbar(barwidth = unit(26, "mm"),
                                barheight = unit(2.6, "mm"),
                                title.position = "top")) +
  labs(x = NULL, y = NULL, title = "B   Which variables carry the diagnosis",
       subtitle = "dark outline: permutation p ≤ 0.05\nagainst that group's own null") +
  base + theme(axis.line = element_blank(), axis.ticks = element_blank(),
               legend.position = "bottom",
               plot.margin = margin(4, 2, 2, 0))

# ---------------- Panel C ----------------------------------------------------
Cd <- subset(B, cohort %in% c("COVID_STEPHENSON","COVID_REN") &
                !is.na(design_auc_frozen))
Cd$cohort <- droplevels(Cd$cohort)
pC <- ggplot(Cd, aes(x = block, y = design_auc_frozen,
                     group = cohort, colour = cohort)) +
  geom_hline(yintercept = .5, linetype = "22", colour = "grey65", linewidth = .3) +
  geom_line(linewidth = .6) +
  geom_point(aes(shape = sig), size = 2.1, fill = "white", stroke = .7) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  scale_colour_manual(values = c(COVID_STEPHENSON = "#D4794E", COVID_REN = "#7A3B22"),
                      labels = c(COVID_STEPHENSON = "Stephenson",
                                 COVID_REN = "Ren"), name = NULL) +
  scale_y_continuous(limits = c(.45, .90), breaks = seq(.5, .9, .1)) +
  labs(x = NULL, y = "design-only AUC",
       title = "C   Same strength, different source",
       subtitle = "solid point: significant in every seed.\nStephenson is affected through sample quality, Ren through hospital") +
  base + theme(legend.position = c(.18, .88),
               legend.background = element_blank(),
               plot.margin = margin(6, 4, 2, 2))

# ---------------- Panel D ----------------------------------------------------
D <- seed |> select(cohort, seed, p_standard, p_collection_preserving) |>
  pivot_longer(c(p_standard, p_collection_preserving),
               names_to = "test", values_to = "p")
D$test <- factor(D$test, levels = c("p_standard","p_collection_preserving"),
                 labels = c("free permutation",
                            "collection-preserving permutation"))
D$cohort <- factor(D$cohort, levels = intersect(ylev, unique(D$cohort)))
D$p <- pmax(D$p, 5e-4)

pD <- ggplot(D, aes(x = p, y = cohort, colour = test)) +
  geom_vline(xintercept = .05, linetype = "22", colour = "grey55", linewidth = .3) +
  geom_point(position = position_dodge(width = .5), size = 1.6, alpha = .85) +
  # 0.001 is the floor of a 1,000-permutation test, so label it rather than a
  # smaller tick no point can reach.
  scale_x_log10(breaks = c(1e-3, .01, .05, .5),
                labels = c("0.001","0.01","0.05","0.5")) +
  coord_cartesian(xlim = c(8e-4, 1)) +
  scale_y_discrete(labels = lab) +
  scale_colour_manual(values = c("grey35", "#B4553C"), name = NULL) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "permutation p (one point per seed)", y = NULL,
       title = "D   Where the two permutation tests disagree",
       subtitle = "only the sepsis-vs-COVID-19 comparison separates them.\n0.001 is the floor of the 1,000-permutation test") +
  base + theme(legend.position = "bottom",
               plot.margin = margin(6, 6, 2, 2))

fig <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1.05, 1), widths = c(1, 1))

ggsave("figures/out/Fig_spectrum.svg", fig, width = 180, height = 160,
       units = "mm", device = svglite::svglite)
ragg::agg_png("figures/out/Fig_spectrum.png", width = 180, height = 160,
              units = "mm", res = 400)
print(fig); invisible(dev.off())
cat("row order (bottom to top):\n"); print(as.character(ylev))
cat("ok\n")
