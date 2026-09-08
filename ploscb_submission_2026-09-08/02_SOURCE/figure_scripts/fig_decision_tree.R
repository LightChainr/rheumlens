# Fig_decision_tree: the decision tree, with four real cohorts traced through it.
# Figure numbers live in the manuscript and in tools/build_html.py, never here.
#
# Every label here must match Section 6 and the gates in scripts/cohorts/
# run_design_screen.py. In particular Q0 is ONLY "no collection stratum holds both
# labels, or the metadata separates the labels exactly": the tuned-versus-frozen
# AUC gap is reported but is not a gate, and calling it a model-selection artefact
# was withdrawn.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
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

BW <- 3.05; BH <- 0.62          # decision box
OW <- 3.35; OH <- 0.72          # outcome box
xd <- 0; xo <- 4.5

nodes <- tribble(
  ~id,   ~x,  ~y,   ~w,  ~h,  ~kind,     ~lab,
  "S",   xd,  7.6,  BW,  .5,  "start",   "Donor-level classifier reports high internal AUC",
  "Q0",  xd,  6.5,  BW,  BH,  "dec",     "Q0  Can this be answered with these data?\na collection stratum holds both labels, AND the recorded\nmetadata does not separate the labels exactly",
  "Q1",  xd,  5.3,  BW,  BH,  "dec",     "Q1  Does the pipeline beat the free permutation?\nlabels permuted freely, whole pipeline refitted",
  "Q2",  xd,  4.1,  BW,  BH,  "dec",     "Q2  Does it beat the collection-preserving permutation?\nlabels permuted within collection strata only",
  "Q3",  xd,  2.9,  BW,  BH,  "dec",     "Q3  Is the metadata-diagnosis link weak?\ndesign-only AUC near chance and not beyond its own null",
  "Q4",  xd,  1.7,  BW,  BH,  "dec",     "Q4  Is a collection stratum large enough to restrict to?",
  "Q5",  xd,  0.5,  BW,  BH,  "dec",     "Q5  Do restriction and residualisation agree?",
  "E",   xd, -0.8,  BW,  .5,  "gate",    "EXIT GATE  ·  external target collected a different way",
  "C0",  xo,  6.5,  OW,  OH,  "stop",    "CANNOT BE ANSWERED\nNo collection stratum holds both labels, or the metadata\nseparates them exactly. → no claim of any kind.",
  "N1",  xo,  5.3,  OW,  OH,  "stop",    "NO SIGNAL\nThe observed AUC lies inside its own null.\n→ nothing to adjust or attribute.",
  "C2",  xo,  4.1,  OW,  OH,  "stop",    "CONSISTENT WITH THE COLLECTION STRATA\nAccuracy does not exceed what the strata alone produce.\n→ separability only.",
  "A",   xo,  2.9,  OW,  OH,  "route",   "ROUTE A · WEAK NUISANCE\nLittle to remove. Adjust only with its cost reported:\na wide design matrix takes accuracy away regardless.",
  "C4",  xo,  1.7,  OW,  OH,  "stop",    "ROUTE C · REDESIGN\nNo stratum is large enough to restrict within. Escalate\nto a crossed-design cohort. → separability only.",
  "B1",  xo,  0.85, OW,  .62, "route",   "ROUTE B1 · CONCORDANT\nThe choice of adjustment is not driving the conclusion.",
  "B2",  xo,  0.05, OW,  .62, "route",   "ROUTE B2 · DISCORDANT\nReport both estimates. They are NOT bounds on a biological effect."
)

edges <- tribble(
  ~from, ~to,  ~lab,  ~side,
  "S","Q0","",        "down",
  "Q0","C0","no",     "right",
  "Q0","Q1","yes",    "down",
  "Q1","N1","no",     "right",
  "Q1","Q2","yes",    "down",
  "Q2","C2","no",     "right",
  "Q2","Q3","yes",    "down",
  "Q3","A","yes",     "right",
  "Q3","Q4","no",     "down",
  "Q4","C4","no",     "right",
  "Q4","Q5","yes",    "down",
  "Q5","B1","yes",    "right",
  "Q5","B2","no",     "right",
  "Q5","E","",        "down"
)
gx <- function(id, f) nodes[[f]][match(id, nodes$id)]
edges <- edges |> mutate(
  x  = gx(from,"x"), y = gx(from,"y"), xe = gx(to,"x"), ye = gx(to,"y"),
  hh = gx(from,"h")/2, hh2 = gx(to,"h")/2,
  x0 = ifelse(side=="down", x, x + gx(from,"w")/2),
  y0 = ifelse(side=="down", y - hh, y),
  x1 = ifelse(side=="down", xe, xe - gx(to,"w")/2),
  y1 = ifelse(side=="down", ye + hh2, ye))

# four real cohorts traced through the tree
traces <- bind_rows(
  tibble(cohort="COMBAT influenza",     step=c("S","Q0","C0")),
  tibble(cohort="CMV HIHA",             step=c("S","Q0","Q1","Q2","Q3","A")),
  tibble(cohort="COMBAT sepsis/COVID",  step=c("S","Q0","Q1","Q2","C2")),
  tibble(cohort="COVID-19 Ren",         step=c("S","Q0","Q1","Q2","Q3","Q4","Q5","B2"))
)
tpal <- c("COMBAT influenza"="#7B6BA8", "CMV HIHA"="#4E8C6E",
          "COMBAT sepsis/COVID"="#B4913C", "COVID-19 Ren"="#C4633E")
off  <- c("COMBAT influenza"=-.09, "CMV HIHA"=-.03,
          "COMBAT sepsis/COVID"=.03, "COVID-19 Ren"=.09)

tseg <- traces |> group_by(cohort) |>
  reframe(from = head(step,-1), to = tail(step,-1)) |>
  left_join(edges |> select(from,to,x0,y0,x1,y1,side), by=c("from","to")) |>
  mutate(d = off[cohort],
         x0 = x0 + ifelse(side=="down", d, 0), x1 = x1 + ifelse(side=="down", d, 0),
         y0 = y0 + ifelse(side=="down", 0, d), y1 = y1 + ifelse(side=="down", 0, d))

npal <- c(start="grey93", dec="#E8EEF4", gate="grey88",
          stop="#F6E7E1", route="#E4EFE6")

p <- ggplot() +
  geom_segment(data=edges, aes(x0,y0,xend=x1,yend=y1), colour="grey72",
               linewidth=.4, arrow=arrow(length=unit(1.5,"mm"), type="closed")) +
  geom_segment(data=tseg, aes(x0,y0,xend=x1,yend=y1, colour=cohort),
               linewidth=.75, lineend="round", alpha=.95) +
  geom_tile(data=nodes, aes(x,y,width=w,height=h, fill=kind),
            colour="grey40", linewidth=.35) +
  geom_text(data=nodes, aes(x,y,label=lab), size=1.85, lineheight=1.05) +
  geom_text(data=filter(edges, lab!=""),
            aes(x = ifelse(side=="down", x0+.10, (x0+x1)/2),
                y = ifelse(side=="down", (y0+y1)/2, (y0+y1)/2+.10), label=lab),
            size=1.85, colour="grey30", hjust=0, fontface="italic") +
  scale_fill_manual(values=npal, guide="none") +
  scale_colour_manual(values=tpal, name="traced cohort") +
  coord_cartesian(xlim=c(xd-BW/2-.15, xo+OW/2+.15), ylim=c(-1.20, 8.00)) +
  guides(colour=guide_legend(nrow=1, override.aes=list(linewidth=1.1))) +
  theme_void(base_size=8) +
  theme(legend.position="bottom", legend.text=element_text(size=6.4),
        legend.title=element_text(size=6.8), legend.key.width=unit(8,"mm"),
        plot.title=element_text(size=8, face="bold"),
        plot.subtitle=element_text(size=6.4, colour="grey35", margin=margin(b=3)),
        plot.margin=margin(4,4,2,4)) +
  labs(title="A decision tree for applying the checks",
       subtitle="Feasibility first, then two permutation questions that are not the same question.\nWhere the two adjustments disagree the tree reports both and claims no bound. Four real comparisons traced.")

ggsave("figures/out/Fig_decision_tree.svg", p, width=180, height=126, units="mm",
       device=svglite::svglite)
ragg::agg_png("figures/out/Fig_decision_tree.png", width=180, height=126,
              units="mm", res=400)
print(p); invisible(dev.off()); cat("ok\n")
