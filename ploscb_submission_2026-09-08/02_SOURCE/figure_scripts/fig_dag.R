# Figure: causal setting for patient-level single-cell classification
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
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

nd <- data.frame(
  id   = c("R","D","Y","X","C","U","P"),
  x    = c(0.00, 1.20, 1.20, 2.75, 2.00, 2.00, 3.85),
  y    = c(1.35, 2.20, 0.50, 1.35, 0.30, 2.55, 1.35),
  lab  = c("R", "D", "Y", "X", "C", "U", "Ŷ"),
  role = c("recruitment","recorded metadata","phenotype label","measured cells",
           "label-derived covariate","unrecorded structure","classifier output"),
  kind = c("process","design","label","data","trap","hidden","data"),
  # +1 puts the role label above the node, -1 below
  side = c(-1, 1, -1, -1, -1, 1, -1))

eg <- data.frame(
  from = c("R","R","D","Y","Y","C","U","X"),
  to   = c("D","Y","X","X","C","X","X","P"),
  type = c("recruit","recruit","technical","biological","mediation","mediation",
           "unrecorded","model"))
eg$x    <- nd$x[match(eg$from, nd$id)]; eg$y    <- nd$y[match(eg$from, nd$id)]
eg$xend <- nd$x[match(eg$to,   nd$id)]; eg$yend <- nd$y[match(eg$to,   nd$id)]
# shorten each arrow so it stops at the node boundary
shrink <- function(d, r1 = 0.20, r2 = 0.235) {
  dx <- d$xend - d$x; dy <- d$yend - d$y; L <- sqrt(dx^2 + dy^2)
  d$x <- d$x + dx/L*r1; d$y <- d$y + dy/L*r1
  d$xend <- d$xend - dx/L*r2; d$yend <- d$yend - dy/L*r2; d }
eg <- shrink(eg)

epal <- c(recruit = "grey45", technical = "#C4633E", biological = "#3C6E9F",
          mediation = "#B4913C", unrecorded = "grey65", model = "grey25")
npal <- c(process = "grey93", design = "#F3DCD3", label = "#D9E4EF",
          data = "grey97", trap = "#F7EDD6", hidden = "white")

p <- ggplot() +
  geom_segment(data = subset(eg, type != "unrecorded"),
               aes(x, y, xend = xend, yend = yend, colour = type),
               linewidth = .55, arrow = arrow(length = unit(1.9, "mm"),
                                              type = "closed")) +
  geom_segment(data = subset(eg, type == "unrecorded"),
               aes(x, y, xend = xend, yend = yend, colour = type),
               linewidth = .55, linetype = "21",
               arrow = arrow(length = unit(1.9, "mm"), type = "closed")) +
  geom_point(data = nd, aes(x, y, fill = kind),
             shape = 21, size = 11, stroke = .5,
             colour = ifelse(nd$kind == "hidden", "grey55", "grey25")) +
  geom_text(data = nd, aes(x, y, label = lab), size = 3.4, fontface = "bold") +
  geom_text(data = nd, aes(x, y + side * .40, label = role), size = 2.15,
            colour = "grey30", lineheight = .95) +
  scale_colour_manual(values = epal, name = NULL,
    breaks = c("recruit","technical","biological","mediation","unrecorded","model"),
    labels = c("recruitment couples D and Y",
               "technical effect  D → X",
               "biological effect  Y → X",
               "mediation  Y → C → X",
               "unrecorded structure (not measurable here)",
               "the fitted model")) +
  scale_fill_manual(values = npal, guide = "none") +
  coord_cartesian(xlim = c(-.50, 4.35), ylim = c(-.20, 3.10)) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(linewidth = .8))) +
  theme_void(base_size = 8) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6.6),
        legend.key.width = unit(9, "mm"), legend.key.height = unit(4, "mm"),
        plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.5, colour = "grey35",
                                     margin = margin(b = 2)),
        plot.margin = margin(4, 4, 2, 4)) +
  labs(title = NULL, subtitle = NULL)

# annotation: the two things the screen can and cannot see
p <- p +
  annotate("segment", x = 1.20, xend = 1.20, y = 1.94, yend = 0.76,
           linetype = "12", colour = "grey55", linewidth = .4) +
  annotate("text", x = 1.32, y = 1.35,
           label = "the association the\nmetadata-only AUC quantifies",
           size = 2.0, colour = "grey35", hjust = 0, lineheight = .95)

ggsave("figures/out/Fig_dag.svg", p, width = 130, height = 92, units = "mm",
       device = svglite::svglite)
ragg::agg_png("figures/out/Fig_dag.png", width = 130, height = 92,
              units = "mm", res = 400)
print(p); invisible(dev.off()); cat("ok\n")
