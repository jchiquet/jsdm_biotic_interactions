## PLNnetwork against the log-GGM on the main design: by abiotic basis
## (smooth_2), where PLNnetwork keeps its advantage under Gompertz (smooth_3),
## partial against marginal readings (smooth_7), and the paired effect of the
## mechanism (smooth_8). Intervals are 95 %.
##   Rscript analysis/figures_bases.R
source("analysis/read_results.R")
source("analysis/figures_common.R")
METRICS <- list(c("AUC", "competition"), c("AUC", "facilitation"),
                c("precision_k", "competition"), c("precision_k", "facilitation"))
MLAB <- unname(METRIC.LAB)

## cell means of every estimator at every basis, both grids of the main design
cells <- function(metric, interaction)
  do.call(rbind, lapply(MAIN.GRIDS, function(g) {
    y <- estimator_cells(g$grid, metric, interaction, by = g$by)
    if (!"dynamics" %in% names(y)) {
      y$dynamics <- "lottery"; y$strength <- paste("intensity", y$competition.intensity)
    } else y$strength <- paste0("κ ", y$dyn.strength, ", σ ", y$dyn.sigma)
    y$cell <- paste(y$dynamics, y$breadth, y$threshold, y$strength)
    y[c("dynamics", "breadth", "threshold", "strength", "cell", "spec", "estimator", "inc")]
  }))
C <- setNames(lapply(METRICS, function(m) cells(m[1], m[2])), MLAB)
best <- function(y, who) { z <- y[y$estimator %in% who, ]; tapply(z$inc, z$cell, max) }

## smooth_2: mean increment by basis, over all worlds and strengths
bs <- do.call(rbind, lapply(names(C), function(nm) {
  y <- C[[nm]]; y <- y[y$estimator %in% c("plnnetwork", "quick"), ]; y$metric.lab <- nm; y }))
bs$mechanism <- mech_factor(bs$dynamics); bs$metric.lab <- factor(bs$metric.lab, MLAB)
bs$base <- factor(BASE.LAB[bs$spec], BASE.LAB)
bs$reading <- factor(EST.LAB[bs$estimator], names(EST.COL))
p <- ggplot(bs, aes(base, inc, colour = reading, shape = reading, group = reading)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", size = 0.45, linewidth = 0.7) +
  facet_grid(metric.lab ~ mechanism, scales = "free_y", labeller = row_wrap) +
  scale_colour_manual(values = EST.COL, name = NULL) + scale_shape_manual(values = EST.SHP, name = NULL) +
  labs(x = "abiotic basis", y = LAB.INC) + th +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(fig_path("smooth_2_bases"), p, width = 12, height = 10, dpi = 300)

## smooth_3: PLNnetwork - log-GGM under Gompertz, by world and strength x noise
fr <- do.call(rbind, lapply(METRIC.LAB[c("AUC_competition", "precision_k_competition")], function(nm) {
  y <- C[[nm]]; y <- y[y$dynamics == "gompertz", ]
  g <- best(y, "plnnetwork") - best(y, "quick")
  k <- unique(y[c("cell", "breadth", "threshold", "strength")]); k$gap <- g[k$cell]; k$metric.lab <- nm; k }))
WT <- "b = %d, threshold %d"
fr$world <- factor(sprintf(WT, fr$breadth, fr$threshold),
                   rev(sprintf(WT, rep(c(2, 5, 10, 20), each = 2), rep(c(10, 50), 4))))
fr$strength <- factor(fr$strength, sort(unique(fr$strength)))
lim <- max(abs(fr$gap))
p <- ggplot(fr, aes(strength, world, fill = gap)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%+.2f", gap), colour = abs(gap) > 0.6 * lim), size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c("FALSE" = "grey15", "TRUE" = "white"), guide = "none") +
  facet_wrap(~ metric.lab) +
  scale_fill_gradient2(low = COLOURS[["quick"]], mid = "grey96", high = COLOURS[["plnnetwork"]],
                       midpoint = 0, limits = c(-lim, lim), name = "PLNnetwork − log-GGM") +
  labs(x = "strength × noise", y = NULL) + th +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank(),
        legend.key.width = unit(1.6, "cm"))
ggsave(fig_path("smooth_3_frontiere"), p, width = 13, height = 6, dpi = 300)

## smooth_7: PLNnetwork - thresholded covariance, same PLN fit, basis and replicate
pm <- partial_marginal()
pm$metric.lab <- metric_factor(pm$metric); pm$mechanism <- mech_factor(pm$mech)
p <- ggplot(pm[is.finite(pm$d), ], aes(factor(breadth), d, colour = mechanism)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey45") +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange",
               position = position_dodge(0.5), size = 0.45, linewidth = 0.7) +
  facet_wrap(~ metric.lab, nrow = 1) +
  scale_colour_manual(values = MECH.COL, name = NULL) +
  labs(x = "niche breadth b", y = "partial − marginal") + th
ggsave(fig_path("smooth_7_partiel_marginal"), p, width = 13, height = 5, dpi = 300)

## smooth_8: PLNnetwork's advantage over the log-GGM under a dynamics, minus
## under the lottery (intensity 2), on the same replicate
pe <- mechanism_effect()
pe$metric.lab <- factor(METRIC.LAB[pe$metric], MLAB)
CT.LAB <- setNames(paste(c("Gompertz", "Ricker"), "−", MECH[["lottery"]]), c("gompertz", "ricker"))
pe$contrast <- factor(CT.LAB[pe$mech], CT.LAB)
p <- ggplot(pe[is.finite(pe$d), ], aes(metric.lab, d, colour = contrast)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey45") +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", position = position_dodge(0.45),
               size = 0.55, linewidth = 0.8) +
  scale_colour_manual(values = setNames(MECH.COL[c("Gompertz", "Ricker")], CT.LAB), name = NULL) +
  scale_x_discrete(labels = function(x) gsub(", ", ",\n", x)) +
  labs(x = NULL, y = "change in the JSDM's advantage") + th
ggsave(fig_path("smooth_8_apparie"), p, width = 10, height = 5.5, dpi = 300)
