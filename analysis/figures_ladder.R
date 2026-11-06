## The four levels of integration (main design, five bases, partial readings):
## ladder_1 (increments), ladder_2 (raw scores), ladder_3 (rank of the latent
## covariance), ladder_4 (what each step brings). Each estimator at its best
## basis in each setting; intervals are 95 % across settings.
##   Rscript analysis/figures_ladder.R
source("analysis/read_results.R")
source("analysis/figures_common.R")

RUNGS <- c(quick = "1. log-GGM", pca10 = "2. PLNPCA (rank 10) + glasso",
           plnglasso = "3. PLN + glasso", plnnetwork = "4. PLNnetwork")
RUNG.COL <- setNames(c(COLOURS[["quick"]], "#E69F00", COLOURS[["plnglasso"]], COLOURS[["plnnetwork"]]), RUNGS)
RUNG.SHP <- setNames(c(17, 15, 18, 16), RUNGS)

C <- ladder_cells()
C$mechanism <- mech_factor(C$mech); C$metric.lab <- metric_factor(C$metric)

## 1-2. the levels, as increments and as raw scores
r <- C[C$estimator %in% names(RUNGS), ]; r$level <- factor(RUNGS[r$estimator], RUNGS)
ladder <- function(y, ylab, hline) {
  ggplot(r, aes(level, .data[[y]], colour = level, shape = level, group = 1)) +
    geom_hline(yintercept = hline, linewidth = 0.3, colour = "grey60") +
    stat_summary(fun = mean, geom = "line", linewidth = 0.6, colour = "grey55") +
    stat_summary(fun.data = mean_cl_normal, geom = "pointrange", size = 0.5, linewidth = 0.8) +
    facet_grid(metric.lab ~ mechanism, scales = "free_y", labeller = row_wrap) +
    scale_colour_manual(values = RUNG.COL, name = NULL) + scale_shape_manual(values = RUNG.SHP, name = NULL) +
    scale_x_discrete(labels = function(x) sub("\\. .*$", "", x)) +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
    labs(x = "level of integration", y = ylab) + th }
ggsave(fig_path("ladder_1_echelle"), ladder("inc", LAB.INC, 0), width = 12, height = 9, dpi = 300)
p <- ladder("raw", "raw score (with interactions)", NULL) +
  geom_hline(data = data.frame(metric.lab = factor(METRIC.LAB[1:2], METRIC.LAB), y = 0.5),
             aes(yintercept = y), linewidth = 0.3, colour = "grey60", linetype = "22")
ggsave(fig_path("ladder_2_brut"), p, width = 12, height = 9, dpi = 300)

## 3. rank of the latent covariance: PLNPCA at ranks 2, 5, 10, PLN (full rank);
## open triangle: gllvm (2 latent variables); dashes: the log-GGM
rk <- C[C$estimator %in% c("pca2", "pca5", "pca10", "plnglasso"), ]
rk$rank <- c(pca2 = 2, pca5 = 5, pca10 = 10, plnglasso = 20)[rk$estimator]
gl <- C[C$estimator == "gllvm", ]; gl$rank <- 2
rf <- aggregate(inc ~ mechanism + metric.lab, C[C$estimator == "quick", ], mean)
p <- ggplot(rk, aes(rank, inc, colour = mechanism, group = mechanism)) +
  geom_hline(data = rf, aes(yintercept = inc, colour = mechanism), linetype = "22", linewidth = 0.5) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", size = 0.4, linewidth = 0.7) +
  stat_summary(data = gl, fun.data = mean_cl_normal, geom = "pointrange", shape = 2, size = 0.5,
               linewidth = 0.6, position = position_nudge(x = -0.6)) +
  facet_wrap(~ metric.lab, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = MECH.COL, name = NULL) +
  scale_x_continuous(breaks = c(2, 5, 10, 20), labels = c("2", "5", "10", "20\n(PLN)")) +
  labs(x = "rank of the latent covariance", y = LAB.INC) + th
ggsave(fig_path("ladder_3_rang"), p, width = 14, height = 4.6, dpi = 300)

## 4. what each step brings: difference of increments between two levels;
## under each interval, the proportion of settings where the step is positive
STEP.LAB <- c(count = "count model (3 − 1)", rank = "full rank (3 − 2)",
              integration = "integration (4 − 3)", total = "overall (4 − 1)")
st <- ladder_steps(C)
st$mechanism <- mech_factor(st$mech); st$metric.lab <- metric_factor(st$metric)
st$step.lab <- factor(STEP.LAB[st$step], STEP.LAB)
sm <- do.call(rbind, lapply(split(st, list(st$mechanism, st$metric.lab, st$step.lab), drop = TRUE), function(g) {
  ci <- mean_cl_normal(g$d)
  data.frame(mechanism = g$mechanism[1], metric.lab = g$metric.lab[1], step.lab = g$step.lab[1],
             y = ci$y, lo = ci$ymin, hi = ci$ymax, label = sprintf("%.0f %%", 100 * mean(g$d > 0))) }))
STEP.COL <- setNames(c("#56B4E9", "#E69F00", "#0072B2", "#999999"), STEP.LAB)
dg <- position_dodge(0.75)
p <- ggplot(sm, aes(mechanism, y, colour = step.lab)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey45") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = dg, size = 0.45, linewidth = 0.7) +
  geom_text(aes(y = lo, label = label), position = dg, vjust = 1.6, size = 2.7, show.legend = FALSE) +
  facet_wrap(~ metric.lab, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = STEP.COL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.05))) +
  labs(x = NULL, y = "difference in increment") + th
ggsave(fig_path("ladder_4_marches"), p, width = 11, height = 8, dpi = 300)
