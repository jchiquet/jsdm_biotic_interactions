## Unequal sampling effort and compositional data (lottery and Gompertz):
## increment against the spread of the effort (compo_1), drift of the
## no-interaction control (compo_2, the spurious signal), raw AUC of both arms
## (compo_2b), and four configurations at a realistic spread (compo_3). Each
## reading at its best basis (poly 3 or spline 10); 95 % intervals across
## replicates.
##   Rscript analysis/figures_compo.R
source("analysis/read_results.R")
source("analysis/figures_common.R")

CORR.LAB <- c(none = "none", unknown = "unknown effort", clr = "unknown effort, CLR", known = "known effort")
CORR.COL <- setNames(c("#999999", "#E69F00", "#CC79A7", "#0072B2"), CORR.LAB)
CORR.SHP <- setNames(c(1, 17, 4, 16), CORR.LAB)
FAM.LAB  <- c(pln = LAB.PLNNET, quick = LAB.QUICK)
INT.LAB  <- c(competition = "competition", facilitation = "facilitation")

D <- compo_cells()
D$metric.lab <- metric_factor(D$metric)
D$mechanism  <- factor(MECH[D$mech], MECH[c("lottery", "gompertz")])
D$family     <- factor(FAM.LAB[D$family], FAM.LAB)
D$correction <- factor(CORR.LAB[D$correction], CORR.LAB)
D$interaction <- factor(INT.LAB[sub("^.*_", "", D$metric)], INT.LAB)
D$panel <- factor(paste(D$family, "·", D$interaction), paste(rep(levels(D$family), each = 2), "·", INT.LAB))
enc <- list(scale_colour_manual(values = CORR.COL, name = "correction"),
            scale_shape_manual(values = CORR.SHP, name = "correction"))
S.BREAKS <- c(0, 0.25, 0.5, 1)
LAB.S <- "spread of the sampling effort s"
b10 <- D$niche.breadth == 10 & D$interaction.threshold == 10
dg <- position_dodge(0.06)
line_plot <- function(d, y, ylab, hline)
  ggplot(d, aes(effort.sd, .data[[y]], colour = correction, shape = correction, group = reading)) +
    geom_hline(yintercept = hline, linewidth = 0.3, colour = "grey60") +
    stat_summary(fun = mean, geom = "line", linewidth = 0.8, position = dg) +
    stat_summary(fun.data = mean_cl_normal, geom = "pointrange", size = 0.4, linewidth = 0.6, position = dg) +
    facet_grid(mechanism ~ panel) + enc +
    scale_x_continuous(breaks = S.BREAKS) +
    labs(x = LAB.S, y = ylab) + th

auc <- D[b10 & D$metric %in% c("AUC_competition", "AUC_facilitation"), ]
ggsave(fig_path("compo_1_effort"), line_plot(auc, "inc", "AUC increment", 0),
       width = 13, height = 6.5, dpi = 300)
ggsave(fig_path("compo_2_faux_signal"), line_plot(auc, "drift", "drift of the raw AUC of the control", 0),
       width = 13, height = 6.5, dpi = 300)

## raw scores of both arms
ARM.LAB <- c("with interactions", "control")
fb <- rbind(transform(auc, y = val, arm = ARM.LAB[1]), transform(auc, y = null, arm = ARM.LAB[2]))
fb$arm <- factor(fb$arm, ARM.LAB)
fb$row <- factor(paste(fb$mechanism, "·", fb$arm), paste(rep(levels(D$mechanism), each = 2), "·", ARM.LAB))
p <- ggplot(fb, aes(effort.sd, y, colour = correction, shape = correction, group = reading)) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, colour = "grey60", linetype = "22") +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8, position = dg) +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", size = 0.4, linewidth = 0.6, position = dg) +
  facet_grid(row ~ panel, labeller = row_wrap) + enc +
  scale_x_continuous(breaks = S.BREAKS) +
  labs(x = LAB.S, y = "raw AUC") + th
ggsave(fig_path("compo_2b_brut"), p, width = 13, height = 10, dpi = 300)

## four configurations at s = 0.5; the last one: the network does not follow
## niche proximity
mw <- D[D$effort.sd == 0.5 & D$metric %in% c("precision_k_competition", "AUC_competition"), ]
W.LAB <- c("b = 5", "b = 10", "b = 20", "b = 10\nthreshold 100")
mw$world <- factor(ifelse(mw$interaction.threshold == 100, W.LAB[4], W.LAB[match(mw$niche.breadth, c(5, 10, 20))]), W.LAB)
mw$metric.lab <- factor(mw$metric.lab, METRIC.LAB[c("precision_k_competition", "AUC_competition")])
mw$column <- factor(paste(mw$mechanism, "·", mw$family),
                    paste(rep(levels(D$mechanism), each = 2), "·", levels(D$family)))
p <- ggplot(mw, aes(world, inc, colour = correction, shape = correction)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_vline(xintercept = 3.5, linewidth = 0.3, colour = "grey80") +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", size = 0.45, linewidth = 0.7,
               position = position_dodge(0.6)) +
  facet_grid(metric.lab ~ column, scales = "free_y", labeller = row_wrap) + enc +
  labs(x = NULL, y = LAB.INC) + th
ggsave(fig_path("compo_3_mondes"), p, width = 13, height = 7, dpi = 300)
