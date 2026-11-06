## Figures of the Supplementary material: the supplementary experiments at
## b = 10, threshold 10 (smooth_4, smooth_4b: noisy covariate and hidden
## gradient; smooth_5a: sites and species; smooth_5b: dense low-rank network),
## and the simulators (annex_1: convergence of the lottery; annex_2, annex_3:
## the network and the counts of one replicate). Intervals are 95 % across
## replicates.
##   Rscript analysis/figures_supplement.R
source("analysis/read_results.R")
source("analysis/figures_common.R")
HID.LAB <- c("without unmeasured gradient", "with unmeasured gradient")
ARM.LAB <- c("with interactions", "no-interaction control")
LAB.NOISE <- "error on the covariate (niche breadths)"
MECH2 <- MECH[c("lottery", "gompertz")]

## replicate increments of PLNnetwork and the log-GGM, each at its best basis
annex_reps <- function(grid_of, by)
  do.call(rbind, lapply(c("lottery", "gompertz"), function(mech) do.call(rbind,
    lapply(list(c("AUC", "competition"), c("precision_k", "competition")), function(m)
      do.call(rbind, lapply(c("plnnetwork", "quick"), function(who) {
        r <- best_reps(grid_of(mech), m[1], m[2], by, who)
        r$reading <- EST.LAB[[who]]
        r$metric.lab <- METRIC.LAB[[paste0(m[1], "_", m[2])]]; r$mechanism <- MECH[[mech]]; r }))))))
lines_plot <- function(dodge) list(
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60"),
  stat_summary(fun = mean, geom = "line", linewidth = 0.8, position = position_dodge(dodge)),
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", size = 0.4, linewidth = 0.6,
               position = position_dodge(dodge)),
  facet_grid(metric.lab ~ mechanism, scales = "free_y", labeller = row_wrap),
  scale_colour_manual(values = EST.COL, name = NULL), scale_shape_manual(values = EST.SHP, name = NULL))

## smooth_4: noisy covariate and hidden gradient
c2 <- annex_reps(function(m) paste0("world_annex_c2p_", m), c("env.observed.rho", "beta.hidden"))
c2$noise <- round(noise_of(c2$env.observed.rho), 2)
c2$hidden <- factor(ifelse(c2$beta.hidden == 0, HID.LAB[1], HID.LAB[2]), HID.LAB)
c2$mechanism <- factor(c2$mechanism, MECH2)
p <- ggplot(c2, aes(noise, inc, colour = reading, shape = reading, linetype = hidden,
                    group = interaction(reading, hidden))) +
  lines_plot(dodge = 0.03) +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  scale_x_continuous(breaks = c(0, 0.1, 0.25, 0.5, 1)) +
  labs(x = LAB.NOISE, y = LAB.INC) + th
ggsave(fig_path("smooth_4_c2p"), p, width = 12, height = 8, dpi = 300)

## smooth_4b: raw AUC of both arms, at the basis retained for smooth_4
c2r <- do.call(rbind, lapply(c("lottery", "gompertz"), function(mech)
  do.call(rbind, lapply(c("plnnetwork", "quick"), function(who) {
    r <- best_raw(paste0("world_annex_c2p_", mech), "AUC", "competition", c("env.observed.rho", "beta.hidden"), who)
    r$reading <- EST.LAB[[who]]; r$mechanism <- MECH[[mech]]; r }))))
c2r$noise <- round(noise_of(c2r$env.observed.rho), 2)
c2r$hidden <- factor(ifelse(c2r$beta.hidden == 0, HID.LAB[1], HID.LAB[2]), HID.LAB)
c2r$arm.lab <- factor(ifelse(c2r$arm == "interactions", ARM.LAB[1], ARM.LAB[2]), ARM.LAB)
c2r$mechanism <- factor(c2r$mechanism, MECH2)
p <- ggplot(c2r[is.finite(c2r$val), ], aes(noise, val, colour = reading, shape = reading, linetype = hidden,
                                           group = interaction(reading, hidden))) +
  geom_hline(yintercept = 0.5, linewidth = 0.3, colour = "grey60", linetype = "22") +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8, position = position_dodge(0.03)) +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", size = 0.4, linewidth = 0.6,
               position = position_dodge(0.03)) +
  facet_grid(arm.lab ~ mechanism) +
  scale_colour_manual(values = EST.COL, name = NULL) + scale_shape_manual(values = EST.SHP, name = NULL) +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  scale_x_continuous(breaks = c(0, 0.1, 0.25, 0.5, 1)) +
  labs(x = LAB.NOISE, y = "raw AUC, competition") + th
ggsave(fig_path("smooth_4b_c2p_brut"), p, width = 12, height = 8, dpi = 300)

## smooth_5a: sites and species
np <- annex_reps(function(m) paste0("world_annex_np_smooth_", m), c("n.env", "n.species"))
np$mechanism <- factor(np$mechanism, MECH2)
np$p <- factor(paste("p =", np$n.species), c("p = 20", "p = 40"))
p <- ggplot(np, aes(n.env, inc, colour = reading, shape = reading, linetype = p,
                    group = interaction(reading, p))) +
  lines_plot(dodge = 4) +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  scale_x_continuous(breaks = c(30, 50, 100)) +
  labs(x = "number of sites n", y = LAB.INC) + th
ggsave(fig_path("smooth_5a_np"), p, width = 12, height = 8, dpi = 300)

## smooth_5b: dense network of low rank
c4 <- annex_reps(function(m) paste0("world_annex_c4_smooth_", m), "network.rank")
c4$mechanism <- factor(c4$mechanism, MECH2)
p <- ggplot(c4, aes(factor(network.rank), inc, colour = reading, shape = reading, group = reading)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8, position = position_dodge(0.15)) +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", size = 0.45, linewidth = 0.7,
               position = position_dodge(0.15)) +
  facet_grid(metric.lab ~ mechanism, scales = "free_y", labeller = row_wrap) +
  scale_colour_manual(values = EST.COL, name = NULL) + scale_shape_manual(values = EST.SHP, name = NULL) +
  labs(x = "network rank", y = LAB.INC) + th
ggsave(fig_path("smooth_5b_c4"), p, width = 11, height = 8, dpi = 300)

## annex_1: convergence of the lottery (describe_simulations.R)
DIAG <- file.path(RESULTS, "diagnostics")
ARM <- c(interactions = "with interactions", control = "no-interaction control")
cv <- read.csv(file.path(DIAG, "lottery_convergence.csv"))
cv$b <- factor(cv$breadth); cv$arm <- factor(ARM[cv$arm], ARM)
cv$D <- factor(paste0("threshold D = ", cv$threshold))
p <- ggplot(cv, aes(year, turnover, colour = b, linetype = arm)) +
  geom_vline(xintercept = 10, colour = "grey55", linetype = "dashed") +
  geom_line(linewidth = 0.7) + facet_wrap(~ D) +
  scale_y_log10() + scale_x_continuous(breaks = c(1, 5, 10, 20, 30, 40)) +
  scale_colour_manual(values = c("2" = "#9ecae1", "5" = "#6baed6", "10" = "#3182bd", "20" = "#08519c"),
                      name = "niche breadth b") +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(x = "year", y = "change in abundance per species and site") + th
ggsave(fig_path("annex_1_convergence"), p, width = 10, height = 4.6, dpi = 300)

## annex_2: the network of one replicate, species ordered by niche optimum
net <- read.csv(file.path(DIAG, "example_network.csv"))
TYPE <- c("competition", "facilitation", "none")
net$type <- factor(ifelse(net$w > 0, TYPE[1], ifelse(net$w < 0, TYPE[2], TYPE[3])), TYPE)
p <- ggplot(net, aes(j, k, fill = type)) + geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = setNames(c("#D55E00", "#0072B2", "grey93"), TYPE), name = "interaction") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20), expand = c(0, 0)) +
  scale_y_reverse(breaks = c(1, 5, 10, 15, 20), expand = c(0, 0)) +
  coord_equal() +
  labs(x = "species (ordered by niche optimum)", y = NULL) + th +
  theme(panel.grid = element_blank())
ggsave(fig_path("annex_2_reseau"), p, width = 5.2, height = 5.2, dpi = 300)

## annex_3: its counts under the three mechanisms, log(1 + count)
ex <- read.csv(file.path(DIAG, "example_counts.csv"))
ex$mech <- mech_factor(ex$mech); ex$arm <- factor(ARM[ex$arm], ARM)
p <- ggplot(ex, aes(species, site, fill = log1p(count))) + geom_raster() +
  facet_grid(arm ~ mech, labeller = labeller(.rows = label_wrap_gen(14))) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", name = "log(1 + Y)") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = "species (ordered by niche optimum)", y = "position of the site on the gradient") + th +
  theme(panel.grid = element_blank())
ggsave(fig_path("annex_3_donnees"), p, width = 10, height = 6, dpi = 300)
