## The simulation campaigns of the paper. From the root of the repository:
##
##   Rscript run_simulations.R <campaign> [grid ...]
##
## writes results/<campaign>/grid_res_<grid>.csv. Environment variables:
## VC_NSIMU (replicates, 40), VC_CORES (default: the number of cores the grid
## was run with for the paper, see README), VC_RESULTS (results folder).
source("R/load.R")

## Settings: the order of the values matters only for the exact reproduction
## of the gllvm fits (see README), and is the order of the published runs.
mech.lottery  <- settings(niche.breadth = c(2, 5, 10, 20), competition.intensity = c(0.5, 2, 5),
                          interaction.threshold = c(10, 50))
mech.dynamics <- settings(niche.breadth = c(2, 5, 10, 20), interaction.threshold = c(10, 50),
                          dynamics = c("gompertz", "ricker"), dyn.strength = c(0.5, 0.9),
                          dyn.sigma = c(0.3, 0.5, 0.8))
gompertz <- list(dynamics = "gompertz", dyn.strength = 0.5, dyn.sigma = 0.5)
## error on the covariate, in niche breadths (b = 10), as a correlation
noise <- c(0, 0.1, 0.25, 0.5, 1)
rho <- sqrt(1 - (noise * 10 / sd(seq(0, 100, length.out = 50)))^2)
one_world <- function(mech, ...) do.call(settings, c(list(...), if (mech == "gompertz") gompertz))
compo <- function(mech) list(
  one_world(mech, niche.breadth = c(5, 10, 20), effort.sd = c(0, 0.25, 0.5, 1)),
  one_world(mech, interaction.threshold = 100, effort.sd = c(0, 0.25, 0.5, 1)))
emission <- function(gllvm, families.lottery = c("gllvm", "empirical"), cores.dynamics = 4) {
  opt <- list(families = c("gllvm", "empirical"), bases = c("cat.2", "poly.3", "spline.2"),
              gllvm = c(gllvm, sd.errors = FALSE))
  list(world_mech_smooth_lottery = list(cores = 4, opt = modifyList(opt, list(families = families.lottery)),
                                        settings = list(settings(niche.breadth = c(5, 10, 20)))),
       world_mech_smooth_dynamics = list(cores = cores.dynamics, opt = opt,
         settings = list(settings(niche.breadth = c(5, 10, 20), dynamics = c("gompertz", "ricker"),
                                  dyn.strength = 0.5, dyn.sigma = 0.5))))
}

## campaign -> grid -> settings (run one table after the other), options, cores
CAMPAIGNS <- list(
  world = list(
    ## main design: three mechanisms, eight configurations, five bases
    world_mech_smooth_lottery  = list(cores = 30, settings = list(mech.lottery)),
    world_mech_smooth_dynamics = list(cores = 30, settings = list(mech.dynamics)),
    ## PLNPCA at ranks 2, 5, 10 on the same worlds
    world_mech_pca_lottery  = list(cores = 30, settings = list(mech.lottery),
                                   opt = list(families = "empirical", pca.ranks = c(2, 5, 10))),
    world_mech_pca_dynamics = list(cores = 30, settings = list(mech.dynamics),
                                   opt = list(families = "empirical", pca.ranks = c(2, 5, 10))),
    ## supplementary experiments, at b = 10, threshold 10
    world_annex_c2p_lottery  = list(cores = 30, settings = list(one_world("lottery", env.observed.rho = rho, beta.hidden = c(0, 1))),
                                    opt = list(bases = c(BASES, "both.quadratic"))),
    world_annex_c2p_gompertz = list(cores = 30, settings = list(one_world("gompertz", env.observed.rho = rho, beta.hidden = c(0, 1))),
                                    opt = list(bases = c(BASES, "both.quadratic"))),
    world_annex_np_smooth_lottery  = list(cores = 30, settings = list(one_world("lottery", n.species = c(20, 40), n.env = c(30, 50, 100)))),
    world_annex_np_smooth_gompertz = list(cores = 30, settings = list(one_world("gompertz", n.species = c(20, 40), n.env = c(30, 50, 100)))),
    world_annex_c4_smooth_lottery  = list(cores = 30, settings = list(one_world("lottery", network.mode = "low.rank", network.rank = c(1, 3, 10)))),
    world_annex_c4_smooth_gompertz = list(cores = 30, settings = list(one_world("gompertz", network.mode = "low.rank", network.rank = c(1, 3, 10)))),
    ## unequal sampling effort and compositional data
    world_compo_lottery  = list(cores = 30, settings = compo("lottery"),
                                opt = list(families = c("pln", "empirical"), bases = c("poly.3", "spline.2"), effort = TRUE)),
    world_compo_gompertz = list(cores = 30, settings = compo("gompertz"),
                                opt = list(families = c("pln", "empirical"), bases = c("poly.3", "spline.2"), effort = TRUE))),
  ## abundance level K under Gompertz
  exploration = list(
    world_abundance_gompertz = list(cores = 14,
      settings = list(settings(K = c(3, 10, 30, 100, 1000), dynamics = "gompertz",
                               dyn.strength = c(0.5, 0.9), dyn.sigma = 0.5)),
      opt = list(families = c("pln", "empirical"), bases = c("cat.2", "quadratic", "poly.3", "spline.2")))),
  ## emission law of gllvm: Poisson or negative binomial, 5 or 10 latent variables
  ## (the lottery grid of emission_pois5 was run with PLN as well)
  emission_pois5 = emission(list(family = "poisson", num.lv = 5),
                            families.lottery = c("pln", "gllvm", "empirical"), cores.dynamics = 14),
  emission_nb5   = emission(list(family = "negative.binomial", num.lv = 5)),
  emission_nb10  = emission(list(family = "negative.binomial", num.lv = 10)))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% names(CAMPAIGNS))
  stop("usage: Rscript run_simulations.R <campaign> [grid ...]; campaigns: ",
       paste(names(CAMPAIGNS), collapse = ", "))
campaign <- args[1]
grids <- CAMPAIGNS[[campaign]]
if (length(args) > 1) grids <- grids[intersect(names(grids), args[-1])]
n.simu <- as.integer(Sys.getenv("VC_NSIMU", "40"))
record_run(campaign, list(grids = names(grids), n.simu = n.simu))
for (name in names(grids)) {
  g <- grids[[name]]
  cores <- as.integer(Sys.getenv("VC_CORES", g$cores))
  for (tab in g$settings) {
    message(name, ": ", nrow(tab), " settings x ", n.simu, " replicates")
    set.seed(20260907)
    res <- run_grid(tab, n.simu, if (is.null(g$opt)) list() else g$opt, cores,
                    checkpoint = file.path(RESULTS, campaign, "checkpoints", name))
    write_grid(res, name, campaign)
  }
}
