## Convergence of Hmsc (Supplementary material): four independent chains on
## each of 30 simulated data sets (lottery, b = 2, threshold 10), at thin 10 and
## with chains five times longer (thin 50). Records the potential scale
## reduction factor of the associations and of the abiotic coefficients, and
## the AUC of each chain. From the root of the repository:
##   VC_CORES=16 Rscript hmsc_convergence.R
## writes results/hmsc_convergence/hmsc_convergence{,_chains}.csv.
source("R/load.R")
suppressMessages({ library(Hmsc); library(coda) })
n.sim   <- as.integer(Sys.getenv("VC_HMSC_SIMS", "30"))
n.chain <- 4
cores   <- as.integer(Sys.getenv("VC_CORES", "16"))
thins   <- c(10, 50)

sim_data <- function(seed) {
  set.seed(seed)
  s <- as.list(settings(niche.breadth = 2)[1, ])
  optima <- draw_niche_optima(20)
  env <- seq(0, 100, length.out = 50)
  W <- build_network(optima, s); diag(W) <- 0
  Y <- simulate_lottery(optima, env, W, 2, K = 100, years = 10)$counts
  present <- which(colSums(Y) > 0)
  counts <- Y[, present]
  drop <- remove_proportional_sites(counts)
  if (length(drop)) { counts <- counts[-drop, ]; env <- env[-drop] }
  list(counts = counts, env = env, W = W[present, present])
}
## the data sets were drawn with R's default generator, the chains with
## parallel streams
RNGkind("Mersenne-Twister")
datasets <- lapply(seq_len(n.sim), function(i) sim_data(20270101 + i))
RNGkind("L'Ecuyer-CMRG")
tasks <- expand.grid(chain = seq_len(n.chain), sim = seq_len(n.sim), thin = thins)
tasks <- tasks[order(-tasks$thin), ]            # long chains first

## AUC of the posterior mean association, thresholded by its posterior support
hmsc_auc <- function(W, model, sign) {
  a <- computeAssociations(model)[[1]]
  est <- lapply(0.01 * 100:1, function(prob) ((a$support > prob) + (a$support < (1 - prob)) > 0) * a$mean)
  roc_auc(roc_points(W, est, sign))
}
fit_chain <- function(i) {
  tk <- tasks[i, ]; d <- datasets[[tk$sim]]
  design <- data.frame(sample = as.factor(seq_len(nrow(d$counts))))
  m <- Hmsc(Y = d$counts, XData = data.frame(V1 = d$env), XFormula = ~ 1 + V1, studyDesign = design,
            ranLevels = list(sample = HmscRandomLevel(units = design$sample)), distr = "lognormal poisson")
  invisible(capture.output(m <- sampleMcmc(m, thin = tk$thin, samples = 1000, transient = 500 * tk$thin,
                                           nChains = 1, verbose = 0), type = "message"))
  post <- convertToCodaObject(m)
  list(sim = tk$sim, thin = tk$thin, chain = tk$chain,
       AUC = hmsc_auc(d$W, m, "all"), AUC_comp = hmsc_auc(d$W, m, "negative"),
       beta = post$Beta[[1]], omega = post$Omega[[1]][[1]])
}
set.seed(20270101)
res <- Filter(is.list, parallel::mclapply(seq_len(nrow(tasks)), fit_chain, mc.cores = cores,
                                          mc.preschedule = FALSE))

chains <- do.call(rbind, lapply(res, function(x) data.frame(sim = x$sim, thin = x$thin, chain = x$chain,
                                                             AUC = x$AUC, AUC_comp = x$AUC_comp)))
summary <- do.call(rbind, lapply(thins, function(th) do.call(rbind, lapply(seq_len(n.sim), function(s) {
  r <- Filter(function(x) x$thin == th && x$sim == s, res)
  if (length(r) < 2) return(NULL)
  rhat <- function(field) tryCatch(gelman.diag(mcmc.list(lapply(r, `[[`, field)), multivariate = FALSE)$psrf[, 1],
                                   error = function(e) NA_real_)
  po <- rhat("omega"); pb <- rhat("beta"); auc <- vapply(r, `[[`, numeric(1), "AUC")
  data.frame(thin = th, sim = s, n.chains = length(r),
             chains.distinct = !any(duplicated(vapply(r, function(x) sum(as.matrix(x$omega)), numeric(1)))),
             rhat.omega.median = median(po, na.rm = TRUE), rhat.omega.max = max(po, na.rm = TRUE),
             rhat.omega.above = mean(po > 1.05, na.rm = TRUE),
             rhat.beta.median = median(pb, na.rm = TRUE), rhat.beta.above = mean(pb > 1.05, na.rm = TRUE),
             auc.mean = mean(auc), auc.range = diff(range(auc)), auc.sd = sd(auc)) }))))
out <- file.path(RESULTS, "hmsc_convergence")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(chains,  file.path(out, "hmsc_convergence_chains.csv"), row.names = FALSE)
write.csv(summary, file.path(out, "hmsc_convergence.csv"), row.names = FALSE)
