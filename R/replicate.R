## A replicate: one world simulated twice from the same random state, with and
## without interactions (the paired no-interaction control), every estimator
## fitted and scored on both. A grid: replicates over a table of settings.

## The abiotic bases, as each estimator takes them. Names are those of the
## result files: cat.2 is categorical with 10 classes, spline.2 a B-spline with
## 10 columns (intercept included), both.quadratic a quadratic on the observed
## and on the hidden gradient.
BASES <- c("no.covar", "cat.2", "quadratic", "poly.3", "spline.2")

abiotic_basis <- function(basis, data) {
  counts <- data$counts; x <- data$env.observed; n <- length(x)
  intercept <- data.frame(Intercept = rep(1, n))
  poly_df <- function(v, degree, prefix) {
    P <- poly(v, degree = degree); colnames(P) <- paste0(prefix, seq_len(degree)); as.data.frame(P) }
  covariates <- function(C)
    list(pln = list(formula = paste("Abundance ~ 1 +", paste(names(C), collapse = " + ")),
                    data = prepare_data(counts = counts, covariates = C)),
         gllvm = list(X = cbind(intercept, C),
                      formula = paste("~ Intercept +", paste(names(C), collapse = " + "))),
         design = cbind("(Intercept)" = 1, as.matrix(C)))
  switch(basis,
    no.covar = list(pln = list(formula = "Abundance ~ 1",
                               data = prepare_data(counts = counts, covariates = matrix(x))),
                    gllvm = list(X = intercept, formula = "~ Intercept"),
                    design = matrix(1, n, 1)),
    cat.2 = {
      classes <- factor(cut(x, 10, labels = 1:10))
      list(pln = list(formula = "Abundance ~ 1 + V1",
                      data = prepare_data(counts = counts, covariates = matrix(classes))),
           gllvm = list(X = cbind(intercept, V1 = classes), formula = "~ Intercept + factor(V1)"),
           design = model.matrix(~ classes)) },
    quadratic = covariates(poly_df(x, 2, "Q")),
    poly.3    = covariates(poly_df(x, 3, "Q")),
    spline.2  = {
      ## orthonormalised B-spline basis: same model, well conditioned
      Q <- qr.Q(qr(cbind(1, splines::bs(x, df = 9))))[, -1, drop = FALSE]
      colnames(Q) <- paste0("S", seq_len(ncol(Q)))
      covariates(as.data.frame(Q)) },
    both.quadratic = {
      H <- poly(data$hidden.observed, degree = 2); colnames(H) <- c("H1", "H2")
      covariates(cbind(poly_df(x, 2, "Q"), as.data.frame(unclass(H))[, c("H1", "H2")])) })
}

## The models fitted on a data set, in the order of the published grids (the
## order matters: see in_fork()). `type` says how a fit is read.
model_specs <- function(data, opt) {
  usable <- if (is.null(data$hidden.observed)) setdiff(opt$bases, "both.quadratic") else opt$bases
  bases <- lapply(setNames(nm = usable), abiotic_basis, data = data)
  main  <- intersect(BASES, opt$bases)
  spec  <- function(method, type, fit, design = NULL) list(method = method, type = type, fit = fit, design = design)
  pln   <- function(b) spec(paste0("PLN.", b), "pln",
                            function() fit_pln(bases[[b]]$pln$formula, bases[[b]]$pln$data, paste0("PLN.", b)),
                            bases[[b]]$design)
  gllvm <- function(b) spec(paste0("gllvm.", b), "gllvm",
                            function() fit_gllvm(data$counts, bases[[b]]$gllvm$X, bases[[b]]$gllvm$formula, opt$gllvm))
  loggm <- function(b) spec(paste0("empirical.", b), "empirical",
                            function() fit_loggm(data$counts, bases[[b]]$design))
  specs <- c(if ("pln" %in% opt$families) lapply(main, pln),
             if ("gllvm" %in% opt$families) lapply(main, gllvm),
             if ("empirical" %in% opt$families) c(
               lapply(main, loggm),
               unlist(lapply(opt$pca.ranks, function(q) lapply(main, function(b)
                 spec(paste0("plnpca.", b, ".q", q), "empirical",
                      function() fit_plnpca(data$counts, bases[[b]]$design, q)))), recursive = FALSE)))
  ## sampling effort: PLN with the offset log(effort) (known) or log(site total)
  ## (unknown); the log-GGM on Y / effort, on counts scaled to the median site
  ## total, or on the centred log-ratio
  if (opt$effort) {
    total <- pmax(rowSums(data$counts), 1)
    effort.specs <- unlist(lapply(main, function(b) { X <- bases[[b]]$design; list(
      spec(paste0("PLN.", b, ".off"), "pln",
           function() fit_pln_offset(data$counts, X, log(data$effort), paste0("PLN.", b, ".off")), X),
      spec(paste0("PLN.", b, ".tss"), "pln",
           function() fit_pln_offset(data$counts, X, log(total), paste0("PLN.", b, ".tss")), X),
      spec(paste0("empirical.", b, ".div"), "empirical", function() fit_loggm(data$counts / data$effort, X)),
      spec(paste0("empirical.", b, ".tss"), "empirical",
           function() fit_loggm(data$counts * stats::median(total) / total, X)),
      spec(paste0("empirical.", b, ".clr"), "empirical", function() fit_loggm(data$counts, X, clr = TRUE))) }),
      recursive = FALSE)
    specs <- c(specs, Filter(function(sp) sp$type %in% c(opt$families), effort.specs))
  }
  ## a hidden gradient given to the models with the right shape
  if ("both.quadratic" %in% usable)
    specs <- c(specs, Filter(function(sp) sp$type %in% opt$families,
                             list(pln("both.quadratic"), gllvm("both.quadratic"), loggm("both.quadratic"))))
  specs
}

## Fits and scores every model on a data set. PLN is read twice: PLNnetwork's
## path, and the unpenalised fit (marginal and PLN + glasso, method ".marginal").
fit_and_score <- function(data, s, opt) {
  rows <- lapply(model_specs(data, opt), function(sp) {
    t0  <- Sys.time()
    fit <- tryCatch(in_fork(sp$fit, opt$fit.timeout), error = function(e) {
      message("  ", sp$method, ": fit failed (", conditionMessage(e), ")"); NULL })
    if (!is.null(fit) && !fit_is_usable(fit, sp$type)) fit <- NULL
    fit.time <- Sys.time() - t0
    score <- function(method, model, type, abiotic = NO.ABIOTIC)
      tryCatch(in_fork(function() evaluate_model(method, model, type, data$nets, fit.time, abiotic),
                       opt$fit.timeout),
               error = function(e) {
                 message("  ", method, ": scoring failed (", conditionMessage(e), ")")
                 data.frame(method = method, time = as.numeric(fit.time, units = "secs")) })
    if (is.null(fit)) {
      empty <- function(m) data.frame(method = m, time = as.numeric(fit.time, units = "secs"))
      return(if (sp$type == "pln") list(empty(sp$method), empty(paste0(sp$method, ".marginal")))
             else list(empty(sp$method)))
    }
    if (sp$type != "pln") return(list(score(sp$method, fit, sp$type)))
    abiotic <- abiotic_quality(fit$inception, sp$design, data$optima, data$env, s$niche.breadth)
    list(score(sp$method, fit$network, "pln", abiotic),
         score(paste0(sp$method, ".marginal"), fit$inception, "pln.marginal", abiotic))
  })
  c(unlist(rows, recursive = FALSE), reference_rows(data, s))
}

## The seed of a setting: a hash of its parameters. The key lists them as the
## code of the published grids did, including parameters that no longer vary
## (fixed values), so that every replicate draws the same world as in those
## grids. The noise on the covariate and the abundance level K are applied
## after the world is drawn, and are keyed at their reference values.
setting_seed <- function(s) {
  key <- list(add.intercept = TRUE, beta.abun = s$beta.abun, beta.comp = s$beta.comp,
              beta.env = 1, beta.hidden = s$beta.hidden, competition = "symmetric",
              competition.density = s$competition.density,
              competition.intensity = s$competition.intensity,
              env.max = 100, env.min = 0, env.observed.rho = 1,
              facilitation.density = s$facilitation.density,
              facilitation.intensity = s$facilitation.intensity, input.network = NA,
              K = 100, n.env = s$n.env, n.env.categories = 5, n.species = s$n.species,
              network.mode = s$network.mode, network.rank = s$network.rank,
              niche.breadth = s$niche.breadth, niche.mean = NA, niche.mode = "pool",
              niche.sd = NA, sp.ratio = 1, sp.ratio.type = "abundance",
              threshold.niche.dist = s$interaction.threshold, years = s$years)
  key <- paste(vapply(key, function(v) as.character(v), character(1)), collapse = "|")
  h <- 0
  for (b in utf8ToInt(key)) h <- (h * 31 + b) %% 2147483647
  h
}

## One replicate. The no-interaction arm replays the random state of the start
## of the replicate with beta.comp = 0, and keeps the species of the other arm.
run_replicate <- function(s, i, opt) {
  seed <- setting_seed(s)
  M <- 2147483647
  covariate.seed <- (seed + i * 7919 + 104729) %% M
  effort.seed    <- if (opt$effort) (seed + i * 7919 + 15485863) %% M
  one_arm <- function(beta.comp, arm, forced = NULL) {
    data <- simulate_data(s, beta.comp, covariate.seed, effort.seed, forced)
    rows <- fit_and_score(data, s, opt)
    info <- c(list(simu = i, arm = arm), s[SETTING.COLUMNS], beta.comp = beta.comp,
              list(n.selected.species = length(data$selected), species.dropped = data$species.dropped,
                   network.signature = data$network.signature, species.signature = data$species.signature))
    list(res = plyr::rbind.fill(lapply(rows, function(r) as.data.frame(c(info, as.list(r))))),
         rai = data$rai, selected = data$selected)
  }
  set.seed((seed + i * 7919) %% M)
  start <- get(".Random.seed", envir = globalenv())
  with  <- one_arm(s$beta.comp, "interactions")
  after <- get(".Random.seed", envir = globalenv())
  assign(".Random.seed", start, envir = globalenv())
  without <- tryCatch(one_arm(0, "no.interaction", with$selected), error = function(e) {
    message("  replicate ", i, ": no-interaction arm failed (", conditionMessage(e), ")"); NULL })
  assign(".Random.seed", after, envir = globalenv())
  list(res = plyr::rbind.fill(with$res, without$res), rai = with$rai)
}

## The data of both arms of a replicate, without fitting anything (the
## no-interaction arm does not depend on the fits). `...` goes to simulate_data().
replicate_data <- function(s, i, ...) {
  seed <- setting_seed(s)
  set.seed((seed + i * 7919) %% 2147483647)
  start <- get(".Random.seed", envir = globalenv())
  with <- simulate_data(s, s$beta.comp, (seed + i * 7919 + 104729) %% 2147483647, ...)
  assign(".Random.seed", start, envir = globalenv())
  without <- simulate_data(s, 0, (seed + i * 7919 + 104729) %% 2147483647,
                           forced.species = with$selected, ...)
  list(interactions = with, control = without)
}

## Columns that describe a setting, in the result files.
SETTING.COLUMNS <- c("dynamics", "dyn.strength", "dyn.sigma", "competition.intensity",
                     "facilitation.intensity", "niche.breadth", "interaction.threshold",
                     "n.env", "n.species", "K", "network.mode", "network.rank",
                     "env.observed.rho", "beta.hidden", "effort.sd")

## A table of settings: every combination of the values given, the other
## parameters at their reference value. Facilitation mirrors competition.
settings <- function(...) {
  g <- expand.grid(..., KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  ref <- list(dynamics = "lottery", competition.intensity = 2, dyn.strength = NA_real_,
              dyn.sigma = NA_real_, niche.breadth = 10, interaction.threshold = 10,
              n.env = 50, n.species = 20, K = 100, years = 10, network.mode = "density.based",
              network.rank = 2, env.observed.rho = 1, beta.hidden = 0, effort.sd = 0,
              competition.density = 0.2, facilitation.density = 0.1, beta.comp = 1, beta.abun = 1)
  for (f in setdiff(names(ref), names(g))) g[[f]] <- ref[[f]]
  g$facilitation.intensity <- -g$competition.intensity
  g
}

## Runs n.simu replicates of every setting, settings in parallel, then
## replicates within a setting on the remaining cores. Each finished setting is
## also written to `checkpoint` (a file prefix), in case the run is stopped.
run_grid <- function(grid, n.simu, opt, cores = 1, checkpoint = NULL) {
  opt <- modifyList(list(families = c("pln", "gllvm", "empirical"), bases = BASES, pca.ranks = NULL,
                         effort = FALSE, gllvm = list(), fit.timeout = 600), opt)
  grid.cores <- max(1L, min(cores, nrow(grid)))
  rep.cores  <- max(1L, cores %/% grid.cores)
  message(sprintf("%d settings on %d core(s), %d replicates on %d core(s) each",
                  nrow(grid), grid.cores, n.simu, rep.cores))
  out <- parallel::mclapply(seq_len(nrow(grid)), function(k) {
    s <- as.list(grid[k, ])
    reps <- parallel::mclapply(seq_len(n.simu), function(i) run_replicate(s, i, opt), mc.cores = rep.cores)
    reps <- Filter(function(r) is.list(r) && is.data.frame(r$res), reps)
    res  <- plyr::rbind.fill(lapply(reps, `[[`, "res"))
    rai  <- cbind(as.data.frame(s[SETTING.COLUMNS]), threshold.niche.dist = s$interaction.threshold,
                  rai_summary(do.call(rbind, lapply(reps, `[[`, "rai"))))
    if (!is.null(checkpoint)) {
      dir.create(dirname(checkpoint), showWarnings = FALSE, recursive = TRUE)
      utils::write.csv(res, sprintf("%s__%s.csv", checkpoint, setting_key(as.data.frame(s), file = TRUE)),
                       row.names = FALSE)
    }
    list(res = res, rai = rai)
  }, mc.cores = grid.cores, mc.preschedule = FALSE)
  failed <- !vapply(out, function(x) is.list(x) && is.data.frame(x$res), logical(1))
  if (any(failed)) warning(sum(failed), " setting(s) failed", call. = FALSE)
  list(res = plyr::rbind.fill(lapply(out[!failed], `[[`, "res")),
       rai = do.call(rbind, lapply(out[!failed], `[[`, "rai")))
}
