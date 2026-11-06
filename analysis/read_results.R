## Reading the grids, shared by paper_numbers.R and the figures, so that a
## number in the text and a figure can never disagree.
##
## Every score is read as an increment: the score of an estimator on a
## replicate with interactions, minus its score on the paired no-interaction
## control of the same replicate. Estimators are compared each at its best
## abiotic basis in each setting (best mean increment).
source("R/results_io.R")

BASES <- c("no.covar", "cat.2", "quadratic", "poly.3", "spline.2")
PAPER.METRICS <- c("AUC_competition", "AUC_facilitation",
                   "precision_k_competition", "precision_k_facilitation")

## basis of a method name: PLN.cat.2.marginal -> cat.2, plnpca.poly.3.q5 -> poly.3
basis_of <- function(method) {
  b <- sub("\\.q[0-9]+$", "", sub("\\.?marginal$", "", sub("^(PLN|gllvm|empirical|plnpca)\\.?", "", method)))
  b[b == ""] <- "continuous.covar"
  b
}

load_grid <- function(grid, campaign = "world") {
  f <- grid_file(campaign, grid)
  if (!file.exists(f)) stop("no such grid: ", f)
  d <- utils::read.csv(f)
  if (!"dynamics" %in% names(d)) d$dynamics <- "lottery"
  d$fam       <- sub("\\..*$", "", d$method)
  d$marginal  <- grepl("\\.marginal$", d$method)
  d$spec      <- basis_of(d$method)
  d$breadth   <- d$niche.breadth
  d$threshold <- d$interaction.threshold
  d$world     <- sprintf("b = %g, threshold %g", d$breadth, d$threshold)
  d
}

## The estimators: rows of a grid, and suffix of the column read ("_alt": the
## graphical-lasso reading of a covariance). quick is the log-GGM.
ESTIMATORS <- list(plnnetwork = list(fam = "PLN",       marginal = FALSE, alt = ""),
                   plnglasso  = list(fam = "PLN",       marginal = TRUE,  alt = "_alt"),
                   gllvm      = list(fam = "gllvm",     marginal = FALSE, alt = "_alt"),
                   quick      = list(fam = "empirical", marginal = FALSE, alt = "_alt"))
COLOURS <- c(plnnetwork = "#0072B2", plnglasso = "#56B4E9", gllvm = "#D55E00", quick = "#CC79A7")

## Paired increment, one row per (setting, replicate, method).
increments <- function(d, value, by = character()) {
  d$.val <- value
  keys <- c("world", "breadth", "threshold", by)
  k <- function(x) do.call(paste, c(x[keys], list(x$simu, x$method), sep = "|"))
  a <- d[d$arm == "interactions", ]
  b <- d[d$arm == "no.interaction", ]
  m <- match(k(a), k(b))
  ok <- !is.na(m) & is.finite(a$.val) & is.finite(b$.val[m])
  out <- a[ok, c(keys, "method", "fam", "spec", "marginal", "simu")]
  out$inc <- a$.val[ok] - b$.val[m[ok]]
  out
}

## Mean increment of every estimator at every basis, in each cell of a grid.
estimator_cells <- function(grid, metric, interaction, by = character(),
                            campaign = "world", bases = BASES) {
  d <- load_grid(grid, campaign)
  d <- d[d$fam %in% c("PLN", "gllvm", "empirical") & d$spec %in% bases, ]
  keys <- c("world", "breadth", "threshold", by)
  do.call(rbind, lapply(names(ESTIMATORS), function(e) {
    s <- ESTIMATORS[[e]]
    x <- d[d$fam == s$fam & d$marginal == s$marginal, ]
    col <- paste0(metric, "_", interaction, s$alt)
    if (!nrow(x) || !col %in% names(x)) return(NULL)
    i <- increments(x, x[[col]], by = by)
    if (!nrow(i)) return(NULL)
    a <- aggregate(inc ~ ., data = i[c(keys, "spec", "inc")], FUN = mean)
    a$n <- aggregate(inc ~ ., data = i[c(keys, "spec", "inc")], FUN = length)$inc
    a$estimator <- e
    a
  }))
}

## Replicate-level increments of an estimator at its best basis in each cell
## of a grid; `by` are the columns that make a cell.
best_reps <- function(grid, metric, interaction, by, who, bases = BASES, campaign = "world") {
  y <- estimator_cells(grid, metric, interaction, by = by, campaign = campaign, bases = bases)
  y$cell <- do.call(paste, c(y[by], sep = "/"))
  d <- load_grid(grid, campaign)
  do.call(rbind, lapply(unique(y$cell), function(cl) {
    z <- y[y$cell == cl & y$estimator %in% who, ]; b <- z[which.max(z$inc), ]
    s <- ESTIMATORS[[b$estimator]]
    x <- d[d$fam == s$fam & d$marginal == s$marginal & d$spec == b$spec, ]
    i <- increments(x, x[[paste0(metric, "_", interaction, s$alt)]], by = by)
    i <- i[do.call(paste, c(i[by], sep = "/")) == cl, ]
    i$estimator <- b$estimator; i$base <- b$spec; i }))
}

## The raw score of both arms at the basis best_reps() retains.
best_raw <- function(grid, metric, interaction, by, who, bases = BASES, campaign = "world") {
  y <- estimator_cells(grid, metric, interaction, by = by, campaign = campaign, bases = bases)
  y$cell <- do.call(paste, c(y[by], sep = "/"))
  d <- load_grid(grid, campaign)
  d$cell <- do.call(paste, c(d[by], sep = "/"))
  do.call(rbind, lapply(unique(y$cell), function(cl) {
    z <- y[y$cell == cl & y$estimator %in% who, ]; b <- z[which.max(z$inc), ]
    s <- ESTIMATORS[[b$estimator]]
    x <- d[d$fam == s$fam & d$marginal == s$marginal & d$spec == b$spec & d$cell == cl, ]
    x$val <- x[[paste0(metric, "_", interaction, s$alt)]]
    x[c(by, "simu", "arm", "val")] }))
}

## The two grids of the main design, and the columns that make a setting.
MAIN.GRIDS <- list(lottery  = list(grid = "world_mech_smooth_lottery",  by = "competition.intensity"),
                   dynamics = list(grid = "world_mech_smooth_dynamics", by = c("dynamics", "dyn.strength", "dyn.sigma")))

## The four levels of integration on the main design: one row per (setting,
## estimator, metric), each estimator at its best basis, with its increment and
## its raw score (simulation with interactions). Estimators: quick (log-GGM),
## pca2, pca5, pca10 (PLNPCA + glasso), gllvm (+ glasso), plnglasso, plnnetwork.
ladder_cells <- function(campaign = "world") {
  long <- function(grp) {
    prep <- function(d) {
      if (grp == "lottery") { d$mech <- "lottery"; d$force <- d$competition.intensity; d$sigma <- 0 }
      else { d$mech <- d$dynamics; d$force <- d$dyn.strength; d$sigma <- d$dyn.sigma }
      d$base <- basis_of(d$method)
      d }
    s <- prep(load_grid(paste0("world_mech_smooth_", grp), campaign))
    p <- prep(load_grid(paste0("world_mech_pca_", grp), campaign))
    est <- list(quick      = list(d = s, sel = s$fam == "empirical", alt = "_alt"),
                pca2       = list(d = p, sel = grepl("^plnpca.*\\.q2$", p$method), alt = "_alt"),
                pca5       = list(d = p, sel = grepl("^plnpca.*\\.q5$", p$method), alt = "_alt"),
                pca10      = list(d = p, sel = grepl("^plnpca.*\\.q10$", p$method), alt = "_alt"),
                gllvm      = list(d = s, sel = s$fam == "gllvm", alt = "_alt"),
                plnglasso  = list(d = s, sel = s$fam == "PLN" & s$marginal, alt = "_alt"),
                plnnetwork = list(d = s, sel = s$fam == "PLN" & !s$marginal, alt = ""))
    do.call(rbind, lapply(names(est), function(e) { x <- est[[e]]$d[est[[e]]$sel, ]
      out <- x[, c("mech", "force", "sigma", "breadth", "threshold", "simu", "arm", "base")]
      for (m in PAPER.METRICS) out[[m]] <- x[[paste0(m, est[[e]]$alt)]]
      out$estimator <- e; out })) }
  L <- rbind(long("lottery"), long("dynamics"))
  L$setting <- paste(L$mech, L$force, L$sigma, L$breadth, L$threshold)
  k <- function(d) paste(d$setting, d$base, d$simu, d$estimator)
  A <- L[L$arm == "interactions", ]; B <- L[L$arm == "no.interaction", ]; mB <- match(k(A), k(B))
  do.call(rbind, lapply(PAPER.METRICS, function(m) {
    x <- data.frame(A[c("setting", "mech", "breadth", "threshold", "force", "sigma", "base", "estimator")],
                    raw = A[[m]], inc = A[[m]] - B[[m]][mB])
    x <- x[is.finite(x$raw) & is.finite(x$inc), ]
    a <- aggregate(cbind(raw, inc) ~ setting + mech + breadth + threshold + force + sigma + base + estimator, x, mean)
    b <- do.call(rbind, lapply(split(a, list(a$setting, a$estimator), drop = TRUE), function(s) s[which.max(s$inc), ]))
    b$metric <- m; b }))
}

## What each step between levels brings: difference of two estimators'
## best-basis increments, per setting and metric.
LADDER.STEPS <- list(count = c("plnglasso", "quick"), rank = c("plnglasso", "pca10"),
                     integration = c("plnnetwork", "plnglasso"), total = c("plnnetwork", "quick"))
ladder_steps <- function(C, steps = LADDER.STEPS) {
  do.call(rbind, lapply(names(steps), function(s) {
    hi <- C[C$estimator == steps[[s]][1], ]; lo <- C[C$estimator == steps[[s]][2], ]
    hi$d <- hi$inc - lo$inc[match(paste(hi$setting, hi$metric), paste(lo$setting, lo$metric))]
    hi$step <- s; hi[c("setting", "mech", "breadth", "threshold", "metric", "step", "d")] }))
}

## Sampling effort and compositional data: each reading at its best basis on
## the raw score of the simulation with interactions, in each cell; raw score
## of both arms, increment, and the drift of the control from its own value at
## effort spread 0 (the spurious signal).
COMPO.READINGS <- list(
  pln.none   = list(re = "^PLN\\.(poly\\.3|spline\\.2)$",             alt = "",     family = "pln",   correction = "none"),
  pln.off    = list(re = "^PLN\\.(poly\\.3|spline\\.2)\\.off$",       alt = "",     family = "pln",   correction = "known"),
  pln.tss    = list(re = "^PLN\\.(poly\\.3|spline\\.2)\\.tss$",       alt = "",     family = "pln",   correction = "unknown"),
  quick.none = list(re = "^empirical\\.(poly\\.3|spline\\.2)$",       alt = "_alt", family = "quick", correction = "none"),
  quick.div  = list(re = "^empirical\\.(poly\\.3|spline\\.2)\\.div$", alt = "_alt", family = "quick", correction = "known"),
  quick.tss  = list(re = "^empirical\\.(poly\\.3|spline\\.2)\\.tss$", alt = "_alt", family = "quick", correction = "unknown"),
  quick.clr  = list(re = "^empirical\\.(poly\\.3|spline\\.2)\\.clr$", alt = "_alt", family = "quick", correction = "clr"))
compo_cells <- function(campaign = "world", mechs = c("lottery", "gompertz")) {
  L <- do.call(rbind, lapply(mechs, function(mech) {
    d <- load_grid(paste0("world_compo_", mech), campaign)
    do.call(rbind, lapply(names(COMPO.READINGS), function(r) {
      s <- COMPO.READINGS[[r]]; x <- d[grepl(s$re, d$method), ]
      x$base <- sub(s$re, "\\1", x$method)
      out <- x[, c("niche.breadth", "interaction.threshold", "effort.sd", "simu", "arm", "base")]
      for (m in PAPER.METRICS) out[[m]] <- x[[paste0(m, s$alt)]]
      out$reading <- r; out$family <- s$family; out$correction <- s$correction; out$mech <- mech
      out })) }))
  L$cell <- paste(L$mech, L$niche.breadth, L$interaction.threshold, L$effort.sd)
  do.call(rbind, lapply(PAPER.METRICS, function(m) {
    x <- L; x$val <- x[[m]]
    a <- x[x$arm == "interactions", ]
    cm <- aggregate(val ~ cell + reading + base, a, mean)
    best <- do.call(rbind, lapply(split(cm, list(cm$cell, cm$reading), drop = TRUE), function(s) s[which.max(s$val), ]))
    x <- x[paste(x$cell, x$reading, x$base) %in% paste(best$cell, best$reading, best$base), ]
    k <- function(z) paste(z$cell, z$reading, z$simu)
    A <- x[x$arm == "interactions", ]; N <- x[x$arm == "no.interaction", ]
    A$null <- N$val[match(k(A), k(N))]; A$inc <- A$val - A$null; A$metric <- m
    z0 <- A[A$effort.sd == 0, ]
    key <- function(z) paste(z$mech, z$niche.breadth, z$interaction.threshold, z$reading, z$simu)
    A$drift <- A$null - z0$null[match(key(A), key(z0))]
    A[c("mech", "niche.breadth", "interaction.threshold", "effort.sd", "simu", "base", "reading",
        "family", "correction", "metric", "val", "null", "inc", "drift")] }))
}

## Error on the covariate in niche breadths, from its correlation with the
## gradient (b = 10, 50 sites on [0, 100]).
noise_of <- function(rho) sqrt(1 - rho^2) * sd(seq(0, 100, length.out = 50)) / 10

## Partial against marginal: the same PLN fit and basis read as PLNnetwork and
## as a thresholded covariance, paired by replicate.
partial_marginal <- function(bases = BASES, campaign = "world") {
  do.call(rbind, lapply(PAPER.METRICS, function(m) do.call(rbind, lapply(MAIN.GRIDS, function(g) {
    d0 <- load_grid(g$grid, campaign)
    do.call(rbind, lapply(bases, function(b) {
      mth <- paste0("PLN.", b)
      d <- d0[d0$method %in% c(mth, paste0(mth, ".marginal")), ]
      i <- increments(d, d[[m]], by = g$by)
      if (!"dynamics" %in% names(i)) i$dynamics <- "lottery"
      k <- function(x) do.call(paste, c(x[c("breadth", "threshold", g$by, "simu")], sep = "|"))
      a <- i[i$method == mth, ]; z <- i[i$method == paste0(mth, ".marginal"), ]
      a$d <- a$inc - z$inc[match(k(a), k(z))]
      data.frame(mech = a$dynamics, breadth = a$breadth, threshold = a$threshold,
                 strength = do.call(paste, c(a[g$by], sep = "/")), base = b, metric = m, d = a$d) })) }))))
}

## Paired effect of the mechanism: PLNnetwork's advantage over the log-GGM,
## each at its best basis in each cell, under a dynamics minus under the
## lottery at intensity 2, on the same world, network and replicate.
mechanism_effect <- function(bases = BASES, campaign = "world", who = c("plnnetwork", "quick")) {
  rep_best <- function(g, metric, interaction, est) {
    y <- estimator_cells(g$grid, metric, interaction, by = g$by, campaign = campaign, bases = bases)
    keys <- c("world", g$by)
    y$cell <- do.call(paste, c(y[keys], sep = "|"))
    y <- y[y$estimator == est, ]
    best <- do.call(rbind, lapply(split(y, y$cell), function(z) z[which.max(z$inc), ]))
    s <- ESTIMATORS[[est]]
    d <- load_grid(g$grid, campaign)
    d <- d[d$fam == s$fam & d$marginal == s$marginal & d$spec %in% bases, ]
    i <- increments(d, d[[paste0(metric, "_", interaction, s$alt)]], by = g$by)
    i$cell <- do.call(paste, c(i[keys], sep = "|"))
    i <- i[paste(i$cell, i$spec) %in% paste(best$cell, best$spec), ]
    if (!"dynamics" %in% names(i)) i$dynamics <- "lottery"
    i }
  do.call(rbind, lapply(PAPER.METRICS, function(m) {
    mi <- strsplit(m, "_(?=[a-z]+$)", perl = TRUE)[[1]]
    gap <- do.call(rbind, lapply(MAIN.GRIDS, function(g) {
      a <- rep_best(g, mi[1], mi[2], who[1]); b <- rep_best(g, mi[1], mi[2], who[2])
      k <- function(x) paste(x$cell, x$simu)
      a$gap <- a$inc - b$inc[match(k(a), k(b))]
      a$strength <- if (is.null(a$competition.intensity)) NA else a$competition.intensity
      a[c("dynamics", "breadth", "threshold", "strength", "simu", "gap")] }))
    lot <- gap[gap$dynamics == "lottery" & gap$strength %in% 2, ]
    dy <- gap[gap$dynamics != "lottery", ]
    dy$d <- dy$gap - lot$gap[match(paste(dy$breadth, dy$threshold, dy$simu), paste(lot$breadth, lot$threshold, lot$simu))]
    data.frame(mech = dy$dynamics, breadth = dy$breadth, threshold = dy$threshold, simu = dy$simu, metric = m, d = dy$d) }))
}
