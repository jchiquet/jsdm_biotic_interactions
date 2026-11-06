## One simulated data set: species pool, network, communities assembled by one
## of the three mechanisms, then the data as an analyst receives them.
##
## `s` is a setting (one row of a grid, as a list). `covariate.seed` and
## `effort.seed` give the noisy covariate and the sampling effort random
## streams of their own, so that both arms of a replicate, and every level of
## error or effort, share them. The order of the random draws is that of the
## code the published grids were computed with: changing it changes the worlds.
simulate_data <- function(s, beta.comp, covariate.seed = NULL, effort.seed = NULL,
                          forced.species = NULL, years = s$years, cv.measure = FALSE) {
  optima <- draw_niche_optima(s$n.species)
  env    <- seq(0, 100, length.out = s$n.env)
  hidden <- if (s$beta.hidden != 0)
    list(env = sample(env), optima = sample(optima), beta = s$beta.hidden)

  ## assembly
  W <- build_network(optima, s)
  diag(W) <- 0
  sim <- switch(s$dynamics,
    lottery  = simulate_lottery(optima, env, W, s$niche.breadth, s$K, years,
                                beta.comp = beta.comp, beta.abun = s$beta.abun,
                                hidden = hidden, cv.measure = cv.measure),
    gompertz = simulate_gompertz(optima, env, W, s$niche.breadth, s$K, beta.comp,
                                 strength = s$dyn.strength, sigma = s$dyn.sigma, hidden = hidden),
    ricker   = simulate_ricker(optima, env, W, s$niche.breadth, s$K, beta.comp,
                               strength = s$dyn.strength, sigma = s$dyn.sigma, hidden = hidden))
  Y <- sim$counts
  network.signature <- sum(abs(W) * seq_along(W))     # equal in both arms of a replicate
  if (s$network.mode == "low.rank") {
    n <- network_edge_counts(optima, s$competition.density, s$facilitation.density,
                             s$interaction.threshold)
    W <- threshold_to_support(W, n[["n.comp"]], n[["n.facil"]])
  }

  ## unequal sampling effort: each individual is kept with the effort of its
  ## site, log(effort) = effort.sd * z, the largest effort being 1
  effort <- rep(1, s$n.env)
  if (!is.null(effort.seed)) {
    thinned <- with_own_stream(effort.seed, function() {
      e <- exp(s$effort.sd * rnorm(s$n.env)); e <- e / max(e)
      list(effort = e, Y = matrix(rbinom(length(Y), round(Y), rep(e, ncol(Y))), nrow(Y), ncol(Y)))
    })
    effort <- thinned$effort; Y <- thinned$Y
  }

  ## the species present, or in the no-interaction arm those of the other arm
  present <- colSums(Y) > 0
  if (is.null(forced.species)) {
    selected <- which(present)
    species.dropped <- 0L
  } else {
    selected <- as.integer(forced.species)[present[as.integer(forced.species)]]
    species.dropped <- length(forced.species) - length(selected)
  }
  ## sites whose counts are proportional to those of another site are removed
  counts  <- Y[, selected]
  removed <- remove_proportional_sites(counts)
  kept    <- setdiff(seq_len(s$n.env), removed)
  if (length(removed)) counts <- counts[-removed, ]

  ## the covariate the models are given: the gradient, or a noisy proxy of it
  env.observed    <- env[kept]
  hidden.observed <- if (!is.null(hidden)) hidden$env[kept]
  if (s$env.observed.rho < 1) {
    proxy <- function(x, offset) with_own_stream(covariate.seed + offset, function()
      observed_gradient(x, s$env.observed.rho, rnorm(length(x))))
    env.observed <- proxy(env, 0)[kept]
    if (!is.null(hidden)) hidden.observed <- proxy(hidden$env, 1)[kept]
  }

  ## the true networks among the species kept
  W <- W[selected, selected]
  facilitation <- W; facilitation[facilitation > 0 & facilitation != 1] <- 0
  competition  <- W; competition[competition < 0] <- 0
  list(counts = counts, env = env[kept], env.observed = env.observed,
       hidden.observed = hidden.observed, effort = effort[kept],
       optima = optima[selected], selected = as.integer(selected),
       nets = list(all = W, competition = competition, facilitation = facilitation),
       network.signature = network.signature,
       species.signature = paste(sort(as.integer(selected)), collapse = "-"),
       species.dropped = species.dropped,
       rai = rai_pairs(Y[, selected], W, optima[selected], s$interaction.threshold),
       turnover = sim$turnover)
}

## Runs `f` on a random stream seeded by `seed`, then puts the main stream back
## where it was.
with_own_stream <- function(seed, f) {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) stats::runif(1)
  state <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", state, envir = globalenv()))
  set.seed(seed %% 2147483647)
  f()
}

## A proxy of the gradient whose correlation with it is exactly rho: the noise
## `e` is made orthogonal to the gradient, and the proxy is put back on the
## scale of the gradient.
observed_gradient <- function(env, rho, e) {
  z <- as.vector(scale(env))
  e <- e - z * sum(z * e) / sum(z * z)
  x <- rho * z + sqrt(1 - rho^2) * (e / sd(e))
  mean(env) + sd(env) * x
}

## Sites whose count profile is perfectly correlated with that of another site.
remove_proportional_sites <- function(Y) {
  r <- cor(t(Y))
  diag(r) <- 0
  dup <- which(r == 1, arr.ind = TRUE)
  if (!nrow(dup)) return(c())
  dup <- dup[dup[, 1] > dup[, 2], ]
  if (is.null(nrow(dup))) dup[[1]] else dup[, 1]
}

## Relative abundance index (Si-Moussi et al. 2020) of every ordered pair of
## species: mean abundance of j where k is present, minus its mean abundance.
## Pairs are typed by the true network, and pairs without interaction split by
## the distance of their niche optima.
rai_pairs <- function(Y, W, optima, threshold) {
  present <- Y > 0
  storage.mode(present) <- "double"
  n.co <- crossprod(present, present)
  rai  <- crossprod(Y, present) / n.co - rep(colMeans(Y), times = ncol(Y))
  rai[n.co == 0] <- NA
  p <- ncol(rai)
  dist <- outer(1:p, 1:p, Vectorize(function(i, j) abs(optima[[i]] - optima[[j]])))
  d <- data.frame(value = as.vector(rai), P = as.vector(W), dist = as.vector(dist),
                  i = rep(1:p, each = p), j = rep(1:p, times = p))
  d <- d[d$i != d$j, ]
  d$status <- factor(ifelse(d$P > 0, "Competition", ifelse(d$P < 0, "Facilitation",
                     ifelse(d$dist <= threshold, "Neutral.local", "Neutral.distant"))),
                     levels = c("Competition", "Facilitation", "Neutral.local", "Neutral.distant"))
  d[c("status", "value")]
}

## Summary of the index by type of pair, over the replicates of a setting.
rai_summary <- function(pairs) {
  types <- levels(pairs$status)
  s <- lapply(types, function(k) summary(pairs$value[pairs$status == k]))
  data.frame(status = types, Min. = sapply(s, `[`, "Min."), `1st Qu.` = sapply(s, `[`, "1st Qu."),
             Median = sapply(s, `[`, "Median"), Mean = sapply(s, `[`, "Mean"),
             `3rd Qu.` = sapply(s, `[`, "3rd Qu."), Max. = sapply(s, `[`, "Max."),
             check.names = FALSE)
}
