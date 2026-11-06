## The three interaction mechanisms. Each turns a niche (optima, breadth), a
## gradient and a network W (zero diagonal) into a sites x species count
## matrix. `hidden`, when given, is a second environmental driver that no model
## is told about: list(env, optima, beta).

## Lottery of VirtualCom (Münkemüller & Gallien 2015): each year, every
## individual of a site is replaced by a recruit drawn with weights
##   environment^beta.env * (1 - competition load / K)^beta.comp * (1 + beta.abun * abundance).
## Rewritten from the original loop: the competition load is updated in O(p)
## per replacement instead of recomputed. `cv.measure` records, per species and
## year, the absolute change in abundance averaged over the sites.
simulate_lottery <- function(optima, env, W, niche.breadth, K, years,
                             beta.env = 1, beta.comp = 1, beta.abun = 1,
                             hidden = NULL, cv.measure = FALSE) {
  n <- length(env); p <- length(optima)
  sp.names <- names(optima)
  Y  <- matrix(0, n, p)
  Wu <- unname(W)
  cv <- if (cv.measure) lapply(1:n, function(l) matrix(0, p, years))
  for (l in 1:n) {
    log.abiotic <- beta.env * unname(dnorm(optima, mean = env[[l]], sd = niche.breadth, log = TRUE) -
                                     dnorm(env[[l]], mean = env[[l]], sd = niche.breadth, log = TRUE))
    if (!is.null(hidden)) {
      h <- hidden$env[[l]]
      log.abiotic <- log.abiotic + hidden$beta *
        unname(dnorm(hidden$optima, mean = h, sd = niche.breadth, log = TRUE) -
               dnorm(h, mean = h, sd = niche.breadth, log = TRUE))
    }
    code <- match(sample(sp.names, K, replace = TRUE), sp.names)
    abundance <- tabulate(code, p)
    previous  <- abundance
    for (year in 1:years) {
      load <- as.vector(crossprod(abundance, Wu))
      for (i in 1:K) {
        weight <- exp(log.abiotic + beta.comp * log(pmax(1 - load / K, 1e-5)) +
                      log(1 + beta.abun * abundance))
        if (anyNA(weight)) weight[is.na(weight)] <- min(weight, na.rm = TRUE)
        rg <- range(weight)
        if (rg[1L] == rg[2L]) weight <- NULL          # all equal: uniform draw
        out <- sample(K, 1)
        new <- sample.int(p, 1, prob = weight)
        old <- code[out]
        if (new != old) {
          load <- load - Wu[old, ] + Wu[new, ]
          abundance[old] <- abundance[old] - 1L
          abundance[new] <- abundance[new] + 1L
          code[out] <- new
        }
      }
      if (cv.measure) {
        cv[[l]][, year] <- abs(abundance - previous)
        previous <- abundance
      }
    }
    Y[l, ] <- abundance
  }
  list(counts = Y, turnover = if (cv.measure) Reduce("+", cv) / n)
}

## Gaussian niche on the log scale, zero at the optimum (sites x species).
log_niche <- function(env, optima, niche.breadth)
  -outer(env, optima, function(x, th) (x - th)^2) / (2 * niche.breadth^2)

hidden_niche <- function(hidden, niche.breadth)
  if (is.null(hidden)) 0 else hidden$beta * log_niche(hidden$env, hidden$optima, niche.breadth)

## Interactions scaled so that the leading eigenvalue of their term is a fixed
## fraction of self-regulation: every network is then stationary.
interaction_scale <- function(W, target) {
  ev <- max(abs(eigen(W, symmetric = TRUE, only.values = TRUE)$values))
  if (!is.finite(ev) || ev <= 0) 0 else target / ev
}

## Multispecies Gompertz, i.e. MAR(1) on the log scale, one community per site:
##   d(t+1) = A d(t) + e(t),  A = (1 - rho) I - gamma W,  counts ~ Poisson(exp(mu + d)),
## gamma W having leading eigenvalue `strength` x rho (kappa in the paper).
simulate_gompertz <- function(optima, env, W, niche.breadth, K, beta.comp = 1,
                              strength = 0.7, sigma = 0.3, rho = 0.5, steps = 200,
                              hidden = NULL) {
  n <- length(env); p <- length(optima)
  mu <- log(K) + log_niche(env, optima, niche.breadth) + hidden_niche(hidden, niche.breadth)
  gamma <- beta.comp * interaction_scale(W, strength * rho)
  A <- (1 - rho) * diag(p) - gamma * unname(W)
  d <- matrix(0, n, p)
  for (t in seq_len(steps))
    d <- d %*% t(A) + matrix(rnorm(n * p, 0, sigma), n, p)
  list(counts = matrix(rpois(n * p, pmin(exp(mu + d), 1e6)), n, p))
}

## Multispecies Ricker (discrete Lotka-Volterra), carrying capacities shaped by
## the niche; the interaction matrix has leading eigenvalue `strength`.
simulate_ricker <- function(optima, env, W, niche.breadth, K, beta.comp = 1,
                            strength = 0.5, sigma = 0.2, r = 0.5, steps = 200,
                            hidden = NULL) {
  n <- length(env); p <- length(optima)
  Kx <- pmax(K * exp(log_niche(env, optima, niche.breadth) + hidden_niche(hidden, niche.breadth)), 1e-3)
  a <- beta.comp * interaction_scale(W, strength) * unname(W)
  N <- Kx / 2
  for (t in seq_len(steps)) {
    load <- (N + N %*% t(a)) / Kx
    ## growth exponent capped at 3: facilitation is unbounded in Lotka-Volterra
    g <- pmin(r * (1 - load) + matrix(rnorm(n * p, 0, sigma), n, p), 3)
    N <- N * exp(g)
    N[!is.finite(N)] <- 0
    N <- pmin(pmax(N, 0), 10 * K)
  }
  list(counts = matrix(rpois(n * p, N), n, p))
}
