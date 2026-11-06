## Interaction networks: species x species matrices, competition positive,
## facilitation negative.

## Numbers of competing and facilitating pairs: proportions of the pairs whose
## niche optima are closer than the threshold.
network_edge_counts <- function(optima, comp.density, facil.density, threshold) {
  d <- as.matrix(dist(optima, diag = TRUE, upper = TRUE))
  n.pairs <- sum((d < threshold) & upper.tri(d))
  c(n.comp = floor(comp.density * n.pairs), n.facil = floor(facil.density * n.pairs))
}

## Sparse network (the main design): the interacting pairs are drawn at random
## among the pairs whose optima are closer than the threshold.
network_density_based <- function(optima, comp.density, facil.density,
                                  comp.intensity, facil.intensity, threshold) {
  p <- length(optima)
  d <- as.matrix(dist(optima, diag = TRUE, upper = TRUE))
  close <- d < threshold
  upper <- which(upper.tri(close), arr.ind = TRUE)
  pairs <- upper[close[upper], , drop = FALSE]
  n.pairs <- nrow(pairs)
  n.comp  <- floor(comp.density * n.pairs)
  n.facil <- floor(facil.density * n.pairs)
  if (n.comp + n.facil == 0)
    warning("no interacting pair (", n.pairs, " pairs closer than ", threshold, ")", call. = FALSE)
  comp  <- sample(seq_len(n.pairs), n.comp)
  facil <- sample(setdiff(seq_len(n.pairs), comp), n.facil)
  W <- diag(0.5, p)
  W[pairs[comp, , drop = FALSE]]  <- comp.intensity
  W[pairs[facil, , drop = FALSE]] <- facil.intensity
  W <- W + t(W)
  dimnames(W) <- list(names(optima), names(optima))
  W
}

## Dense network of low rank (supplementary experiment): a smooth niche-overlap
## matrix truncated to `rank` eigenvectors (the first, all positive, left out),
## scaled to the interaction energy of the sparse network and shifted so that
## the proportion of competing pairs is the same. Draws no random number.
network_low_rank <- function(optima, rank, comp.density, facil.density,
                             comp.intensity, facil.intensity, threshold) {
  p <- length(optima)
  d <- as.matrix(dist(optima, diag = TRUE, upper = TRUE))
  S <- exp(-(d / threshold)^2)
  k <- max(1L, min(as.integer(rank), p - 1L))
  e <- eigen(S, symmetric = TRUE)
  idx <- seq.int(2L, length.out = k)
  idx <- idx[idx <= p]
  L <- e$vectors[, idx, drop = FALSE] %*% diag(sqrt(pmax(e$values[idx], 0)), length(idx))
  W <- L %*% t(L)
  diag(W) <- 0
  counts <- network_edge_counts(optima, comp.density, facil.density, threshold)
  n.comp <- counts[["n.comp"]]; n.facil <- counts[["n.facil"]]
  target.rms  <- sqrt((n.comp * comp.intensity^2 + n.facil * facil.intensity^2) /
                      max(p * (p - 1) / 2, 1))
  current.rms <- sqrt(mean(W[upper.tri(W)]^2))
  if (is.finite(current.rms) && current.rms > 0) W <- W * (target.rms / current.rms)
  share.comp <- if (n.comp + n.facil > 0) n.comp / (n.comp + n.facil) else 0.5
  W <- W - stats::quantile(W[upper.tri(W)], probs = 1 - share.comp, names = FALSE)
  diag(W) <- 0
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  diag(W) <- 1
  dimnames(W) <- list(names(optima), names(optima))
  W
}

## The network of a setting.
build_network <- function(optima, s) {
  if (s$network.mode == "low.rank")
    network_low_rank(optima, s$network.rank, s$competition.density, s$facilitation.density,
                     s$competition.intensity, s$facilitation.intensity, s$interaction.threshold)
  else
    network_density_based(optima, s$competition.density, s$facilitation.density,
                          s$competition.intensity, s$facilitation.intensity,
                          s$interaction.threshold)
}

## A dense network cannot be scored (it has no true negatives): its scoring
## target keeps the n.comp most competitive and n.facil most facilitative pairs.
threshold_to_support <- function(W, n.comp, n.facil) {
  M <- W; diag(M) <- 0
  ut  <- upper.tri(M)
  val <- M[ut]
  keep <- rep(FALSE, length(val))
  if (n.comp  > 0) keep[order(val, decreasing = TRUE)[seq_len(min(n.comp,  sum(val > 0)))]] <- TRUE
  if (n.facil > 0) keep[order(val)[seq_len(min(n.facil, sum(val < 0)))]] <- TRUE
  val[!keep] <- 0
  out <- matrix(0, nrow(M), ncol(M), dimnames = dimnames(M))
  out[ut] <- val
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  diag(out) <- 1
  out
}

## Reference ranking without data: a sequence of random networks that adds one
## pair at a time, each a competition or, with `facil.probability`, a
## facilitation.
random_networks_sequence <- function(n.species, facil.probability) {
  networks <- list()
  network  <- diag(0.5, n.species)
  upper    <- which(upper.tri(network), arr.ind = TRUE)
  n.pairs  <- 0.5 * n.species * (n.species - 1)
  for (k in 1:(n.pairs - 1)) {
    edge <- sample(1:nrow(upper), 1)
    type <- rbinom(1, 1, facil.probability)
    network[upper[edge, ][[1]], upper[edge, ][[2]]] <- 1 - 2 * type
    upper <- upper[-edge, ]
    networks[[k]] <- network
  }
  type <- rbinom(1, 1, facil.probability)
  network[upper[[1]], upper[[2]]] <- 1 - 2 * type
  networks[[n.pairs]] <- network
  lapply(networks, function(M) M + t(M))
}
