## Scoring an estimated network against the true one.
##
## Every estimator ranks the pairs of species, for competition and for
## facilitation separately. Two metrics are read off one ranking: the AUC, and
## the top-k precision, k being the number of true interactions of that sign.
## Readings: PLNnetwork ranks the pairs by the order in which they enter its
## penalty path; the other estimators report a covariance, ranked either as a
## thresholded correlation ("marginal") or through a graphical-lasso path built
## on it (columns with suffix "_alt").

## Sign of an interaction in an estimated matrix. True networks: competition
## positive. A covariance shows competition as a negative entry; a precision
## matrix, whose sign is that of minus the partial correlation, as a positive one.
sign_for <- function(interaction, kind) {
  if (kind == "precision") if (interaction == "competition") "positive" else "negative"
  else if (interaction == "competition") "negative" else "positive"
}

## Graphical-lasso path on a covariance, 25 penalties a decade, from the
## largest |Sigma_ij| down to depth x that value.
sparse_precision_path <- function(Sigma, depth = 1e-4, per_decade = 25) {
  n_penalties <- round(-log10(depth) * per_decade) + 1
  hi <- max(abs(Sigma[upper.tri(Sigma)]))
  ## a collapsed covariance: no residual association
  if (!is.finite(hi) || hi <= 0 || max(abs(diag(Sigma))) < 1e-8)
    return(rep(list(diag(1 / pmax(diag(Sigma), 1e-8))), n_penalties))
  ## a tiny covariance is rescaled, or the solver does not converge
  scale <- max(abs(diag(Sigma)))
  if (scale < 1e-6) { Sigma <- Sigma / scale; hi <- hi / scale }
  rhos <- hi * 10^(-(seq_len(n_penalties) - 1) / per_decade)
  lapply(rhos, function(rho) PLNmodels::graphical_lasso(Sigma, rho = rho)$wi)
}

## The covariance an estimator reports (PLN: the Sigma of its unpenalised fit).
family_covariance <- function(model, type)
  switch(type,
    pln.marginal = {
      s <- model$model_par$Sigma
      if (!is.matrix(s)) stop("diagonal Sigma: no marginal association to rank")
      s },
    gllvm     = getResidualCov(model)$cov,
    empirical = model$Sigma)

## Ranking from a path of sparse matrices: order of entry with the right sign,
## sparsest first, ties broken by the size of the entry in the densest matrix.
path_ranking <- function(mats, orient) {
  n.edges <- vapply(mats, function(M) sum(M[upper.tri(M)] != 0), numeric(1))
  mats <- mats[order(n.edges)]
  n <- length(mats); p <- nrow(mats[[1]])
  first <- matrix(Inf, p, p)
  for (i in rev(seq_len(n))) first[is.finite(orient(mats[[i]]))] <- i
  dense <- orient(mats[[n]]); dense[!is.finite(dense)] <- 0
  tie <- if (max(dense) > 0) dense / (2 * max(dense)) else dense
  S <- (n + 1 - first) + tie
  S[!is.finite(first)] <- -Inf
  diag(S) <- -Inf
  S
}

## Score matrix of the pairs (larger = stronger evidence; wrong sign = -Inf).
ranking_score <- function(model, type, sign, reading = c("native", "glasso"), path = NULL) {
  reading <- match.arg(reading)
  orient <- function(M) {
    S <- if (sign == "positive") M else -M
    S[!(S > 0)] <- -Inf
    diag(S) <- -Inf
    S
  }
  if (type == "pln")
    return(path_ranking(lapply(model$models, function(m) m$model_par$Omega), orient))
  if (reading == "glasso")
    return(path_ranking(path, orient))
  Sigma <- family_covariance(model, type)
  orient(if (type == "gllvm") Sigma else cov2cor(Sigma))
}

## AUC of a ranking, exactly (mid-ranks for ties); NA without both classes.
ranking_auc <- function(score, truth) {
  up <- upper.tri(truth)
  s  <- score[up]; y <- truth[up] != 0
  if (!any(y) || all(y)) return(NA_real_)
  lev <- sort(unique(s), decreasing = TRUE)
  g   <- match(s, lev)
  f <- c(0, cumsum(tabulate(g[!y], length(lev))) / sum(!y))
  t <- c(0, cumsum(tabulate(g[y],  length(lev))) / sum(y))
  sum(diff(f) * (head(t, -1) + tail(t, -1)) / 2)
}

## Top-k precision: proportion of true pairs among the k best ranked, k the
## number of true pairs.
precision_k <- function(score, truth) {
  up <- upper.tri(truth)
  y  <- truth[up] != 0
  k  <- sum(y)
  s  <- score[up]
  cand <- which(is.finite(s))
  sel  <- utils::head(cand[order(-s[cand])], k)
  if (k > 0) sum(y[sel]) / k else NA_real_
}

## One row of measures for one fitted model.
evaluate_model <- function(method, model, type, nets, fit.time, abiotic = NO.ABIOTIC) {
  if (type == "pln.marginal" && !is.matrix(model$model_par$Sigma))   # diagonal fallback
    return(data.frame(method = method, time = as.numeric(fit.time, units = "secs")))
  kind <- if (type == "pln") "precision" else "covariance"
  has.glasso <- type != "pln"
  path <- if (has.glasso)
    tryCatch(sparse_precision_path(family_covariance(model, type)), error = function(e) NULL)
  measures <- function(interaction, reading) {
    S <- if (reading == "glasso" && is.null(path)) NULL else
      tryCatch(ranking_score(model, type, sign_for(interaction, if (reading == "glasso") "precision" else kind),
                             reading, path), error = function(e) NULL)
    net <- nets[[interaction]]
    sfx <- paste0("_", interaction, if (reading == "glasso") "_alt" else "")
    setNames(list(if (is.null(S)) NA_real_ else ranking_auc(S, net),
                  if (is.null(S)) NA_real_ else tryCatch(precision_k(S, net), error = function(e) NA_real_)),
             paste0(c("AUC", "precision_k"), sfx))
  }
  readings <- if (has.glasso) c("native", "glasso") else "native"
  as.data.frame(c(list(method = method, time = as.numeric(fit.time, units = "secs")),
                  as.list(abiotic),
                  unlist(lapply(readings, function(r)
                    c(measures("competition", r), measures("facilitation", r))), recursive = FALSE)))
}

NO.ABIOTIC <- c(abiotic_cor = NA_real_, niche_optimum_cor = NA_real_, abiotic_share = NA_real_)

## How well a PLN fit recovers the niche: mean correlation of the fitted and
## true log-responses over species, correlation of the fitted and true optima,
## and the proportion of latent variance taken by the abiotic part.
abiotic_quality <- function(model, design, optima, env, niche.breadth) {
  B <- tryCatch(model$model_par$B, error = function(e) NULL)
  if (is.null(B) || nrow(B) != ncol(design) || nrow(design) != length(env)) return(NO.ABIOTIC)
  fitted <- design %*% B
  truth  <- outer(env, optima, function(e, o) dnorm(o, mean = e, sd = niche.breadth, log = TRUE))
  per.species <- vapply(seq_len(ncol(fitted)), function(j)
    if (sd(fitted[, j]) == 0 || sd(truth[, j]) == 0) NA_real_ else cor(fitted[, j], truth[, j]), numeric(1))
  optimum.hat <- env[apply(fitted, 2, which.max)]
  Sigma <- model$model_par$Sigma
  residual <- if (is.matrix(Sigma)) diag(Sigma) else as.numeric(Sigma)
  env.var  <- apply(fitted, 2, var)
  c(abiotic_cor = mean(per.species, na.rm = TRUE),
    niche_optimum_cor = if (sd(optimum.hat) == 0 || sd(optima) == 0) NA_real_ else cor(optimum.hat, optima),
    abiotic_share = if (length(residual) != length(env.var)) NA_real_
                    else mean(env.var / (env.var + residual), na.rm = TRUE))
}

## ROC of a sequence of networks (recall, fallout), and the AUC of the curve:
## used for the random reference networks only.
roc_points <- function(truth, estimates, sign) {
  diag(truth) <- 0; p <- nrow(truth)
  true.edge <- truth != 0
  do.call(rbind, lapply(estimates, function(M) {
    if (sign == "positive") M[M < 0] <- 0
    if (sign == "negative") M[M > 0] <- 0
    diag(M) <- 0
    est <- M != 0
    TP <- sum(est & true.edge); FN <- sum(!est & true.edge)
    FP <- sum(est & !true.edge); TN <- sum(!est & !true.edge) - p
    recall  <- if (TP + FN == 0) NA else round(TP / (TP + FN), 3)
    fallout <- if (TN + FP == 0) NA else round(FP / (FP + TN), 3)
    c(fallout = fallout, recall = recall)
  }))
}
roc_auc <- function(roc) {
  ok <- !is.na(roc[, "fallout"]) & !is.na(roc[, "recall"])
  f <- c(0, roc[ok, "fallout"], 1); t <- c(0, roc[ok, "recall"], 1)
  o <- order(f, t); f <- f[o]; t <- t[o]
  sum(diff(f) * (head(t, -1) + tail(t, -1)) / 2)
}

## Rankings that use no count data: random networks (with the density and the
## niche-proximity rule of the simulation, or fully random) and the proximity
## of the niche optima. Their random draws come after the fits, as in the
## published grids.
reference_rows <- function(data, s) {
  nets <- data$nets
  W <- nets$all; facil <- nets$facilitation; comp <- nets$competition
  d <- as.matrix(dist(data$optima, diag = TRUE, upper = TRUE))
  upper <- which(upper.tri(d, diag = FALSE), arr.ind = TRUE)
  interacting <- upper[W[upper] != 0, ]
  threshold <- max(d[interacting])
  ratio <- s$facilitation.density / s$competition.density
  informed <- lapply(seq(0.01, s$competition.density, length.out = 50), function(x)
    network_density_based(data$optima, x, ratio * x, 1, -1, threshold))
  informed.row <- list(method = "informed.random",
    AUC_competition  = roc_auc(roc_points(comp, informed, "positive")),
    AUC_facilitation = roc_auc(roc_points(facil, informed, "negative")))

  p <- length(data$optima)
  n.facil <- nrow(upper[facil[upper] != 0, ]); if (is.null(n.facil)) n.facil <- 1
  n.comp  <- nrow(upper[comp[upper] != 0, ]);  if (is.null(n.comp))  n.comp  <- 1
  n.pairs <- 0.5 * p * (p - 1)
  r <- (n.facil / n.pairs) / (n.comp / n.pairs)
  fully <- random_networks_sequence(p, facil.probability = r / (1 + r))
  fully.row <- list(method = "fully.random",
    AUC_competition  = roc_auc(roc_points(comp, fully, "positive")),
    AUC_facilitation = if (n.facil > 0) roc_auc(roc_points(facil, fully, "negative")) else NA)

  proximity <- -d; diag(proximity) <- -Inf
  proximity.row <- list(method = "niche.proximity",
    AUC_competition = ranking_auc(proximity, comp), AUC_facilitation = ranking_auc(proximity, facil),
    precision_k_competition = precision_k(proximity, comp),
    precision_k_facilitation = precision_k(proximity, facil))
  lapply(list(informed.row, fully.row, proximity.row), as.data.frame)
}
