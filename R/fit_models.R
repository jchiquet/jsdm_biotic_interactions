## The estimators. Each fit_*() returns a fitted model; scoring.R reads it.

## PLN (unpenalised, full covariance) and PLNnetwork (its sparse precision path,
## seeded by the PLN fit). When the default optimiser fails, the PLN fit is
## tried with the second backend, then without profiling the covariance out,
## then with a diagonal covariance. The PLN fit is kept for the marginal
## reading and PLN + glasso.
fit_pln <- function(formula, data, label = "PLN") {
  formula <- as.formula(formula)
  post <- list(rsquared = FALSE)       # the R2 post-treatment diverges on rare species
  inception <- NULL
  for (control in list(PLN_param(config_post = post),
                       PLN_param(backend = "builtin", config_post = post),
                       PLN_param(config_post = post, config_optim = list(profiled = FALSE)))) {
    if (!is.null(inception)) break
    tryCatch(inception <- PLN(formula, data, control = control),
             error = function(e) message(label, ": PLN failed (", conditionMessage(e), ")"))
  }
  if (is.null(inception))
    inception <- PLN(formula, data, control = PLN_param(covariance = "diagonal", config_post = post))

  ## 76 penalties down to 1e-3 of the largest, 25 a decade
  model <- NULL
  tryCatch(model <- PLNnetwork(formula, data,
                               control = PLNnetwork_param(inception = inception, penalize_diagonal = FALSE,
                                                          min_ratio = 1e-3, n_penalties = 76)),
           error = function(e) message(label, ": PLNnetwork failed (", conditionMessage(e), ")"))
  if (is.null(model)) {
    diagonal <- PLN(formula, data, control = PLN_param(covariance = "diagonal", config_post = post))
    model <- PLNnetwork(formula, data, penalties = 10^seq(log10(30), log10(0.001), length.out = 100),
                        control = PLNnetwork_param(inception = diagonal, penalize_diagonal = FALSE))
  }
  list(network = model, inception = inception)
}

## PLN on a design matrix with an offset (sampling effort).
fit_pln_offset <- function(counts, design, offset, label = "PLN") {
  dat <- data.frame(site = seq_len(nrow(counts)))
  dat$Abundance <- as.matrix(counts); dat$X <- as.matrix(design); dat$O <- offset
  fit_pln("Abundance ~ 0 + X + offset(O)", dat, label = label)
}

## gllvm, variational approximation. `control`: family, num.lv (NULL: the
## package default, 2) and sd.errors.
fit_gllvm <- function(counts, X, formula, control = list()) {
  family <- if (is.null(control$family)) "poisson" else control$family
  se     <- if (is.null(control$sd.errors)) TRUE else control$sd.errors
  if (is.null(control$num.lv))
    gllvm(y = counts, X, family = family, formula = as.formula(formula), sd.errors = se)
  else
    gllvm(y = counts, X, family = family, formula = as.formula(formula),
          num.lv = control$num.lv, sd.errors = se)
}

## PLNPCA at one rank, returned as a covariance so that it is read like the
## log-GGM.
fit_plnpca <- function(counts, design, rank) {
  dat <- data.frame(site = seq_len(nrow(counts)))
  dat$Abundance <- as.matrix(counts); dat$X <- as.matrix(design)
  fam <- PLNPCA(Abundance ~ 0 + X, data = dat, ranks = rank,
                control = PLNPCA_param(trace = 0, config_post = list(rsquared = FALSE)))
  m <- getModel(fam, rank)
  list(Sigma = m$model_par$Sigma, loglik = m$loglik)
}

## log-GGM, no JSDM: log(1 + Y) regressed on the abiotic design; the network is
## read from the residual covariance. `clr` centres log(1 + Y) per site first.
fit_loggm <- function(counts, design, clr = FALSE) {
  L <- log1p(as.matrix(counts))
  if (clr) L <- L - rowMeans(L)
  X <- as.matrix(design)
  res <- stats::lm.fit(X, L)$residuals
  df  <- max(nrow(X) - qr(X)$rank, 1L)
  list(Sigma = crossprod(res) / df, df = df)
}

## A fit is discarded when its log-likelihood (or, for the log-GGM, its
## covariance) is not finite: an interrupted gllvm fit returns -Inf.
fit_is_usable <- function(fit, type) {
  if (type == "empirical")
    return(is.matrix(fit$Sigma) && all(is.finite(fit$Sigma)))
  ll <- switch(type,
    pln   = tryCatch(fit$inception$loglik, error = function(e) NULL),
    gllvm = tryCatch(as.numeric(stats::logLik(fit)), error = function(e) NULL))
  if (is.null(ll) || !length(ll)) return(TRUE)
  all(is.finite(ll))
}

## Runs f() in a forked process killed after `timeout` seconds (PLNmodels and
## the graphical lasso run compiled code that R cannot interrupt). Every call
## also moves the main random stream to its next substream: the number and
## order of calls are part of what makes a replicate reproducible.
in_fork <- function(f, timeout) {
  job <- parallel::mcparallel(R.utils::withTimeout(f(), timeout = timeout, onTimeout = "error"),
                              detached = FALSE)
  got <- parallel::mccollect(job, wait = FALSE, timeout = timeout + 30)
  if (is.null(got)) {
    tools::pskill(job$pid, tools::SIGKILL)
    parallel::mccollect(job, wait = FALSE)
    stop("killed after ", timeout, " s")
  }
  out <- got[[1]]
  if (inherits(out, "try-error")) stop(attr(out, "condition")$message)
  out
}
