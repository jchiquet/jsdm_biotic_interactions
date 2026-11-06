## Descriptions of the simulated data, without fitting any model (Supplementary
## material): what the data look like (data_description.csv), one replicate
## under the three mechanisms (example_counts.csv, example_network.csv), and the
## convergence of the lottery over 40 years (lottery_convergence.csv). Written
## to results/diagnostics/. From the root of the repository:
##   Rscript describe_simulations.R
source("R/load.R")
OUT <- file.path(RESULTS, "diagnostics")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cores <- as.integer(Sys.getenv("VC_CORES", "8"))
## the reference strengths of the main design
reference <- rbind(settings(niche.breadth = c(2, 5, 10, 20), interaction.threshold = c(10, 50)),
                   settings(niche.breadth = c(2, 5, 10, 20), interaction.threshold = c(10, 50),
                            dynamics = c("gompertz", "ricker"), dyn.strength = 0.5, dyn.sigma = 0.5))
upper_edges <- function(W) sum(W[upper.tri(W)] != 0)

## 1. 24 settings x 40 replicates, both arms
D <- do.call(rbind, parallel::mclapply(seq_len(nrow(reference)), function(k) {
  s <- as.list(reference[k, ])
  do.call(rbind, lapply(1:40, function(i) {
    r <- replicate_data(s, i)
    do.call(rbind, lapply(names(r), function(arm) { Y <- as.matrix(r[[arm]]$counts); nets <- r[[arm]]$nets
      data.frame(mech = s$dynamics, breadth = s$niche.breadth, threshold = s$interaction.threshold,
                 simu = i, arm = arm, sites = nrow(Y), species = ncol(Y), zeros = mean(Y == 0),
                 individuals.mean = mean(rowSums(Y)), individuals.min = min(rowSums(Y)),
                 individuals.max = max(rowSums(Y)), richness.mean = mean(rowSums(Y > 0)),
                 prevalence.median = median(colMeans(Y > 0)), count.max = max(Y),
                 edges.competition = upper_edges(nets$competition),
                 edges.facilitation = upper_edges(nets$facilitation)) })) })) }, mc.cores = cores))
D <- D[order(D$mech, D$breadth, D$threshold, D$simu, D$arm), ]
num <- sapply(D, is.numeric); D[num] <- lapply(D[num], signif, 6)
write.csv(D, file.path(OUT, "data_description.csv"), row.names = FALSE)

## 2. replicate 1 of b = 10, threshold 10, under the three mechanisms; species
## ordered by niche optimum
ex <- reference[reference$niche.breadth == 10 & reference$interaction.threshold == 10, ]
L <- unlist(lapply(seq_len(nrow(ex)), function(k) {
  s <- as.list(ex[k, ]); r <- replicate_data(s, 1)
  lapply(names(r), function(arm) c(r[[arm]], mech = s$dynamics, arm = arm)) }), recursive = FALSE)
counts <- do.call(rbind, lapply(L, function(x) { o <- order(x$optima); Y <- as.matrix(x$counts)
  data.frame(mech = x$mech, arm = x$arm, site = rep(round(x$env, 3), ncol(Y)),
             species = rep(seq_along(o), each = nrow(Y)), count = as.vector(Y[, o])) }))
x <- L[[1]]; o <- order(x$optima)
W <- (x$nets$competition - abs(x$nets$facilitation))[o, o]
network <- data.frame(j = rep(seq_along(o), each = length(o)), k = rep(seq_along(o), length(o)),
                      w = as.vector(W), optimum = rep(round(sort(x$optima), 3), each = length(o)))
write.csv(counts, file.path(OUT, "example_counts.csv"), row.names = FALSE)
write.csv(network, file.path(OUT, "example_network.csv"), row.names = FALSE)

## 3. the lottery run for 40 years instead of 10, 10 replicates: change in
## abundance of a species at a site from one year to the next
lottery <- reference[reference$dynamics == "lottery", ]
CV <- do.call(rbind, parallel::mclapply(seq_len(nrow(lottery)), function(k) {
  s <- as.list(lottery[k, ])
  do.call(rbind, lapply(1:10, function(i) {
    r <- replicate_data(s, i, years = 40, cv.measure = TRUE)
    do.call(rbind, lapply(names(r), function(arm)
      data.frame(breadth = s$niche.breadth, threshold = s$interaction.threshold, arm = arm,
                 year = 1:40, turnover = colMeans(r[[arm]]$turnover)))) })) }, mc.cores = cores))
CV <- aggregate(turnover ~ breadth + threshold + arm + year, CV, mean)
CV$turnover <- signif(CV$turnover, 6)
write.csv(CV[order(CV$arm, CV$threshold, CV$breadth, CV$year), ], file.path(OUT, "lottery_convergence.csv"),
          row.names = FALSE)
