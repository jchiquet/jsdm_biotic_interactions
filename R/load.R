## Packages and functions of the simulation study. Every script starts with
## source("R/load.R"), from the root of the repository.
suppressMessages({
  library(PLNmodels)
  library(gllvm)
})
## one BLAS thread per process: the runs are parallel over forked processes
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1)
}
for (f in c("simulate_species.R", "simulate_networks.R", "simulate_mechanisms.R",
            "simulate_data.R", "fit_models.R", "scoring.R", "replicate.R", "results_io.R"))
  source(file.path("R", f))
## parallel and reproducible random streams
RNGkind("L'Ecuyer-CMRG")
