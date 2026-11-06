## Where results live: results/<campaign>/grid_res_<grid>.csv (one row per
## setting x replicate x arm x estimator) and RAI.stats_<grid>.csv. A grid may
## be stored gzipped. Rerunning some settings of a grid replaces them in the
## file and keeps the others.

RESULTS <- Sys.getenv("VC_RESULTS", "results")

grid_file <- function(campaign, grid, root = RESULTS) {
  plain <- file.path(root, campaign, sprintf("grid_res_%s.csv", grid))
  gz <- paste0(plain, ".gz")
  if (!file.exists(plain) && file.exists(gz)) gz else plain
}

setting_key <- function(d, file = FALSE) {
  cols <- intersect(c(SETTING.COLUMNS, "threshold.niche.dist"), names(d))
  if (file) return(gsub("[^A-Za-z0-9._-]", "", do.call(paste, c(lapply(d[cols], as.character), sep = "_"))))
  do.call(paste, c(lapply(d[cols], as.character), sep = "|"))
}

write_grid <- function(res, grid, campaign, root = RESULTS) {
  dir <- file.path(root, campaign)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (part in list(list(x = res$res, file = grid_file(campaign, grid, root)),
                    list(x = res$rai, file = file.path(dir, sprintf("RAI.stats_%s.csv", grid))))) {
    if (is.null(part$x) || !nrow(part$x)) next
    new <- as.data.frame(part$x)
    if (file.exists(part$file)) {
      old <- utils::read.csv(part$file)
      new <- plyr::rbind.fill(old[!setting_key(old) %in% setting_key(new), , drop = FALSE], new)
    }
    utils::write.csv(new, sub("\\.gz$", "", part$file), row.names = FALSE)
  }
  message("  written: ", grid_file(campaign, grid, root))
}

## Rebuilds a grid from the checkpoints of a run that was stopped.
assemble_checkpoints <- function(campaign, grid, root = RESULTS) {
  files <- list.files(file.path(root, campaign, "checkpoints"),
                      pattern = sprintf("^%s__.*\\.csv$", grid), full.names = TRUE)
  if (!length(files)) stop("no checkpoint for ", grid)
  write_grid(list(res = plyr::rbind.fill(lapply(files, utils::read.csv))), grid, campaign, root)
}

## One entry per run in results/<campaign>/run_info.json: date, machine, commit,
## R and package versions.
record_run <- function(campaign, info = list(), root = RESULTS) {
  path <- file.path(root, campaign, "run_info.json")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  runs <- if (file.exists(path)) jsonlite::read_json(path) else list()
  git <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE),
                  error = function(e) "", warning = function(w) "")
  pkgs <- c("PLNmodels", "gllvm", "Hmsc")
  runs[[length(runs) + 1]] <- c(list(
    started = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), host = unname(Sys.info()[["nodename"]]),
    commit = paste(git, collapse = ""), R = R.version.string,
    packages = lapply(setNames(nm = pkgs), function(p)
      if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA)), info)
  jsonlite::write_json(runs, path, auto_unbox = TRUE, pretty = TRUE, na = "null", null = "null")
}
