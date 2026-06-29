## Main replication script for the HRR comment.
##
## Run this file from the R project root.  By default, it reuses the existing
## nominal-Q estimates and regenerates the tables/figures used in the comment.
## Set REESTIMATE <- TRUE to rerun the date-by-date estimation first.

REESTIMATE <- FALSE
AREAS <- c("US", "EZ")
N_CORES <- 8
NB_LOOP <- 3
NLMINB_ITER <- 20
NM_ITER <- 1000

run_side_effect_script <- function(file) {
  sys.source(file.path("scripts", file), envir = new.env(parent = globalenv()))
}

source_script <- function(file) {
  source(file.path("scripts", file), local = globalenv())
}

run_script_main <- function(file) {
  source_script(file)
  script_main <- get("main", envir = globalenv())
  script_main()
}

set_env <- function(...) {
  values <- list(...)
  keys <- names(values)
  old <- Sys.getenv(keys, unset = NA_character_)
  do.call(Sys.setenv, values)
  invisible(old)
}

restore_env <- function(old) {
  for (key in names(old)) {
    if (is.na(old[[key]])) {
      Sys.unsetenv(key)
    } else {
      do.call(Sys.setenv, setNames(list(old[[key]]), key))
    }
  }
}

estimate_area <- function(area) {
  theta_center_file <- file.path(
    "tempo", "smooth_row_poly3_shrinkage",
    paste0(area, "_theta_center_mean_poly3_scaled.csv")
  )
  old <- set_env(
    AREA = area,
    N_CORES = as.character(N_CORES),
    NB_LOOP = as.character(NB_LOOP),
    NLMINB_ITER = as.character(NLMINB_ITER),
    NM_ITER = as.character(NM_ITER),
    MODEL_VARIANT = "smooth_row_poly3",
    TAIL_WEIGHT_SIDE = "high",
    TAIL_BIN_WEIGHT = "2",
    TAIL_MOMENT_WEIGHT = "2",
    INCLUDE_5Y5Y_TARGET = "0",
    FWD_BIN_WEIGHT = "0",
    FWD_TAIL_BIN_WEIGHT = "0",
    FWD_BOUNDS_WEIGHT = "0",
    LR_MEAN_TARGET_PCT = "2.5",
    LR_MEAN_WEIGHT = "10",
    LR_EXTREME_CAP = "0.20",
    LR_EXTREME_WEIGHT = "5",
    PERSISTENCE_WEIGHT = "0",
    THETA_CENTER_FILE = theta_center_file,
    THETA_SHRINK_WEIGHT = "0.05"
  )
  on.exit(restore_env(old), add = TRUE)

  run_script_main("estimate_nominal_q.R")
}

required_estimation_outputs <- function() {
  suffix <- paste0(
    "_nominal_Q_refit_smooth_row_poly3_finegrid_targets_",
    "hightailbin2_tailmoment2_",
    "spot_only_",
    "lrMean2p5W10_lrExtremeCap20W5_thetaShrinkW0p05"
  )
  files <- unlist(lapply(AREAS, function(area) {
    out_dir <- file.path("outputs", paste0(area, suffix))
    file.path(out_dir, c("all_results.rds", "ngmmr_Q_matrices.rds", "diagnostics.csv"))
  }))
  files
}

check_estimation_outputs <- function() {
  missing <- required_estimation_outputs()[!file.exists(required_estimation_outputs())]
  if (length(missing) > 0) {
    stop(
      "Existing estimates are missing. Set REESTIMATE <- TRUE or restore:\n",
      paste0(" - ", missing, collapse = "\n"),
      call. = FALSE
    )
  }
}

run_all <- function() {
  message("\n1. Three-date Epstein-Zin example")
  run_side_effect_script("three_date_ez_example.R")

  if (REESTIMATE) {
    message("\n2. Estimating nominal risk-neutral matrices")
    for (area in AREAS) {
      message(sprintf("\n--- %s ---", area))
      estimate_area(area)
    }
  } else {
    message("\n2. Reusing existing nominal risk-neutral matrix estimates")
    check_estimation_outputs()
  }

  message("\n3. Ratio and physical-probability figures")
  run_script_main("make_ratio_and_probability_figures.R")

  message("\n4. Fit diagnostics, HRR-consistent-P figures, and RMSE table")
  run_script_main("make_fit_diagnostics.R")

  message("\n5. Frechet bounds for the 5y5y high-inflation probability")
  run_script_main("make_5y5y_bounds_diagnostic.R")

  message("\nDone. Main outputs have been regenerated.")
}

run_all()
