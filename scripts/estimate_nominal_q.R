## Estimate nominal RN transition matrices from 5y, 10y, and 5y5y targets.

suppressPackageStartupMessages({
  library(haven)
})

source(file.path("scripts", "estimation_helpers.R"))

make_tail_weight_spec <- function(tail_bin_weight = 3,
                                  tail_moment_weight = 5,
                                  fwd_bin_weight = 1,
                                  tail_weight_side = c("high", "both")) {
  tail_weight_side <- match.arg(tail_weight_side)
  bin_weights <- rep(1, 24)
  ## In each 8-bin block, bins 7 and 8 are the >4% inflation tail. Optionally
  ## bins 1 and 2 (<0%) can also be overweighted.  The 5y5y block is kept as
  ## a loose auxiliary target, so it receives its own lower uniform weight.
  tail_bins_one_block <- if (tail_weight_side == "both") c(1:2, 7:8) else 7:8
  tail_bins <- c(tail_bins_one_block,
                 8 + tail_bins_one_block)
  bin_weights[tail_bins] <- tail_bin_weight
  bin_weights[17:24] <- fwd_bin_weight

  list(
    bin_weights = bin_weights,
    tail_moment_weight = tail_moment_weight,
    tail_bin_weight = tail_bin_weight,
    fwd_bin_weight = fwd_bin_weight,
    tail_weight_side = tail_weight_side
  )
}

weighted_finegrid_loss <- function(model_g5, model_g10, model_fwd,
                                   target_g5, target_g10, target_fwd,
                                   weight_spec,
                                   include_fwd = TRUE,
                                   fwd_tail4_proxy = NA_real_,
                                   fwd_tail4_proxy_weight = 0,
                                   fwd_tail4_bounds = c(NA_real_, NA_real_),
                                   fwd_tail4_bounds_weight = 0) {
  if (include_fwd) {
    bin_errors <- c(model_g5 - target_g5,
                    model_g10 - target_g10,
                    model_fwd - target_fwd)
    bin_weights <- weight_spec$bin_weights
  } else {
    bin_errors <- c(model_g5 - target_g5,
                    model_g10 - target_g10)
    bin_weights <- weight_spec$bin_weights[1:16]
  }
  bin_loss <- mean(bin_weights * bin_errors^2)

  tail_errors <- c(
    sum(model_g5[7:8]) - sum(target_g5[7:8]),
    sum(model_g10[7:8]) - sum(target_g10[7:8])
  )
  if (include_fwd) {
    ## The 5y5y target, when used, is a loose auxiliary Gaussian proxy.  Do
    ## not also impose the explicit tail-moment overweight on it.
  }
  if (weight_spec$tail_weight_side == "both") {
    low_tail_errors <- c(
      sum(model_g5[1:2]) - sum(target_g5[1:2]),
      sum(model_g10[1:2]) - sum(target_g10[1:2])
    )
    if (include_fwd) {
      ## Same logic as above: low-tail moment overweight applies only to the
      ## directly option-implied 5y and 10y targets.
    }
    tail_errors <- as.vector(rbind(low_tail_errors, tail_errors))
  }
  tail_loss <- weight_spec$tail_moment_weight * mean(tail_errors^2)

  proxy_loss <- 0
  if (is.finite(fwd_tail4_proxy) && fwd_tail4_proxy_weight > 0) {
    proxy_loss <- fwd_tail4_proxy_weight *
      (sum(model_fwd[7:8]) - fwd_tail4_proxy)^2
  }

  bounds_loss <- 0
  if (fwd_tail4_bounds_weight > 0 &&
      length(fwd_tail4_bounds) == 2 &&
      all(is.finite(fwd_tail4_bounds))) {
    model_tail4 <- sum(model_fwd[7:8])
    below <- max(fwd_tail4_bounds[1] - model_tail4, 0)
    above <- max(model_tail4 - fwd_tail4_bounds[2], 0)
    bounds_loss <- fwd_tail4_bounds_weight * (below^2 + above^2)
  }

  bin_loss + tail_loss + proxy_loss + bounds_loss
}

finegrid_zc_bins <- function(dists, year_id, month_id, horizon,
                             pi_lim_pct = c(seq(-1, 5, by = 1), Inf)) {
  support_pct <- 100 * as.numeric(dists$sprt.Z[, 1])
  pmf <- as.numeric(dists$data.ZC.N[, horizon, year_id, month_id])
  ## HRR's reported tail probabilities include grid points that sit on the
  ## threshold, e.g. the 4.0% support point enters the ">4%" tail because the
  ## stored floating-point value is just above 0.04.  Assign threshold points
  ## to the upper bin so bins 7:8 match those reported-tail conventions.
  bin <- length(pi_lim_pct) + 1L -
    rowSums(outer(support_pct, pi_lim_pct, function(x, lim) lim - x > 1e-7))
  out <- numeric(8)
  for (b in seq_len(8)) {
    out[b] <- sum(pmf[bin == b], na.rm = TRUE)
  }
  out / sum(out)
}

adjust_forward_proxy_bins <- function(base_bins, tail_row) {
  out <- pmax(as.numeric(base_bins), 0)

  fixed <- numeric(8)
  fixed[1] <- max(tail_row$tailm1_5y5y, 0)
  fixed[2] <- max(tail_row$tail0_5y5y - tail_row$tailm1_5y5y, 0)
  fixed[7] <- max(tail_row$tail4_5y5y - tail_row$tail5_5y5y, 0)
  fixed[8] <- max(tail_row$tail5_5y5y, 0)

  rem <- 1 - sum(fixed[c(1, 2, 7, 8)])
  if (!is.finite(rem) || rem < 0) {
    stop("Invalid 5y5y tail proxies: fixed tails exceed one.")
  }

  middle <- 3:6
  middle_base <- pmax(out[middle], 0)
  if (sum(middle_base) <= 0) {
    out[middle] <- rem / length(middle)
  } else {
    out[middle] <- rem * middle_base / sum(middle_base)
  }
  out[c(1, 2, 7, 8)] <- fixed[c(1, 2, 7, 8)]
  out / sum(out)
}

moment_proxy_5y5y_tail4 <- function(g5, g10,
                                    bin_values_pct = c(-2, -0.5, 0.5, 1.5,
                                                       2.5, 3.5, 4.5, 6),
                                    rho = 0,
                                    min_sd = 0.25) {
  mean5 <- sum(g5 * bin_values_pct)
  mean10 <- sum(g10 * bin_values_pct)
  var5 <- sum(g5 * (bin_values_pct - mean5)^2)
  var10 <- sum(g10 * (bin_values_pct - mean10)^2)

  mean_fwd <- 2 * mean10 - mean5
  disc <- 4 * var10 - var5 * (1 - rho^2)
  if (!is.finite(disc) || disc < 0) {
    return(NA_real_)
  }
  sd_fwd <- -rho * sqrt(var5) + sqrt(disc)
  sd_fwd <- max(sd_fwd, min_sd)
  out <- 1 - pnorm(4, mean = mean_fwd, sd = sd_fwd)
  min(max(out, 0), 1)
}

gaussian_5y5y_bins <- function(g5, g10,
                               rho,
                               bin_values_pct = c(-2, -0.5, 0.5, 1.5,
                                                  2.5, 3.5, 4.5, 6),
                               pi_lim_pct = c(seq(-1, 5, by = 1), Inf),
                               min_sd = 0.25) {
  mean5 <- sum(g5 * bin_values_pct)
  mean10 <- sum(g10 * bin_values_pct)
  var5 <- sum(g5 * (bin_values_pct - mean5)^2)
  var10 <- sum(g10 * (bin_values_pct - mean10)^2)

  mean_fwd <- 2 * mean10 - mean5
  disc <- 4 * var10 - var5 * (1 - rho^2)
  if (!is.finite(disc) || disc < 0) {
    return(rep(NA_real_, length(pi_lim_pct)))
  }
  sd_fwd <- -rho * sqrt(var5) + sqrt(disc)
  sd_fwd <- max(sd_fwd, min_sd)

  lower <- c(-Inf, pi_lim_pct[-length(pi_lim_pct)])
  upper <- pi_lim_pct
  out <- pnorm(upper, mean = mean_fwd, sd = sd_fwd) -
    pnorm(lower, mean = mean_fwd, sd = sd_fwd)
  out <- pmax(out, 0)
  out / sum(out)
}

estimate_realized_forward_avg_rho <- function(infl_pct, window = 60L) {
  infl_raw <- as.numeric(infl_pct)
  infl_pct <- infl_raw
  infl_pct[infl_pct == 0] <- NA_real_
  n <- length(infl_pct)
  if (n < 2 * window) {
    return(NA_real_)
  }

  compute_rho <- function(x) {
    starts <- seq_len(length(x) - 2 * window + 1L)
    first_avg <- vapply(starts, function(i) {
      mean(x[i:(i + window - 1L)], na.rm = FALSE)
    }, numeric(1))
    second_avg <- vapply(starts, function(i) {
      mean(x[(i + window):(i + 2L * window - 1L)], na.rm = FALSE)
    }, numeric(1))
    ok <- is.finite(first_avg) & is.finite(second_avg)
    if (sum(ok) < 5) {
      return(NA_real_)
    }
    cor(first_avg[ok], second_avg[ok])
  }

  rho <- compute_rho(infl_pct)
  if (!is.finite(rho) && any(infl_raw == 0, na.rm = TRUE)) {
    ## Useful for the EA file, where the zero-as-missing convention is too strict.
    rho <- compute_rho(infl_raw)
  }
  rho
}

resolve_fwd_proxy_rho <- function(value, dists) {
  if (identical(tolower(value), "data")) {
    rho <- estimate_realized_forward_avg_rho(as.numeric(dists$infl))
  } else {
    rho <- suppressWarnings(as.numeric(value))
  }
  if (!is.finite(rho)) {
    stop("FWD_PROXY_RHO must be numeric or 'data', and the data estimate must be finite.")
  }
  max(min(rho, 0.999), -0.999)
}

make_year_month_index <- function(dists) {
  years <- as.integer(dists$years[, 1])
  last_month <- as.integer(dists$last.month[1, 1])
  do.call(
    rbind,
    lapply(seq_along(years), function(year_id) {
      month_max <- if (year_id == length(years)) last_month else 12L
      data.frame(
        year_id = year_id,
        year = years[year_id],
        month = seq_len(month_max)
      )
    })
  )
}

make_finegrid_targets <- function(res, dists,
                                  area = Sys.getenv("AREA", unset = "US"),
                                  fwd_target_source = c("hrr", "gaussian"),
                                  fwd_proxy_rho = 0) {
  area <- normalize_hrr_area(area)
  fwd_target_source <- match.arg(fwd_target_source)
  fs_data <- res$fs.data
  n_month <- dim(fs_data)[2]
  ym <- make_year_month_index(dists)
  tails_55 <- read_dta(file.path("input", sprintf("%s_55tails_monthly.dta", area)))

  g5 <- matrix(NA_real_, 8, n_month)
  g10 <- matrix(NA_real_, 8, n_month)
  fwd <- matrix(NA_real_, 8, n_month)
  fwd_moment_proxy_tail4 <- rep(NA_real_, n_month)
  fwd_tail4_lower <- rep(NA_real_, n_month)
  fwd_tail4_upper <- rep(NA_real_, n_month)

  bounds_file <- file.path("outputs", "hrr_nominal_real_q_diagnostics",
                           sprintf("%s_5y5y_tail_bounds_from_5y_10y.csv", area))
  if (!file.exists(bounds_file) && area == "US") {
    bounds_file <- file.path("outputs", "hrr_nominal_real_q_diagnostics",
                             "5y5y_tail_bounds_from_5y_10y.csv")
  }
  if (file.exists(bounds_file)) {
    bounds_df <- read.csv(bounds_file)
    if (all(c("month", "lower", "upper") %in% names(bounds_df))) {
      bounds_df <- bounds_df[match(seq_len(n_month), bounds_df$month), ]
      fwd_tail4_lower <- bounds_df$lower
      fwd_tail4_upper <- bounds_df$upper
    }
  }

  for (m in seq_len(n_month)) {
    g5[, m] <- finegrid_zc_bins(dists, ym$year_id[m], ym$month[m], horizon = 5)
    g10[, m] <- finegrid_zc_bins(dists, ym$year_id[m], ym$month[m], horizon = 10)
    fwd_moment_proxy_tail4[m] <- moment_proxy_5y5y_tail4(
      g5[, m], g10[, m], rho = fwd_proxy_rho
    )
    tail_row <- tails_55[tails_55$year == ym$year[m] & tails_55$month == ym$month[m], ]
    if (nrow(tail_row) != 1) {
      stop(sprintf("Could not find unique 5y5y tail row for month %03d.", m))
    }
    if (fwd_target_source == "gaussian") {
      fwd[, m] <- gaussian_5y5y_bins(g5[, m], g10[, m], rho = fwd_proxy_rho)
    } else {
      fwd[, m] <- adjust_forward_proxy_bins(fs_data[, m], tail_row)
    }
  }

  list(g5 = g5, g10 = g10, fwd = fwd,
       fwd_moment_proxy_tail4 = fwd_moment_proxy_tail4,
       fwd_tail4_lower = fwd_tail4_lower,
       fwd_tail4_upper = fwd_tail4_upper,
       fwd_target_source = fwd_target_source,
       fwd_proxy_rho = fwd_proxy_rho,
       year_month = ym)
}

fit_one_month_finegrid <- function(month_id,
                                   targets,
                                   gs_model,
                                   fs_model,
                                   bts,
                                   nb_loop,
                                   nlminb_iter,
                                   nm_iter,
                                   variant,
                                   mask,
                                   eta0,
                                   weight_spec,
                                   include_fwd_target,
                                   fwd_proxy_weight,
                                   fwd_bounds_weight) {
  target_g5 <- targets$g5[, month_id]
  target_g10 <- targets$g10[, month_id]
  target_fwd <- targets$fwd[, month_id]
  fwd_tail4_proxy <- targets$fwd_moment_proxy_tail4[month_id]
  fwd_tail4_bounds <- c(targets$fwd_tail4_lower[month_id],
                        targets$fwd_tail4_upper[month_id])
  bt <- bts[month_id]

  obj <- function(eta) {
    A <- eta_to_structured_mlogit_A(eta, mask = mask, variant = variant)
    mod <- avg_inf_510_fast(A, bt)
    fwd <- avg_yoy_6to10_model(A, bt)
    weighted_finegrid_loss(
      model_g5 = mod$g5,
      model_g10 = mod$g10,
      model_fwd = fwd,
      target_g5 = target_g5,
      target_g10 = target_g10,
      target_fwd = target_fwd,
      weight_spec = weight_spec,
      include_fwd = include_fwd_target,
      fwd_tail4_proxy = fwd_tail4_proxy,
      fwd_tail4_proxy_weight = fwd_proxy_weight,
      fwd_tail4_bounds = fwd_tail4_bounds,
      fwd_tail4_bounds_weight = fwd_bounds_weight
    )
  }

  fit <- fit_with_nlminb_nm_loops(
    z0 = eta0,
    obj = obj,
    nb_loop = nb_loop,
    nlminb_iter = nlminb_iter,
    nm_iter = nm_iter,
    trace = FALSE
  )

  A <- eta_to_structured_mlogit_A(fit$par, mask = mask, variant = variant)
  mod <- avg_inf_510_fast(A, bt)
  fwd <- avg_yoy_6to10_model(A, bt)

  hrr_g5 <- gs_model[, 1, month_id]
  hrr_g10 <- gs_model[, 2, month_id]
  hrr_fwd <- fs_model[, month_id]

  diagnostics <- data.frame(
    month = month_id,
    objective_ngmmr_24 = fit$value,
    objective_hrr_24 = weighted_finegrid_loss(
      model_g5 = hrr_g5,
      model_g10 = hrr_g10,
      model_fwd = hrr_fwd,
      target_g5 = target_g5,
      target_g10 = target_g10,
      target_fwd = target_fwd,
      weight_spec = weight_spec,
      include_fwd = include_fwd_target,
      fwd_tail4_proxy = fwd_tail4_proxy,
      fwd_tail4_proxy_weight = fwd_proxy_weight,
      fwd_tail4_bounds = fwd_tail4_bounds,
      fwd_tail4_bounds_weight = fwd_bounds_weight
    ),
    rmse24_hrr_pp = rmse_pp(c(hrr_g5 - target_g5,
                              hrr_g10 - target_g10,
                              hrr_fwd - target_fwd)),
    rmse24_ngmmr_pp = rmse_pp(c(mod$g5 - target_g5,
                                mod$g10 - target_g10,
                                fwd - target_fwd)),
    rmse5_hrr_pp = rmse_pp(hrr_g5 - target_g5),
    rmse5_ngmmr_pp = rmse_pp(mod$g5 - target_g5),
    rmse10_hrr_pp = rmse_pp(hrr_g10 - target_g10),
    rmse10_ngmmr_pp = rmse_pp(mod$g10 - target_g10),
    rmse_fwd_hrr_pp = rmse_pp(hrr_fwd - target_fwd),
    rmse_fwd_ngmmr_pp = rmse_pp(fwd - target_fwd),
    qtail0_5y_target = sum(target_g5[1:2]),
    qtail0_5y_hrr = sum(hrr_g5[1:2]),
    qtail0_5y_ngmmr = sum(mod$g5[1:2]),
    qtail4_5y_target = sum(target_g5[7:8]),
    qtail4_5y_hrr = sum(hrr_g5[7:8]),
    qtail4_5y_ngmmr = sum(mod$g5[7:8]),
    qtail0_10y_target = sum(target_g10[1:2]),
    qtail0_10y_hrr = sum(hrr_g10[1:2]),
    qtail0_10y_ngmmr = sum(mod$g10[1:2]),
    qtail4_10y_target = sum(target_g10[7:8]),
    qtail4_10y_hrr = sum(hrr_g10[7:8]),
    qtail4_10y_ngmmr = sum(mod$g10[7:8]),
    qtail0_5y5y_target = sum(target_fwd[1:2]),
    qtail0_5y5y_hrr = sum(hrr_fwd[1:2]),
    qtail0_5y5y_ngmmr = sum(fwd[1:2]),
    qtail4_5y5y_target = sum(target_fwd[7:8]),
    qtail4_5y5y_hrr = sum(hrr_fwd[7:8]),
    qtail4_5y5y_ngmmr = sum(fwd[7:8]),
    qtail4_5y5y_moment_proxy = fwd_tail4_proxy,
    qtail4_5y5y_lower_bound = fwd_tail4_bounds[1],
    qtail4_5y5y_upper_bound = fwd_tail4_bounds[2],
    qtail4_5y5y_bound_violation_hrr =
      max(fwd_tail4_bounds[1] - sum(hrr_fwd[7:8]), 0) +
      max(sum(hrr_fwd[7:8]) - fwd_tail4_bounds[2], 0),
    qtail4_5y5y_bound_violation_ngmmr =
      max(fwd_tail4_bounds[1] - sum(fwd[7:8]), 0) +
      max(sum(fwd[7:8]) - fwd_tail4_bounds[2], 0),
    fwd_proxy_weight = fwd_proxy_weight,
    fwd_bounds_weight = fwd_bounds_weight,
    tail_bin_weight = weight_spec$tail_bin_weight,
    tail_moment_weight = weight_spec$tail_moment_weight,
    tail_weight_side = weight_spec$tail_weight_side,
    valid_A = is_valid_transition_matrix(A)
  )

  list(
    diagnostics = diagnostics,
    eta = fit$par,
    A = A,
    history = fit$history,
    target = list(g5 = target_g5, g10 = target_g10, fwd = target_fwd),
    hrr = list(g5 = hrr_g5, g10 = hrr_g10, fwd = hrr_fwd),
    ngmmr = list(g5 = mod$g5, g10 = mod$g10, fwd = fwd)
  )
}

plot_finegrid_rmse <- function(diagnostics, output_file) {
  pdf(output_file, width = 11, height = 7)
  oldpar <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  plot_one <- function(hrr, ngmmr, title) {
    ymax <- max(c(hrr, ngmmr), na.rm = TRUE) * 1.08
    plot(diagnostics$date, hrr, type = "l", col = "gray55", lwd = 2,
         ylim = c(0, ymax), xlab = "", ylab = "RMSE, pp",
         main = title, las = 1)
    grid()
    lines(diagnostics$date, ngmmr, col = "black", lwd = 2)
    legend("topleft", legend = c("HRR", "NGMMR"),
           col = c("gray55", "black"), lwd = 2, bty = "n")
  }

  plot_one(diagnostics$rmse24_hrr_pp, diagnostics$rmse24_ngmmr_pp,
           "5y + 10y + 5y5y")
  plot_one(diagnostics$rmse5_hrr_pp, diagnostics$rmse5_ngmmr_pp, "5y")
  plot_one(diagnostics$rmse10_hrr_pp, diagnostics$rmse10_ngmmr_pp, "10y")
  plot_one(diagnostics$rmse_fwd_hrr_pp, diagnostics$rmse_fwd_ngmmr_pp, "5y5y")
  mtext("Fine-grid/proxy 24-target fit: HRR matrix vs. NGMMR",
        outer = TRUE, font = 2)
}

main <- function() {
  area <- normalize_hrr_area(Sys.getenv("AREA", unset = "US"))
  n_cores <- as.integer(Sys.getenv("N_CORES", unset = "8"))
  nb_loop <- as.integer(Sys.getenv("NB_LOOP", unset = "3"))
  nlminb_iter <- as.integer(Sys.getenv("NLMINB_ITER", unset = "20"))
  nm_iter <- as.integer(Sys.getenv("NM_ITER", unset = "1000"))
  tail_bin_weight <- as.numeric(Sys.getenv("TAIL_BIN_WEIGHT", unset = "1.5"))
  tail_moment_weight <- as.numeric(Sys.getenv("TAIL_MOMENT_WEIGHT", unset = "2"))
  fwd_bin_weight <- as.numeric(Sys.getenv("FWD_BIN_WEIGHT", unset = "1"))
  tail_weight_side <- Sys.getenv("TAIL_WEIGHT_SIDE", unset = "high")
  include_fwd_target <- as.logical(as.integer(Sys.getenv("INCLUDE_5Y5Y_TARGET", unset = "1")))
  fwd_target_source <- Sys.getenv("FWD_TARGET_SOURCE", unset = "gaussian")
  fwd_proxy_rho_input <- Sys.getenv("FWD_PROXY_RHO", unset = "data")
  fwd_proxy_weight <- as.numeric(Sys.getenv("FWD_PROXY_WEIGHT", unset = "0.2"))
  fwd_bounds_weight <- as.numeric(Sys.getenv("FWD_BOUNDS_WEIGHT", unset = "0"))
  eta_start_file <- Sys.getenv("ETA_START_FILE", unset = "")
  variant <- Sys.getenv("MODEL_VARIANT", unset = "smooth_row_poly2")
  weight_spec <- make_tail_weight_spec(
    tail_bin_weight = tail_bin_weight,
    tail_moment_weight = tail_moment_weight,
    fwd_bin_weight = fwd_bin_weight,
    tail_weight_side = tail_weight_side
  )

  mask <- hrr_zero_pattern()
  eta_length <- qn_variant_eta_length(variant, mask = mask)
  eta0 <- rep(0, eta_length)

  weight_tag <- sprintf("%stailbin%s_tailmoment%s",
                        tail_weight_side,
                        gsub("\\.", "p", format(tail_bin_weight, trim = TRUE)),
                        gsub("\\.", "p", format(tail_moment_weight, trim = TRUE)))
  out_name <- paste0("nominal_Q_refit_", variant, "_finegrid_targets_", weight_tag)
  if (!include_fwd_target) {
    out_name <- paste0(out_name, "_spot_only")
  }
  if (fwd_proxy_weight > 0) {
    proxy_tag <- gsub("\\.", "p", format(fwd_proxy_weight, trim = TRUE))
    out_name <- paste0(out_name, "_moment5y5yW", proxy_tag)
  }
  if (include_fwd_target && fwd_target_source == "gaussian") {
    rho_tag <- gsub("\\.", "p", gsub("-", "m", fwd_proxy_rho_input))
    fwd_weight_tag <- gsub("\\.", "p", format(fwd_bin_weight, trim = TRUE))
    out_name <- paste0(out_name, "_gauss5y5yRho", rho_tag,
                       "_fwdBinW", fwd_weight_tag)
  }
  if (fwd_bounds_weight > 0) {
    bounds_tag <- gsub("\\.", "p", format(fwd_bounds_weight, trim = TRUE))
    out_name <- paste0(out_name, "_frechet5y5yW", bounds_tag)
  }
  out_name <- paste0(area, "_", out_name)
  out_dir <- file.path(getwd(), "outputs", out_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  res <- load_hrr_monthquart(area)
  dists <- load_hrr_dists(area)
  fwd_proxy_rho <- resolve_fwd_proxy_rho(fwd_proxy_rho_input, dists)
  targets <- make_finegrid_targets(
    res, dists,
    area = area,
    fwd_target_source = fwd_target_source,
    fwd_proxy_rho = fwd_proxy_rho
  )

  gs_model <- unwrap_matlab_cell(res$gs.model)
  fs_model <- unwrap_matlab_cell(res$fs.model)
  infl <- as.numeric(dists$infl)
  n_month <- dim(res$gs.data)[3]
  bts <- compute_bts(infl)[seq_len(n_month)]
  month_ids <- parse_month_ids(n_month)
  month_labels <- make_month_labels(n_month)

  eta_start_by_month <- NULL
  if (nzchar(eta_start_file)) {
    eta_start_df <- read.csv(eta_start_file, check.names = FALSE)
    eta_cols <- grep("^eta[0-9]+$", names(eta_start_df), value = TRUE)
    if (!("month" %in% names(eta_start_df)) || length(eta_cols) != length(eta0)) {
      stop(sprintf("ETA_START_FILE must contain month and eta1,...,eta%d columns.",
                   length(eta0)))
    }
    eta_start_by_month <- as.matrix(eta_start_df[, eta_cols])
    rownames(eta_start_by_month) <- as.character(eta_start_df$month)
  }

  message(sprintf("Running %s fine-grid/proxy 24-target fit for %d dates with %d cores.",
                  area, length(month_ids), n_cores))
  message(sprintf("Variant: %s.", variant))
  message(sprintf("Weights: side=%s, tail bins=%g, explicit tail moments=%g.",
                  tail_weight_side, tail_bin_weight, tail_moment_weight))
  message(sprintf("5y5y bin weight: %g.", fwd_bin_weight))
  message(sprintf("Objective includes 5y5y target: %s.",
                  ifelse(include_fwd_target, "yes", "no")))
  message(sprintf("5y5y target source: %s; rho=%.4f.",
                  fwd_target_source, fwd_proxy_rho))
  if (fwd_proxy_weight > 0) {
    message(sprintf("Objective includes 5y/10y moment-implied 5y5y >4 proxy with weight %g.",
                    fwd_proxy_weight))
  }
  if (fwd_bounds_weight > 0) {
    message(sprintf("Objective penalizes 5y5y >4 Frechet-bound violations with weight %g.",
                    fwd_bounds_weight))
  }
  if (!is.null(eta_start_by_month)) {
    message(sprintf("Warm starts: %s", eta_start_file))
  }

  worker <- function(month_id, verbose = n_cores <= 1) {
    if (verbose) {
      message(sprintf("Starting month %03d (%s)", month_id, month_labels[month_id]))
    }
    out <- fit_one_month_finegrid(
      month_id = month_id,
      targets = targets,
      gs_model = gs_model,
      fs_model = fs_model,
      bts = bts,
      nb_loop = nb_loop,
      nlminb_iter = nlminb_iter,
      nm_iter = nm_iter,
      variant = variant,
      mask = mask,
      eta0 = if (!is.null(eta_start_by_month) &&
                 as.character(month_id) %in% rownames(eta_start_by_month)) {
        unname(eta_start_by_month[as.character(month_id), ])
      } else {
        eta0
      },
      weight_spec = weight_spec,
      include_fwd_target = include_fwd_target,
      fwd_proxy_weight = fwd_proxy_weight,
      fwd_bounds_weight = fwd_bounds_weight
    )
    if (verbose) {
      message(sprintf("Finished month %03d: HRR %.3f pp, NGMMR %.3f pp.",
                      month_id,
                      out$diagnostics$rmse24_hrr_pp,
                      out$diagnostics$rmse24_ngmmr_pp))
    }
    out
  }

  progress_file <- file.path(out_dir, "diagnostics_in_progress.csv")
  write_progress <- function(results_so_far) {
    diagnostics_so_far <- do.call(rbind, lapply(results_so_far, `[[`, "diagnostics"))
    diagnostics_so_far$label <- month_labels[diagnostics_so_far$month]
    diagnostics_so_far$date <- as.Date(paste0(diagnostics_so_far$label, "-01"))
    diagnostics_so_far <- diagnostics_so_far[order(diagnostics_so_far$month), ]
    write.csv(diagnostics_so_far, progress_file, row.names = FALSE)
    invisible(diagnostics_so_far)
  }

  run_batch <- function(batch_ids) {
    if (.Platform$OS.type == "unix" && n_cores > 1) {
      parallel::mclapply(batch_ids, worker,
                         mc.cores = min(n_cores, length(batch_ids)),
                         mc.preschedule = FALSE)
    } else {
      lapply(batch_ids, worker)
    }
  }

  start_time <- proc.time()[["elapsed"]]
  batches <- split(month_ids, ceiling(seq_along(month_ids) / max(1, n_cores)))
  results <- vector("list", 0)
  message(sprintf("Progress will be written to %s.", progress_file))

  if (.Platform$OS.type == "unix" && n_cores > 1) {
    for (batch_id in seq_along(batches)) {
      batch_ids <- batches[[batch_id]]
      message(sprintf("Batch %03d/%03d: months %03d-%03d (%s to %s).",
                      batch_id, length(batches),
                      min(batch_ids), max(batch_ids),
                      month_labels[min(batch_ids)], month_labels[max(batch_ids)]))
      batch_results <- run_batch(batch_ids)
      results <- c(results, batch_results)
      diagnostics_so_far <- write_progress(results)
      elapsed_min <- (proc.time()[["elapsed"]] - start_time) / 60
      message(sprintf("Completed %d/%d dates in %.1f minutes. Mean NGMMR RMSE: %.3f pp.",
                      nrow(diagnostics_so_far), length(month_ids), elapsed_min,
                      mean(diagnostics_so_far$rmse24_ngmmr_pp)))
    }
  } else {
    for (batch_id in seq_along(batches)) {
      batch_ids <- batches[[batch_id]]
      message(sprintf("Batch %03d/%03d: month %03d (%s).",
                      batch_id, length(batches),
                      batch_ids, month_labels[batch_ids]))
      batch_results <- run_batch(batch_ids)
      results <- c(results, batch_results)
      diagnostics_so_far <- write_progress(results)
      elapsed_min <- (proc.time()[["elapsed"]] - start_time) / 60
      message(sprintf("Completed %d/%d dates in %.1f minutes. Mean NGMMR RMSE: %.3f pp.",
                      nrow(diagnostics_so_far), length(month_ids), elapsed_min,
                      mean(diagnostics_so_far$rmse24_ngmmr_pp)))
    }
  }

  diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))
  diagnostics$label <- month_labels[diagnostics$month]
  diagnostics$date <- as.Date(paste0(diagnostics$label, "-01"))
  diagnostics <- diagnostics[order(diagnostics$month), ]

  saveRDS(results, file.path(out_dir, "all_results.rds"))
  save_matrix_stack(results, "A", file.path(out_dir, "ngmmr_Q_matrices.rds"))
  write.csv(diagnostics, file.path(out_dir, "diagnostics.csv"), row.names = FALSE)

  summary_stats <- data.frame(
    metric = c("RMSE 5y+10y+5y5y", "RMSE 5y", "RMSE 10y", "RMSE 5y5y",
               "Tail <0 5y", "Tail >4 5y",
               "Tail <0 10y", "Tail >4 10y",
               "Tail <0 5y5y", "Tail >4 5y5y",
               "Tail >4 5y5y moment proxy"),
    mean_hrr_pp = c(mean(diagnostics$rmse24_hrr_pp),
                    mean(diagnostics$rmse5_hrr_pp),
                    mean(diagnostics$rmse10_hrr_pp),
                    mean(diagnostics$rmse_fwd_hrr_pp),
                    rmse_pp(diagnostics$qtail0_5y_hrr -
                              diagnostics$qtail0_5y_target),
                    rmse_pp(diagnostics$qtail4_5y_hrr -
                              diagnostics$qtail4_5y_target),
                    rmse_pp(diagnostics$qtail0_10y_hrr -
                              diagnostics$qtail0_10y_target),
                    rmse_pp(diagnostics$qtail4_10y_hrr -
                              diagnostics$qtail4_10y_target),
                    rmse_pp(diagnostics$qtail0_5y5y_hrr -
                              diagnostics$qtail0_5y5y_target),
                    rmse_pp(diagnostics$qtail4_5y5y_hrr -
                              diagnostics$qtail4_5y5y_target),
                    rmse_pp(diagnostics$qtail4_5y5y_hrr -
                              diagnostics$qtail4_5y5y_moment_proxy)),
    mean_ngmmr_pp = c(mean(diagnostics$rmse24_ngmmr_pp),
                      mean(diagnostics$rmse5_ngmmr_pp),
                      mean(diagnostics$rmse10_ngmmr_pp),
                      mean(diagnostics$rmse_fwd_ngmmr_pp),
                      rmse_pp(diagnostics$qtail0_5y_ngmmr -
                                diagnostics$qtail0_5y_target),
                      rmse_pp(diagnostics$qtail4_5y_ngmmr -
                                diagnostics$qtail4_5y_target),
                      rmse_pp(diagnostics$qtail0_10y_ngmmr -
                                diagnostics$qtail0_10y_target),
                      rmse_pp(diagnostics$qtail4_10y_ngmmr -
                                diagnostics$qtail4_10y_target),
                      rmse_pp(diagnostics$qtail0_5y5y_ngmmr -
                                diagnostics$qtail0_5y5y_target),
                      rmse_pp(diagnostics$qtail4_5y5y_ngmmr -
                                diagnostics$qtail4_5y5y_target),
                      rmse_pp(diagnostics$qtail4_5y5y_ngmmr -
                                diagnostics$qtail4_5y5y_moment_proxy))
  )
  message("\nDone.")
  message(sprintf("Diagnostics: %s", file.path(out_dir, "diagnostics.csv")))
  print(summary_stats, digits = 4)
}

if (sys.nframe() == 0) {
  main()
}
