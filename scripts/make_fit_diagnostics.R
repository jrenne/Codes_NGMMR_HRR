## Fit diagnostics:
## 1. HRR reported P vs model-consistent P using HRR's own matrices.
## 2. Nominal-Q high-tail fit: option-implied target, HRR fit, NGMMR fit.
## 3. Compact RMSE table for 5y and 10y nominal-Q bin probabilities.

suppressPackageStartupMessages({
  library(haven)
  library(optimx)
})

source(file.path("scripts", "make_ratio_and_probability_figures.R"))
source(file.path("scripts", "estimate_nominal_q.R"))

rmse_pp <- function(x) {
  100 * sqrt(mean(x^2, na.rm = TRUE))
}

tail_idx <- function(threshold) {
  if (threshold == 4) {
    7:8
  } else if (threshold == 5) {
    8
  } else {
    stop("Only thresholds 4 and 5 are implemented.")
  }
}

make_hrr_consistent_p_series <- function(area, psi = 1.5) {
  area <- normalize_hrr_area(area)
  hrr_p <- read_hrr_p(area)

  series <- compute_area_ies_series(area, psi_values = c(1 / 3, 1.5),
                                    matrix_source = "hrr")
  series <- series[abs(series$psi - psi) < 1e-10, ]
  months <- sort(unique(series$month))
  rows <- lapply(months, function(k) {
    sk <- series[series$month == k, ]
    out <- data.frame(
      area = area,
      month = k,
      date = as.Date(sk$date[1]),
      psi = psi,
      ez_converged = all(sk$ez_converged)
    )
    for (threshold in c(4, 5)) {
      st <- sk[sk$threshold == threshold, ]
      out[[paste0("p_model_5y_gt", threshold)]] <- st$p_5y
      out[[paste0("p_model_10y_gt", threshold)]] <- st$p_10y
      out[[paste0("p_model_5y5y_gt", threshold)]] <-
        st$p_5y5y
      out[[paste0("p_hrr_5y_gt", threshold)]] <-
        as.numeric(hrr_p[[hrr_p_column("5y", threshold)]][k])
      out[[paste0("p_hrr_10y_gt", threshold)]] <-
        as.numeric(hrr_p[[hrr_p_column("10y", threshold)]][k])
      out[[paste0("p_hrr_5y5y_gt", threshold)]] <-
        as.numeric(hrr_p[[hrr_p_column("5y5y", threshold)]][k])
    }
    out
  })

  do.call(rbind, rows)
}

ez_disaster_inputs <- function(n_state) {
  beta <- 0.99
  gamma <- 3
  g <- 0.02
  z0_h <- 1.03
  alpha_h <- 5.45
  z0_l <- 1.06
  alpha_l <- 15.18
  p_bar_h <- 0.356
  p_bar_l <- 0.085

  kappa <- gamma - 1
  g_tilde_l <- exp(-kappa * g) *
    (p_bar_l * alpha_l / (alpha_l - kappa) * z0_l^kappa + 1 - p_bar_l)
  g_tilde_h <- exp(-kappa * g) *
    (p_bar_h * alpha_h / (alpha_h - kappa) * z0_h^kappa + 1 - p_bar_h)
  kappa <- gamma
  g_star_l <- exp(-kappa * g) *
    (p_bar_l * alpha_l / (alpha_l - kappa) * z0_l^kappa + 1 - p_bar_l)
  g_star_h <- exp(-kappa * g) *
    (p_bar_h * alpha_h / (alpha_h - kappa) * z0_h^kappa + 1 - p_bar_h)

  list(
    beta = beta,
    gamma = gamma,
    g_tilde = c(g_tilde_l, rep(exp((1 - gamma) * g), n_state - 2L), g_tilde_h),
    g_star = c(g_star_l, rep(exp(-gamma * g), n_state - 2L), g_star_h)
  )
}

consistent_p_tails_from_qn <- function(Qn, bt, psi = 1.5,
                                       pi_values = extended_pi_grid(nrow(Qn))) {
  ez <- ez_disaster_inputs(nrow(Qn))
  Qr <- nominal_to_real_q(Qn, pi_values)
  sol <- solve_ez_fixed_point_from_q(
    Qr, ez$g_tilde, ez$g_star,
    beta = ez$beta, gamma = ez$gamma, psi = psi
  )
  P <- sol$P
  rows <- list()
  row_id <- 1L
  for (threshold in c(4, 5)) {
    tp <- tail_prob_5_10_threshold(P, bt, threshold, pi_bar = pi_values)
    rows[[row_id]] <- data.frame(
      threshold = threshold,
      p_5y = tp[["y5"]],
      p_10y = tp[["y10"]],
      p_5y5y = tail_prob_5y5y_threshold(P, bt, threshold, pi_bar = pi_values),
      ez_converged = sol$converged
    )
    row_id <- row_id + 1L
  }
  do.call(rbind, rows)
}

make_nominal_q_fit_series <- function(area) {
  area <- normalize_hrr_area(area)
  results <- readRDS(file.path(default_out_dir(area), "all_results.rds"))
  diagnostics <- read.csv(file.path(default_out_dir(area), "diagnostics.csv"))
  diagnostics$date <- as.Date(diagnostics$date)

  rows <- list()
  row_id <- 1L
  for (i in seq_along(results)) {
    for (threshold in c(4, 5)) {
      idx <- tail_idx(threshold)
      for (h in c("5y", "10y")) {
        slot <- if (h == "5y") "g5" else "g10"
        rows[[row_id]] <- data.frame(
          area = area,
          month = diagnostics$month[i],
          date = diagnostics$date[i],
          threshold = threshold,
          horizon = h,
          q_data = sum(results[[i]]$target[[slot]][idx]),
          q_hrr_model = sum(results[[i]]$hrr[[slot]][idx]),
          q_ngmmr_model = sum(results[[i]]$ngmmr[[slot]][idx])
        )
        row_id <- row_id + 1L
      }
    }
  }

  do.call(rbind, rows)
}

fit_hrr_form_one_month <- function(result, bt, eta0,
                                   weight_spec,
                                   nlminb_iter = 60) {
  obj <- function(eta) {
    A <- amatrix_hrr_101(z_to_hrr_theta(eta))
    mod <- avg_inf_510_fast(A, bt)
    fwd <- forward_average_bins_model(A, bt)
    weighted_finegrid_loss(
      model_g5 = mod$g5,
      model_g10 = mod$g10,
      model_fwd = fwd,
      target_g5 = result$target$g5,
      target_g10 = result$target$g10,
      target_fwd = result$target$fwd,
      weight_spec = weight_spec,
      include_fwd = FALSE,
      fwd_tail4_proxy = NA_real_,
      fwd_tail4_proxy_weight = 0
    )
  }

  fit <- nlminb(
    start = eta0,
    objective = obj,
    control = list(iter.max = nlminb_iter,
                   eval.max = max(200, 20 * nlminb_iter))
  )
  A <- amatrix_hrr_101(z_to_hrr_theta(fit$par))
  mod <- avg_inf_510_fast(A, bt)
  fwd <- forward_average_bins_model(A, bt)
  list(par = fit$par, value = fit$objective, g5 = mod$g5, g10 = mod$g10,
       fwd = fwd, A = A)
}

make_hrr_form_consistent_p_series <- function(area, output_dir,
                                              psi = 1.5,
                                              force = FALSE) {
  area <- normalize_hrr_area(area)
  cache_dir <- file.path(getwd(), "cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  model_tag <- basename(default_out_dir(area))
  cache_file <- file.path(
    cache_dir,
    paste0(area, "_hrr_form_consistent_p_",
           model_tag, "_psi", gsub("[.]", "p", as.character(psi)), ".csv")
  )
  if (file.exists(cache_file) && !force) {
    out <- read.csv(cache_file)
    out$date <- as.Date(out$date)
    return(out)
  }

  results <- readRDS(file.path(default_out_dir(area), "all_results.rds"))
  diagnostics <- read.csv(file.path(default_out_dir(area), "diagnostics.csv"))
  diagnostics$date <- as.Date(diagnostics$date)

  res <- load_hrr_monthquart(area)
  dists <- load_hrr_dists(area)
  pinew <- unwrap_matlab_cell(res$pinew)
  bts <- compute_bts(as.numeric(dists$infl))[seq_len(ncol(pinew))]
  weight_spec <- make_tail_weight_spec(
    tail_bin_weight = 2,
    tail_moment_weight = 2,
    fwd_bin_weight = 0,
    fwd_tail_bin_weight = 0,
    tail_weight_side = "high"
  )

  n_core <- min(8L, max(1L, parallel::detectCores() - 1L))
  month_rows <- parallel::mclapply(seq_along(results), function(i) {
    eta0 <- hrr_theta_to_z(pinew[, i])
    fit <- fit_hrr_form_one_month(results[[i]], bts[i], eta0, weight_spec)
    p_rows <- consistent_p_tails_from_qn(fit$A, bts[i], psi = psi)
    p_rows$area <- area
    p_rows$month <- diagnostics$month[i]
    p_rows$date <- diagnostics$date[i]
    p_rows$psi <- psi
    p_rows
  }, mc.cores = n_core)

  out <- do.call(rbind, month_rows)
  out <- out[order(out$area, out$threshold, out$date), ]
  write.csv(out, cache_file, row.names = FALSE)
  out
}

make_ngmmr_consistent_p_series <- function(area, psi = 1.5) {
  s <- compute_area_ies_series(area, psi_values = c(1 / 3, 1.5),
                               matrix_source = "ngmmr")
  s <- s[abs(s$psi - psi) < 1e-10, ]
  s <- s[, c("area", "month", "date", "psi", "threshold",
             "p_5y", "p_10y", "p_5y5y", "ez_converged")]
  s$date <- as.Date(s$date)
  s
}

plot_consistent_p_comparison <- function(hrr_series, hrr_form_series,
                                         ngmmr_series, threshold,
                                         output_file) {
  areas <- c("US", "EZ")
  area_names <- c(US = "U.S.", EZ = "Euro area")
  horizons <- c("5y", "10y", "5y5y")

  hrr_form_series <- hrr_form_series[hrr_form_series$threshold == threshold, ]
  ngmmr_series <- ngmmr_series[ngmmr_series$threshold == threshold, ]

  pdf(output_file, pointsize = 15, width = 12.5, height = 7.4)
  oldpar <- par(mfrow = c(2, 3), mar = c(3.8, 4.9, 2.8, 1.0),
                oma = c(3.2, 0, 0, 0),
                cex.axis = 1.15, cex.lab = 1.22, cex.main = 1.2)
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  for (area in areas) {
    s_hrr <- hrr_series[hrr_series$area == area, ]
    s_hrr <- s_hrr[order(s_hrr$date), ]
    s_form <- hrr_form_series[hrr_form_series$area == area, ]
    s_form <- s_form[order(s_form$date), ]
    s_ngmmr <- ngmmr_series[ngmmr_series$area == area, ]
    s_ngmmr <- s_ngmmr[order(s_ngmmr$date), ]

    for (h in horizons) {
      model_col <- paste0("p_model_", h, "_gt", threshold)
      reported_col <- paste0("p_hrr_", h, "_gt", threshold)
      p_col <- paste0("p_", h)
      ymax <- 1.15 * max(100 * c(s_hrr[[model_col]],
                                  s_hrr[[reported_col]],
                                  s_form[[p_col]],
                                  s_ngmmr[[p_col]]),
                         na.rm = TRUE)
      plot(s_hrr$date, 100 * s_hrr[[model_col]], type = "l",
           col = "black", lwd = 2.2, ylim = c(0, ymax), xlab = "",
           ylab = "Probability (%)", las = 1,
           main = sprintf("%s, %s", area_names[[area]], h))
      grid()
      lines(s_form$date, 100 * s_form[[p_col]],
            col = "gray45", lwd = 2.0, lty = 2)
      lines(s_ngmmr$date, 100 * s_ngmmr[[p_col]],
            col = "black", lwd = 2.0, lty = 3)
      lines(s_hrr$date, 100 * s_hrr[[reported_col]],
            col = "gray65", lwd = 2.4, lty = 1)
    }
  }
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new()
  legend("bottom",
         legend = c("HRR reported P",
                    "HRR matrices, consistent P",
                    "HRR-form refit, consistent P",
                    "NGMMR refit, consistent P"),
         col = c("gray65", "black", "gray45", "black"),
         lty = c(1, 1, 2, 3), lwd = c(2.4, 2.2, 2.0, 2.0),
         bg = "white", box.col = "gray80", cex = 0.95,
         horiz = TRUE, inset = 0.01, xpd = NA)
}

make_hrr_form_refit_series <- function(area, output_dir,
                                       force = FALSE) {
  area <- normalize_hrr_area(area)
  cache_dir <- file.path(getwd(), "cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  model_tag <- basename(default_out_dir(area))
  cache_file <- file.path(cache_dir,
                          paste0(area, "_hrr_form_refit_series_",
                                 model_tag, ".csv"))
  if (file.exists(cache_file) && !force) {
    out <- read.csv(cache_file)
    out$date <- as.Date(out$date)
    return(out)
  }

  results <- readRDS(file.path(default_out_dir(area), "all_results.rds"))
  diagnostics <- read.csv(file.path(default_out_dir(area), "diagnostics.csv"))
  diagnostics$date <- as.Date(diagnostics$date)

  res <- load_hrr_monthquart(area)
  dists <- load_hrr_dists(area)
  pinew <- unwrap_matlab_cell(res$pinew)
  bts <- compute_bts(as.numeric(dists$infl))[seq_len(ncol(pinew))]
  weight_spec <- make_tail_weight_spec(
    tail_bin_weight = 2,
    tail_moment_weight = 2,
    fwd_bin_weight = 0,
    fwd_tail_bin_weight = 0,
    tail_weight_side = "high"
  )

  n_core <- min(8L, max(1L, parallel::detectCores() - 1L))
  month_rows <- parallel::mclapply(seq_along(results), function(i) {
    eta0 <- hrr_theta_to_z(pinew[, i])
    fit <- fit_hrr_form_one_month(results[[i]], bts[i], eta0, weight_spec)
    rows_i <- list()
    row_id <- 1L
    for (threshold in c(4, 5)) {
      idx <- tail_idx(threshold)
      for (h in c("5y", "10y", "5y5y")) {
        slot <- switch(h, "5y" = "g5", "10y" = "g10", "5y5y" = "fwd")
        rows_i[[row_id]] <- data.frame(
          area = area,
          month = diagnostics$month[i],
          date = diagnostics$date[i],
          threshold = threshold,
          horizon = h,
          q_hrr_form_refit = sum(fit[[slot]][idx]),
          objective_hrr_form_refit = fit$value
        )
        row_id <- row_id + 1L
      }
    }
    do.call(rbind, rows_i)
  }, mc.cores = n_core)

  out <- do.call(rbind, month_rows)
  write.csv(out, cache_file, row.names = FALSE)
  out
}

make_hrr_form_refit_bins <- function(area, output_dir,
                                     force = FALSE) {
  area <- normalize_hrr_area(area)
  cache_dir <- file.path(getwd(), "cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  model_tag <- basename(default_out_dir(area))
  cache_file <- file.path(cache_dir,
                          paste0(area, "_hrr_form_refit_bins_",
                                 model_tag, ".csv"))
  if (file.exists(cache_file) && !force) {
    out <- read.csv(cache_file)
    out$date <- as.Date(out$date)
    return(out)
  }

  results <- readRDS(file.path(default_out_dir(area), "all_results.rds"))
  diagnostics <- read.csv(file.path(default_out_dir(area), "diagnostics.csv"))
  diagnostics$date <- as.Date(diagnostics$date)

  res <- load_hrr_monthquart(area)
  dists <- load_hrr_dists(area)
  pinew <- unwrap_matlab_cell(res$pinew)
  bts <- compute_bts(as.numeric(dists$infl))[seq_len(ncol(pinew))]
  weight_spec <- make_tail_weight_spec(
    tail_bin_weight = 2,
    tail_moment_weight = 2,
    fwd_bin_weight = 0,
    fwd_tail_bin_weight = 0,
    tail_weight_side = "high"
  )

  n_core <- min(8L, max(1L, parallel::detectCores() - 1L))
  month_rows <- parallel::mclapply(seq_along(results), function(i) {
    eta0 <- hrr_theta_to_z(pinew[, i])
    fit <- fit_hrr_form_one_month(results[[i]], bts[i], eta0, weight_spec)
    rows_i <- list()
    row_id <- 1L
    for (h in c("5y", "10y", "5y5y")) {
      slot <- switch(h, "5y" = "g5", "10y" = "g10", "5y5y" = "fwd")
      for (b in seq_along(fit[[slot]])) {
        rows_i[[row_id]] <- data.frame(
          area = area,
          month = diagnostics$month[i],
          date = diagnostics$date[i],
          horizon = h,
          bin = b,
          q_hrr_form_refit_bin = fit[[slot]][b]
        )
        row_id <- row_id + 1L
      }
    }
    do.call(rbind, rows_i)
  }, mc.cores = n_core)

  out <- do.call(rbind, month_rows)
  write.csv(out, cache_file, row.names = FALSE)
  out
}

plot_hrr_consistent_p <- function(series, threshold, output_file) {
  areas <- c("US", "EZ")
  area_names <- c(US = "U.S.", EZ = "Euro area")
  horizons <- c("5y", "10y", "5y5y")

  pdf(output_file, pointsize = 15, width = 12.5, height = 7.4)
  oldpar <- par(mfrow = c(2, 3), mar = c(3.8, 4.9, 2.8, 1.0),
                oma = c(3.2, 0, 0, 0),
                cex.axis = 1.15, cex.lab = 1.22, cex.main = 1.2)
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  for (area in areas) {
    s_area <- series[series$area == area, ]
    s_area <- s_area[order(s_area$date), ]
    for (h in horizons) {
      model_col <- paste0("p_model_", h, "_gt", threshold)
      hrr_col <- paste0("p_hrr_", h, "_gt", threshold)
      ymax <- 1.15 * max(100 * c(s_area[[model_col]], s_area[[hrr_col]]),
                         na.rm = TRUE)
      plot(s_area$date, 100 * s_area[[model_col]], type = "l",
           col = "black", lwd = 2, ylim = c(0, ymax), xlab = "",
           ylab = "Probability (%)", las = 1,
           main = sprintf("%s, %s", area_names[[area]], h))
      grid()
      lines(s_area$date, 100 * s_area[[hrr_col]],
            col = "gray45", lwd = 2, lty = 3)
    }
  }
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new()
  legend("bottom",
         legend = c("Model-consistent P", "HRR reported P"),
         col = c("black", "gray45"), lty = c(1, 3), lwd = 2,
         bg = "white", box.col = "gray80", cex = 1.08,
         horiz = TRUE, inset = 0.01, xpd = NA)
}

plot_nominal_q_fit <- function(series, threshold, output_file) {
  areas <- c("US", "EZ")
  area_names <- c(US = "U.S.", EZ = "Euro area")
  horizons <- c("5y", "10y")

  pdf(output_file, pointsize = 15, width = 10.8, height = 7.6)
  oldpar <- par(mfrow = c(2, 2), mar = c(3.8, 5.0, 2.8, 1.0),
                oma = c(3.2, 0, 0, 0),
                cex.axis = 1.15, cex.lab = 1.22, cex.main = 1.2)
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  for (area in areas) {
    for (h in horizons) {
      s <- series[series$area == area &
                    series$horizon == h &
                    series$threshold == threshold, ]
      s <- s[order(s$date), ]
      ymax <- 1.15 * max(100 * c(s$q_data, s$q_hrr_model,
                                  s$q_hrr_form_refit, s$q_ngmmr_model),
                         na.rm = TRUE)
      plot(s$date, 100 * s$q_ngmmr_model, type = "n",
           ylim = c(0, ymax), xlab = "", ylab = "Probability (%)",
           las = 1, main = sprintf("%s, %s", area_names[[area]], h))
      grid()
      lines(s$date, 100 * s$q_hrr_model, col = "gray70", lwd = 2, lty = 1)
      lines(s$date, 100 * s$q_hrr_form_refit,
            col = "gray45", lwd = 1.8, lty = 3)
      lines(s$date, 100 * s$q_ngmmr_model, col = "black", lwd = 2.4, lty = 1)
      points(s$date, 100 * s$q_data, col = "gray20", pch = 4, cex = 0.75,
             lwd = 1.2)
    }
  }
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new()
  legend("bottom",
         legend = c("Option-implied Q", "HRR fitted Q",
                    "HRR-form refit", "NGMMR fitted Q"),
         col = c("gray20", "gray70", "gray45", "black"),
         lty = c(NA, 1, 3, 1), lwd = c(NA, 2, 1.8, 2.4),
         pch = c(4, NA, NA, NA), pt.cex = c(0.85, NA, NA, NA),
         bg = "white", box.col = "gray80", cex = 1.08,
         horiz = TRUE, inset = 0.01, xpd = NA)
}

make_full_bin_rmse_rows <- function(output_dir) {
  rows <- list()
  row_id <- 1L
  for (area in c("US", "EZ")) {
    results <- readRDS(file.path(default_out_dir(area), "all_results.rds"))
    hrr_form_bins <- make_hrr_form_refit_bins(area, output_dir)
    for (h in c("5y", "10y", "All")) {
      hrr_err <- numeric(0)
      hrr_form_err <- numeric(0)
      ngmmr_err <- numeric(0)
      for (i in seq_along(results)) {
        result <- results[[i]]
        slots <- switch(h,
          "5y" = "g5",
          "10y" = "g10",
          "All" = c("g5", "g10")
        )
        for (slot in slots) {
          h_label <- switch(slot, "g5" = "5y", "g10" = "10y")
          hrr_form_i <- hrr_form_bins[
            hrr_form_bins$month == i & hrr_form_bins$horizon == h_label,
          ]
          hrr_form_i <- hrr_form_i[order(hrr_form_i$bin), ]
          hrr_err <- c(hrr_err, result$hrr[[slot]] - result$target[[slot]])
          hrr_form_err <- c(hrr_form_err,
                            hrr_form_i$q_hrr_form_refit_bin -
                              result$target[[slot]])
          ngmmr_err <- c(ngmmr_err, result$ngmmr[[slot]] - result$target[[slot]])
        }
      }
      rows[[row_id]] <- data.frame(
        area = area,
        horizon = h,
        object = "8-bin distribution",
        hrr_rmse_pp = rmse_pp(hrr_err),
        hrr_form_refit_rmse_pp = rmse_pp(hrr_form_err),
        ngmmr_rmse_pp = rmse_pp(ngmmr_err)
      )
      row_id <- row_id + 1L
    }
  }
  do.call(rbind, rows)
}

make_rmse_table <- function(qfit_series, output_dir) {
  rows <- list()
  row_id <- 1L
  for (area in c("US", "EZ")) {
    for (h in c("5y", "10y")) {
      for (threshold in c(4, 5)) {
        s <- qfit_series[qfit_series$area == area &
                           qfit_series$horizon == h &
                           qfit_series$threshold == threshold, ]
        rows[[row_id]] <- data.frame(
          area = area,
          horizon = h,
          object = paste0("$>", threshold, "\\%$ tail"),
          hrr_rmse_pp = rmse_pp(s$q_hrr_model - s$q_data),
          hrr_form_refit_rmse_pp =
            rmse_pp(s$q_hrr_form_refit - s$q_data),
          ngmmr_rmse_pp = rmse_pp(s$q_ngmmr_model - s$q_data)
        )
        row_id <- row_id + 1L
      }
    }
  }
  full_table <- make_full_bin_rmse_rows(output_dir)
  tail_table <- do.call(rbind, rows)

  fmt <- function(x) {
    ifelse(is.na(x), "--", sprintf("%.2f", x))
  }
  full_tex_rows <- sprintf("%s & %s & %s & %.2f & %s & %.2f \\\\",
                           full_table$area, full_table$horizon,
                           full_table$object, full_table$hrr_rmse_pp,
                           fmt(full_table$hrr_form_refit_rmse_pp),
                           full_table$ngmmr_rmse_pp)
  tail_tex_rows <- sprintf("%s & %s & %s & %.2f & %.2f & %.2f \\\\",
                           tail_table$area, tail_table$horizon,
                           tail_table$object, tail_table$hrr_rmse_pp,
                           tail_table$hrr_form_refit_rmse_pp,
                           tail_table$ngmmr_rmse_pp)
  tex <- c(
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}llcccc@{}}",
    "\\toprule",
    "Area & Horizon & Object & HRR & HRR-form & NGMMR \\\\",
    "\\midrule",
    "\\multicolumn{6}{l}{\\emph{Panel A: full binned distributions}}\\\\",
    full_tex_rows,
    "\\addlinespace",
    "\\multicolumn{6}{l}{\\emph{Panel B: displayed high-inflation tails}}\\\\",
    tail_tex_rows,
    "\\bottomrule",
    "\\end{tabular*}"
  )
  writeLines(tex, file.path(output_dir, "nominal_Q_tail_rmse_table.tex"))
  rbind(full_table, tail_table)
}

main <- function() {
  out_dir <- file.path(getwd(), "outputs", "comment_figures")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  hrr_p_series <- do.call(rbind, lapply(c("US", "EZ"),
                                        make_hrr_consistent_p_series))
  qfit_series <- do.call(rbind, lapply(c("US", "EZ"),
                                       make_nominal_q_fit_series))
  hrr_form_refit_series <- do.call(
    rbind,
    lapply(c("US", "EZ"), make_hrr_form_refit_series, output_dir = out_dir)
  )
  hrr_form_p_series <- do.call(
    rbind,
    lapply(c("US", "EZ"), make_hrr_form_consistent_p_series,
           output_dir = out_dir)
  )
  ngmmr_p_series <- do.call(
    rbind,
    lapply(c("US", "EZ"), make_ngmmr_consistent_p_series)
  )
  qfit_series <- merge(
    qfit_series,
    hrr_form_refit_series[
      c("area", "month", "date", "threshold", "horizon",
        "q_hrr_form_refit")
    ],
    by = c("area", "month", "date", "threshold", "horizon"),
    all.x = TRUE
  )
  qfit_series <- qfit_series[
    order(qfit_series$area, qfit_series$threshold,
          qfit_series$horizon, qfit_series$date),
  ]

  plot_consistent_p_comparison(
    hrr_p_series, hrr_form_p_series, ngmmr_p_series, 4,
    file.path(out_dir, "model_consistent_P_comparison_gt4.pdf")
  )
  plot_nominal_q_fit(
    qfit_series, 4,
    file.path(out_dir, "nominal_Q_tail_fit_HRR_NGMMR_gt4.pdf")
  )

  rmse_table <- make_rmse_table(qfit_series, out_dir)
  print(rmse_table)

  message("Done.")
  message(sprintf("Outputs written to: %s", out_dir))
}

if (sys.nframe() == 0) {
  main()
}
