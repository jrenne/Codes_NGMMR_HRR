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
  res <- load_hrr_monthquart(area)
  dists <- load_hrr_dists(area)
  hrr_p <- read_hrr_p(area)

  pinew <- unwrap_matlab_cell(res$pinew)
  infl <- as.numeric(dists$infl)
  bts <- compute_bts(infl)[seq_len(ncol(pinew))]

  pi_values <- c(-2, seq(-0.5, 4.5, by = 1), 6)
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

  g_tilde <- c(g_tilde_l, rep(exp((1 - gamma) * g), 6), g_tilde_h)
  g_star <- c(g_star_l, rep(exp(-gamma * g), 6), g_star_h)

  rows <- lapply(seq_len(ncol(pinew)), function(k) {
    Qn <- amatrix_hrr_101(pinew[, k])
    Qr <- nominal_to_real_q(Qn, pi_values)
    sol <- solve_ez_fixed_point_from_q(Qr, g_tilde, g_star,
                                       beta = beta, gamma = gamma, psi = psi)
    P <- sol$P
    bt <- bts[k]

    out <- data.frame(
      area = area,
      month = k,
      date = as.Date(hrr_p$date_stata[k]),
      psi = psi,
      ez_converged = sol$converged
    )
    for (threshold in c(4, 5)) {
      tp <- tail_prob_5_10_threshold(P, bt, threshold)
      out[[paste0("p_model_5y_gt", threshold)]] <- tp[["y5"]]
      out[[paste0("p_model_10y_gt", threshold)]] <- tp[["y10"]]
      out[[paste0("p_model_5y5y_gt", threshold)]] <-
        tail_prob_5y5y_threshold(P, bt, threshold)
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
    fwd <- avg_yoy_6to10_model(A, bt)
    weighted_finegrid_loss(
      model_g5 = mod$g5,
      model_g10 = mod$g10,
      model_fwd = fwd,
      target_g5 = result$target$g5,
      target_g10 = result$target$g10,
      target_fwd = result$target$fwd,
      weight_spec = weight_spec,
      include_fwd = TRUE,
      fwd_tail4_proxy = sum(result$target$fwd[7:8]),
      fwd_tail4_proxy_weight = 0.2
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
  list(par = fit$par, value = fit$objective, g5 = mod$g5, g10 = mod$g10)
}

make_hrr_form_refit_series <- function(area, output_dir,
                                       force = FALSE) {
  area <- normalize_hrr_area(area)
  cache_dir <- file.path(getwd(), "cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(cache_dir, paste0(area, "_hrr_form_refit_series.csv"))
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
    tail_bin_weight = 1.5,
    tail_moment_weight = 2,
    fwd_bin_weight = 1,
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
      for (h in c("5y", "10y")) {
        slot <- if (h == "5y") "g5" else "g10"
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
          threshold = paste0(">", threshold, "\\%"),
          hrr_rmse_pp = rmse_pp(s$q_hrr_model - s$q_data),
          hrr_form_refit_rmse_pp =
            rmse_pp(s$q_hrr_form_refit - s$q_data),
          ngmmr_rmse_pp = rmse_pp(s$q_ngmmr_model - s$q_data)
        )
        row_id <- row_id + 1L
      }
    }
  }
  table <- do.call(rbind, rows)

  tex <- c(
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}llcccc@{}}",
    "\\toprule",
    "Area & Horizon & Threshold & HRR & HRR-form & NGMMR \\\\",
    "\\midrule",
    sprintf("%s & %s & %s & %.2f & %.2f & %.2f \\\\",
            table$area, table$horizon, table$threshold,
            table$hrr_rmse_pp, table$hrr_form_refit_rmse_pp,
            table$ngmmr_rmse_pp),
    "\\bottomrule",
    "\\end{tabular*}"
  )
  writeLines(tex, file.path(output_dir, "nominal_Q_tail_rmse_table.tex"))
  table
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

  plot_hrr_consistent_p(
    hrr_p_series, 4,
    file.path(out_dir, "HRR_model_consistent_P_gt4.pdf")
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
