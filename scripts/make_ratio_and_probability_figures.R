## Ratio and physical-probability figures used in the comment.

suppressPackageStartupMessages({
  library(haven)
})

source(file.path("scripts", "pq_mapping_helpers.R"))

figure_n_cores <- function() {
  if (exists("N_CORES", inherits = TRUE)) {
    max(1L, as.integer(get("N_CORES", inherits = TRUE)))
  } else {
    max(1L, as.integer(Sys.getenv("N_CORES", "1")))
  }
}

default_out_dir <- function(area) {
  file.path(getwd(), "outputs",
            paste0(area, "_nominal_Q_refit_smooth_row_poly3_finegrid_targets_",
                   "hightailbin2_tailmoment2_",
                   "spot_only_",
                   "lrMean2p5W10_lrExtremeCap20W5_thetaShrinkW0p05"))
}

ies_series_cache_file <- function(area, matrix_source, psi_values) {
  cache_dir <- file.path(getwd(), "cache")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  source_tag <- if (matrix_source == "ngmmr") {
    basename(default_out_dir(area))
  } else {
    "hrr_matrices"
  }
  psi_tag <- paste0("psi", paste(gsub("[.]", "p", as.character(psi_values)),
                                 collapse = "_"))
  file.path(cache_dir, paste0(area, "_ies_series_", source_tag, "_",
                              psi_tag, ".csv"))
}

tail_bins <- function(threshold) {
  if (threshold == 4) {
    7:8
  } else if (threshold == 5) {
    8
  } else {
    stop("Only thresholds 4 and 5 are implemented.")
  }
}

tail_prob_5_10_threshold <- function(A, bt, threshold,
                                     pi_bar = extended_pi_grid(nrow(A))) {
  bins <- avg_inf_510_fast(A, bt, pi_bar = pi_bar)
  idx <- tail_bins(threshold)
  c(
    y5 = sum(bins$g5[idx]),
    y10 = sum(bins$g10[idx])
  )
}

tail_prob_5y5y_threshold <- function(A, bt, threshold,
                                     pi_bar = extended_pi_grid(nrow(A))) {
  bins <- forward_average_bins_model(A, bt, pi_bar = pi_bar)
  sum(bins[tail_bins(threshold)])
}

compute_area_ies_series <- function(area, psi_values = c(1 / 3, 1.5),
                                    matrix_source = c("ngmmr", "hrr"),
                                    force = FALSE) {
  area <- normalize_hrr_area(area)
  matrix_source <- match.arg(matrix_source)
  cache_file <- ies_series_cache_file(area, matrix_source, psi_values)
  if (file.exists(cache_file) && !force) {
    out <- read.csv(cache_file)
    out$date <- as.Date(out$date)
    return(out)
  }

  res <- load_hrr_monthquart(area)
  dists <- load_hrr_dists(area)

  if (matrix_source == "ngmmr") {
    out_dir <- default_out_dir(area)
    q_file <- file.path(out_dir, "ngmmr_Q_matrices.rds")
    Qn_array <- readRDS(q_file)
    diagnostics <- read.csv(file.path(out_dir, "diagnostics.csv"))
    diagnostics$date <- as.Date(diagnostics$date)
  } else {
    pinew <- unwrap_matlab_cell(res$pinew)
    Qn_array <- array(NA_real_, dim = c(8, 8, ncol(pinew)))
    for (k in seq_len(ncol(pinew))) {
      Qn_array[, , k] <- amatrix_hrr_101(pinew[, k])
    }
    hrr_p <- read_hrr_p(area)
    diagnostics <- data.frame(
      month = seq_len(ncol(pinew)),
      date = as.Date(hrr_p$date_stata[seq_len(ncol(pinew))])
    )
  }

  infl <- as.numeric(dists$infl)
  pi_values <- extended_pi_grid(dim(Qn_array)[1])
  bts <- compute_bts(infl)[seq_len(dim(res$gs.data)[3])]
  sprt_z <- as.numeric(dists$sprt.Z[, 1])
  years <- as.integer(dists$years[, 1])
  last_month <- as.integer(dists$last.month[1, 1])
  year_month_index <- do.call(
    rbind,
    lapply(seq_along(years), function(year_id) {
      month_max <- if (year_id == length(years)) last_month else 12L
      data.frame(year_id = year_id, month = seq_len(month_max))
    })
  )

  pi_values_dec <- pi_values / 100
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

  g_tilde <- c(g_tilde_l, rep(exp((1 - gamma) * g), length(pi_values) - 2L), g_tilde_h)
  g_star <- c(g_star_l, rep(exp(-gamma * g), length(pi_values) - 2L), g_star_h)

  rows <- list()
  row_id <- 1L
  for (k in seq_len(dim(Qn_array)[3])) {
    Qn <- Qn_array[, , k]
    Qr <- nominal_to_real_q(Qn, pi_values)
    bt <- bts[diagnostics$month[k]]

    for (psi in psi_values) {
      sol <- solve_ez_fixed_point_from_q(Qr, g_tilde, g_star,
                                         beta = beta, gamma = gamma, psi = psi)
      P <- sol$P
      for (threshold in c(4, 5)) {
        q_spot <- tail_prob_5_10_threshold(Qn, bt, threshold, pi_bar = pi_values)
        qr_spot <- tail_prob_5_10_threshold(Qr, bt, threshold, pi_bar = pi_values)
        p_spot <- tail_prob_5_10_threshold(P, bt, threshold, pi_bar = pi_values)
        q_bins <- avg_inf_510_fast(Qn, bt, pi_bar = pi_values)
        report_pi_values_dec <- c(-2, seq(-0.5, 4.5, by = 1), 6) / 100
        w5 <- exp(5 * report_pi_values_dec)
        w10 <- exp(10 * report_pi_values_dec)
        qrh_5y_bins <- q_bins$g5 * w5 / sum(q_bins$g5 * w5)
        qrh_10y_bins <- q_bins$g10 * w10 / sum(q_bins$g10 * w10)
        idx_tail <- tail_bins(threshold)

        ym <- year_month_index[diagnostics$month[k], ]
        ## HRR's reported P series matches 0.662 times HRR real-Q tails with
        ## slightly different effective 5% conventions: US includes the 5%
        ## support point, while EZ excludes it. The 4% convention is common.
        support_idx <- if (threshold == 4 || area == "US") sprt_z >= threshold / 100 else sprt_z > 0.05

        rows[[row_id]] <- data.frame(
          area = area,
          month = diagnostics$month[k],
          date = diagnostics$date[k],
          psi = psi,
          threshold = threshold,
          q_5y = q_spot[["y5"]],
          q_10y = q_spot[["y10"]],
          q_5y5y = tail_prob_5y5y_threshold(Qn, bt, threshold, pi_bar = pi_values),
          qr_5y = qr_spot[["y5"]],
          qr_10y = qr_spot[["y10"]],
          qr_5y5y = tail_prob_5y5y_threshold(Qr, bt, threshold, pi_bar = pi_values),
          qrh_5y = sum(qrh_5y_bins[idx_tail]),
          qrh_10y = sum(qrh_10y_bins[idx_tail]),
          q_data_n_5y = sum(dists$data.ZC.N[support_idx, 5, ym$year_id, ym$month]),
          q_data_n_10y = sum(dists$data.ZC.N[support_idx, 10, ym$year_id, ym$month]),
          q_data_r_5y = sum(dists$data.ZC.Q[support_idx, 5, ym$year_id, ym$month]),
          q_data_r_10y = sum(dists$data.ZC.Q[support_idx, 10, ym$year_id, ym$month]),
          p_5y = p_spot[["y5"]],
          p_10y = p_spot[["y10"]],
          p_5y5y = tail_prob_5y5y_threshold(P, bt, threshold, pi_bar = pi_values),
          ez_converged = sol$converged
        )
        row_id <- row_id + 1L
      }
    }
  }

  out <- do.call(rbind, rows)
  write.csv(out, cache_file, row.names = FALSE)
  out
}

read_hrr_p <- function(area) {
  area <- normalize_hrr_area(area)
  read_dta(file.path("input", sprintf("%swestimates.dta", area)))
}

hrr_p_column <- function(horizon, threshold) {
  if (horizon == "5y5y") {
    paste0("higher", threshold, "_5y5y")
  } else {
    paste0("zc_higher", threshold, "_", horizon)
  }
}

plot_probability_threshold_figure <- function(all_series, threshold, output_file) {
  hrr_list <- list(US = read_hrr_p("US"), EZ = read_hrr_p("EZ"))
  areas <- c("US", "EZ")
  area_names <- c(US = "U.S.", EZ = "Euro area")
  horizons <- c("5y", "10y", "5y5y")
  horizon_labels <- c("5y", "10y", "5y5y")
  horizon_expr <- list(
    y5 = expression(pi[t*","*t+5]),
    y10 = expression(pi[t*","*t+10]),
    y5y5 = expression(pi[t+5*","*t+10])
  )
  psi_values <- sort(unique(all_series$psi))

  pdf(output_file, pointsize = 15, width = 11.2, height = 7.4)
  oldpar <- par(mfrow = c(2, 3), mar = c(3.8, 4.6, 2.8, 1.0),
                oma = c(3.2, 0, 0, 0),
                cex.axis = 1.15, cex.lab = 1.22, cex.main = 1.2)
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  for (area in areas) {
    series <- all_series[all_series$area == area &
                           all_series$threshold == threshold, ]
    series <- series[order(series$date, series$psi), ]
    hrr <- hrr_list[[area]]
    for (h in horizons) {
      p_col <- paste0("p_", h)
      hrr_col <- hrr_p_column(h, threshold)
      hrr_y <- as.numeric(hrr[[hrr_col]])

      y_for_range <- c(series[[p_col]], hrr_y)
      ymax <- 1.18 * max(100 * y_for_range, na.rm = TRUE)
      plot(series$date[series$psi == psi_values[1] & series$threshold == threshold],
           100 * series[[p_col]][series$psi == psi_values[1] & series$threshold == threshold],
           type = "n", ylim = c(0, ymax), xlab = "", ylab = "Probability (%)",
           las = 1,
           main = sprintf("%s, %s", area_names[[area]], horizon_labels[match(h, horizons)]))
      grid()

      for (i in seq_along(psi_values)) {
        idx <- series$psi == psi_values[i] & series$threshold == threshold
        lines(series$date[idx], 100 * series[[p_col]][idx],
              col = "black", lwd = 1 + 2 * (i - 1), lty = 1)
      }
      lines(hrr$date_stata, 100 * hrr_y, col = "gray35", lwd = 2, lty = 3)

    }
  }
  legend_text <- c(lapply(psi_values, function(x) bquote(psi == .(round(x, 2)))),
                   "HRR reported P")
  legend_col <- c(rep("black", length(psi_values)), "gray35")
  legend_lty <- c(rep(1, length(psi_values)), 3)
  legend_lwd <- c(1 + 2 * (seq_along(psi_values) - 1), 2)
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new()
  legend("bottom",
         legend = legend_text,
         col = legend_col,
         lty = legend_lty,
         lwd = legend_lwd,
         bg = "white", box.col = "gray80", ncol = 3, seg.len = 2,
         cex = 1.08, inset = 0.01, xpd = NA)
}

plot_ratio_figure <- function(all_series, output_file) {
  areas <- c("US", "EZ")
  area_names <- c(US = "U.S.", EZ = "Euro area")
  horizons <- c("5y", "10y")
  horizon_labels <- c("5y", "10y")
  psi_values <- sort(unique(all_series$psi))
  threshold <- 4

  pdf(output_file, pointsize = 15, width = 10.6, height = 5.0)
  oldpar <- par(mfrow = c(1, 2), mar = c(3.8, 5.4, 2.8, 1.0),
                oma = c(4.0, 0, 0, 0),
                cex.axis = 1.15, cex.lab = 1.22, cex.main = 1.2)
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  for (area in areas) {
    s_area <- all_series[all_series$area == area & all_series$threshold == threshold, ]
    s_area <- s_area[order(s_area$date, s_area$psi), ]
    ratios <- unlist(lapply(horizons, function(h) {
      s_area[[paste0("p_", h)]] / s_area[[paste0("qr_", h)]]
    }))
    ymax <- min(1.4, 1.15 * max(ratios, 0.66, na.rm = TRUE))

    plot(s_area$date[s_area$psi == psi_values[1]],
         rep(NA_real_, sum(s_area$psi == psi_values[1])),
         type = "n", ylim = c(0, ymax), xlab = "",
         ylab = expression(P / Q[horizon]^r),
         las = 1,
         main = bquote(.(area_names[[area]]) * ", " * pi > .(threshold) * "%"))
    grid()
    abline(h = 0.66, col = "gray65", lwd = 2.5, lty = 3)

    for (j in seq_along(horizons)) {
      for (i in seq_along(psi_values)) {
        idx <- s_area$psi == psi_values[i]
        ratio <- s_area[[paste0("p_", horizons[j])]][idx] /
          s_area[[paste0("qr_", horizons[j])]][idx]
        lines(s_area$date[idx], ratio,
              col = "black",
              lwd = 1 + 2 * (i - 1),
              lty = j)
      }
    }
  }

  legend_text <- c(
    horizon_labels,
    lapply(psi_values, function(x) bquote(psi == .(round(x, 2)))),
    "HRR 0.66"
  )
  legend_col <- c(rep("black", length(horizons)),
                  rep("black", length(psi_values)), "gray65")
  legend_lty <- c(seq_along(horizons), rep(1, length(psi_values)), 3)
  legend_lwd <- c(2, 2, 1 + 2 * (seq_along(psi_values) - 1), 2.5)
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new()
  legend("bottom",
         legend = legend_text,
         col = legend_col,
         lty = legend_lty,
         lwd = legend_lwd,
         bg = "white", box.col = "gray80", cex = 1.0,
         ncol = 3, inset = 0.01, xpd = NA)
}

compute_series_tasks <- function(tasks, psi_values, n_cores) {
  worker <- function(task) {
    message(sprintf("  Computing %s series with %s matrices",
                    task$area, toupper(task$matrix_source)))
    series <- compute_area_ies_series(task$area, psi_values,
                                      matrix_source = task$matrix_source)
    message(sprintf("  Finished %s series with %s matrices",
                    task$area, toupper(task$matrix_source)))
    series
  }

  if (.Platform$OS.type == "unix" && n_cores > 1L && length(tasks) > 1L) {
    parallel::mclapply(tasks, worker,
                       mc.cores = min(n_cores, length(tasks)),
                       mc.preschedule = FALSE)
  } else {
    lapply(tasks, worker)
  }
}

main <- function() {
  out_dir <- file.path(getwd(), "outputs", "comment_figures")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  psi_values <- c(1 / 3, 1.5)
  n_cores <- figure_n_cores()
  message(sprintf("Using %d core(s) for series computations.", n_cores))

  hrr_tasks <- lapply(c("EZ", "US"), function(area) {
    list(area = area, matrix_source = "hrr")
  })
  hrr_series_list <- compute_series_tasks(hrr_tasks, psi_values, n_cores)
  names(hrr_series_list) <- c("EZ", "US")
  ratio_series <- do.call(rbind, hrr_series_list)
  message("Drawing P/Q-ratio figure.")
  plot_ratio_figure(ratio_series, file.path(out_dir, "figure_ratios_HRR.pdf"))

  ngmmr_tasks <- lapply(c("EZ", "US"), function(area) {
    list(area = area, matrix_source = "ngmmr")
  })
  ngmmr_series_list <- compute_series_tasks(ngmmr_tasks, psi_values, n_cores)
  names(ngmmr_series_list) <- c("EZ", "US")

  all_series <- do.call(rbind, ngmmr_series_list)
  message("Drawing NGMMR P/Q-ratio figure.")
  plot_ratio_figure(all_series, file.path(out_dir, "figure_ratios_NGMMR.pdf"))
  message("Drawing physical-probability figure for the 4% threshold.")
  plot_probability_threshold_figure(
    all_series, 4,
    file.path(out_dir, "figure_IES_6plots_HRR_gt4.pdf")
  )

  message("Done.")
  message("Figure outputs have been written.")
}

if (sys.nframe() == 0) {
  main()
}
