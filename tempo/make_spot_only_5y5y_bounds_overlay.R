## Temporary diagnostic: reproduce the 5y5y bounds figure for pi>4% and
## pi>5%, and add the spot-only 5y/10y refit as a red line.

source(file.path("scripts", "make_5y5y_bounds_diagnostic.R"))

tempo_dir <- file.path(getwd(), "tempo", "spot_only_5y10y")

default_main_estimation_out_dir <- function(area) {
  suffix <- paste0(
    "_nominal_Q_refit_smooth_row_poly2_finegrid_targets_",
    "hightailbin1p5_tailmoment2_moment5y5yW0p2_",
    "gauss5y5yRhodata_fwdBinW1"
  )
  file.path("outputs", paste0(area, suffix))
}

tail_idx_for_threshold <- function(threshold) {
  if (threshold == 4) {
    7:8
  } else if (threshold == 5) {
    8
  } else {
    stop("Only thresholds 4 and 5 are implemented.")
  }
}

forward_average_bins_model <- function(A, bt,
                                       pi_bar = c(-2, seq(-0.5, 4.5, by = 1), 6),
                                       pi_lim = c(seq(-1, 5, by = 1), Inf)) {
  n_state <- length(pi_bar)
  step <- 2L
  x_int <- as.integer(round(step * pi_bar))
  min_sum <- 5L * min(x_int)
  max_sum <- 5L * max(x_int)
  offset <- 1L - min_sum
  n_sum <- max_sum - min_sum + 1L

  dist_state <- numeric(n_state)
  dist_state[bt] <- 1
  dist_fwd <- matrix(0, n_state, n_sum)

  for (h in seq_len(10)) {
    if (h <= 5) {
      dist_state <- as.numeric(dist_state %*% A)
    } else {
      next_fwd <- matrix(0, n_state, n_sum)
      for (j in seq_len(n_state)) {
        if (h == 6) {
          mass_by_sum <- sum(dist_state * A[, j])
          idx_old <- offset
          idx_new <- idx_old + x_int[j]
          next_fwd[j, idx_new] <- next_fwd[j, idx_new] + mass_by_sum
        } else {
          mass_by_sum <- as.numeric(crossprod(A[, j], dist_fwd))
          idx_old <- which(mass_by_sum != 0)
          idx_new <- idx_old + x_int[j]
          keep <- idx_new >= 1L & idx_new <= n_sum
          next_fwd[j, idx_new[keep]] <- next_fwd[j, idx_new[keep]] +
            mass_by_sum[idx_old[keep]]
        }
      }
      dist_fwd <- next_fwd
    }
  }

  cum_prob <- colSums(dist_fwd)
  used <- which(cum_prob != 0)
  cum_sum <- used - offset
  avg <- (cum_sum / step) / 5

  bin <- length(pi_bar) + 1L -
    rowSums(outer(avg, pi_lim, function(a, lim) lim - a > -1e-7))
  out <- numeric(n_state)
  for (b in seq_len(n_state)) {
    out[b] <- sum(cum_prob[used[bin == b]])
  }
  out / sum(out)
}

gaussian_forward_tail_proxy_threshold <- function(p5, p10, rho,
                                                  support_pct,
                                                  threshold,
                                                  min_sd = 0.25) {
  mean5 <- sum(p5 * support_pct)
  mean10 <- sum(p10 * support_pct)
  var5 <- sum(p5 * (support_pct - mean5)^2)
  var10 <- sum(p10 * (support_pct - mean10)^2)

  mean_fwd <- 2 * mean10 - mean5
  disc <- 4 * var10 - var5 * (1 - rho^2)
  if (!is.finite(disc) || disc < 0) {
    return(NA_real_)
  }

  sd_fwd <- -rho * sqrt(var5) + sqrt(disc)
  sd_fwd <- max(sd_fwd, min_sd)
  min(max(1 - pnorm(threshold, mean = mean_fwd, sd = sd_fwd), 0), 1)
}

make_bounds_series_with_spot_only <- function(area = "US", threshold = 4) {
  area <- normalize_hrr_area(area)
  threshold <- as.numeric(threshold)
  dists <- load_hrr_dists(area)
  ym <- make_year_month_index(dists)
  tails_55 <- read_dta(file.path("input", sprintf("%s_55tails_monthly.dta", area)))
  main_results <- readRDS(file.path(default_main_estimation_out_dir(area),
                                    "all_results.rds"))
  spot_results <- readRDS(file.path(tempo_dir, area, "all_results.rds"))
  infl <- as.numeric(dists$infl)
  bts <- compute_bts(infl)[seq_len(nrow(ym))]

  spot_file <- file.path(tempo_dir, area, "diagnostics.csv")
  if (!file.exists(spot_file)) {
    stop(sprintf("Missing spot-only diagnostics: %s", spot_file))
  }
  spot_diag <- read.csv(spot_file)
  spot_diag$date <- as.Date(spot_diag$date)

  support <- as.numeric(dists$sprt.Z[, 1])
  event <- outer(support, support, function(x, z) {
    2 * z - x > threshold / 100
  })
  tail_idx <- tail_idx_for_threshold(threshold)
  hrr_col <- paste0("tail", threshold, "_5y5y")

  rows <- lapply(seq_len(nrow(ym)), function(m) {
    p5 <- as.numeric(dists$data.ZC.N[, 5, ym$year_id[m], ym$month[m]])
    p10 <- as.numeric(dists$data.ZC.N[, 10, ym$year_id[m], ym$month[m]])
    b <- bound_event_probability(p5, p10, event)
    tail_row <- tails_55[tails_55$year == ym$year[m] &
                           tails_55$month == ym$month[m], ]
    if (nrow(tail_row) != 1) {
      stop(sprintf("Could not find unique HRR 5y5y row for %s month %03d.",
                   area, m))
    }

    data.frame(
      area = area,
      month = m,
      date = as.Date(sprintf("%04d-%02d-01", ym$year[m], ym$month[m])),
      threshold = threshold,
      lower = b[["lower"]],
      upper = b[["upper"]],
      gaussian_rho_min = gaussian_forward_tail_proxy_threshold(
        p5, p10, rho = -0.7, support_pct = 100 * support, threshold = threshold
      ),
      gaussian_rho_max = gaussian_forward_tail_proxy_threshold(
        p5, p10, rho = 0.7, support_pct = 100 * support, threshold = threshold
      ),
      hrr_proxy = tail_row[[hrr_col]],
      ngmmr_5y5y = sum(forward_average_bins_model(
        main_results[[m]]$A, bts[m]
      )[tail_idx]),
      ngmmr_spot_only_5y5y = sum(forward_average_bins_model(
        spot_results[[m]]$A, bts[m]
      )[tail_idx])
    )
  })

  do.call(rbind, rows)
}

plot_bounds_series_with_spot_only <- function(bounds, output_file) {
  areas <- c("US", "EZ")
  area_names <- c(US = "U.S.", EZ = "Euro area")

  pdf(output_file, pointsize = 15, width = 12.2, height = 8.8)
  oldpar <- par(mfrow = c(2, 2), mar = c(3.5, 4.4, 2.2, 1),
                oma = c(1.45, 0, 0, 0))
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  for (area in areas) {
    for (threshold in c(4, 5)) {
    b <- bounds[bounds$area == area & bounds$threshold == threshold, ]
    ymax <- max(100 * c(b$lower, b$upper,
                        b$gaussian_rho_min, b$gaussian_rho_max,
                        b$hrr_proxy, b$ngmmr_5y5y,
                        b$ngmmr_spot_only_5y5y),
                na.rm = TRUE) * 1.08
    plot(b$date, 100 * b$upper, type = "n",
         ylim = c(0, ymax), las = 1, xlab = "", ylab = "Probability (%)",
         main = sprintf("%s, pi > %d%%", area_names[[area]], threshold))
    grid()
    polygon(c(b$date, rev(b$date)),
            100 * c(b$lower, rev(b$upper)),
            col = "gray88", border = NA)
    lines(b$date, 100 * b$gaussian_rho_min,
          col = "gray45", lwd = 1.7, lty = 3)
    lines(b$date, 100 * b$gaussian_rho_max,
          col = "gray45", lwd = 1.7, lty = 3)
    lines(b$date, 100 * b$ngmmr_spot_only_5y5y,
          col = "gray62", lwd = 2.2, lty = 1)
    lines(b$date, 100 * b$ngmmr_5y5y,
          col = "black", lwd = 2.6, lty = 1)
    lines(b$date, 100 * b$hrr_proxy,
          col = "black", lwd = 2.6, lty = 2)
    }
  }

  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new()
  legend("bottom",
         legend = c("Bounds",
                    "Gaussian",
                    "NGMMR",
                    "HRR proxy",
                    "Spot-only"),
         col = c("gray88", "gray45", "black", "black", "gray62"),
         lty = c(1, 3, 1, 2, 1),
         lwd = c(8, 1.7, 2.6, 2.6, 2.2),
         bg = "white",
         box.col = "gray80",
         cex = 0.76,
         horiz = TRUE,
         xpd = NA,
         inset = -0.01)

  invisible(output_file)
}

bounds <- do.call(rbind, lapply(c("US", "EZ"), function(area) {
  do.call(rbind, lapply(c(4, 5), function(threshold) {
    make_bounds_series_with_spot_only(area, threshold)
  }))
}))
out_file <- file.path(tempo_dir, "5y5y_tail_bounds_from_5y_10y_spot_only_overlay.pdf")
csv_file <- file.path(tempo_dir, "5y5y_tail_bounds_from_5y_10y_spot_only_overlay.csv")

write.csv(bounds, csv_file, row.names = FALSE)
plot_bounds_series_with_spot_only(bounds, out_file)

message(sprintf("Wrote %s", csv_file))
message(sprintf("Wrote %s", out_file))
