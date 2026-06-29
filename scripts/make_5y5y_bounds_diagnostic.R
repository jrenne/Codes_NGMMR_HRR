## Frechet bounds for the 5y5y high-inflation probability.

suppressPackageStartupMessages({
  library(haven)
})

source(file.path("scripts", "estimate_nominal_q.R"))

max_flow <- function(cap, source, sink, tol = 1e-12) {
  n <- nrow(cap)
  flow <- matrix(0, n, n)
  total <- 0

  repeat {
    parent <- rep(NA_integer_, n)
    parent[source] <- source
    path_cap <- rep(0, n)
    path_cap[source] <- Inf
    queue <- source

    while (length(queue) > 0 && is.na(parent[sink])) {
      u <- queue[1]
      queue <- queue[-1]
      residual <- cap[u, ] - flow[u, ]
      candidates <- which(is.na(parent) & residual > tol)
      for (v in candidates) {
        parent[v] <- u
        path_cap[v] <- min(path_cap[u], residual[v])
        queue <- c(queue, v)
        if (v == sink) {
          break
        }
      }
    }

    if (is.na(parent[sink])) {
      break
    }

    add <- path_cap[sink]
    v <- sink
    while (v != source) {
      u <- parent[v]
      flow[u, v] <- flow[u, v] + add
      flow[v, u] <- flow[v, u] - add
      v <- u
    }
    total <- total + add
  }

  total
}

max_event_probability <- function(p_x, p_z, event_matrix) {
  p_x <- pmax(as.numeric(p_x), 0)
  p_z <- pmax(as.numeric(p_z), 0)
  p_x <- p_x / sum(p_x)
  p_z <- p_z / sum(p_z)

  nx <- length(p_x)
  nz <- length(p_z)
  source <- 1L
  x_nodes <- 1L + seq_len(nx)
  z_nodes <- 1L + nx + seq_len(nz)
  sink <- 2L + nx + nz

  cap <- matrix(0, sink, sink)
  cap[source, x_nodes] <- p_x
  cap[z_nodes, sink] <- p_z

  for (i in seq_len(nx)) {
    for (j in seq_len(nz)) {
      if (event_matrix[i, j]) {
        cap[x_nodes[i], z_nodes[j]] <- 1
      }
    }
  }

  max_flow(cap, source, sink)
}

bound_event_probability <- function(p_x, p_z, event_matrix) {
  upper <- max_event_probability(p_x, p_z, event_matrix)
  lower <- 1 - max_event_probability(p_x, p_z, !event_matrix)
  c(lower = max(lower, 0), upper = min(upper, 1))
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
          idx_new <- offset + x_int[j]
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

default_estimation_out_dir <- function(area) {
  suffix <- paste0(
    "_nominal_Q_refit_smooth_row_poly3_finegrid_targets_",
    "hightailbin2_tailmoment2_",
    "spot_only_",
    "lrMean2p5W10_lrExtremeCap20W5_thetaShrinkW0p05"
  )
  file.path("outputs", paste0(area, suffix))
}

make_bounds_series <- function(area = "US", threshold = 4) {
  area <- normalize_hrr_area(area)
  threshold <- as.numeric(threshold)
  dists <- load_hrr_dists(area)
  ym <- make_year_month_index(dists)
  tails_55 <- read_dta(file.path("input", sprintf("%s_55tails_monthly.dta", area)))
  results <- readRDS(file.path(default_estimation_out_dir(area), "all_results.rds"))
  bts <- compute_bts(as.numeric(dists$infl))[seq_len(nrow(ym))]

  support <- as.numeric(dists$sprt.Z[, 1])
  event <- outer(support, support, function(x, z) 2 * z - x > threshold / 100)
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
      hrr_series = tail_row[[hrr_col]],
      ngmmr_5y5y = sum(forward_average_bins_model(
        results[[m]]$A, bts[m]
      )[tail_idx])
    )
  })

  do.call(rbind, rows)
}

plot_bounds_series <- function(bounds, output_file) {
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
    ymax <- max(100 * c(b$lower, b$upper, b$hrr_series, b$ngmmr_5y5y),
                na.rm = TRUE) * 1.08
    plot(b$date, 100 * b$upper, type = "n",
         ylim = c(0, ymax), las = 1, xlab = "", ylab = "Probability (%)",
         main = sprintf("%s, pi > %d%%", area_names[[area]], threshold))
    grid()
    polygon(c(b$date, rev(b$date)),
            100 * c(b$lower, rev(b$upper)),
            col = "gray88", border = NA)
    lines(b$date, 100 * b$ngmmr_5y5y,
          col = "black", lwd = 2.6, lty = 1)
    lines(b$date, 100 * b$hrr_series,
          col = "black", lwd = 2.6, lty = 2)
    }
  }

  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new()
  legend("bottom",
         legend = c("Bounds",
                    "NGMMR",
                    "HRR 5y5y"),
         col = c("gray88", "black", "black"),
         lty = c(1, 1, 2),
         lwd = c(8, 2.6, 2.6),
         bg = "white",
         box.col = "gray80",
         cex = 0.76,
         horiz = TRUE,
         xpd = NA,
         inset = -0.01)

  invisible(output_file)
}

main <- function() {
  out_dir <- file.path(getwd(), "outputs", "hrr_nominal_real_q_diagnostics")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  bounds <- do.call(rbind, lapply(c("US", "EZ"), function(area) {
    do.call(rbind, lapply(c(4, 5), function(threshold) {
      make_bounds_series(area, threshold)
    }))
  }))
  bounds$below_lower <- bounds$hrr_series < bounds$lower - 1e-8
  bounds$above_upper <- bounds$hrr_series > bounds$upper + 1e-8

  csv_file <- file.path(out_dir, "5y5y_tail_bounds_from_5y_10y.csv")
  summary_file <- file.path(out_dir, "5y5y_tail_bounds_from_5y_10y_summary.csv")
  pdf_file <- file.path(out_dir, "5y5y_tail_bounds_from_5y_10y.pdf")

  summary_rows <- lapply(c("US", "EZ"), function(area) {
    do.call(rbind, lapply(c(4, 5), function(threshold) {
      b <- bounds[bounds$area == area & bounds$threshold == threshold, ]
      data.frame(
        area = area,
        threshold = threshold,
        mean_lower = mean(b$lower),
        mean_upper = mean(b$upper),
        mean_hrr_series = mean(b$hrr_series),
        mean_ngmmr = mean(b$ngmmr_5y5y),
        share_hrr_below_lower = mean(b$below_lower),
        share_hrr_above_upper = mean(b$above_upper),
        min_hrr_minus_lower_pp = 100 * min(b$hrr_series - b$lower),
        max_hrr_minus_upper_pp = 100 * max(b$hrr_series - b$upper)
      )
    }))
  })
  summary <- do.call(rbind, summary_rows)

  write.csv(bounds, csv_file, row.names = FALSE)
  write.csv(summary, summary_file, row.names = FALSE)
  plot_bounds_series(bounds, pdf_file)

  message(sprintf("Wrote %s", csv_file))
  message(sprintf("Wrote %s", summary_file))
  message(sprintf("Wrote %s", pdf_file))
  print(summary, digits = 4)
  invisible(bounds)
}

if (sys.nframe() == 0) {
  main()
}
