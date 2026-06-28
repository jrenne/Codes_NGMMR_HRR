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

gaussian_forward_tail_proxy <- function(p5, p10, rho,
                                        support_pct,
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
  min(max(1 - pnorm(4, mean = mean_fwd, sd = sd_fwd), 0), 1)
}

default_estimation_out_dir <- function(area) {
  suffix <- paste0(
    "_nominal_Q_refit_smooth_row_poly2_finegrid_targets_",
    "hightailbin1p5_tailmoment2_moment5y5yW0p2_",
    "gauss5y5yRhodata_fwdBinW1"
  )
  file.path("outputs", paste0(area, suffix))
}

make_bounds_series <- function(area = "US") {
  area <- normalize_hrr_area(area)
  dists <- load_hrr_dists(area)
  ym <- make_year_month_index(dists)
  tails_55 <- read_dta(file.path("input", sprintf("%s_55tails_monthly.dta", area)))
  diag <- read.csv(file.path(default_estimation_out_dir(area), "diagnostics.csv"))

  support <- as.numeric(dists$sprt.Z[, 1])
  event <- outer(support, support, function(x, z) 2 * z - x > 0.04)

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
    diag_row <- diag[diag$month == m, ]

    data.frame(
      area = area,
      month = m,
      date = as.Date(sprintf("%04d-%02d-01", ym$year[m], ym$month[m])),
      lower = b[["lower"]],
      upper = b[["upper"]],
      gaussian_rho_min = gaussian_forward_tail_proxy(p5, p10, rho = -1,
                                                     support_pct = 100 * support),
      gaussian_rho_max = gaussian_forward_tail_proxy(p5, p10, rho = 1,
                                                     support_pct = 100 * support),
      hrr_proxy = tail_row$tail4_5y5y,
      ngmmr_5y5y = diag_row$qtail4_5y5y_ngmmr
    )
  })

  do.call(rbind, rows)
}

plot_bounds_series <- function(bounds, output_file) {
  pdf(output_file, pointsize = 15, width = 10, height = 5.8)
  oldpar <- par(mar = c(4, 4.4, 1.2, 1))
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  ymax <- max(100 * c(bounds$lower, bounds$upper,
                      bounds$gaussian_rho_min, bounds$gaussian_rho_max,
                      bounds$hrr_proxy, bounds$ngmmr_5y5y),
              na.rm = TRUE) * 1.08
  plot(bounds$date, 100 * bounds$upper, type = "n",
       ylim = c(0, ymax), las = 1, xlab = "", ylab = "Probability (%)")
  grid()
  polygon(c(bounds$date, rev(bounds$date)),
          100 * c(bounds$lower, rev(bounds$upper)),
          col = "gray88", border = NA)
  lines(bounds$date, 100 * bounds$gaussian_rho_min,
        col = "gray45", lwd = 1.7, lty = 3)
  lines(bounds$date, 100 * bounds$gaussian_rho_max,
        col = "gray45", lwd = 1.7, lty = 3)
  lines(bounds$date, 100 * bounds$ngmmr_5y5y,
        col = "black", lwd = 2.6, lty = 1)
  lines(bounds$date, 100 * bounds$hrr_proxy,
        col = "black", lwd = 2.6, lty = 2)
  legend("topright",
         legend = c("Frechet bounds",
                    "Gaussian proxy, rho=-1",
                    "Gaussian proxy, rho=1",
                    "NGMMR Q 5y5y",
                    "HRR 5y5y proxy"),
         col = c("gray88", "gray45", "gray45", "black", "black"),
         lty = c(1, 3, 3, 1, 2),
         lwd = c(8, 1.7, 1.7, 2.6, 2.6),
         bg = "white",
         box.col = "gray80")

  invisible(output_file)
}

main <- function() {
  out_dir <- file.path(getwd(), "outputs", "hrr_nominal_real_q_diagnostics")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  bounds <- make_bounds_series("US")
  bounds$below_lower <- bounds$hrr_proxy < bounds$lower - 1e-8
  bounds$above_upper <- bounds$hrr_proxy > bounds$upper + 1e-8

  csv_file <- file.path(out_dir, "5y5y_tail_bounds_from_5y_10y.csv")
  summary_file <- file.path(out_dir, "5y5y_tail_bounds_from_5y_10y_summary.csv")
  pdf_file <- file.path(out_dir, "5y5y_tail_bounds_from_5y_10y.pdf")

  summary <- data.frame(
    mean_lower = mean(bounds$lower),
    mean_upper = mean(bounds$upper),
    mean_gaussian_rho_min = mean(bounds$gaussian_rho_min, na.rm = TRUE),
    mean_gaussian_rho_max = mean(bounds$gaussian_rho_max, na.rm = TRUE),
    mean_hrr_proxy = mean(bounds$hrr_proxy),
    mean_ngmmr = mean(bounds$ngmmr_5y5y),
    share_hrr_below_lower = mean(bounds$below_lower),
    share_hrr_above_upper = mean(bounds$above_upper),
    min_hrr_minus_lower_pp = 100 * min(bounds$hrr_proxy - bounds$lower),
    max_hrr_minus_upper_pp = 100 * max(bounds$hrr_proxy - bounds$upper)
  )

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
