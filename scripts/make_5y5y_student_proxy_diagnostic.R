## Diagnostic comparison of Gaussian and Student-t 5y5y tail proxies.

source(file.path("scripts", "estimate_nominal_q.R"))

student_5y5y_tail4 <- function(g5, g10,
                               rho,
                               df = 5,
                               bin_values_pct = c(-2, -0.5, 0.5, 1.5,
                                                  2.5, 3.5, 4.5, 6),
                               min_sd = 0.25) {
  mean5 <- sum(g5 * bin_values_pct)
  mean10 <- sum(g10 * bin_values_pct)
  var5 <- sum(g5 * (bin_values_pct - mean5)^2)
  var10 <- sum(g10 * (bin_values_pct - mean10)^2)

  mean_fwd <- 2 * mean10 - mean5
  disc <- 4 * var10 - var5 * (1 - rho^2)
  if (!is.finite(disc) || disc < 0 || df <= 2) {
    return(NA_real_)
  }

  sd_fwd <- -rho * sqrt(var5) + sqrt(disc)
  sd_fwd <- max(sd_fwd, min_sd)
  scale_fwd <- sd_fwd * sqrt((df - 2) / df)

  out <- 1 - pt((4 - mean_fwd) / scale_fwd, df = df)
  min(max(out, 0), 1)
}

make_area_proxy_series <- function(area, student_df = 5) {
  area <- normalize_hrr_area(area)
  dists <- load_hrr_dists(area)
  ym <- make_year_month_index(dists)
  rho <- resolve_fwd_proxy_rho("data", dists)
  tails_55 <- read_dta(file.path("input", sprintf("%s_55tails_monthly.dta", area)))

  rows <- lapply(seq_len(nrow(ym)), function(m) {
    g5 <- finegrid_zc_bins(dists, ym$year_id[m], ym$month[m], horizon = 5)
    g10 <- finegrid_zc_bins(dists, ym$year_id[m], ym$month[m], horizon = 10)
    tail_row <- tails_55[tails_55$year == ym$year[m] &
                           tails_55$month == ym$month[m], ]
    if (nrow(tail_row) != 1) {
      stop(sprintf("Could not find unique 5y5y tail row for %s month %03d.",
                   area, m))
    }

    data.frame(
      area = area,
      month = m,
      date = as.Date(sprintf("%04d-%02d-01", ym$year[m], ym$month[m])),
      rho = rho,
      gaussian = moment_proxy_5y5y_tail4(g5, g10, rho = rho),
      student = student_5y5y_tail4(g5, g10, rho = rho, df = student_df),
      hrr_proxy = tail_row$tail4_5y5y
    )
  })

  do.call(rbind, rows)
}

plot_student_proxy_diagnostic <- function(series, output_file, student_df = 5) {
  areas <- c("US", "EZ")
  ylim <- c(0, 1.08 * max(100 * c(series$gaussian,
                                  series$student,
                                  series$hrr_proxy),
                          na.rm = TRUE))

  pdf(output_file, pointsize = 15, width = 10.8, height = 5.8)
  oldpar <- par(mfrow = c(1, 2), mar = c(3.4, 4.6, 1.2, 0.8),
                oma = c(2.0, 0, 0, 0))
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  for (area in areas) {
    s <- series[series$area == area, ]
    plot(s$date, 100 * s$gaussian, type = "n", ylim = ylim,
         xlab = "", ylab = "Probability (%)", las = 1)
    grid()
    lines(s$date, 100 * s$gaussian, col = "black", lwd = 2.4, lty = 1)
    lines(s$date, 100 * s$student, col = "black", lwd = 2.4, lty = 2)
    lines(s$date, 100 * s$hrr_proxy, col = "gray55", lwd = 2.2, lty = 3)
    mtext(area, side = 3, line = 0.2, cex = 1.05)
  }

  par(fig = c(0, 1, 0, 0.20), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
      new = TRUE)
  plot.new()
  legend("center",
         legend = c("Gaussian proxy",
                    sprintf("Student-t proxy, df=%d", student_df),
                    "HRR 5y5y proxy"),
         col = c("black", "black", "gray55"),
         lty = c(1, 2, 3),
         lwd = c(2.4, 2.4, 2.2),
         horiz = TRUE,
         bg = "white",
         box.col = "gray80",
         cex = 0.9)

  invisible(output_file)
}

main <- function() {
  student_df <- as.integer(Sys.getenv("STUDENT_DF", unset = "5"))
  out_dir <- file.path(getwd(), "outputs", "hrr_nominal_real_q_diagnostics")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  series <- do.call(rbind, lapply(c("US", "EZ"), make_area_proxy_series,
                                  student_df = student_df))

  csv_file <- file.path(out_dir, sprintf("5y5y_gaussian_student_proxy_df%d.csv",
                                         student_df))
  pdf_file <- file.path(out_dir, sprintf("5y5y_gaussian_student_proxy_df%d.pdf",
                                         student_df))

  write.csv(series, csv_file, row.names = FALSE)
  plot_student_proxy_diagnostic(series, pdf_file, student_df = student_df)

  summary <- aggregate(cbind(gaussian, student, hrr_proxy) ~ area,
                       data = series, FUN = mean, na.rm = TRUE)
  message(sprintf("Wrote %s", csv_file))
  message(sprintf("Wrote %s", pdf_file))
  print(summary, digits = 4)
  invisible(series)
}

if (sys.nframe() == 0) {
  main()
}
