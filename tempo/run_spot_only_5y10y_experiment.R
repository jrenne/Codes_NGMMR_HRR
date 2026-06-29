## Temporary diagnostic: fit the 6-parameter nominal-Q model using only
## the 5y and 10y distributions, then inspect the implied 5y5y tail.

REESTIMATE_SPOT_ONLY <- TRUE
AREAS <- c("US", "EZ")
N_CORES <- 8

spot_only_suffix <- paste0(
  "_nominal_Q_refit_smooth_row_poly2_finegrid_targets_",
  "hightailbin1p5_tailmoment2_spot_only"
)

tempo_dir <- file.path(getwd(), "tempo", "spot_only_5y10y")
dir.create(tempo_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path("scripts", "estimate_nominal_q.R"))

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

copy_dir_contents <- function(from, to) {
  dir.create(to, showWarnings = FALSE, recursive = TRUE)
  files <- list.files(from, full.names = TRUE)
  file.copy(files, to, overwrite = TRUE, recursive = TRUE)
}

run_spot_only_estimation <- function(area) {
  out_dir <- file.path("outputs", paste0(area, spot_only_suffix))
  done <- file.exists(file.path(out_dir, "diagnostics.csv")) &&
    file.exists(file.path(out_dir, "all_results.rds")) &&
    file.exists(file.path(out_dir, "ngmmr_Q_matrices.rds"))

  if (!done || REESTIMATE_SPOT_ONLY) {
    message(sprintf("\nEstimating spot-only model for %s.", area))
    old <- set_env(
      AREA = area,
      N_CORES = as.character(N_CORES),
      MODEL_VARIANT = "smooth_row_poly2",
      TAIL_WEIGHT_SIDE = "high",
      TAIL_BIN_WEIGHT = "1.5",
      TAIL_MOMENT_WEIGHT = "2",
      INCLUDE_5Y5Y_TARGET = "0",
      FWD_TARGET_SOURCE = "gaussian",
      FWD_PROXY_RHO = "data",
      FWD_PROXY_WEIGHT = "0",
      FWD_BIN_WEIGHT = "1",
      FWD_BOUNDS_WEIGHT = "0"
    )
    on.exit(restore_env(old), add = TRUE)
    main()
  } else {
    message(sprintf("\nReusing existing spot-only estimates for %s.", area))
  }

  copy_dir_contents(out_dir, file.path(tempo_dir, area))
  out_dir
}

make_area_series <- function(area, out_dir) {
  dists <- load_hrr_dists(area)
  ym <- make_year_month_index(dists)
  tails_55 <- read_dta(file.path("input", sprintf("%s_55tails_monthly.dta", area)))
  diag <- read.csv(file.path(out_dir, "diagnostics.csv"))

  rows <- lapply(seq_len(nrow(diag)), function(i) {
    m <- diag$month[i]
    tail_row <- tails_55[tails_55$year == ym$year[m] &
                           tails_55$month == ym$month[m], ]
    data.frame(
      area = area,
      month = m,
      date = as.Date(diag$date[i]),
      q5_target = diag$qtail4_5y_target[i],
      q5_model = diag$qtail4_5y_ngmmr[i],
      q10_target = diag$qtail4_10y_target[i],
      q10_model = diag$qtail4_10y_ngmmr[i],
      q5y5y_ngmmr = diag$qtail4_5y5y_ngmmr[i],
      q5y5y_gaussian_proxy = diag$qtail4_5y5y_moment_proxy[i],
      q5y5y_hrr_proxy = tail_row$tail4_5y5y
    )
  })
  do.call(rbind, rows)
}

rmse_pp <- function(x) 100 * sqrt(mean(x^2, na.rm = TRUE))

plot_spot_only_diagnostic <- function(series, output_file) {
  pdf(output_file, pointsize = 15, width = 10.8, height = 7.4)
  oldpar <- par(mfrow = c(2, 2), mar = c(3.7, 4.8, 2.5, 1.0),
                oma = c(3.0, 0, 0, 0),
                cex.axis = 1.12, cex.lab = 1.18, cex.main = 1.16)
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  area_names <- c(US = "U.S.", EZ = "Euro area")
  for (area in AREAS) {
    s <- series[series$area == area, ]
    s <- s[order(s$date), ]

    ymax_spot <- 1.12 * max(100 * c(s$q5_target, s$q5_model,
                                     s$q10_target, s$q10_model), na.rm = TRUE)
    plot(s$date, 100 * s$q5_model, type = "l",
         col = "black", lwd = 2.2, ylim = c(0, ymax_spot),
         las = 1, xlab = "", ylab = "Probability (%)",
         main = paste(area_names[[area]], "5y and 10y fit"))
    grid()
    points(s$date, 100 * s$q5_target, pch = 4, col = "gray20",
           cex = 0.65, lwd = 1.1)
    lines(s$date, 100 * s$q10_model, col = "gray45", lwd = 2.2, lty = 2)
    points(s$date, 100 * s$q10_target, pch = 1, col = "gray45",
           cex = 0.65, lwd = 1.1)

    ymax_fwd <- 1.12 * max(100 * c(s$q5y5y_ngmmr, s$q5y5y_gaussian_proxy,
                                    s$q5y5y_hrr_proxy), na.rm = TRUE)
    plot(s$date, 100 * s$q5y5y_ngmmr, type = "l",
         col = "black", lwd = 2.4, ylim = c(0, ymax_fwd),
         las = 1, xlab = "", ylab = "Probability (%)",
         main = paste(area_names[[area]], "implied 5y5y"))
    grid()
    lines(s$date, 100 * s$q5y5y_gaussian_proxy,
          col = "gray45", lwd = 2.0, lty = 3)
    lines(s$date, 100 * s$q5y5y_hrr_proxy,
          col = "gray45", lwd = 2.0, lty = 2)
  }

  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new()
  legend("bottom",
         legend = c("5y model / 5y5y NGMMR",
                    "5y target",
                    "10y model",
                    "10y target",
                    "Gaussian 5y5y proxy",
                    "HRR 5y5y proxy"),
         col = c("black", "gray20", "gray45", "gray45", "gray45", "gray45"),
         lty = c(1, NA, 2, NA, 3, 2),
         lwd = c(2.4, NA, 2.2, NA, 2.0, 2.0),
         pch = c(NA, 4, NA, 1, NA, NA),
         bg = "white", box.col = "gray80",
         horiz = TRUE, cex = 0.82, xpd = NA, inset = 0.005)
}

out_dirs <- setNames(lapply(AREAS, run_spot_only_estimation), AREAS)
series <- do.call(rbind, lapply(AREAS, function(area) {
  make_area_series(area, out_dirs[[area]])
}))

summary <- do.call(rbind, lapply(AREAS, function(area) {
  s <- series[series$area == area, ]
  data.frame(
    area = area,
    rmse_5y_gt4_pp = rmse_pp(s$q5_model - s$q5_target),
    rmse_10y_gt4_pp = rmse_pp(s$q10_model - s$q10_target),
    rmse_5y5y_vs_gaussian_pp = rmse_pp(s$q5y5y_ngmmr - s$q5y5y_gaussian_proxy),
    rmse_5y5y_vs_hrr_proxy_pp = rmse_pp(s$q5y5y_ngmmr - s$q5y5y_hrr_proxy),
    mean_5y5y_ngmmr_pct = 100 * mean(s$q5y5y_ngmmr, na.rm = TRUE),
    mean_5y5y_gaussian_pct = 100 * mean(s$q5y5y_gaussian_proxy, na.rm = TRUE),
    mean_5y5y_hrr_proxy_pct = 100 * mean(s$q5y5y_hrr_proxy, na.rm = TRUE)
  )
}))

write.csv(series, file.path(tempo_dir, "spot_only_5y10y_series.csv"), row.names = FALSE)
write.csv(summary, file.path(tempo_dir, "spot_only_5y10y_summary.csv"), row.names = FALSE)
plot_spot_only_diagnostic(series,
                          file.path(tempo_dir, "spot_only_5y10y_diagnostic.pdf"))

message("\nSpot-only diagnostic written to:")
message(sprintf("  %s", tempo_dir))
print(summary, digits = 4)
