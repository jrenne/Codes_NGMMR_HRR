## Temporary diagnostic: refit the HRR transition-matrix form using the
## true model-implied 5y5y forward-average distribution.

source(file.path("scripts", "make_fit_diagnostics.R"))

out_dir <- file.path(getwd(), "tempo", "hrr_form_true_5y5y_refit")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

make_target_series <- function(area) {
  results <- readRDS(file.path(default_out_dir(area), "all_results.rds"))
  diagnostics <- read.csv(file.path(default_out_dir(area), "diagnostics.csv"))
  diagnostics$date <- as.Date(diagnostics$date)

  rows <- list()
  row_id <- 1L
  for (i in seq_along(results)) {
    for (threshold in c(4, 5)) {
      idx <- tail_idx(threshold)
      for (h in c("5y", "10y", "5y5y")) {
        slot <- switch(h, "5y" = "g5", "10y" = "g10", "5y5y" = "fwd")
        rows[[row_id]] <- data.frame(
          area = area,
          month = diagnostics$month[i],
          date = diagnostics$date[i],
          threshold = threshold,
          horizon = h,
          q_target = sum(results[[i]]$target[[slot]][idx]),
          q_ngmmr = sum(results[[i]]$ngmmr[[slot]][idx])
        )
        row_id <- row_id + 1L
      }
    }
  }
  do.call(rbind, rows)
}

message("Refitting HRR-form matrices with the true 5y5y forward-average object.")
hrr_form <- do.call(rbind, lapply(c("US", "EZ"), function(area) {
  make_hrr_form_refit_series(area, output_dir = out_dir, force = TRUE)
}))
targets <- do.call(rbind, lapply(c("US", "EZ"), make_target_series))

series <- merge(
  targets,
  hrr_form[c("area", "month", "date", "threshold", "horizon",
             "q_hrr_form_refit")],
  by = c("area", "month", "date", "threshold", "horizon"),
  all.x = TRUE
)
series <- series[order(series$area, series$threshold,
                       series$horizon, series$date), ]

rmse_pp <- function(x) 100 * sqrt(mean(x^2, na.rm = TRUE))
rmse_rows <- list()
row_id <- 1L
for (area in c("US", "EZ")) {
  for (threshold in c(4, 5)) {
    for (h in c("5y", "10y", "5y5y")) {
      s <- series[series$area == area &
                    series$threshold == threshold &
                    series$horizon == h, ]
      rmse_rows[[row_id]] <- data.frame(
        area = area,
        threshold = threshold,
        horizon = h,
        hrr_form_rmse_pp = rmse_pp(s$q_hrr_form_refit - s$q_target),
        ngmmr_rmse_pp = rmse_pp(s$q_ngmmr - s$q_target)
      )
      row_id <- row_id + 1L
    }
  }
}
rmse_table <- do.call(rbind, rmse_rows)

plot_fit <- function(series, threshold, output_file) {
  areas <- c("US", "EZ")
  area_names <- c(US = "U.S.", EZ = "Euro area")
  horizons <- c("5y", "10y", "5y5y")

  pdf(output_file, pointsize = 15, width = 12.5, height = 7.4)
  oldpar <- par(mfrow = c(2, 3), mar = c(3.8, 4.9, 2.8, 1.0),
                oma = c(3.1, 0, 0, 0),
                cex.axis = 1.12, cex.lab = 1.18, cex.main = 1.16)
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  for (area in areas) {
    for (h in horizons) {
      s <- series[series$area == area &
                    series$threshold == threshold &
                    series$horizon == h, ]
      s <- s[order(s$date), ]
      ymax <- 1.12 * max(100 * c(s$q_target, s$q_hrr_form_refit, s$q_ngmmr),
                         na.rm = TRUE)
      plot(s$date, 100 * s$q_target, type = "n",
           ylim = c(0, ymax), xlab = "", ylab = "Probability (%)",
           las = 1, main = sprintf("%s, %s", area_names[[area]], h))
      grid()
      lines(s$date, 100 * s$q_ngmmr, col = "gray60", lwd = 2.2, lty = 1)
      lines(s$date, 100 * s$q_hrr_form_refit,
            col = "black", lwd = 2.2, lty = 1)
      points(s$date, 100 * s$q_target, col = "black", pch = 4,
             cex = 0.62, lwd = 1.0)
    }
  }

  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  plot.new()
  legend("bottom",
         legend = c("Target", "HRR-form refit", "NGMMR baseline"),
         col = c("black", "black", "gray60"),
         lty = c(NA, 1, 1), lwd = c(NA, 2.2, 2.2),
         pch = c(4, NA, NA), pt.cex = c(0.75, NA, NA),
         bg = "white", box.col = "gray80",
         horiz = TRUE, cex = 0.95, xpd = NA, inset = -0.005)
}

write.csv(series, file.path(out_dir, "hrr_form_true_5y5y_refit_series.csv"),
          row.names = FALSE)
write.csv(rmse_table, file.path(out_dir, "hrr_form_true_5y5y_refit_rmse.csv"),
          row.names = FALSE)
plot_fit(series, 4, file.path(out_dir, "hrr_form_true_5y5y_refit_fit_gt4.pdf"))
plot_fit(series, 5, file.path(out_dir, "hrr_form_true_5y5y_refit_fit_gt5.pdf"))

message("Wrote diagnostics to:")
message(sprintf("  %s", out_dir))
print(rmse_table, digits = 4)
