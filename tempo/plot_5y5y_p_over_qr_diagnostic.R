## Temporary diagnostic: P/Q^r ratios for true 5y5y forward-average
## inflation above 4%, with the same smoothing convention as figure_ratios.pdf.

source(file.path("scripts", "make_ratio_and_probability_figures.R"))

out_dir <- file.path(getwd(), "tempo", "p_over_q_5y5y")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

psi_values <- c(1 / 3, 1.5)
series <- do.call(rbind, lapply(c("US", "EZ"), function(area) {
  compute_area_ies_series(area, psi_values = psi_values, matrix_source = "ngmmr")
}))

series <- series[series$threshold == 4, ]
series$ratio_p_over_qn_5y5y <- series$p_5y5y / series$q_5y5y
series$ratio_p_over_qr_5y5y <- series$p_5y5y / series$qr_5y5y

csv_file <- file.path(out_dir, "NGMMR_5y5y_P_over_Qr_gt4.csv")
pdf_file <- file.path(out_dir, "NGMMR_5y5y_P_over_Qr_gt4_smoothed.pdf")
write.csv(series, csv_file, row.names = FALSE)

pdf(pdf_file, pointsize = 15, width = 10.8, height = 4.8)
oldpar <- par(mfrow = c(1, 2), mar = c(3.8, 4.8, 2.2, 1.0),
              oma = c(1.4, 0, 0, 0),
              cex.axis = 1.12, cex.lab = 1.18, cex.main = 1.16)
on.exit({
  par(oldpar)
  dev.off()
}, add = TRUE)

area_names <- c(US = "U.S.", EZ = "Euro area")
for (area in c("US", "EZ")) {
  s_area <- series[series$area == area, ]
  s_area <- s_area[order(s_area$date, s_area$psi), ]
  ratios <- unlist(lapply(psi_values, function(psi) {
    smooth_na(s_area$ratio_p_over_qr_5y5y[s_area$psi == psi])
  }))
  ymax <- 1.12 * max(ratios, na.rm = TRUE)
  plot(s_area$date[s_area$psi == psi_values[1]],
       rep(NA_real_, sum(s_area$psi == psi_values[1])),
       type = "n", ylim = c(0, ymax), xlab = "",
       ylab = expression(P / Q^r),
       las = 1, main = area_names[[area]])
  grid()
  for (i in seq_along(psi_values)) {
    idx <- s_area$psi == psi_values[i]
    lines(s_area$date[idx], smooth_na(s_area$ratio_p_over_qr_5y5y[idx]),
          col = "black", lwd = 1 + 2 * (i - 1), lty = 1)
  }
}

par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
plot.new()
legend("bottom",
       legend = lapply(psi_values, function(x) bquote(psi == .(round(x, 2)))),
       col = "black", lty = 1, lwd = 1 + 2 * (seq_along(psi_values) - 1),
       bg = "white", box.col = "gray80", horiz = TRUE,
       cex = 0.95, inset = -0.01, xpd = NA)

message(sprintf("Wrote %s", csv_file))
message(sprintf("Wrote %s", pdf_file))
