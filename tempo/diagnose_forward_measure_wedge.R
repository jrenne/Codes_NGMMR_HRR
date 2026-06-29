## Temporary diagnostic: common nominal Q versus H-forward nominal Q.
##
## The current NGMMR solution gives a physical transition matrix P and a
## nominal risk-neutral transition matrix Qn.  This script reconstructs the
## nominal one-period Arrow-price matrix implied by the Epstein-Zin mapping,
## computes nominal bond prices, and compares annual-inflation marginals under
## the common nominal measure and under the corresponding H-forward measure.

rm(list = ls())

source(file.path("scripts", "pq_mapping_helpers.R"))

out_dir <- file.path("tempo", "forward_measure_wedge")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

result_suffix <- paste0(
  "_nominal_Q_refit_smooth_row_poly2_finegrid_targets_",
  "hightailbin1p5_tailmoment2_moment5y5yW0p2_",
  "gauss5y5yRhodata_fwdBinW1"
)

model_primitives <- function() {
  beta <- 0.99
  gamma <- 3
  psi <- 1.5
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
    psi = psi,
    g_tilde = c(g_tilde_l, rep(exp((1 - gamma) * g), 6), g_tilde_h),
    g_star = c(g_star_l, rep(exp(-gamma * g), 6), g_star_h)
  )
}

nominal_arrow_prices <- function(Qn, pi_values_pct, prim) {
  Qr <- nominal_to_real_q(Qn, pi_values_pct)
  sol <- solve_ez_fixed_point_from_q(Qr, prim$g_tilde, prim$g_star,
                                     beta = prim$beta,
                                     gamma = prim$gamma,
                                     psi = prim$psi)
  P <- sol$P
  u <- sol$u

  m <- as.numeric(P %*% (prim$g_tilde * u^(1 - prim$gamma)))
  row_factor <- prim$beta * m^((prim$gamma - 1 / prim$psi) / (1 - prim$gamma))
  col_factor <- prim$g_star * u^(1 / prim$psi - prim$gamma)
  Mr <- outer(row_factor, col_factor)
  Mn <- sweep(Mr, 2, exp(pi_values_pct / 100), "/")
  qn_arrow <- P * Mn

  check <- sweep(qn_arrow, 1, rowSums(qn_arrow), "/")
  list(
    qn_arrow = qn_arrow,
    P = P,
    Qr = Qr,
    max_qn_error = max(abs(check - Qn)),
    ez_converged = sol$converged
  )
}

bond_prices <- function(qn_arrow, max_h = 10L) {
  n_state <- nrow(qn_arrow)
  B <- matrix(NA_real_, n_state, max_h)
  B[, 1] <- rowSums(qn_arrow)
  if (max_h >= 2L) {
    for (h in 2:max_h) {
      B[, h] <- as.numeric(qn_arrow %*% B[, h - 1L])
    }
  }
  B
}

forward_marginal <- function(qn_arrow, B, bt, H) {
  dist <- numeric(nrow(qn_arrow))
  dist[bt] <- 1
  for (remaining in H:1) {
    if (remaining == 1L) {
      trans <- sweep(qn_arrow, 1, B[, 1], "/")
    } else {
      trans <- qn_arrow * matrix(B[, remaining - 1L],
                                 nrow(qn_arrow), ncol(qn_arrow),
                                 byrow = TRUE)
      trans <- sweep(trans, 1, B[, remaining], "/")
    }
    dist <- as.numeric(dist %*% trans)
  }
  dist / sum(dist)
}

common_marginal <- function(Qn, bt, H) {
  dist <- numeric(nrow(Qn))
  dist[bt] <- 1
  for (h in seq_len(H)) {
    dist <- as.numeric(dist %*% Qn)
  }
  dist
}

make_area_series <- function(area) {
  area <- normalize_hrr_area(area)
  out_est <- file.path("outputs", paste0(area, result_suffix))
  Qn_array <- readRDS(file.path(out_est, "ngmmr_Q_matrices.rds"))
  diagnostics <- read.csv(file.path(out_est, "diagnostics.csv"))
  diagnostics$date <- as.Date(diagnostics$date)

  dists <- load_hrr_dists(area)
  res <- load_hrr_monthquart(area)
  bts <- compute_bts(as.numeric(dists$infl))[seq_len(dim(res$gs.data)[3])]
  pi_values <- c(-2, seq(-0.5, 4.5, by = 1), 6)
  prim <- model_primitives()

  rows <- list()
  row_id <- 1L
  for (k in seq_len(dim(Qn_array)[3])) {
    Qn <- Qn_array[, , k]
    bt <- bts[diagnostics$month[k]]
    ap <- nominal_arrow_prices(Qn, pi_values, prim)
    B <- bond_prices(ap$qn_arrow, max_h = 10L)

    for (H in 6:10) {
      qc <- common_marginal(Qn, bt, H)
      qf <- forward_marginal(ap$qn_arrow, B, bt, H)
      for (threshold in c(4, 5)) {
        idx <- which(pi_values > threshold)
        rows[[row_id]] <- data.frame(
          area = area,
          month = diagnostics$month[k],
          date = diagnostics$date[k],
          H = H,
          threshold = threshold,
          q_common_tail = sum(qc[idx]),
          q_forward_tail = sum(qf[idx]),
          wedge = sum(qf[idx]) - sum(qc[idx]),
          max_qn_error = ap$max_qn_error,
          ez_converged = ap$ez_converged
        )
        row_id <- row_id + 1L
      }
    }
  }

  do.call(rbind, rows)
}

series <- do.call(rbind, lapply(c("US", "EZ"), make_area_series))
write.csv(series, file.path(out_dir, "forward_measure_wedge_annual.csv"),
          row.names = FALSE)

avg_series <- aggregate(
  cbind(q_common_tail, q_forward_tail, wedge) ~ area + month + date + threshold,
  data = series,
  FUN = mean
)
write.csv(avg_series, file.path(out_dir, "forward_measure_wedge_avg6to10.csv"),
          row.names = FALSE)

summary_tab <- aggregate(
  abs(wedge) ~ area + threshold,
  data = avg_series,
  FUN = function(x) c(mean = mean(x), p90 = quantile(x, 0.9), max = max(x))
)
summary_out <- do.call(data.frame, summary_tab)
names(summary_out) <- c("area", "threshold",
                        "mean_abs_wedge", "p90_abs_wedge", "max_abs_wedge")
write.csv(summary_out, file.path(out_dir, "forward_measure_wedge_summary.csv"),
          row.names = FALSE)

pdf(file.path(out_dir, "forward_measure_wedge_avg6to10.pdf"),
    width = 10.5, height = 6.8)
oldpar <- par(mfrow = c(2, 2), mar = c(3.2, 4.4, 2.0, 0.8),
              oma = c(0, 0, 0, 0))
on.exit({
  par(oldpar)
  dev.off()
}, add = TRUE)

area_names <- c(US = "U.S.", EZ = "Euro area")
for (area in c("US", "EZ")) {
  for (threshold in c(4, 5)) {
    sub <- avg_series[avg_series$area == area & avg_series$threshold == threshold, ]
    ymax <- max(100 * c(sub$q_common_tail, sub$q_forward_tail), na.rm = TRUE) * 1.08
    plot(sub$date, 100 * sub$q_common_tail, type = "l", lwd = 2,
         col = "black", ylim = c(0, ymax), las = 1,
         xlab = "", ylab = "Probability (%)",
         main = sprintf("%s, annual forward avg., >%d%%", area_names[[area]], threshold))
    lines(sub$date, 100 * sub$q_forward_tail, lwd = 2, col = "gray60", lty = 2)
    grid()
    legend("topleft",
           legend = c(expression(Q^n), expression(Q^{n * "," * H})),
           col = c("black", "gray60"), lty = c(1, 2), lwd = 2,
           bty = "o", bg = "white", cex = 0.95)
  }
}

message("Wrote diagnostics to ", out_dir)
print(summary_out)
