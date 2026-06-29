## Plot NGMMR P and Q probabilities of average inflation in both tails.
##
## Inputs:
## - nominal risk-neutral Q matrices estimated by the refit scripts.
## - HRR reported physical probabilities from AREAwestimates.dta.
##
## Convention:
## - Q series are computed from the fitted nominal risk-neutral matrices.
## - For 5y and 10y, the data-Q series is the nominal fine-grid RN object,
##   because this is the object targeted by the headline-probability refit.
## - P series are computed by first converting nominal Q to real Q and then
##   applying the Epstein-Zin Q-to-P mapping used in the comment.

suppressPackageStartupMessages({
  library(haven)
})

source(file.path("scripts", "qn_parameterization.R"))

q_to_p <- function(Q, w, tol = 1e-12) {
  Q <- as.matrix(Q)
  w <- as.numeric(w)
  Q_div_w <- sweep(Q, 2, w, "/")
  s <- rowSums(Q_div_w)
  if (any(s <= tol)) {
    stop("Some rows have nearly zero normalization in Q-to-P mapping.")
  }
  Q_div_w / s
}

solve_ez_fixed_point_from_q <- function(Q, g_tilde, g_star,
                                        beta, gamma, psi,
                                        tol = 1e-10,
                                        max_iter = 100000) {
  J <- nrow(Q)
  u <- rep(1, J)
  a <- (1 - 1 / psi) / (1 - gamma)
  b <- 1 / (1 - 1 / psi)

  for (it in seq_len(max_iter)) {
    w <- as.numeric(g_star) * u^(1 / psi - gamma)
    P <- q_to_p(Q, w)
    m <- as.numeric(P %*% (as.numeric(g_tilde) * u^(1 - gamma)))
    if (any(m <= 0) || any(!is.finite(m))) {
      stop("Non-positive/non-finite fixed-point intermediate value.")
    }
    u_new <- ((1 - beta) + beta * (m^a))^b
    err <- max(abs(u_new - u))
    if (err < tol) {
      w_new <- as.numeric(g_star) * u_new^(1 / psi - gamma)
      return(list(u = u_new, w = w_new, P = q_to_p(Q, w_new),
                  converged = TRUE, iter = it, err = err))
    }
    u <- u_new
  }

  w <- as.numeric(g_star) * u^(1 / psi - gamma)
  list(u = u, w = w, P = q_to_p(Q, w),
       converged = FALSE, iter = max_iter, err = err)
}

nominal_to_real_q <- function(Qn, pi_values_pct) {
  v <- exp(pi_values_pct / 100)
  diag(as.vector(1 / (Qn %*% v))) %*% Qn %*% diag(as.vector(v))
}

tail_prob_5_10 <- function(A, bt, side = c("high", "low")) {
  side <- match.arg(side)
  bins <- avg_inf_510_fast(A, bt)
  idx <- if (side == "high") 7:8 else 1:2
  c(
    p5 = sum(bins$g5[idx]),
    p10 = sum(bins$g10[idx])
  )
}

tail_prob_5y5y_proxy <- function(A, bt, side = c("high", "low")) {
  side <- match.arg(side)
  fwd_bins <- forward_average_bins_model(A, bt)
  idx <- if (side == "high") 7:8 else 1:2
  sum(fwd_bins[idx])
}

make_tail_probability_series <- function(psi = 1.5, out_dir,
                                         area = Sys.getenv("AREA", unset = "US")) {
  area <- normalize_hrr_area(area)
  q_file <- file.path(out_dir, "ngmmr_Q_matrices.rds")
  diagnostics_file <- file.path(out_dir, "diagnostics.csv")

  Qn_array <- readRDS(q_file)
  diagnostics <- read.csv(diagnostics_file)
  diagnostics$date <- as.Date(diagnostics$date)

  res <- load_hrr_monthquart(area)
  dists <- load_hrr_dists(area)
  infl <- as.numeric(dists$infl)
  n_month <- dim(res$gs.data)[3]
  bts <- compute_bts(infl)[seq_len(n_month)]

  pi_values <- c(-2, seq(-0.5, 4.5, by = 1), 6)
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
  beta <- 0.99
  gamma <- 3
  psi <- psi
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

  rows <- lapply(seq_len(dim(Qn_array)[3]), function(k) {
    month_id <- diagnostics$month[k]
    Qn <- Qn_array[, , k]
    Qr <- nominal_to_real_q(Qn, pi_values)
    sol <- solve_ez_fixed_point_from_q(Qr, g_tilde, g_star,
                                       beta = beta, gamma = gamma, psi = psi)
    P <- sol$P
    bt <- bts[month_id]

    q_tail_high <- tail_prob_5_10(Qn, bt, side = "high")
    q_tail_low <- tail_prob_5_10(Qn, bt, side = "low")
    qr_tail_high <- tail_prob_5_10(Qr, bt, side = "high")
    qr_tail_low <- tail_prob_5_10(Qr, bt, side = "low")
    p_tail_high <- tail_prob_5_10(P, bt, side = "high")
    p_tail_low <- tail_prob_5_10(P, bt, side = "low")
    ym <- year_month_index[month_id, ]

    data.frame(
      month = month_id,
      date = diagnostics$date[k],
      label = diagnostics$label[k],
      q_5y_gt4 = q_tail_high[["p5"]],
      q_10y_gt4 = q_tail_high[["p10"]],
      q_5y5y_gt4 = tail_prob_5y5y_proxy(Qn, bt, side = "high"),
      q_data_5y_gt4 = diagnostics$qtail4_5y_target[k],
      q_data_10y_gt4 = diagnostics$qtail4_10y_target[k],
      q_data_5y5y_gt4 = diagnostics$qtail4_5y5y_target[k],
      q_real_5y_gt4 = qr_tail_high[["p5"]],
      q_real_10y_gt4 = qr_tail_high[["p10"]],
      q_real_5y5y_gt4 = tail_prob_5y5y_proxy(Qr, bt, side = "high"),
      q_data_real_5y_gt4 = sum(dists$data.ZC.Q[sprt_z > 0.04, 5, ym$year_id, ym$month]),
      q_data_real_10y_gt4 = sum(dists$data.ZC.Q[sprt_z > 0.04, 10, ym$year_id, ym$month]),
      q_data_real_5y5y_gt4 = diagnostics$qtail4_5y5y_target[k],
      p_5y_gt4 = p_tail_high[["p5"]],
      p_10y_gt4 = p_tail_high[["p10"]],
      p_5y5y_gt4 = tail_prob_5y5y_proxy(P, bt, side = "high"),
      q_5y_lt0 = q_tail_low[["p5"]],
      q_10y_lt0 = q_tail_low[["p10"]],
      q_5y5y_lt0 = tail_prob_5y5y_proxy(Qn, bt, side = "low"),
      q_data_5y_lt0 = diagnostics$qtail0_5y_target[k],
      q_data_10y_lt0 = diagnostics$qtail0_10y_target[k],
      q_data_5y5y_lt0 = diagnostics$qtail0_5y5y_target[k],
      q_real_5y_lt0 = qr_tail_low[["p5"]],
      q_real_10y_lt0 = qr_tail_low[["p10"]],
      q_real_5y5y_lt0 = tail_prob_5y5y_proxy(Qr, bt, side = "low"),
      q_data_real_5y_lt0 = sum(dists$data.ZC.Q[sprt_z < 0, 5, ym$year_id, ym$month]),
      q_data_real_10y_lt0 = sum(dists$data.ZC.Q[sprt_z < 0, 10, ym$year_id, ym$month]),
      q_data_real_5y5y_lt0 = diagnostics$qtail0_5y5y_target[k],
      p_5y_lt0 = p_tail_low[["p5"]],
      p_10y_lt0 = p_tail_low[["p10"]],
      p_5y5y_lt0 = tail_prob_5y5y_proxy(P, bt, side = "low"),
      psi = psi,
      ez_converged = sol$converged,
      ez_iter = sol$iter,
      ez_err = sol$err
    )
  })

  do.call(rbind, rows)
}

use_real_q_series <- function(series) {
  out <- series
  out$q_5y_gt4 <- series$q_real_5y_gt4
  out$q_10y_gt4 <- series$q_real_10y_gt4
  out$q_5y5y_gt4 <- series$q_real_5y5y_gt4
  out$q_data_5y_gt4 <- series$q_data_real_5y_gt4
  out$q_data_10y_gt4 <- series$q_data_real_10y_gt4
  out$q_data_5y5y_gt4 <- series$q_data_real_5y5y_gt4
  out$q_5y_lt0 <- series$q_real_5y_lt0
  out$q_10y_lt0 <- series$q_real_10y_lt0
  out$q_5y5y_lt0 <- series$q_real_5y5y_lt0
  out$q_data_5y_lt0 <- series$q_data_real_5y_lt0
  out$q_data_10y_lt0 <- series$q_data_real_10y_lt0
  out$q_data_5y5y_lt0 <- series$q_data_real_5y5y_lt0
  out
}

plot_tail_probability_series <- function(series, hrr, output_file,
                                         side = c("high", "low")) {
  side <- match.arg(side)
  pdf(output_file, width = 11, height = 7)
  oldpar <- par(mfrow = c(3, 1), mar = c(3.4, 4.4, 2.4, 1), oma = c(0, 0, 2, 0))
  on.exit({
    par(oldpar)
    dev.off()
  }, add = TRUE)

  panel <- function(q, q_data, p, hrr_p, main) {
    ymax <- max(100 * c(q, q_data, p, hrr_p), na.rm = TRUE) * 1.08
    col_p <- "black"
    col_q <- "gray70"
    plot(series$date, 100 * p, type = "l", lwd = 2, col = col_p,
         ylim = c(0, ymax), xlab = "", ylab = "Probability (%)",
         main = main, las = 1)
    grid()
    lines(series$date, 100 * q, lwd = 2, col = col_q)
    lines(series$date, 100 * q_data, lwd = 2, col = col_q, lty = 2)
    lines(hrr$date_stata, 100 * hrr_p, lwd = 2, col = col_p, lty = 3)
    legend("topleft",
           legend = c("NGMMR Q", "data/proxy Q", "NGMMR P", "HRR reported P"),
           col = c(col_q, col_q, col_p, col_p),
           lty = c(1, 2, 1, 3),
           lwd = 2,
           bg = "white",
           box.col = "gray80")
  }

  if (side == "high") {
    panel(series$q_5y_gt4, series$q_data_5y_gt4,
          series$p_5y_gt4, hrr$zc_higher4_5y,
          expression(paste(pi[t*","*t+5] > 4, "%")))
    panel(series$q_10y_gt4, series$q_data_10y_gt4,
          series$p_10y_gt4, hrr$zc_higher4_10y,
          expression(paste(pi[t*","*t+10] > 4, "%")))
    panel(series$q_5y5y_gt4, series$q_data_5y5y_gt4,
          series$p_5y5y_gt4, hrr$higher4_5y5y,
          expression(paste(pi[t+5*","*t+10] > 4, "%")))

    mtext("Inflation tail probabilities above 4%: NGMMR P/Q and HRR reported P",
          outer = TRUE, font = 2)
  } else {
    panel(series$q_5y_lt0, series$q_data_5y_lt0,
          series$p_5y_lt0, hrr$zc_lower0_5y,
          expression(paste(pi[t*","*t+5] < 0, "%")))
    panel(series$q_10y_lt0, series$q_data_10y_lt0,
          series$p_10y_lt0, hrr$zc_lower0_10y,
          expression(paste(pi[t*","*t+10] < 0, "%")))
    panel(series$q_5y5y_lt0, series$q_data_5y5y_lt0,
          series$p_5y5y_lt0, hrr$lower0_5y5y,
          expression(paste(pi[t+5*","*t+10] < 0, "%")))

    mtext("Inflation tail probabilities below 0%: NGMMR P/Q and HRR reported P",
          outer = TRUE, font = 2)
  }
}
