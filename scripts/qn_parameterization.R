## Nominal RN transition-matrix parameterizations used in the comment.

suppressPackageStartupMessages({
  library(R.matlab)
})

source(file.path("scripts", "hrr_model_bins.R"))

hrr_zero_pattern <- function() {
  amatrix_hrr_101(rep(0.05, 6)) > 0
}

extended_pi_grid <- function(n_state = 8L) {
  if (n_state == 8L) {
    return(c(-2, seq(-0.5, 4.5, by = 1), 6))
  }
  stop("Only the 8-state HRR grid is active.")
}

stable_softmax <- function(x) {
  ex <- exp(x - max(x))
  ex / sum(ex)
}

scale_state <- function(k, n_state = 8) {
  (k - (n_state + 1) / 2) / ((n_state - 1) / 2)
}

eta_to_smooth_row_poly2_A <- function(eta, mask = hrr_zero_pattern()) {
  ## Six-parameter NGMMR specification.  Each row is a quadratic logit in the
  ## destination state, and the row-specific slope/curvature vary smoothly with
  ## the current state:
  ##   beta1(i) = a0 + a1 s_i + a2 s_i^2,
  ##   beta2(i) = b0 + b1 s_i + b2 s_i^2.
  stopifnot(length(eta) == 6L)

  n_state <- nrow(mask)
  A <- matrix(0, n_state, n_state)

  for (i in seq_len(n_state)) {
    allowed <- which(mask[i, ])
    si <- scale_state(i, n_state)
    sj <- scale_state(allowed, n_state)
    beta1 <- eta[1] + eta[2] * si + eta[3] * si^2
    beta2 <- eta[4] + eta[5] * si + eta[6] * si^2
    A[i, allowed] <- stable_softmax(beta1 * sj + beta2 * sj^2)
  }

  A
}

eta_to_smooth_row_poly3_A <- function(eta, mask = hrr_zero_pattern()) {
  ## Nine-parameter smooth cubic logit.  Each row is a cubic logit in the
  ## destination state; the three destination-polynomial coefficients vary
  ## quadratically with the current state.
  stopifnot(length(eta) == 9L)

  n_state <- nrow(mask)
  A <- matrix(0, n_state, n_state)

  for (i in seq_len(n_state)) {
    allowed <- which(mask[i, ])
    si <- scale_state(i, n_state)
    sj <- scale_state(allowed, n_state)
    beta1 <- eta[1] + eta[2] * si + eta[3] * si^2
    beta2 <- eta[4] + eta[5] * si + eta[6] * si^2
    beta3 <- eta[7] + eta[8] * si + eta[9] * si^2
    A[i, allowed] <- stable_softmax(beta1 * sj + beta2 * sj^2 + beta3 * sj^3)
  }

  A
}

eta_to_smooth_row_poly3_extreme_A <- function(eta, mask = hrr_zero_pattern()) {
  ## Sixteen-parameter smooth cubic logit with separate extreme-destination
  ## scores.  Middle destinations use the smooth cubic specification; states 1
  ## and 8 have their own quadratic functions of the current state, plus one
  ## common tilt for persistence in the two extreme rows.
  stopifnot(length(eta) == 16L)

  n_state <- nrow(mask)
  A <- matrix(0, n_state, n_state)

  for (i in seq_len(n_state)) {
    allowed <- which(mask[i, ])
    si <- scale_state(i, n_state)
    sj <- scale_state(allowed, n_state)
    beta1 <- eta[1] + eta[2] * si + eta[3] * si^2
    beta2 <- eta[4] + eta[5] * si + eta[6] * si^2
    beta3 <- eta[7] + eta[8] * si + eta[9] * si^2
    scores <- beta1 * sj + beta2 * sj^2 + beta3 * sj^3

    low_pos <- which(allowed == 1L)
    if (length(low_pos) == 1L) {
      scores[low_pos] <- eta[10] + eta[11] * si + eta[12] * si^2
    }

    high_pos <- which(allowed == n_state)
    if (length(high_pos) == 1L) {
      scores[high_pos] <- eta[13] + eta[14] * si + eta[15] * si^2
    }

    diag_extreme_pos <- which(allowed == i & i %in% c(1L, n_state))
    if (length(diag_extreme_pos) == 1L) {
      scores[diag_extreme_pos] <- scores[diag_extreme_pos] + eta[16]
    }

    A[i, allowed] <- stable_softmax(scores)
  }

  A
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

eta_to_monotone_tail_common_middle_A <- function(eta, mask = hrr_zero_pattern()) {
  ## Eight-parameter exploratory specification.
  ## - Low-tail and high-tail masses are monotone in the current state.
  ## - Remaining mass is allocated over admissible middle states using a common
  ##   quadratic destination-state shape.
  stopifnot(length(eta) == 8L)

  n_state <- nrow(mask)
  A <- matrix(0, n_state, n_state)
  s <- scale_state(seq_len(n_state), n_state)

  low_shape <- eta[1] - exp(eta[2]) * s
  high_shape <- eta[3] + exp(eta[4]) * s
  total_tail <- 0.90 * inv_logit(eta[5])
  tail_split <- exp(cbind(low_shape, high_shape))
  low_mass <- total_tail * tail_split[, 1] / rowSums(tail_split)
  high_mass <- total_tail * tail_split[, 2] / rowSums(tail_split)

  for (i in seq_len(n_state)) {
    allowed <- which(mask[i, ])
    middle <- setdiff(allowed, c(1L, n_state))
    if (1L %in% allowed) {
      A[i, 1L] <- low_mass[i]
    }
    if (n_state %in% allowed) {
      A[i, n_state] <- high_mass[i]
    }

    rem <- 1 - sum(A[i, ])
    if (length(middle) > 0) {
      sj <- scale_state(middle, n_state)
      scores <- eta[6] * sj + eta[7] * sj^2 + eta[8] * s[i] * sj
      A[i, middle] <- rem * stable_softmax(scores)
    } else {
      A[i, allowed] <- A[i, allowed] / sum(A[i, allowed])
    }
  }

  A
}

eta_to_common_shape_tail_shift_A <- function(eta, mask = hrr_zero_pattern()) {
  ## Seven-parameter exploratory specification.
  ## - Middle destinations follow one common quadratic shape, shifted by the
  ##   current state.
  ## - Low- and high-tail scores move monotonically with the current state.
  stopifnot(length(eta) == 7L)

  n_state <- nrow(mask)
  A <- matrix(0, n_state, n_state)
  s <- scale_state(seq_len(n_state), n_state)

  for (i in seq_len(n_state)) {
    allowed <- which(mask[i, ])
    sj <- scale_state(allowed, n_state)
    scores <- eta[1] * sj + eta[2] * sj^2 + eta[3] * s[i] * sj

    low_pos <- which(allowed == 1L)
    high_pos <- which(allowed == n_state)
    if (length(low_pos) == 1L) {
      scores[low_pos] <- scores[low_pos] + eta[4] - exp(eta[5]) * s[i]
    }
    if (length(high_pos) == 1L) {
      scores[high_pos] <- scores[high_pos] + eta[6] + exp(eta[7]) * s[i]
    }

    A[i, allowed] <- stable_softmax(scores)
  }

  A
}

eta_to_mean_reverting_shape_A <- function(eta, mask = hrr_zero_pattern()) {
  ## Six-parameter ordered-state specification.  Each row is centered around a
  ## mean-reverting conditional mean, with a common smooth shape over ordered
  ## destination states.
  stopifnot(length(eta) == 6L)

  n_state <- nrow(mask)
  pi_bar <- extended_pi_grid(n_state)
  A <- matrix(0, n_state, n_state)
  s <- scale_state(seq_len(n_state), n_state)

  rho <- inv_logit(eta[1])
  pi_star <- -0.5 + 6.5 * inv_logit(eta[2])
  sigma_level <- 0.40 + exp(eta[3])
  sigma_curv <- 1.5 * tanh(eta[4])
  skew <- 1.5 * tanh(eta[5])
  tail_shape <- 0.4 * tanh(eta[6])

  for (i in seq_len(n_state)) {
    allowed <- which(mask[i, ])
    mu_i <- (1 - rho) * pi_star + rho * pi_bar[i]
    sigma_i <- sigma_level * exp(sigma_curv * s[i]^2)
    x <- (pi_bar[allowed] - mu_i) / sigma_i
    scores <- -0.5 * x^2 + skew * x + tail_shape * log1p(x^2)
    probs <- pmax(stable_softmax(scores), 1e-8)
    A[i, allowed] <- probs / sum(probs)
  }

  A
}

eta_to_hrr_endpoints_smooth_middle_A <- function(eta, mask = hrr_zero_pattern()) {
  ## Eight-parameter hybrid.  Rows 1 and 8 use HRR's endpoint structure, with
  ## persistence capped at 0.90 by construction.  Rows 2--7 use the same smooth
  ## quadratic logit shape as eta_to_smooth_row_poly2_A.
  stopifnot(length(eta) == 8L, nrow(mask) == 8L)

  n_state <- nrow(mask)
  A <- matrix(0, n_state, n_state)
  x_min <- 0.02
  x1 <- x_min + (0.2 - x_min) * inv_logit(eta[1])
  x2 <- x_min + (0.2 - x_min) * inv_logit(eta[2])

  A[1, 1] <- 1 - 5 * x1
  A[1, 2:6] <- x1
  A[8, 3:7] <- x2
  A[8, 8] <- 1 - 5 * x2

  middle_eta <- eta[3:8]
  for (i in 2:7) {
    allowed <- which(mask[i, ])
    si <- scale_state(i, n_state)
    sj <- scale_state(allowed, n_state)
    beta1 <- middle_eta[1] + middle_eta[2] * si + middle_eta[3] * si^2
    beta2 <- middle_eta[4] + middle_eta[5] * si + middle_eta[6] * si^2
    A[i, allowed] <- stable_softmax(beta1 * sj + beta2 * sj^2)
  }

  A
}

eta_to_structured_mlogit_A <- function(eta,
                                       mask = hrr_zero_pattern(),
                                       variant = "smooth_row_poly2") {
  if (variant == "smooth_row_poly2") {
    return(eta_to_smooth_row_poly2_A(eta, mask = mask))
  }
  if (variant == "smooth_row_poly3") {
    return(eta_to_smooth_row_poly3_A(eta, mask = mask))
  }
  if (variant == "smooth_row_poly3_extreme") {
    return(eta_to_smooth_row_poly3_extreme_A(eta, mask = mask))
  }
  if (variant == "monotone_tail_common_middle") {
    return(eta_to_monotone_tail_common_middle_A(eta, mask = mask))
  }
  if (variant == "common_shape_tail_shift") {
    return(eta_to_common_shape_tail_shift_A(eta, mask = mask))
  }
  if (variant == "mean_reverting_shape") {
    return(eta_to_mean_reverting_shape_A(eta, mask = mask))
  }
  if (variant == "hrr_endpoints_smooth_middle") {
    return(eta_to_hrr_endpoints_smooth_middle_A(eta, mask = mask))
  }
  stop("Unknown nominal-Q parameterization: ", variant)
}

qn_variant_eta_length <- function(variant = "smooth_row_poly2",
                                  mask = hrr_zero_pattern()) {
  if (variant == "smooth_row_poly2") {
    return(6L)
  }
  if (variant == "smooth_row_poly3") {
    return(9L)
  }
  if (variant == "smooth_row_poly3_extreme") {
    return(16L)
  }
  if (variant == "monotone_tail_common_middle") {
    return(8L)
  }
  if (variant == "common_shape_tail_shift") {
    return(7L)
  }
  if (variant == "mean_reverting_shape") {
    return(6L)
  }
  if (variant == "hrr_endpoints_smooth_middle") {
    return(8L)
  }
  stop("Unknown nominal-Q parameterization: ", variant)
}

matrix_power <- function(A, n) {
  stopifnot(n >= 0, n == as.integer(n), nrow(A) == ncol(A))
  if (n == 0) {
    return(diag(nrow(A)))
  }

  out <- diag(nrow(A))
  for (k in seq_len(n)) {
    out <- out %*% A
  }
  out
}

avg_yoy_6to10_model <- function(A, bt) {
  powers <- lapply(5:9, function(h) matrix_power(A, h))
  as.numeric(A[bt, ] %*% Reduce("+", powers) / 5)
}

forward_average_bins_model <- function(A, bt,
                                       pi_bar = c(-2, seq(-0.5, 4.5, by = 1), 6),
                                       pi_lim = c(seq(-1, 5, by = 1), Inf)) {
  n_state <- length(pi_bar)
  n_bin <- length(pi_lim)
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

  bin <- n_bin + 1L -
    rowSums(outer(avg, pi_lim, function(a, lim) lim - a > -1e-7))
  out <- numeric(n_bin)
  for (b in seq_len(n_bin)) {
    out[b] <- sum(cum_prob[used[bin == b]])
  }
  out / sum(out)
}
