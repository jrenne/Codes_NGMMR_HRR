## Nominal RN transition-matrix parameterizations used in the comment.

suppressPackageStartupMessages({
  library(R.matlab)
})

source(file.path("scripts", "hrr_model_bins.R"))

hrr_zero_pattern <- function() {
  amatrix_hrr_101(rep(0.05, 6)) > 0
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

eta_to_structured_mlogit_A <- function(eta,
                                       mask = hrr_zero_pattern(),
                                       variant = "smooth_row_poly2") {
  if (variant != "smooth_row_poly2") {
    stop("Only the active 'smooth_row_poly2' specification is kept here. ",
         "Exploratory variants are not needed to reproduce the comment outputs.")
  }
  eta_to_smooth_row_poly2_A(eta, mask = mask)
}

qn_variant_eta_length <- function(variant = "smooth_row_poly2",
                                  mask = hrr_zero_pattern()) {
  if (variant != "smooth_row_poly2") {
    stop("Only the active 'smooth_row_poly2' specification is kept here.")
  }
  6L
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
