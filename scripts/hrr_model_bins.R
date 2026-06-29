## Fast model-implied 5y/10y bin probabilities for HRR model 101.
##
## This script implements a pure-R version of HRR's avg_inf_510.m, but
## stores cumulative inflation on an integer half-percentage-point grid.
## It then validates the output against HRR's saved model-implied bins.

suppressPackageStartupMessages({
  library(R.matlab)
})

input_dir <- file.path(getwd(), "input")

unwrap_matlab_cell <- function(x) {
  while (is.list(x) && length(x) == 1) {
    x <- x[[1]]
  }
  x
}

logistic <- function(z) {
  1 / (1 + exp(-z))
}

logit <- function(p) {
  log(p / (1 - p))
}

clip01 <- function(p, eps = 1e-10) {
  pmin(pmax(p, eps), 1 - eps)
}

z_to_hrr_theta <- function(z) {
  stopifnot(length(z) == 6)
  L <- logistic(z)

  x1 <- 0.2 * L[1]
  x2 <- 0.2 * L[2]

  S <- L[3]
  x3 <- S * L[4]
  x4 <- S * (1 - L[4])

  x5 <- 0.5 * (1 - S) * L[5]
  x6 <- (1 - S - x5) * L[6]

  c(x1, x2, x3, x4, x5, x6)
}

hrr_theta_to_z <- function(theta, eps = 1e-10) {
  stopifnot(length(theta) == 6)
  x <- pmax(as.numeric(theta), 0)

  S <- x[3] + x[4]
  denom_56 <- 1 - S
  denom_6 <- 1 - S - x[5]

  c(
    logit(clip01(5 * x[1], eps)),
    logit(clip01(5 * x[2], eps)),
    logit(clip01(S, eps)),
    logit(clip01(x[3] / S, eps)),
    logit(clip01(2 * x[5] / denom_56, eps)),
    logit(clip01(x[6] / denom_6, eps))
  )
}

amatrix_hrr_101 <- function(x) {
  stopifnot(length(x) == 6)

  dm <- c(
    -x[3] - x[4] - x[5] - x[6],
    -x[3] - x[4] - 2 * x[5]
  )
  dm <- c(dm, rev(dm))

  d <- c(
    -5 * x[1],
    -x[3] - x[5] - x[6],
    dm,
    -x[4] - x[5] - x[6],
    -5 * x[2]
  )

  d1 <- c(x[1], x[6], x[6], x[5], x[5], x[5], x[4] + x[5])
  d_1 <- c(x[3] + x[5], x[5], x[5], x[5], x[6], x[6], x[2])

  A <- diag(8)
  A <- A + diag(d, nrow = 8)
  A[cbind(1:7, 2:8)] <- A[cbind(1:7, 2:8)] + d1
  A[cbind(2:8, 1:7)] <- A[cbind(2:8, 1:7)] + d_1
  A[3:6, 1] <- x[3]
  A[3:6, 8] <- x[4]
  A[1, 3:6] <- x[1]
  A[8, 3:6] <- x[2]
  A
}

is_valid_transition_matrix <- function(A, tol = 1e-10) {
  is.matrix(A) &&
    all(is.finite(A)) &&
    min(A) >= -tol &&
    max(abs(rowSums(A) - 1)) <= tol
}

avg_inf_510_fast <- function(A, bt,
                             pi_bar = c(-2, seq(-0.5, 4.5, by = 1), 6),
                             pi_lim = c(seq(-1, 5, by = 1), Inf)) {
  n_state <- length(pi_bar)
  n_bin <- length(pi_lim)
  stopifnot(
    all(dim(A) == c(n_state, n_state)),
    bt >= 1,
    bt <= n_state
  )

  ## HRR's pi_bar values are half-integers, so cumulative sums can be
  ## stored exactly as integer multiples of 0.5 percentage points.
  step <- 2L
  x_int <- as.integer(round(step * pi_bar))
  min_step <- min(x_int)
  max_step <- max(x_int)
  h2 <- 10L

  min_sum <- h2 * min_step
  max_sum <- h2 * max_step
  offset <- 1L - min_sum
  n_sum <- max_sum - min_sum + 1L

  F <- matrix(0, nrow = n_state, ncol = n_sum)
  F[bt, offset] <- 1

  out <- list()
  for (h in seq_len(h2)) {
    F_new <- matrix(0, nrow = n_state, ncol = n_sum)

    for (j in seq_len(n_state)) {
      ## Probability mass over previous cumulative sums after moving to j.
      mass_by_sum <- as.numeric(crossprod(A[, j], F))
      idx_old <- which(mass_by_sum != 0)
      idx_new <- idx_old + x_int[j]
      keep <- idx_new >= 1L & idx_new <= n_sum
      F_new[j, idx_new[keep]] <- F_new[j, idx_new[keep]] + mass_by_sum[idx_old[keep]]
    }

    F <- F_new

    if (h %in% c(5L, 10L)) {
      cum_prob <- colSums(F)
      used <- which(cum_prob != 0)
      cum_sum <- used - offset
      avg <- (cum_sum / step) / h

      bin <- n_bin + 1L - rowSums(outer(avg, pi_lim, function(a, lim) lim - a > -1e-7))
      gg <- numeric(n_bin)
      for (b in seq_len(n_bin)) {
        gg[b] <- sum(cum_prob[used[bin == b]])
      }
      out[[as.character(h)]] <- gg
    }
  }

  list(g5 = out[["5"]], g10 = out[["10"]])
}

normalize_hrr_area <- function(area = Sys.getenv("AREA", unset = "US")) {
  area <- toupper(area)
  if (area %in% c("EA", "EU", "EURO", "EUROAREA", "EURO_AREA")) {
    area <- "EZ"
  }
  if (!(area %in% c("US", "EZ"))) {
    stop("AREA must be US or EZ.")
  }
  area
}

load_hrr_monthquart <- function(area = Sys.getenv("AREA", unset = "US")) {
  area <- normalize_hrr_area(area)
  readMat(file.path(input_dir, sprintf("results_est_101_%s_monthquart.mat", area)))
}

load_hrr_dists <- function(area = Sys.getenv("AREA", unset = "US")) {
  area <- normalize_hrr_area(area)
  readMat(file.path(input_dir, sprintf("dists_%s.mat", area)))
}

load_hrr_monthquart_us <- function() {
  load_hrr_monthquart("US")
}

load_hrr_dists_us <- function() {
  load_hrr_dists("US")
}

compute_bts <- function(infl, pi_lim = c(seq(-1, 5, by = 1), Inf)) {
  vapply(as.numeric(infl), function(z) which(pi_lim >= z)[1], integer(1))
}

validate_against_hrr_saved <- function(test_months = c(1, 40, 80, 121, 122, 123, 166)) {
  res <- load_hrr_monthquart_us()
  dists <- load_hrr_dists_us()

  pinew <- unwrap_matlab_cell(res$pinew)
  gs_model <- unwrap_matlab_cell(res$gs.model)
  infl <- as.numeric(dists$infl)
  bts <- compute_bts(infl)[seq_len(ncol(pinew))]

  rows <- lapply(test_months, function(i) {
    x <- pinew[, i]
    A <- amatrix_hrr_101(x)
    fast <- avg_inf_510_fast(A, bts[i])
    hrr_g5 <- gs_model[, 1, i]
    hrr_g10 <- gs_model[, 2, i]

    data.frame(
      month = i,
      max_abs_g5 = max(abs(fast$g5 - hrr_g5)),
      max_abs_g10 = max(abs(fast$g10 - hrr_g10)),
      sum_g5 = sum(fast$g5),
      sum_g10 = sum(fast$g10),
      row_min = min(rowSums(A)),
      row_max = max(rowSums(A))
    )
  })

  do.call(rbind, rows)
}

validate_parameterization <- function(n_draws = 1000, seed = 123) {
  set.seed(seed)

  random_ok <- replicate(n_draws, {
    z <- rnorm(6, sd = 3)
    theta <- z_to_hrr_theta(z)
    A <- amatrix_hrr_101(theta)
    is_valid_transition_matrix(A)
  })

  res <- load_hrr_monthquart_us()
  pinew <- unwrap_matlab_cell(res$pinew)
  month_id <- seq_len(ncol(pinew))

  roundtrip <- lapply(month_id, function(i) {
    theta0 <- pinew[, i]
    z <- hrr_theta_to_z(theta0)
    theta1 <- z_to_hrr_theta(z)
    A <- amatrix_hrr_101(theta1)
    data.frame(
      month = i,
      max_abs_theta = max(abs(theta1 - theta0)),
      valid_A = is_valid_transition_matrix(A)
    )
  })
  roundtrip <- do.call(rbind, roundtrip)

  list(
    random_valid_share = mean(random_ok),
    roundtrip_max_abs_theta = max(roundtrip$max_abs_theta),
    roundtrip_all_valid = all(roundtrip$valid_A),
    roundtrip = roundtrip
  )
}
