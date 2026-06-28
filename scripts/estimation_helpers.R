## Helpers for date-by-date nominal RN matrix estimation.

suppressPackageStartupMessages({
  library(parallel)
  library(optimx)
})

source(file.path("scripts", "qn_parameterization.R"))

parse_month_ids <- function(n_month) {
  month_ids_env <- Sys.getenv("MONTH_IDS", unset = "")
  if (!nzchar(month_ids_env)) {
    return(seq_len(n_month))
  }
  ids <- as.integer(strsplit(month_ids_env, ",", fixed = TRUE)[[1]])
  ids[!is.na(ids)]
}

make_month_labels <- function(n_month, start_date = as.Date("2011-01-01")) {
  format(seq.Date(start_date, by = "month", length.out = n_month), "%Y-%m")
}

extract_optimx_par <- function(fit, n_par) {
  as.numeric(fit[1, paste0("p", seq_len(n_par))])
}

fit_with_nlminb_nm_loops <- function(z0, obj,
                                     nb_loop = 3,
                                     nlminb_iter = 20,
                                     nm_iter = 1000,
                                     trace = TRUE) {
  z_cur <- z0
  z_best <- z0
  f_best <- obj(z_best)

  history <- data.frame(loop = 0L, method = "start",
                        value = f_best, convergence = 0L)

  for (loop_id in seq_len(nb_loop)) {
    fit_nl <- nlminb(
      start = z_cur,
      objective = obj,
      control = list(iter.max = nlminb_iter,
                     eval.max = max(100, 20 * nlminb_iter))
    )
    if (is.finite(fit_nl$objective) && fit_nl$objective < f_best) {
      z_best <- fit_nl$par
      f_best <- fit_nl$objective
    }

    history <- rbind(
      history,
      data.frame(loop = loop_id, method = "nlminb",
                 value = fit_nl$objective,
                 convergence = fit_nl$convergence)
    )

    fit_nm <- optimx(
      par = fit_nl$par,
      fn = obj,
      method = "Nelder-Mead",
      control = list(maxit = nm_iter)
    )
    z_nm <- extract_optimx_par(fit_nm, length(z0))
    f_nm <- as.numeric(fit_nm$value[1])
    if (is.finite(f_nm) && f_nm < f_best) {
      z_best <- z_nm
      f_best <- f_nm
    }
    z_cur <- z_best

    history <- rbind(
      history,
      data.frame(loop = loop_id, method = "Nelder-Mead",
                 value = f_nm,
                 convergence = as.integer(fit_nm$convcode[1]))
    )

    if (trace) {
      message(sprintf("Loop %d: nlminb %.8g, NM %.8g, best %.8g",
                      loop_id, fit_nl$objective, f_nm, f_best))
    }
  }

  list(par = z_best, value = f_best, history = history)
}

rmse_pp <- function(x) {
  100 * sqrt(mean(x^2))
}

save_matrix_stack <- function(results, slot, output_file) {
  month_ids <- vapply(results, function(x) x$diagnostics$month, integer(1))
  arr <- array(NA_real_, dim = c(8, 8, length(results)),
               dimnames = list(paste0("from_", seq_len(8)),
                               paste0("to_", seq_len(8)),
                               paste0("month_", month_ids)))
  for (i in seq_along(results)) {
    arr[, , i] <- results[[i]][[slot]]
  }
  saveRDS(arr, output_file)
  invisible(arr)
}
