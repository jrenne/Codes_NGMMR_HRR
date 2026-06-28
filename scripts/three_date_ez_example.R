# Three-date Epstein-Zin example used in the comment.
#
# The script computes the physical, real risk-neutral, and nominal
# risk-neutral probabilities reported in the simple two-state illustration. It
# writes a CSV table and a LaTeX table fragment to
# outputs/three_date_ez_example/.

rm(list = ls())

out_dir <- file.path("outputs", "three_date_ez_example")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

beta <- 0.95
gamma <- 3
psi <- 1.5
g <- c(L = 1, H = 0.75)
infl_gross <- c(L = 1.10, H = 1.25)
p1 <- c(L = 0.90, H = 0.10)

ez_mu <- function(values, probs) {
  sum(probs * values^(1 - gamma))^(1 / (1 - gamma))
}

ez_value <- function(cons, next_values, probs) {
  cont <- ez_mu(next_values, probs)
  ((1 - beta) * cons^(1 - 1 / psi) +
     beta * cont^(1 - 1 / psi))^(1 / (1 - 1 / psi))
}

ez_sdf <- function(cons_growth, next_value, mu_next) {
  beta * cons_growth^(-1 / psi) * (next_value / mu_next)^(1 / psi - gamma)
}

compute_case <- function(trans, case_name) {
  states <- c("L", "H")
  paths <- expand.grid(s1 = states, s2 = states, stringsAsFactors = FALSE)
  paths$path <- paste0(paths$s1, paths$s2)

  c0 <- 1
  c1 <- c(L = c0 * unname(g["L"]), H = c0 * unname(g["H"]))
  c2 <- setNames(numeric(nrow(paths)), paths$path)
  for (ii in seq_len(nrow(paths))) {
    c2[paths$path[ii]] <- c1[paths$s1[ii]] * g[paths$s2[ii]]
  }

  v2 <- c2
  v1 <- setNames(numeric(2), states)
  mu1 <- setNames(numeric(2), states)
  for (s1 in states) {
    idx <- paths$s1 == s1
    probs <- trans[s1, paths$s2[idx]]
    vals <- v2[paths$path[idx]]
    mu1[s1] <- ez_mu(vals, probs)
    v1[s1] <- ez_value(c1[s1], vals, probs)
  }

  mu0 <- ez_mu(v1, p1)
  m01 <- setNames(numeric(2), states)
  for (s1 in states) {
    m01[s1] <- ez_sdf(c1[s1] / c0, v1[s1], mu0)
  }

  qr1 <- p1 * m01
  qr1 <- qr1 / sum(qr1)
  qn1 <- qr1 / infl_gross
  qn1 <- qn1 / sum(qn1)

  p_path <- setNames(numeric(nrow(paths)), paths$path)
  m12 <- setNames(numeric(nrow(paths)), paths$path)
  for (ii in seq_len(nrow(paths))) {
    s1 <- paths$s1[ii]
    s2 <- paths$s2[ii]
    path <- paths$path[ii]
    idx <- paths$s1 == s1
    p_path[path] <- p1[s1] * trans[s1, s2]
    m12[path] <- ez_sdf(as.numeric(c2[path]) / as.numeric(c1[s1]),
                        as.numeric(v2[path]), as.numeric(mu1[s1]))
  }

  qr2 <- matrix(NA_real_, nrow = 2, ncol = 2, dimnames = list(states, states))
  qn2 <- matrix(NA_real_, nrow = 2, ncol = 2, dimnames = list(states, states))
  for (s1 in states) {
    idx <- paths$s1 == s1
    children <- paths$s2[idx]
    qr2[s1, children] <- trans[s1, children] * as.numeric(m12[paths$path[idx]])
    qr2[s1, ] <- qr2[s1, ] / sum(qr2[s1, ])
    qn2[s1, ] <- qr2[s1, ] / infl_gross
    qn2[s1, ] <- qn2[s1, ] / sum(qn2[s1, ])
  }

  ad <- as.numeric(p_path) * as.numeric(m01[paths$s1]) * as.numeric(m12)
  names(ad) <- paths$path
  q_path <- ad / sum(ad)
  qn_path <- q_path / (infl_gross[paths$s1] * infl_gross[paths$s2])
  qn_path <- qn_path / sum(qn_path)
  names(qn_path) <- paths$path

  event_prob <- function(path_prob, event) {
    sum(path_prob[event])
  }
  cond_prob <- function(path_prob, event, cond) {
    sum(path_prob[event & cond]) / sum(path_prob[cond])
  }

  event_defs <- list(
    "s1_H" = list(label = "$\\{s_1=H\\}$",
                  p = event_prob(p_path, paths$s1 == "H"),
                  q = event_prob(q_path, paths$s1 == "H"),
                  qn = event_prob(qn_path, paths$s1 == "H")),
    "s2_H_given_s1_L" = list(label = "$\\{s_2=H\\mid s_1=L\\}$",
                             p = cond_prob(p_path, paths$s2 == "H", paths$s1 == "L"),
                             q = cond_prob(q_path, paths$s2 == "H", paths$s1 == "L"),
                             qn = cond_prob(qn_path, paths$s2 == "H", paths$s1 == "L")),
    "s2_H_given_s1_H" = list(label = "$\\{s_2=H\\mid s_1=H\\}$",
                             p = cond_prob(p_path, paths$s2 == "H", paths$s1 == "H"),
                             q = cond_prob(q_path, paths$s2 == "H", paths$s1 == "H"),
                             qn = cond_prob(qn_path, paths$s2 == "H", paths$s1 == "H")),
    "any_H" = list(label = "$\\{s_1=H \\text{ or } s_2=H\\}$",
                   p = event_prob(p_path, paths$s1 == "H" | paths$s2 == "H"),
                   q = event_prob(q_path, paths$s1 == "H" | paths$s2 == "H"),
                   qn = event_prob(qn_path, paths$s1 == "H" | paths$s2 == "H"))
  )

  rows <- lapply(names(event_defs), function(name) {
    ev <- event_defs[[name]]
    data.frame(
      case = case_name,
      event_id = name,
      event = ev$label,
      P = ev$p,
      Q_real = ev$q,
      P_over_Q_real = ev$p / ev$q,
      Q_nominal = ev$qn,
      P_over_Q_nominal = ev$p / ev$qn,
      stringsAsFactors = FALSE
    )
  })

  list(
    rows = do.call(rbind, rows),
    paths = data.frame(
      case = case_name,
      path = paths$path,
      P = as.numeric(p_path[paths$path]),
      AD = as.numeric(ad[paths$path]),
      Q_real = as.numeric(q_path[paths$path]),
      Q_nominal_T2 = as.numeric(qn_path[paths$path]),
      stringsAsFactors = FALSE
    )
  )
}

trans_iid <- matrix(
  c(0.90, 0.10,
    0.90, 0.10),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("L", "H"), c("L", "H"))
)

trans_markov <- matrix(
  c(0.90, 0.10,
    0.30, 0.70),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("L", "H"), c("L", "H"))
)

res_iid <- compute_case(trans_iid, "iid")
res_markov <- compute_case(trans_markov, "persistent_markov")

event_table <- rbind(res_iid$rows, res_markov$rows)
path_table <- rbind(res_iid$paths, res_markov$paths)

wide <- merge(
  event_table[event_table$case == "iid",
              c("event_id", "event", "P", "Q_real", "P_over_Q_real",
                "P_over_Q_nominal")],
  event_table[event_table$case == "persistent_markov",
              c("event_id", "P", "Q_real", "P_over_Q_real",
                "P_over_Q_nominal")],
  by = "event_id",
  suffixes = c("_iid", "_markov")
)
wide <- wide[match(c("s1_H", "s2_H_given_s1_L", "s2_H_given_s1_H",
                     "any_H"),
                   wide$event_id), ]

fmt <- function(x) sprintf("%.3f", x)
latex_rows <- apply(wide, 1, function(z) {
  paste0(
    z[["event"]], " & ",
    fmt(as.numeric(z[["P_iid"]])), " & ",
    fmt(as.numeric(z[["Q_real_iid"]])), " & ",
    fmt(as.numeric(z[["P_over_Q_real_iid"]])), " & ",
    fmt(as.numeric(z[["P_over_Q_nominal_iid"]])), " & ",
    fmt(as.numeric(z[["P_markov"]])), " & ",
    fmt(as.numeric(z[["Q_real_markov"]])), " & ",
    fmt(as.numeric(z[["P_over_Q_real_markov"]])), " & ",
    fmt(as.numeric(z[["P_over_Q_nominal_markov"]])), " \\\\"
  )
})

latex_file <- file.path(out_dir, "three_date_ez_table_rows.tex")
writeLines(latex_rows, latex_file)

print(wide[, c("event", "P_iid", "Q_real_iid", "P_over_Q_real_iid",
               "P_over_Q_nominal_iid", "P_markov", "Q_real_markov",
               "P_over_Q_real_markov", "P_over_Q_nominal_markov")])
message("\nWrote outputs to: ", normalizePath(out_dir))
