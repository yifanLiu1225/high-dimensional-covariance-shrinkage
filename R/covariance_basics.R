# 组员 A：Wishart 基础与样本协方差风险模拟
#
# 本脚本对应 project_plan.Rmd 中的组员 A 任务：
# 1. 生成多元正态样本并计算样本协方差；
# 2. 数值验证 Wishart 关系、无偏性和样本均值/协方差独立性；
# 3. 比较不同 p/n 下的元素误差、Frobenius 风险和谱稳定性；
# 4. 输出本地表格和图形，不涉及 GitHub 上传。

required_packages <- c(
  "MASS", "mvtnorm", "Matrix", "ggplot2",
  "dplyr", "tidyr", "purrr", "readr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

make_covariance <- function(p, structure = c("identity", "ar1", "block", "spiked"),
                            rho = 0.6, spike = 6) {
  structure <- match.arg(structure)

  if (structure == "identity") {
    return(diag(p))
  }

  if (structure == "ar1") {
    return(rho ^ abs(outer(seq_len(p), seq_len(p), "-")))
  }

  if (structure == "block") {
    groups <- rep(seq_len(ceiling(p / 10)), each = 10, length.out = p)
    Sigma <- matrix(0.1, nrow = p, ncol = p)
    Sigma[groups == groups[col(Sigma)]] <- 0.7
    diag(Sigma) <- 1
    return(Sigma)
  }

  if (structure == "spiked") {
    set.seed(1000 + p)
    Q <- qr.Q(qr(matrix(rnorm(p * p), nrow = p)))
    eig <- c(spike, 3, rep(1, max(0, p - 2)))
    return(Q %*% diag(eig, nrow = p) %*% t(Q))
  }
}

theoretical_frobenius_risk <- function(Sigma, n) {
  nu <- n - 1
  ((sum(diag(Sigma))^2) + sum(Sigma^2)) / nu
}

draw_normal_sample <- function(n, Sigma, mu = NULL) {
  p <- nrow(Sigma)
  if (is.null(mu)) mu <- rep(0, p)
  mvtnorm::rmvnorm(n = n, mean = mu, sigma = Sigma)
}

safe_spectral_metrics <- function(S) {
  eig <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  eig <- sort(Re(eig), decreasing = TRUE)
  tolerance <- max(1, abs(eig[1])) * 1e-8
  positive_eig <- eig[eig > tolerance]

  list(
    largest_eigenvalue = eig[1],
    smallest_eigenvalue = eig[length(eig)],
    effective_smallest_eigenvalue = if (length(positive_eig)) min(positive_eig) else NA_real_,
    condition_number = if (length(positive_eig) == nrow(S)) {
      max(positive_eig) / min(positive_eig)
    } else {
      Inf
    },
    invertible = length(positive_eig) == nrow(S)
  )
}

simulate_covariance_once <- function(n, p, structure = "identity", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- make_covariance(p, structure)
  X <- draw_normal_sample(n, Sigma)
  S <- stats::cov(X)
  error_matrix <- S - Sigma
  spectral <- safe_spectral_metrics(S)

  data.frame(
    n = n,
    p = p,
    p_over_n = p / n,
    structure = structure,
    frobenius_error = sum(error_matrix^2),
    element_mse = mean(error_matrix^2),
    theoretical_frobenius_risk = theoretical_frobenius_risk(Sigma, n),
    largest_eigenvalue = spectral$largest_eigenvalue,
    smallest_eigenvalue = spectral$smallest_eigenvalue,
    effective_smallest_eigenvalue = spectral$effective_smallest_eigenvalue,
    condition_number = spectral$condition_number,
    invertible = spectral$invertible
  )
}

run_covariance_risk_simulation <- function(
    regimes = data.frame(
      regime = c("low", "moderate", "high", "near_singular", "singular"),
      n = c(500, 200, 200, 150, 100),
      p = c(10, 50, 150, 140, 200)
    ),
    structures = c("identity", "ar1"),
    reps = 100,
    seed = 20261030
) {
  set.seed(seed)
  purrr::map_dfr(
    seq_len(nrow(regimes)),
    function(i) {
      purrr::map_dfr(
        structures,
        function(structure) {
          purrr::map_dfr(
            seq_len(reps),
            function(rep_id) {
              simulate_covariance_once(
                n = regimes$n[i],
                p = regimes$p[i],
                structure = structure,
                seed = sample.int(1e8, 1)
              ) |>
                dplyr::mutate(
                  regime = regimes$regime[i],
                  rep = rep_id
                )
            }
          )
        }
      )
    }
  )
}

run_wishart_checks <- function(reps = 2000, n = 120, seed = 20261030) {
  set.seed(seed)
  Sigma <- matrix(
    c(1.0, 0.55, 0.25,
      0.55, 1.4, 0.35,
      0.25, 0.35, 0.8),
    nrow = 3,
    byrow = TRUE
  )

  samples <- purrr::map_dfr(
    seq_len(reps),
    function(rep_id) {
      X <- draw_normal_sample(n, Sigma)
      S <- stats::cov(X)
      xbar <- colMeans(X)
      data.frame(
        rep = rep_id,
        mean_1 = xbar[1],
        mean_2 = xbar[2],
        mean_3 = xbar[3],
        S11 = S[1, 1],
        S12 = S[1, 2],
        S22 = S[2, 2]
      )
    }
  )

  empirical_covariance <- c(
    S11 = mean(samples$S11),
    S12 = mean(samples$S12),
    S22 = mean(samples$S22)
  )
  target_covariance <- c(
    S11 = Sigma[1, 1],
    S12 = Sigma[1, 2],
    S22 = Sigma[2, 2]
  )

  mean_covariance <- c(
    mean_1_S12 = cor(samples$mean_1, samples$S12),
    mean_1_S11 = cor(samples$mean_1, samples$S11),
    mean_2_S12 = cor(samples$mean_2, samples$S12)
  )

  list(
    target_covariance = target_covariance,
    empirical_covariance = empirical_covariance,
    bias = empirical_covariance - target_covariance,
    mean_covariance_correlations = mean_covariance
  )
}

summarise_risk_results <- function(results) {
  results |>
    dplyr::group_by(regime, n, p, p_over_n, structure) |>
    dplyr::summarise(
      reps = dplyr::n(),
      frobenius_risk_mean = mean(frobenius_error),
      frobenius_risk_sd = stats::sd(frobenius_error),
      theoretical_risk = mean(theoretical_frobenius_risk),
      element_mse_mean = mean(element_mse),
      largest_eigenvalue_mean = mean(largest_eigenvalue),
      smallest_eigenvalue_mean = mean(smallest_eigenvalue),
      condition_number_median = if (all(is.infinite(condition_number))) {
        Inf
      } else {
        stats::median(condition_number[is.finite(condition_number)])
      },
      singular_rate = mean(!invertible),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      risk_ratio = frobenius_risk_mean / theoretical_risk,
      risk_mc_se = frobenius_risk_sd / sqrt(reps)
    )
}

plot_frobenius_risk <- function(summary_results, structure_name = "identity") {
  plot_data <- summary_results |>
    dplyr::filter(structure == structure_name) |>
    dplyr::mutate(
      lower = pmax(0, frobenius_risk_mean - 1.96 * risk_mc_se),
      upper = frobenius_risk_mean + 1.96 * risk_mc_se
    )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = p_over_n, y = frobenius_risk_mean, color = regime, group = 1)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      width = 0.015,
      alpha = 0.6
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = theoretical_risk),
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title = "Sample covariance Frobenius risk",
      subtitle = "Dashed line: exact theoretical risk under the same covariance structure",
      x = "p / n",
      y = "Frobenius risk (log scale)",
      color = "Regime"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

plot_singularity_rate <- function(summary_results) {
  ggplot2::ggplot(
    summary_results,
    ggplot2::aes(x = p_over_n, y = singular_rate, color = structure)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70") +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Singularity rate of the sample covariance",
      x = "p / n",
      y = "Empirical singularity rate",
      color = "Covariance structure"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

save_member_A_outputs <- function(
    output_table_dir = "tables/member_A",
    output_figure_dir = "figures/member_A",
    reps = 100,
    seed = 20261030
) {
  dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)

  risk_results <- run_covariance_risk_simulation(reps = reps, seed = seed)
  risk_summary <- summarise_risk_results(risk_results)
  wishart_checks <- run_wishart_checks(seed = seed)

  readr::write_csv(risk_results, file.path(output_table_dir, "covariance_risk_simulations.csv"))
  readr::write_csv(risk_summary, file.path(output_table_dir, "covariance_risk_summary.csv"))
  readr::write_csv(
    data.frame(
      quantity = names(wishart_checks$target_covariance),
      target = unname(wishart_checks$target_covariance),
      empirical = unname(wishart_checks$empirical_covariance),
      bias = unname(wishart_checks$bias)
    ),
    file.path(output_table_dir, "wishart_unbiasedness_check.csv")
  )
  readr::write_csv(
    data.frame(
      quantity = names(wishart_checks$mean_covariance_correlations),
      correlation = unname(wishart_checks$mean_covariance_correlations)
    ),
    file.path(output_table_dir, "mean_covariance_independence_check.csv")
  )

  risk_plot <- plot_frobenius_risk(risk_summary, structure = "identity")
  singularity_plot <- plot_singularity_rate(risk_summary)

  ggplot2::ggsave(
    file.path(output_figure_dir, "frobenius_risk_vs_p_over_n.png"),
    risk_plot,
    width = 9,
    height = 6,
    dpi = 160
  )
  ggplot2::ggsave(
    file.path(output_figure_dir, "sample_covariance_singularity_rate.png"),
    singularity_plot,
    width = 9,
    height = 6,
    dpi = 160
  )

  invisible(list(
    risk_results = risk_results,
    risk_summary = risk_summary,
    wishart_checks = wishart_checks
  ))
}

