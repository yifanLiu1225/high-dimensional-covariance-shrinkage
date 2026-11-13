# 组员 B：线性收缩与最优强度
# 1. 比较球形、对角和常相关目标；
# 2. 实现固定强度网格搜索和数据驱动强度；
# 3. 与 corpcor 比较，输出风险曲线和偏差—方差分解。
# 在项目根目录运行：
# source("R/02_linear_shrinkage.R", encoding = "UTF-8")
# results_B <- save_member_B_outputs()

#install.packages(c("ggplot2", "corpcor"))
required_packages_B <- c("ggplot2", "corpcor")
missing_packages_B <- required_packages_B[
  !vapply(required_packages_B, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages_B) > 0) {
  stop("Please install: ", paste(missing_packages_B, collapse = ", "))
}

# 一、协方差结构：参数与 A 保持一致
make_covariance_B <- function(p, structure) {
  if (structure == "identity") return(diag(p))
  if (structure == "ar1") return(0.6 ^ abs(outer(seq_len(p), seq_len(p), "-")))
  if (structure == "block") {
    groups <- rep(seq_len(ceiling(p / 10)), each = 10, length.out = p)
    Sigma <- matrix(0.1, p, p)
    Sigma[outer(groups, groups, "==")] <- 0.7
    diag(Sigma) <- 1
    return(Sigma)
  }
  if (structure == "spiked") {
    set.seed(1000 + p)
    Q <- qr.Q(qr(matrix(rnorm(p * p), p)))
    return(Q %*% diag(c(6, 3, rep(1, p - 2))) %*% t(Q))
  }
  stop("Unknown covariance structure.")
}

# 二、目标矩阵与线性收缩 
make_shrinkage_target <- function(S, target) {
  p <- ncol(S)
  v <- diag(S)
  if (target == "spherical") return(diag(mean(v), p))
  if (target == "diagonal") return(diag(v, p))
  if (target != "constant_correlation") stop("Unknown target.")
  R <- stats::cov2cor(S)
  rho <- mean(R[upper.tri(R)])
  T <- rho * sqrt(outer(v, v))
  diag(T) <- v
  T
}

linear_shrinkage <- function(S, T, lambda) {
  (1 - lambda) * S + lambda * T
}

# 三、数据驱动强度：不使用真实 Sigma 
estimate_shrinkage_intensity <- function(S, n, target) {
  # 球形、对角目标：正态模型下的矩估计，S 的分母为 n-1。
  # T 由样本估计，分子已考虑 S 与 T 的相关项。
  m <- n - 1
  den <- (m - 1) * (m + 2)
  V <- ((m - 2) * S^2 + m * outer(diag(S), diag(S))) / den
  if (target == "spherical") {
    trace_square <- (m^2 * sum(S^2) - m * sum(diag(S))^2) / den
    numerator <- sum(V) - 2 * trace_square / (m * ncol(S))
  } else {
    numerator <- sum(V) - sum(diag(V))
  }
  distance <- sum((S - make_shrinkage_target(S, target))^2)
  if (distance <= .Machine$double.eps * sum(S^2)) return(1)
  min(1, max(0, numerator / distance))
}

select_constant_correlation_intensity <- function(X, lambda_grid, folds = 5) {
  # 每折只用训练集构造目标，用验证集协方差计算损失。
  # 选择后在全样本上重估，CV 不保证得到全样本的精确最优强度。
  fold_id <- sample(rep(seq_len(folds), length.out = nrow(X)))
  scores <- matrix(0, folds, length(lambda_grid))
  for (k in seq_len(folds)) {
    S <- stats::cov(X[fold_id != k, , drop = FALSE])
    V <- stats::cov(X[fold_id == k, , drop = FALSE])
    D <- make_shrinkage_target(S, "constant_correlation") - S
    scores[k, ] <- sum((S - V)^2) + 2 * lambda_grid * sum((S - V) * D) +
      lambda_grid^2 * sum(D^2)
  }
  lambda_grid[which.min(colMeans(scores))]
}

# 四、固定目标的理论分解 
# 固定 T：R(lambda)=(1-lambda)^2*a+lambda^2*b；令导数为0得 lambda*=a/(a+b)。
# a 是 S 的方差风险，b=||T-Sigma||_F^2 是目标偏差；lambda 越大，偏差代价越大。
# 若 T=T(S)，固定 lambda 的风险还含 2*lambda*(1-lambda)*c，
# c=E<S-Sigma,T-Sigma>，最优强度为截断到[0,1]的 (a-c)/(a+b-2c)。
fixed_target_decomposition <- function(Sigma, n, targets, lambda_grid) {
  a <- (sum(diag(Sigma))^2 + sum(Sigma^2)) / (n - 1)
  do.call(rbind, lapply(targets, function(target) {
    # 用总体构造后保持固定的 T0，只用于理论演示。
    T0 <- make_shrinkage_target(Sigma, target)
    b <- sum((T0 - Sigma)^2)
    data.frame(target = target, lambda = lambda_grid, a = a, b = b,
               variance = (1 - lambda_grid)^2 * a, bias_squared = lambda_grid^2 * b,
               total_risk = (1 - lambda_grid)^2 * a + lambda_grid^2 * b,
               lambda_star = a / (a + b))
  }))
}

# 五、模拟：同一份样本比较不同方法 
run_shrinkage_simulation <- function(reps = 500, seed = 20261103) {
  regimes <- data.frame(
    regime = c("low", "moderate", "high", "near_singular", "singular"),
    n = c(500, 200, 200, 150, 100), p = c(10, 50, 150, 140, 200)
  )
  structures <- c("identity", "ar1", "block", "spiked")
  targets <- c("spherical", "diagonal", "constant_correlation")
  lambda_grid <- seq(0, 1, by = 0.01)
  risk_results <- grid_results <- theory_results <- list()
  index <- 0
  for (i in seq_len(nrow(regimes))) {
    for (structure in structures) {
      index <- index + 1
      n <- regimes$n[i]
      p <- regimes$p[i]
      Sigma <- make_covariance_B(p, structure)
      L <- chol(Sigma)
      # 在构造 Sigma 后设置抽样种子，避免尖峰构造的固定种子导致重复样本。
      set.seed(seed + index)
      message(sprintf("%s: %s, n=%d, p=%d", regimes$regime[i], structure, n, p))
      meta <- data.frame(n = n, p = p, p_over_n = p / n,
                         structure = structure, regime = regimes$regime[i])
      rows <- vector("list", reps)
      grid_loss <- array(0, c(reps, length(lambda_grid), length(targets)))
      for (r in seq_len(reps)) {
        X <- matrix(rnorm(n * p), nrow = n) %*% L
        S <- stats::cov(X)
        T <- lapply(targets, function(target) make_shrinkage_target(S, target))
        lambda <- c(estimate_shrinkage_intensity(S, n, "spherical"),
                    estimate_shrinkage_intensity(S, n, "diagonal"),
                    select_constant_correlation_intensity(X, lambda_grid))
        estimates <- lapply(seq_along(targets), function(j) linear_shrinkage(S, T[[j]], lambda[j]))
        corp <- corpcor::cov.shrink(X, verbose = FALSE)
        estimates <- c(list(S), estimates, list(corp))
        error <- vapply(estimates, function(A) sum((A - Sigma)^2), numeric(1))
        # corpcor 默认同时收缩相关与方差，两种强度分别记录。
        rows[[r]] <- data.frame(
          meta, rep = r,
          method = c("sample_covariance", "spherical_gaussian", "diagonal_gaussian",
                     "constant_correlation_cv", "corpcor_default"),
          lambda = c(0, lambda, attr(corp, "lambda")),
          lambda_variance = c(rep(0, 4), attr(corp, "lambda.var")),
          frobenius_error = error, element_mse = error / p^2,
          error_difference_vs_sample = error - error[1]
        )
        for (j in seq_along(targets)) {
          D <- T[[j]] - S
          grid_loss[r, , j] <- sum((S - Sigma)^2) +
            2 * lambda_grid * sum((S - Sigma) * D) + lambda_grid^2 * sum(D^2)
        }
      }
      risk_results[[index]] <- do.call(rbind, rows)
      grid_results[[index]] <- do.call(rbind, lapply(seq_along(targets), function(j) {
        risk <- colMeans(grid_loss[, , j])
        # 已知 Sigma 下的平均风险网格最小点，仅为模拟参考，不是独立测试成绩。
        data.frame(meta, target = targets[j], lambda = lambda_grid,
                   frobenius_risk_mean = risk,
                   risk_mc_se = apply(grid_loss[, , j], 2, stats::sd) / sqrt(reps),
                   lambda_grid_best = lambda_grid[which.min(risk)])
      }))
      theory_results[[index]] <- cbind(meta, fixed_target_decomposition(Sigma, n, targets, lambda_grid))
    }
  }
  list(risk_results = do.call(rbind, risk_results),
       grid_results = do.call(rbind, grid_results),
       bias_variance = do.call(rbind, theory_results))
}

summarise_shrinkage_results <- function(results) {
  groups <- split(results, interaction(results$regime, results$structure, results$method, drop = TRUE))
  summary <- lapply(groups, function(d) {
    risk <- mean(d$frobenius_error)
    baseline <- mean(d$frobenius_error - d$error_difference_vs_sample)
    data.frame(d[1, c("n", "p", "p_over_n", "structure", "regime", "method")],
               reps = nrow(d), lambda_mean = mean(d$lambda),
               lambda_variance_mean = mean(d$lambda_variance),
               frobenius_risk_mean = risk, frobenius_risk_sd = stats::sd(d$frobenius_error),
               element_mse_mean = mean(d$element_mse),
               risk_mc_se = stats::sd(d$frobenius_error) / sqrt(nrow(d)),
               risk_reduction_percent = 100 * (1 - risk / baseline))
  })
  summary <- do.call(rbind, summary)
  rownames(summary) <- NULL
  summary
}

# 六、作图
plot_shrinkage_risk_curve <- function(grid_results) {
  d <- subset(grid_results, regime == "high")
  baseline <- subset(d, lambda == 0 & target == "spherical")
  ggplot2::ggplot(d, ggplot2::aes(lambda, frobenius_risk_mean, color = target)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_hline(data = baseline, ggplot2::aes(yintercept = frobenius_risk_mean),
                       linetype = "dashed", color = "grey40") +
    ggplot2::facet_wrap(~ structure, scales = "free_y") +
    ggplot2::labs(title = "Frobenius risk versus shrinkage intensity",
                  subtitle = "n = 200, p = 150; dashed line: sample covariance risk",
                  x = "Shrinkage intensity (lambda)", y = "Frobenius risk", color = "Target") +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "bottom")
}

plot_shrinkage_comparison <- function(summary_results) {
  ggplot2::ggplot(summary_results, ggplot2::aes(p_over_n, frobenius_risk_mean, color = method)) +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~ structure, scales = "free_y") + ggplot2::scale_y_log10() +
    ggplot2::labs(title = "Shrinkage methods and corpcor: risk comparison",
                  x = "p / n", y = "Frobenius risk (log scale)", color = "Method") +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "bottom") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2))
}

plot_shrinkage_bias_variance <- function(bias_variance) {
  d <- subset(bias_variance, regime == "moderate" & structure == "ar1" & target == "spherical")
  plot_data <- rbind(data.frame(lambda = d$lambda, component = "Variance", risk = d$variance),
                     data.frame(lambda = d$lambda, component = "Squared bias", risk = d$bias_squared),
                     data.frame(lambda = d$lambda, component = "Total risk", risk = d$total_risk))
  ggplot2::ggplot(plot_data, ggplot2::aes(lambda, risk, color = component)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = d$lambda_star[1], linetype = "dashed") +
    ggplot2::labs(title = "Bias-variance tradeoff for a fixed spherical target",
                  subtitle = sprintf("AR(1), n = 200, p = 50; lambda* = %.3f", d$lambda_star[1]),
                  x = "Shrinkage intensity (lambda)", y = "Risk component", color = "Component") +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(legend.position = "bottom")
}

# 七、保存入口
save_member_B_outputs <- function(output_table_dir = "tables/member_B",
                                  output_figure_dir = "figures/member_B",
                                  reps = 500, seed = 20261103) {
  stopifnot(reps >= 2)
  dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_figure_dir, recursive = TRUE, showWarnings = FALSE)
  results <- run_shrinkage_simulation(reps, seed)
  summary <- summarise_shrinkage_results(results$risk_results)
  tables <- list(shrinkage_risk_simulations = results$risk_results,
                 shrinkage_risk_summary = summary,
                 shrinkage_intensity_grid = results$grid_results,
                 shrinkage_bias_variance = results$bias_variance)
  for (name in names(tables)) {
    utils::write.csv(tables[[name]], file.path(output_table_dir, paste0(name, ".csv")),
                     row.names = FALSE, fileEncoding = "UTF-8")
  }
  plots <- list(frobenius_risk_vs_shrinkage_intensity = plot_shrinkage_risk_curve(results$grid_results),
                frobenius_risk_vs_p_over_n = plot_shrinkage_comparison(summary),
                shrinkage_bias_variance_tradeoff = plot_shrinkage_bias_variance(results$bias_variance))
  for (name in names(plots)) {
    filename <- file.path(output_figure_dir, paste0(name, ".png"))
    ggplot2::ggsave(filename, plots[[name]], width = 10, height = 7, dpi = 160, bg = "white")
    if (file.info(filename)$size == 0) stop("Image export failed: ", filename)
  }
  invisible(tables)
}
