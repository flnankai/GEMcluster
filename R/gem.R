.elliptical_family_levels <- c("gaussian", "t", "laplace", "slash")

.rcpp_fun_exists <- function(name) {
  exists(name, mode = "function", inherits = TRUE)
}

weighted_mean <- function(x, w) {
  sum(w * x) / sum(w)
}

weighted_var <- function(x, w) {
  mu <- weighted_mean(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

trapz <- function(x, y) {
  idx <- 2:length(x)
  sum((x[idx] - x[idx - 1L]) * (y[idx] + y[idx - 1L]) / 2)
}

project_pd <- function(M, eps = 1e-6) {
  if (.rcpp_fun_exists("project_pd_cpp")) {
    return(project_pd_cpp(M, eps))
  }
  M <- (M + t(M)) / 2
  eig <- eigen(M, symmetric = TRUE)
  eig$values[eig$values < eps] <- eps
  out <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  (out + t(out)) / 2
}

normalize_trace <- function(Sigma) {
  p <- ncol(Sigma)
  Sigma <- Sigma * (p / sum(diag(Sigma)))
  (Sigma + t(Sigma)) / 2
}

#' AR(1) covariance matrix
#'
#' @param p Dimension.
#' @param rho AR correlation.
#' @export
make_ar_cov <- function(p, rho = 0.5) {
  outer(seq_len(p), seq_len(p), function(i, j) rho^abs(i - j))
}

#' Dense block-constant mean matrix
#'
#' @param K Number of classes. The current template supports K = 3.
#' @param p Dimension.
#' @param signal Signal size multiplying the block pattern.
#' @export
dense_block_means <- function(K = 3, p = 200, signal = 0.5) {
  if (K != 3L) {
    stop("dense_block_means() currently supports K = 3.")
  }
  block_id <- cut(seq_len(p), breaks = 4L, labels = FALSE, include.lowest = TRUE)
  block_levels <- rbind(
    c(1.5, 0.5, -0.5, -1.5),
    c(-0.5, 1.5, 0.5, -1.5),
    c(0.5, -1.5, 1.5, -0.5)
  )
  Mu <- matrix(0, nrow = K, ncol = p)
  for (k in seq_len(K)) {
    Mu[k, ] <- signal * block_levels[k, block_id]
  }
  Mu
}

all_permutations <- function(x) {
  if (length(x) <= 1L) {
    return(list(x))
  }
  out <- vector("list", 0L)
  for (i in seq_along(x)) {
    rest <- all_permutations(x[-i])
    out <- c(out, lapply(rest, function(r) c(x[i], r)))
  }
  out
}

#' Clustering accuracy up to label permutation
#'
#' @param pred Predicted labels.
#' @param truth True labels.
#' @export
cluster_accuracy <- function(pred, truth) {
  pred <- as.integer(pred)
  truth <- as.integer(truth)
  pred_levels <- sort(unique(pred))
  true_levels <- sort(unique(truth))

  if (length(pred_levels) != length(true_levels) || length(pred_levels) > 8L) {
    mapped <- integer(length(pred))
    for (lev in pred_levels) {
      idx <- which(pred == lev)
      tab <- table(truth[idx])
      mapped[idx] <- as.integer(names(tab)[which.max(tab)])
    }
    return(list(accuracy = mean(mapped == truth), mapped = mapped))
  }

  perms <- all_permutations(true_levels)
  best_acc <- -Inf
  best_map <- pred
  for (perm in perms) {
    map <- setNames(perm, pred_levels)
    mapped <- as.integer(map[as.character(pred)])
    acc <- mean(mapped == truth)
    if (acc > best_acc) {
      best_acc <- acc
      best_map <- mapped
    }
  }
  list(accuracy = best_acc, mapped = best_map)
}

#' Relative Frobenius error
#' @export
relative_frobenius <- function(A, B) {
  norm(A - B, type = "F") / norm(B, type = "F")
}

softmax_rows <- function(log_mat) {
  if (.rcpp_fun_exists("softmax_rows_cpp")) {
    return(softmax_rows_cpp(log_mat))
  }
  row_max <- apply(log_mat, 1, max)
  x <- exp(log_mat - row_max)
  x / rowSums(x)
}

row_logsumexp <- function(log_mat) {
  if (.rcpp_fun_exists("row_logsumexp_cpp")) {
    return(row_logsumexp_cpp(log_mat))
  }
  row_max <- apply(log_mat, 1, max)
  row_max + log(rowSums(exp(log_mat - row_max)))
}

l1_distance_matrix <- function(X, centers, use_dims = NULL) {
  if (is.null(use_dims)) {
    use_dims <- integer(0)
  } else {
    use_dims <- as.integer(use_dims)
  }
  if (.rcpp_fun_exists("l1_distance_matrix_cpp")) {
    return(l1_distance_matrix_cpp(X, centers, use_dims))
  }
  if (length(use_dims) == 0L) {
    use_dims <- seq_len(ncol(X))
  }
  sapply(seq_len(nrow(centers)), function(k) {
    rowSums(abs(sweep(X[, use_dims, drop = FALSE], 2, centers[k, use_dims], "-")))
  })
}

calc_delta <- function(X, Mu, Omega, eps = 1e-8) {
  if (.rcpp_fun_exists("calc_delta_cpp")) {
    return(calc_delta_cpp(X, Mu, Omega, eps))
  }
  n <- nrow(X)
  K <- nrow(Mu)
  Delta <- matrix(0, nrow = n, ncol = K)
  for (g in seq_len(K)) {
    Rg <- sweep(X, 2, Mu[g, ], "-")
    Delta[, g] <- rowSums((Rg %*% Omega) * Rg)
  }
  pmax(Delta, eps)
}

col_medians <- function(X) {
  if (.rcpp_fun_exists("col_medians_cpp")) {
    return(as.numeric(col_medians_cpp(X)))
  }
  apply(X, 2, stats::median)
}

cluster_centers_median <- function(X, cluster, K) {
  if (.rcpp_fun_exists("cluster_centers_median_cpp")) {
    centers <- cluster_centers_median_cpp(X, as.integer(cluster), as.integer(K))
    centers <- as.matrix(centers)
    empty <- which(!stats::complete.cases(centers))
    if (length(empty) == 0L) {
      return(centers)
    }
  }
  p <- ncol(X)
  centers <- matrix(0, nrow = K, ncol = p)
  for (k in seq_len(K)) {
    idx <- which(cluster == k)
    if (length(idx) == 0L) {
      centers[k, ] <- X[sample.int(nrow(X), 1L), ]
    } else {
      centers[k, ] <- col_medians(X[idx, , drop = FALSE])
    }
  }
  centers
}

#' Simulate an elliptical mixture sample
#'
#' @param n Total sample size.
#' @param Mu K by p mean matrix.
#' @param Sigma Common p by p scatter/covariance matrix.
#' @param family One of gaussian, t, laplace, slash.
#' @param t_df Degrees of freedom for t.
#' @param slash_df Tail parameter for slash; must be greater than 2.
#' @param pi Optional class proportions.
#' @param seed Optional random seed.
#' @export
simulate_elliptical_mixture <- function(n, Mu, Sigma,
                                        family = .elliptical_family_levels,
                                        t_df = 5, slash_df = 4,
                                        pi = NULL, seed = NULL) {
  family <- match.arg(family)
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required for simulation.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  K <- nrow(Mu)
  p <- ncol(Mu)
  if (is.null(pi)) {
    counts <- rep(floor(n / K), K)
    if (n %% K > 0L) {
      counts[seq_len(n %% K)] <- counts[seq_len(n %% K)] + 1L
    }
  } else {
    counts <- as.vector(stats::rmultinom(1, size = n, prob = pi))
  }
  z <- rep(seq_len(K), counts)
  X <- matrix(NA_real_, nrow = length(z), ncol = p)
  start <- 1L
  for (k in seq_len(K)) {
    nk <- counts[k]
    if (nk == 0L) {
      next
    }
    idx <- start:(start + nk - 1L)
    X[idx, ] <- simulate_elliptical_block(
      nk, Mu[k, ], Sigma, family = family, t_df = t_df, slash_df = slash_df
    )
    start <- start + nk
  }
  ord <- sample.int(length(z))
  list(X = X[ord, , drop = FALSE], label = z[ord], family = family)
}

simulate_elliptical_block <- function(n, mu, Sigma, family, t_df = 5, slash_df = 4) {
  p <- length(mu)
  if (family == "gaussian") {
    return(mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma))
  }
  if (family == "t") {
    if (t_df <= 2) {
      stop("t_df must be > 2 so that covariance is standardized.")
    }
    Sigma_scaled <- ((t_df - 2) / t_df) * Sigma
    return(mvtnorm::rmvt(n, sigma = Sigma_scaled, df = t_df, delta = mu, type = "shifted"))
  }
  Z <- mvtnorm::rmvnorm(n, sigma = Sigma)
  if (family == "laplace") {
    scale_vec <- sqrt(stats::rexp(n, rate = 1))
  } else if (family == "slash") {
    if (slash_df <= 2) {
      stop("slash_df must be > 2 so that covariance is standardized.")
    }
    u <- stats::runif(n)
    scale_vec <- sqrt((slash_df - 2) / slash_df) * u^(-1 / slash_df)
  } else {
    stop("Unsupported family.")
  }
  sweep(Z, 1, scale_vec, "*") + matrix(mu, nrow = n, ncol = p, byrow = TRUE)
}

kmedian_single <- function(X, K, maxit = 60, cluster_init = NULL) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(cluster_init)) {
    centers <- X[sample.int(n, K, replace = FALSE), , drop = FALSE]
    cluster <- rep(NA_integer_, n)
  } else {
    cluster <- as.integer(cluster_init)
    centers <- cluster_centers_median(X, cluster, K)
  }
  objective <- Inf
  for (iter in seq_len(maxit)) {
    dist_mat <- l1_distance_matrix(X, centers)
    cluster_new <- max.col(-dist_mat, ties.method = "first")
    for (k in seq_len(K)) {
      if (!any(cluster_new == k)) {
        cluster_new[sample.int(n, 1L)] <- k
      }
    }
    centers_new <- cluster_centers_median(X, cluster_new, K)
    objective_new <- sum(dist_mat[cbind(seq_len(n), cluster_new)])
    if (!any(is.na(cluster)) && all(cluster_new == cluster)) {
      cluster <- cluster_new
      centers <- centers_new
      objective <- objective_new
      break
    }
    cluster <- cluster_new
    centers <- centers_new
    objective <- objective_new
  }
  list(cluster = cluster, centers = centers, objective = objective)
}

kmedian_multistart <- function(X, K, nstart = 12L, maxit = 60) {
  best_fit <- NULL
  for (s in seq_len(nstart)) {
    fit <- kmedian_single(X, K = K, maxit = maxit)
    if (is.null(best_fit) || fit$objective < best_fit$objective) {
      best_fit <- fit
    }
  }
  best_fit
}

permute_columns <- function(X) {
  out <- X
  for (j in seq_len(ncol(X))) {
    out[, j] <- sample(X[, j], size = nrow(X), replace = FALSE)
  }
  out
}

compute_between_ss_l1 <- function(X, cluster) {
  K <- length(unique(cluster))
  if (.rcpp_fun_exists("compute_between_ss_l1_cpp")) {
    return(compute_between_ss_l1_cpp(X, as.integer(cluster), as.integer(K)))
  }
  overall_med <- col_medians(X)
  val <- 0
  for (k in sort(unique(cluster))) {
    idx <- which(cluster == k)
    if (length(idx) == 0L) {
      next
    }
    med_k <- col_medians(X[idx, , drop = FALSE])
    val <- val + length(idx) * sum(abs(med_k - overall_med))
  }
  val
}

make_sparse_kmedian_tau_grid <- function(X, K, nstart = 4L, max_iter = 50L) {
  init <- kmedian_multistart(X, K = K, nstart = nstart, maxit = max_iter)
  center_mean <- colMeans(init$centers)
  disp <- colSums(abs(init$centers - matrix(center_mean, K, ncol(X), byrow = TRUE)))
  tau_grid <- unique(as.numeric(stats::quantile(disp, probs = seq(0.2, 0.8, length.out = 6), names = FALSE)))
  tau_grid <- tau_grid[is.finite(tau_grid) & tau_grid > 0]
  if (length(tau_grid) == 0L) {
    tau_grid <- mean(disp)
  }
  sort(unique(tau_grid))
}

sparse_kmedian_single <- function(X, K, tau, max_iter = 50, cluster_init = NULL) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  if (is.null(cluster_init)) {
    cluster <- sample.int(K, size = n, replace = TRUE)
  } else {
    cluster <- as.integer(cluster_init)
  }
  centers <- matrix(0, nrow = K, ncol = p)
  use_dims <- seq_len(p)
  objective <- Inf
  for (iter in seq_len(max_iter)) {
    centers <- cluster_centers_median(X, cluster, K)
    center_mean <- colMeans(centers)
    diff_sum <- colSums(abs(centers - matrix(center_mean, K, p, byrow = TRUE)))
    use_dims <- which(diff_sum >= tau)
    if (length(use_dims) == 0L) {
      use_dims <- seq_len(p)
    }
    dist_mat <- l1_distance_matrix(X, centers, use_dims = use_dims)
    cluster_new <- max.col(-dist_mat, ties.method = "first")
    for (k in seq_len(K)) {
      if (!any(cluster_new == k)) {
        cluster_new[sample.int(n, 1L)] <- k
      }
    }
    objective_new <- sum(dist_mat[cbind(seq_len(n), cluster_new)])
    if (all(cluster_new == cluster)) {
      cluster <- cluster_new
      objective <- objective_new
      break
    }
    cluster <- cluster_new
    objective <- objective_new
  }
  list(cluster = cluster, centers = centers, selected_dims = use_dims, objective = objective)
}

sparse_kmedian_multistart <- function(X, K, tau, nstart = 8L, max_iter = 50) {
  best_fit <- NULL
  for (s in seq_len(nstart)) {
    fit <- sparse_kmedian_single(X, K = K, tau = tau, max_iter = max_iter)
    if (is.null(best_fit) || fit$objective < best_fit$objective) {
      best_fit <- fit
    }
  }
  best_fit
}

select_tau_sparse_kmedian <- function(X, K, tau_grid, B = 5L, nstart = 5L) {
  gap <- rep(-Inf, length(tau_grid))
  for (i in seq_along(tau_grid)) {
    tau <- tau_grid[i]
    fit <- sparse_kmedian_multistart(X, K = K, tau = tau, nstart = nstart)
    stat_obs <- compute_between_ss_l1(X[, fit$selected_dims, drop = FALSE], fit$cluster)
    boot_stats <- numeric(B)
    for (b in seq_len(B)) {
      Xb <- permute_columns(X)
      fit_b <- sparse_kmedian_multistart(
        Xb, K = K, tau = tau, nstart = max(2L, ceiling(nstart / 2))
      )
      boot_stats[b] <- compute_between_ss_l1(Xb[, fit_b$selected_dims, drop = FALSE], fit_b$cluster)
    }
    gap[i] <- log(pmax(stat_obs, 1e-8)) - mean(log(pmax(boot_stats, 1e-8)))
  }
  list(best_tau = tau_grid[which.max(gap)], gap = gap, tau_grid = tau_grid)
}

weighted_kde_eval <- function(y_eval, y_obs, weights, bandwidth) {
  if (.rcpp_fun_exists("weighted_kde_eval_cpp")) {
    return(as.numeric(weighted_kde_eval_cpp(y_eval, y_obs, weights, bandwidth)))
  }
  z <- outer(y_eval, y_obs, "-") / bandwidth
  drop(stats::dnorm(z) %*% weights) / bandwidth
}

build_generator <- function(delta, tau, p, bandwidth = NULL, grid_size = 250L,
                            eps_u = 1e-4, eps_density = 1e-10,
                            omega_clip = c(1e-3, 50)) {
  y_obs <- log1p(as.vector(delta))
  weights <- as.vector(tau)
  weights <- weights / sum(weights)
  n_eff <- (sum(weights)^2) / sum(weights^2)
  if (is.null(bandwidth)) {
    sd_y <- sqrt(max(weighted_var(y_obs, weights), 1e-8))
    bandwidth <- 1.06 * sd_y * n_eff^(-1 / 5)
  }
  bandwidth <- max(bandwidth, 0.05)
  y_min <- max(min(y_obs) - 3 * bandwidth, log1p(eps_u))
  y_max <- max(y_obs) + 3 * bandwidth
  if (y_max <= y_min) {
    y_max <- y_min + 1
  }
  y_grid <- seq(y_min, y_max, length.out = grid_size)
  u_grid <- pmax(expm1(y_grid), eps_u)
  if (.rcpp_fun_exists("weighted_kde_eval_grid_cpp")) {
    fy <- as.numeric(weighted_kde_eval_grid_cpp(y_grid[1L], y_grid[2L] - y_grid[1L], grid_size, y_obs, weights, bandwidth))
  } else {
    fy <- weighted_kde_eval(y_grid, y_obs, weights, bandwidth)
  }
  log_q <- log(pmax(fy, eps_density)) - log1p(u_grid)
  log_tilde_g <- (1 - p / 2) * log(u_grid) + log_q
  normalizer <- trapz(u_grid, exp(log_q))
  log_g_raw <- log_tilde_g - log(normalizer)
  spline_fit <- try(stats::smooth.spline(y_grid, log_g_raw, spar = 0.55), silent = TRUE)
  if (inherits(spline_fit, "try-error")) {
    log_g_grid <- log_g_raw
    slope <- diff(log_g_grid) / diff(y_grid)
    dlogg_dy <- c(slope, tail(slope, 1L))
  } else {
    log_g_grid <- stats::predict(spline_fit, x = y_grid, deriv = 0)$y
    dlogg_dy <- stats::predict(spline_fit, x = y_grid, deriv = 1)$y
  }
  score_grid <- -dlogg_dy / (1 + u_grid)
  score_grid <- pmin(pmax(score_grid, omega_clip[1]), omega_clip[2])
  list(
    log_g = function(u) {
      y <- log1p(pmax(u, eps_u))
      stats::approx(y_grid, log_g_grid, xout = y, rule = 2)$y
    },
    omega = function(u) {
      y <- log1p(pmax(u, eps_u))
      vals <- stats::approx(y_grid, score_grid, xout = y, rule = 2)$y
      pmin(pmax(vals, omega_clip[1]), omega_clip[2])
    },
    bandwidth = bandwidth,
    n_eff = n_eff,
    y_grid = y_grid
  )
}

soft_threshold_offdiag <- function(A, lambda) {
  out <- A
  off_idx <- row(out) != col(out)
  out[off_idx] <- sign(out[off_idx]) * pmax(abs(out[off_idx]) - lambda, 0)
  out
}

select_factor_count_gr <- function(A, M = 8L, eps = 1e-8) {
  A <- project_pd(A)
  vals <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  vals <- sort(pmax(vals, eps), decreasing = TRUE)
  q <- length(vals)
  if (q < 4L) {
    return(0L)
  }
  M_eff <- min(as.integer(M), q - 2L)
  if (M_eff < 1L) {
    return(0L)
  }
  q1 <- q - 1L
  V <- numeric(q1)
  for (j in 0:(q1 - 1L)) {
    if (j == 0L) {
      V[j + 1L] <- sum(vals[1:q1])
    } else {
      V[j + 1L] <- sum(vals[(j + 1L):q1])
    }
  }
  gr <- rep(-Inf, M_eff)
  for (j in seq_len(M_eff)) {
    num <- log1p(vals[j] / max(V[j], eps))
    den <- log1p(vals[j + 1L] / max(V[j + 1L], eps))
    gr[j] <- num / max(den, eps)
  }
  as.integer(which.max(gr))
}

poet_map <- function(A, lambda_u, m = 0L) {
  if (.rcpp_fun_exists("poet_map_cpp")) {
    return(poet_map_cpp(A, lambda_u, as.integer(m)))
  }
  A <- project_pd(A)
  p <- ncol(A)
  eig <- eigen(A, symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)
  vals <- pmax(eig$values[ord], 1e-8)
  vecs <- eig$vectors[, ord, drop = FALSE]
  if (m > 0L) {
    spike_vals <- vals[seq_len(m)]
    spike_vecs <- vecs[, seq_len(m), drop = FALSE]
    spike <- spike_vecs %*% diag(spike_vals, nrow = m) %*% t(spike_vecs)
  } else {
    spike <- matrix(0, nrow = p, ncol = p)
  }
  remainder <- A - spike
  remainder_thr <- soft_threshold_offdiag(remainder, lambda_u)
  diag(remainder_thr) <- diag(remainder)
  project_pd(spike + remainder_thr)
}

weighted_spatial_sign_scatter <- function(residuals, weights, eps_r = 1e-6) {
  if (.rcpp_fun_exists("weighted_spatial_sign_scatter_cpp")) {
    return(weighted_spatial_sign_scatter_cpp(residuals, weights, eps_r))
  }
  norms2 <- rowSums(residuals^2)
  coeff <- weights / pmax(norms2, eps_r)
  crossprod(sweep(residuals, 1, sqrt(coeff), "*")) / sum(weights)
}

weighted_tyler <- function(residuals, weights, Sigma_init, maxit = 40L,
                           tol = 1e-4, ridge = 0.02, eps_r = 1e-6) {
  p <- ncol(residuals)
  Sigma <- normalize_trace(project_pd(Sigma_init))
  weight_sum <- sum(weights)
  for (iter in seq_len(maxit)) {
    Omega <- solve(Sigma)
    if (.rcpp_fun_exists("quadform_rows_cpp")) {
      quad <- as.numeric(quadform_rows_cpp(residuals, Omega))
    } else {
      quad <- rowSums((residuals %*% Omega) * residuals)
    }
    coeff <- p * weights / pmax(quad, eps_r) / weight_sum
    if (.rcpp_fun_exists("weighted_crossprod_cpp")) {
      Sigma_new <- weighted_crossprod_cpp(residuals, coeff)
    } else {
      Sigma_new <- crossprod(sweep(residuals, 1, sqrt(coeff), "*"))
    }
    Sigma_new <- normalize_trace((1 - ridge) * Sigma_new + ridge * diag(p))
    Sigma_new <- project_pd(Sigma_new)
    rel_change <- norm(Sigma_new - Sigma, type = "F") / max(1, norm(Sigma, type = "F"))
    Sigma <- Sigma_new
    if (rel_change < tol) {
      break
    }
  }
  Sigma
}

select_glasso_cov <- function(Sigma_pt, lambda_base, n_eff,
                              ebic_gamma = 0.5,
                              lambda_grid = NULL) {
  if (!requireNamespace("huge", quietly = TRUE)) {
    stop("Package 'huge' is required for GEM covariance estimation.")
  }
  p <- ncol(Sigma_pt)
  Sigma_pt <- normalize_trace(project_pd(Sigma_pt))
  if (is.null(lambda_grid)) {
    lambda_grid <- lambda_base * exp(seq(log(1.8), log(0.45), length.out = 8))
  }
  lambda_grid <- sort(unique(lambda_grid), decreasing = TRUE)
  fit <- try(
    huge::huge.glasso(Sigma_pt, lambda = lambda_grid, cov.output = TRUE, verbose = FALSE),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    Sigma_fallback <- normalize_trace(project_pd(Sigma_pt + lambda_base * diag(p)))
    Omega_fallback <- project_pd(solve(Sigma_fallback))
    return(list(Omega = Omega_fallback, Sigma = Sigma_fallback, lambda = lambda_base, ebic = NA_real_))
  }
  ebic <- -n_eff * fit$loglik + log(n_eff) * fit$df + 4 * ebic_gamma * log(p) * fit$df
  best_idx <- which.min(ebic)
  Sigma_best <- normalize_trace(project_pd(fit$cov[[best_idx]]))
  Omega_best <- project_pd(solve(Sigma_best))
  list(Omega = Omega_best, Sigma = Sigma_best, lambda = lambda_grid[best_idx], ebic = ebic[best_idx])
}

gem_fit_single <- function(X, K,
                           max_outer = 25L,
                           tol = 2e-3,
                           init_tau = NULL,
                           init_n_boot = 5L,
                           init_nstart = 6L,
                           init_max_iter = 50L,
                           eta_mu = 0.7,
                           eta_omega = 0.7,
                           lambda_u = NULL,
                           lambda_omega = NULL,
                           lambda_u_scale = 0.55,
                           lambda_omega_scale = 0.28,
                           m_gr_upper = 8L,
                           bandwidth = NULL,
                           verbose = FALSE) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  pick_m <- function(A) select_factor_count_gr(A, M = m_gr_upper)
  if (is.null(init_tau) || !is.finite(init_tau) || init_tau <= 0) {
    tau_grid <- make_sparse_kmedian_tau_grid(
      X, K = K, nstart = max(4L, ceiling(init_nstart / 2)), max_iter = init_max_iter
    )
    tau_sel <- select_tau_sparse_kmedian(
      X, K = K, tau_grid = tau_grid, B = init_n_boot, nstart = max(3L, ceiling(init_nstart / 2))
    )
    init_tau <- tau_sel$best_tau
  }
  init <- sparse_kmedian_multistart(
    X, K = K, tau = init_tau, nstart = init_nstart, max_iter = init_max_iter
  )
  Mu <- init$centers
  z0 <- init$cluster
  Tau <- matrix(0, nrow = n, ncol = K)
  Tau[cbind(seq_len(n), z0)] <- 1
  pi_vec <- colMeans(Tau)
  hard_residuals <- X - Mu[z0, , drop = FALSE]
  n_eff0 <- n
  lambda_u_now <- if (is.null(lambda_u)) lambda_u_scale * sqrt(log(p) / n_eff0) else lambda_u
  ss0 <- weighted_spatial_sign_scatter(hard_residuals, rep(1, n))
  m0 <- pick_m(ss0)
  m_history <- integer(max_outer + 1L)
  m_history[1L] <- m0
  Sigma_pss0 <- poet_map(ss0, lambda_u_now, m = m0)
  Sigma_ty0 <- weighted_tyler(hard_residuals, rep(1, n), Sigma_pss0)
  Sigma_pt0 <- poet_map(Sigma_ty0, lambda_u_now, m = m0)
  lambda_omega_now <- if (is.null(lambda_omega)) lambda_omega_scale * sqrt(log(p) / n_eff0) else lambda_omega
  glass0 <- select_glasso_cov(Sigma_pt0, lambda_base = lambda_omega_now, n_eff = n_eff0)
  Sigma <- glass0$Sigma
  Omega <- glass0$Omega
  Delta0 <- calc_delta(X, Mu, Omega)
  generator <- build_generator(Delta0, Tau, p = p, bandwidth = bandwidth)

  for (iter in seq_len(max_outer)) {
    Delta <- calc_delta(X, Mu, Omega)
    log_g_mat <- matrix(generator$log_g(as.vector(Delta)), nrow = n, ncol = K)
    log_tau <- sweep(log_g_mat, 2, log(pmax(pi_vec, 1e-10)), "+")
    Tau_new <- softmax_rows(log_tau)
    pi_new <- pmax(colMeans(Tau_new), 1e-8)
    pi_new <- pi_new / sum(pi_new)
    generator_new <- build_generator(Delta, Tau_new, p = p, bandwidth = bandwidth)
    score_mat <- matrix(generator_new$omega(as.vector(Delta)), nrow = n, ncol = K)
    Mu_prop <- Mu
    for (g in seq_len(K)) {
      wg <- Tau_new[, g] * score_mat[, g]
      if (sum(wg) > 1e-8) {
        Mu_prop[g, ] <- colSums(sweep(X, 1, wg, "*")) / sum(wg)
      }
    }
    Mu_new <- (1 - eta_mu) * Mu + eta_mu * Mu_prop
    residual_blocks <- lapply(seq_len(K), function(g) sweep(X, 2, Mu_new[g, ], "-"))
    residuals <- do.call(rbind, residual_blocks)
    weights <- as.vector(Tau_new)
    n_eff <- (sum(weights)^2) / sum(weights^2)
    lambda_u_now <- if (is.null(lambda_u)) lambda_u_scale * sqrt(log(p) / n_eff) else lambda_u
    ss <- weighted_spatial_sign_scatter(residuals, weights)
    m_now <- pick_m(ss)
    m_history[iter + 1L] <- m_now
    Sigma_pss <- poet_map(ss, lambda_u_now, m = m_now)
    Sigma_ty <- weighted_tyler(residuals, weights, Sigma_pss)
    Sigma_pt <- poet_map(Sigma_ty, lambda_u_now, m = m_now)
    lambda_omega_now <- if (is.null(lambda_omega)) lambda_omega_scale * sqrt(log(p) / n_eff) else lambda_omega
    glass_fit <- select_glasso_cov(Sigma_pt, lambda_base = lambda_omega_now, n_eff = n_eff)
    Omega_prop <- glass_fit$Omega
    Omega_new <- project_pd((1 - eta_omega) * Omega + eta_omega * Omega_prop)
    Sigma_new <- normalize_trace(project_pd(solve(Omega_new)))
    Omega_new <- project_pd(solve(Sigma_new))
    delta_mu <- max(abs(Mu_new - Mu))
    delta_omega <- norm(Omega_new - Omega, type = "F") / max(1, norm(Omega, type = "F"))
    delta_pi <- max(abs(pi_new - pi_vec))
    Mu <- Mu_new
    Sigma <- Sigma_new
    Omega <- Omega_new
    Tau <- Tau_new
    pi_vec <- pi_new
    generator <- generator_new
    if (isTRUE(verbose)) {
      cat(sprintf(
        "GEM iter=%d, m=%d, dMu=%.3e, dOmega=%.3e, dPi=%.3e\n",
        iter, m_now, delta_mu, delta_omega, delta_pi
      ))
    }
    if (max(delta_mu, delta_omega, delta_pi) < tol) {
      break
    }
  }

  Delta_final <- calc_delta(X, Mu, Omega)
  generator_final <- build_generator(Delta_final, Tau, p = p, bandwidth = bandwidth)
  log_g_final <- matrix(generator_final$log_g(as.vector(Delta_final)), nrow = n, ncol = K)
  log_tau_final <- sweep(log_g_final, 2, log(pmax(pi_vec, 1e-10)), "+")
  Tau_final <- softmax_rows(log_tau_final)
  cluster <- max.col(Tau_final, ties.method = "first")
  objective <- sum(row_logsumexp(log_tau_final))
  structure(list(
    cluster = cluster,
    Tau = Tau_final,
    pi = pi_vec,
    Mu = Mu,
    Sigma = Sigma,
    Omega = Omega,
    generator = generator_final,
    init = init,
    init_tau = init_tau,
    init_selected_dims = init$selected_dims,
    iterations = iter,
    m_initial = m0,
    m_final = m_history[iter + 1L],
    m_history = m_history[seq_len(iter + 1L)],
    objective = objective
  ), class = c("gem_fit", "list"))
}

#' Fit GEM clustering
#'
#' @param X n by p data matrix.
#' @param K Number of clusters.
#' @param max_outer Maximum number of GEM outer iterations.
#' @param tol Convergence tolerance.
#' @param init_n_boot Number of permutation bootstraps for sparse K-median tau selection.
#' @param init_nstart Number of sparse K-median starts.
#' @param init_max_iter Maximum sparse K-median iterations.
#' @param outer_nstart Number of outer GEM starts; the largest final objective is retained.
#' @param eta_mu Mean update damping.
#' @param eta_omega Precision update damping.
#' @param lambda_u Optional POET threshold.
#' @param lambda_omega Optional graphical-lasso threshold base.
#' @param lambda_u_scale Default POET threshold scale.
#' @param lambda_omega_scale Default graphical-lasso threshold scale.
#' @param m_gr_upper Maximum factor count in the GR factor-number estimator.
#' @param bandwidth Optional kernel bandwidth for the nonparametric generator.
#' @param verbose Print progress.
#' @export
gem_fit <- function(X, K,
                    max_outer = 25L,
                    tol = 2e-3,
                    init_n_boot = 5L,
                    init_nstart = 6L,
                    init_max_iter = 50L,
                    outer_nstart = 3L,
                    eta_mu = 0.7,
                    eta_omega = 0.7,
                    lambda_u = NULL,
                    lambda_omega = NULL,
                    lambda_u_scale = 0.55,
                    lambda_omega_scale = 0.28,
                    m_gr_upper = 8L,
                    bandwidth = NULL,
                    verbose = FALSE) {
  X <- as.matrix(X)
  outer_nstart <- max(1L, as.integer(outer_nstart))
  init_tau_grid <- make_sparse_kmedian_tau_grid(
    X, K = K, nstart = max(4L, ceiling(init_nstart / 2)), max_iter = init_max_iter
  )
  tau_sel <- select_tau_sparse_kmedian(
    X, K = K, tau_grid = init_tau_grid, B = init_n_boot, nstart = max(3L, ceiling(init_nstart / 2))
  )
  init_tau <- tau_sel$best_tau
  init_tau_gap <- tau_sel$gap
  best_fit <- NULL
  start_summaries <- vector("list", outer_nstart)
  for (s in seq_len(outer_nstart)) {
    fit_s <- gem_fit_single(
      X = X, K = K, max_outer = max_outer, tol = tol, init_tau = init_tau,
      init_n_boot = init_n_boot, init_nstart = init_nstart, init_max_iter = init_max_iter,
      eta_mu = eta_mu, eta_omega = eta_omega, lambda_u = lambda_u,
      lambda_omega = lambda_omega, lambda_u_scale = lambda_u_scale,
      lambda_omega_scale = lambda_omega_scale, m_gr_upper = m_gr_upper,
      bandwidth = bandwidth, verbose = FALSE
    )
    start_summaries[[s]] <- data.frame(
      start = s,
      init_tau = fit_s$init_tau,
      objective = fit_s$objective,
      iterations = fit_s$iterations,
      m_final = fit_s$m_final
    )
    if (isTRUE(verbose)) {
      cat(sprintf(
        "GEM outer start=%d, objective=%.4f, iter=%d, m_final=%d\n",
        s, fit_s$objective, fit_s$iterations, fit_s$m_final
      ))
    }
    if (is.null(best_fit) || fit_s$objective > best_fit$objective) {
      best_fit <- fit_s
      best_fit$best_start <- s
    }
  }
  best_fit$outer_nstart <- outer_nstart
  best_fit$init_tau <- init_tau
  best_fit$init_tau_grid <- init_tau_grid
  best_fit$init_tau_gap <- init_tau_gap
  best_fit$start_summaries <- do.call(rbind, start_summaries)
  class(best_fit) <- c("gem_fit", "list")
  best_fit
}

#' Predict clusters from a fitted GEM model
#'
#' @param object A fitted object returned by [gem_fit()].
#' @param newdata Optional new data matrix. If NULL, returns fitted labels or posterior probabilities.
#' @param type Either "class" or "posterior".
#' @param ... Ignored.
#' @export
predict.gem_fit <- function(object, newdata = NULL, type = c("class", "posterior"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    Tau <- object$Tau
  } else {
    X <- as.matrix(newdata)
    Delta <- calc_delta(X, object$Mu, object$Omega)
    log_g <- matrix(object$generator$log_g(as.vector(Delta)), nrow = nrow(X), ncol = nrow(object$Mu))
    log_tau <- sweep(log_g, 2, log(pmax(object$pi, 1e-10)), "+")
    Tau <- softmax_rows(log_tau)
  }
  if (type == "posterior") {
    return(Tau)
  }
  max.col(Tau, ties.method = "first")
}

#' Wrapper for GEM prediction
#' @export
gem_predict <- function(object, newdata = NULL, type = c("class", "posterior"), ...) {
  stats::predict(object, newdata = newdata, type = match.arg(type), ...)
}
