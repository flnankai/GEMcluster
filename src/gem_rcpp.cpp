// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

constexpr double INV_SQRT_2PI = 0.39894228040143267794;

arma::mat symmetrize_arma(const arma::mat& M) {
  return 0.5 * (M + M.t());
}

arma::mat project_pd_impl(const arma::mat& M, const double eps) {
  arma::mat Ms = symmetrize_arma(M);
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Ms);
  eigval.transform([&](double x) { return x < eps ? eps : x; });
  arma::mat out = eigvec * arma::diagmat(eigval) * eigvec.t();
  return symmetrize_arma(out);
}

int select_factor_count_from_desc(
    const arma::vec& vals_desc,
    const int M,
    const double eps) {
  const int q = static_cast<int>(vals_desc.n_elem);
  if (q < 4) {
    return 0;
  }

  const int M_eff = std::min(M, q - 2);
  if (M_eff < 1) {
    return 0;
  }

  const int q1 = q - 1;
  arma::vec V(q1, arma::fill::zeros);
  double tail_sum = 0.0;
  for (int j = q1 - 1; j >= 0; --j) {
    tail_sum += vals_desc[j];
    V[j] = tail_sum;
  }

  arma::vec gr(M_eff, arma::fill::zeros);
  gr.fill(-arma::datum::inf);
  for (int j = 0; j < M_eff; ++j) {
    const double num = std::log1p(vals_desc[j] / std::max(V[j], eps));
    const double den = std::log1p(vals_desc[j + 1] / std::max(V[j + 1], eps));
    gr[j] = num / std::max(den, eps);
  }

  const arma::uword best_idx = gr.index_max();
  return static_cast<int>(best_idx) + 1;
}

arma::mat soft_threshold_offdiag_impl(const arma::mat& A, const double lambda) {
  arma::mat out = A;
  const arma::uword p = out.n_cols;

  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword i = 0; i < p; ++i) {
      if (i == j) {
        continue;
      }
      const double a = out(i, j);
      if (a > lambda) {
        out(i, j) = a - lambda;
      } else if (a < -lambda) {
        out(i, j) = a + lambda;
      } else {
        out(i, j) = 0.0;
      }
    }
  }

  return out;
}

arma::mat poet_from_projected_impl(
    const arma::mat& A_pd,
    const arma::vec& vals_desc,
    const arma::mat& vecs_desc,
    const double lambda_u,
    const int m,
    const double pd_eps) {
  const arma::uword p = A_pd.n_cols;
  arma::mat spike(p, p, arma::fill::zeros);

  if (m > 0) {
    const arma::mat spike_vecs = vecs_desc.cols(0, m - 1);
    const arma::vec spike_vals = vals_desc.subvec(0, m - 1);
    spike = spike_vecs * arma::diagmat(spike_vals) * spike_vecs.t();
  }

  const arma::mat remainder = A_pd - spike;
  arma::mat remainder_thr = soft_threshold_offdiag_impl(remainder, lambda_u);
  remainder_thr.diag() = remainder.diag();

  return project_pd_impl(spike + remainder_thr, pd_eps);
}

double median_from_vec(std::vector<double>& vals) {
  const int n = static_cast<int>(vals.size());
  if (n == 0) {
    return NA_REAL;
  }

  const int mid = n / 2;
  std::nth_element(vals.begin(), vals.begin() + mid, vals.end());
  const double upper = vals[mid];

  if (n % 2 == 1) {
    return upper;
  }

  std::nth_element(vals.begin(), vals.begin() + mid - 1, vals.begin() + mid);
  const double lower = vals[mid - 1];
  return 0.5 * (lower + upper);
}

// [[Rcpp::export]]
arma::mat project_pd_cpp(
    const arma::mat& M,
    const double eps = 1e-6) {
  return project_pd_impl(M, eps);
}

// [[Rcpp::export]]
List factor_poet_map_cpp(
    const arma::mat& A,
    const double lambda_u,
    const int M = 8,
    const double eig_eps = 1e-8,
    const double pd_eps = 1e-6) {
  arma::mat A_sym = symmetrize_arma(A);
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, A_sym);
  eigval.transform([&](double x) { return x < pd_eps ? pd_eps : x; });

  const arma::uvec ord = arma::sort_index(eigval, "descend");
  const arma::vec vals_desc = eigval.elem(ord);
  const arma::mat vecs_desc = eigvec.cols(ord);
  const arma::mat A_pd = eigvec * arma::diagmat(eigval) * eigvec.t();

  const int m = select_factor_count_from_desc(vals_desc, M, eig_eps);
  const arma::mat Sigma = poet_from_projected_impl(
    symmetrize_arma(A_pd),
    vals_desc,
    vecs_desc,
    lambda_u,
    m,
    pd_eps
  );

  return List::create(
    _["m"] = m,
    _["Sigma"] = Sigma
  );
}

// [[Rcpp::export]]
arma::mat poet_map_cpp(
    const arma::mat& A,
    const double lambda_u,
    const int m = 0,
    const double pd_eps = 1e-6) {
  arma::mat A_sym = symmetrize_arma(A);
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, A_sym);
  eigval.transform([&](double x) { return x < pd_eps ? pd_eps : x; });

  const arma::uvec ord = arma::sort_index(eigval, "descend");
  const arma::vec vals_desc = eigval.elem(ord);
  const arma::mat vecs_desc = eigvec.cols(ord);
  const arma::mat A_pd = eigvec * arma::diagmat(eigval) * eigvec.t();

  return poet_from_projected_impl(
    symmetrize_arma(A_pd),
    vals_desc,
    vecs_desc,
    lambda_u,
    m,
    pd_eps
  );
}

// [[Rcpp::export]]
NumericMatrix l1_distance_matrix_cpp(
    const NumericMatrix& X,
    const NumericMatrix& centers,
    const IntegerVector& use_dims) {
  const int n = X.nrow();
  const int p = X.ncol();
  const int K = centers.nrow();
  const int d = use_dims.size();
  NumericMatrix out(n, K);

  if (d == 0) {
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < K; ++k) {
        double acc = 0.0;
        for (int j = 0; j < p; ++j) {
          acc += std::abs(X(i, j) - centers(k, j));
        }
        out(i, k) = acc;
      }
    }
    return out;
  }

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < K; ++k) {
      double acc = 0.0;
      for (int idx = 0; idx < d; ++idx) {
        const int j = use_dims[idx] - 1;
        acc += std::abs(X(i, j) - centers(k, j));
      }
      out(i, k) = acc;
    }
  }

  return out;
}

// [[Rcpp::export]]
arma::mat calc_delta_cpp(
    const arma::mat& X,
    const arma::mat& Mu,
    const arma::mat& Omega,
    const double eps) {
  const arma::uword n = X.n_rows;
  const arma::uword K = Mu.n_rows;
  arma::mat out(n, K, arma::fill::zeros);

  for (arma::uword g = 0; g < K; ++g) {
    arma::mat centered = X.each_row() - Mu.row(g);
    arma::vec quad = arma::sum((centered * Omega) % centered, 1);
    quad.transform([&](double x) { return x < eps ? eps : x; });
    out.col(g) = quad;
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix softmax_rows_cpp(const NumericMatrix& log_mat) {
  const int n = log_mat.nrow();
  const int K = log_mat.ncol();
  NumericMatrix out(n, K);

  for (int i = 0; i < n; ++i) {
    double row_max = log_mat(i, 0);
    for (int k = 1; k < K; ++k) {
      if (log_mat(i, k) > row_max) {
        row_max = log_mat(i, k);
      }
    }

    double row_sum = 0.0;
    for (int k = 0; k < K; ++k) {
      const double val = std::exp(log_mat(i, k) - row_max);
      out(i, k) = val;
      row_sum += val;
    }

    if (!R_finite(row_sum) || row_sum <= 0.0) {
      const double uniform_prob = 1.0 / static_cast<double>(K);
      for (int k = 0; k < K; ++k) {
        out(i, k) = uniform_prob;
      }
    } else {
      for (int k = 0; k < K; ++k) {
        out(i, k) /= row_sum;
      }
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericVector row_logsumexp_cpp(const NumericMatrix& log_mat) {
  const int n = log_mat.nrow();
  const int K = log_mat.ncol();
  NumericVector out(n);

  for (int i = 0; i < n; ++i) {
    double row_max = log_mat(i, 0);
    for (int k = 1; k < K; ++k) {
      if (log_mat(i, k) > row_max) {
        row_max = log_mat(i, k);
      }
    }

    double row_sum = 0.0;
    for (int k = 0; k < K; ++k) {
      row_sum += std::exp(log_mat(i, k) - row_max);
    }

    out[i] = row_max + std::log(row_sum);
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix weighted_spatial_sign_scatter_cpp(
    const NumericMatrix& residuals,
    const NumericVector& weights,
    const double eps_r) {
  const int n = residuals.nrow();
  const int p = residuals.ncol();
  NumericMatrix out(p, p);
  double weight_sum = 0.0;

  for (int i = 0; i < n; ++i) {
    double norm2 = 0.0;
    for (int j = 0; j < p; ++j) {
      const double val = residuals(i, j);
      norm2 += val * val;
    }

    const double coeff = weights[i] / std::max(norm2, eps_r);
    weight_sum += weights[i];

    if (coeff == 0.0) {
      continue;
    }

    for (int a = 0; a < p; ++a) {
      const double ra = residuals(i, a);
      for (int b = a; b < p; ++b) {
        out(a, b) += coeff * ra * residuals(i, b);
      }
    }
  }

  if (weight_sum > 0.0) {
    for (int a = 0; a < p; ++a) {
      for (int b = a; b < p; ++b) {
        out(a, b) /= weight_sum;
        if (a != b) {
          out(b, a) = out(a, b);
        }
      }
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericVector quadform_rows_cpp(
    const NumericMatrix& residuals,
    const NumericMatrix& Omega) {
  const int n = residuals.nrow();
  const int p = residuals.ncol();
  NumericVector out(n);
  NumericVector tmp(p);

  for (int i = 0; i < n; ++i) {
    for (int a = 0; a < p; ++a) {
      double acc = 0.0;
      for (int b = 0; b < p; ++b) {
        acc += residuals(i, b) * Omega(b, a);
      }
      tmp[a] = acc;
    }

    double quad = 0.0;
    for (int j = 0; j < p; ++j) {
      quad += tmp[j] * residuals(i, j);
    }
    out[i] = quad;
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix weighted_crossprod_cpp(
    const NumericMatrix& residuals,
    const NumericVector& coeff) {
  const int n = residuals.nrow();
  const int p = residuals.ncol();
  NumericMatrix out(p, p);

  for (int i = 0; i < n; ++i) {
    const double ci = coeff[i];
    if (ci == 0.0) {
      continue;
    }

    for (int a = 0; a < p; ++a) {
      const double ra = residuals(i, a);
      for (int b = a; b < p; ++b) {
        out(a, b) += ci * ra * residuals(i, b);
      }
    }
  }

  for (int a = 0; a < p; ++a) {
    for (int b = a + 1; b < p; ++b) {
      out(b, a) = out(a, b);
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericVector weighted_kde_eval_cpp(
    const NumericVector& y_eval,
    const NumericVector& y_obs,
    const NumericVector& weights,
    const double bandwidth) {
  const int m = y_eval.size();
  const int n = y_obs.size();
  NumericVector out(m);

  for (int i = 0; i < m; ++i) {
    double acc = 0.0;
    const double yi = y_eval[i];
    for (int j = 0; j < n; ++j) {
      const double z = (yi - y_obs[j]) / bandwidth;
      acc += R::dnorm4(z, 0.0, 1.0, 0) * weights[j];
    }
    out[i] = acc / bandwidth;
  }

  return out;
}

// [[Rcpp::export]]
NumericVector weighted_kde_eval_grid_cpp(
    const double y_start,
    const double y_step,
    const int grid_size,
    const NumericVector& y_obs,
    const NumericVector& weights,
    const double bandwidth) {
  const int n = y_obs.size();
  NumericVector out(grid_size);

  if (grid_size <= 0) {
    return out;
  }

  const double step_scaled = y_step / bandwidth;
  const double ratio_ratio = std::exp(-(step_scaled * step_scaled));

  for (int j = 0; j < n; ++j) {
    const double z0 = (y_start - y_obs[j]) / bandwidth;
    double kernel_val = std::exp(-0.5 * z0 * z0);
    double ratio = std::exp(-(z0 * step_scaled) - 0.5 * step_scaled * step_scaled);
    const double wj = weights[j];

    for (int i = 0; i < grid_size; ++i) {
      out[i] += wj * kernel_val;
      kernel_val *= ratio;
      ratio *= ratio_ratio;
    }
  }

  const double scale = INV_SQRT_2PI / bandwidth;
  for (int i = 0; i < grid_size; ++i) {
    out[i] *= scale;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector col_medians_cpp(const NumericMatrix& X) {
  const int n = X.nrow();
  const int p = X.ncol();
  NumericVector out(p);
  std::vector<double> vals(n);

  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < n; ++i) {
      vals[i] = X(i, j);
    }
    out[j] = median_from_vec(vals);
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix cluster_centers_median_cpp(
    const NumericMatrix& X,
    const IntegerVector& cluster,
    const int K) {
  const int n = X.nrow();
  const int p = X.ncol();
  NumericMatrix centers(K, p);
  IntegerVector counts(K);

  for (int i = 0; i < n; ++i) {
    const int k = cluster[i] - 1;
    if (k >= 0 && k < K) {
      counts[k] += 1;
    }
  }

  std::vector<std::vector<double>> buffers(K);
  for (int k = 0; k < K; ++k) {
    if (counts[k] > 0) {
      buffers[k].reserve(counts[k]);
    }
  }

  for (int k = 0; k < K; ++k) {
    if (counts[k] == 0) {
      for (int j = 0; j < p; ++j) {
        centers(k, j) = NA_REAL;
      }
      continue;
    }

    for (int j = 0; j < p; ++j) {
      buffers[k].clear();
      for (int i = 0; i < n; ++i) {
        if (cluster[i] == (k + 1)) {
          buffers[k].push_back(X(i, j));
        }
      }
      centers(k, j) = median_from_vec(buffers[k]);
    }
  }

  return centers;
}

// [[Rcpp::export]]
double compute_between_ss_l1_cpp(
    const NumericMatrix& X,
    const IntegerVector& cluster,
    const int K) {
  const int n = X.nrow();
  const int p = X.ncol();
  NumericVector overall_med = col_medians_cpp(X);
  NumericMatrix cluster_meds = cluster_centers_median_cpp(X, cluster, K);
  IntegerVector counts(K);
  double out = 0.0;

  for (int i = 0; i < n; ++i) {
    const int k = cluster[i] - 1;
    if (k >= 0 && k < K) {
      counts[k] += 1;
    }
  }

  for (int k = 0; k < K; ++k) {
    if (counts[k] == 0) {
      continue;
    }
    double acc = 0.0;
    for (int j = 0; j < p; ++j) {
      acc += std::abs(cluster_meds(k, j) - overall_med[j]);
    }
    out += static_cast<double>(counts[k]) * acc;
  }

  return out;
}
