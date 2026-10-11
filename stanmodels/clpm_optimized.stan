functions {
  real exchangeable_normal_kernel(
      vector y, vector mu, real diagonal_variance,
      real off_diagonal_covariance) {
    int K = rows(y);
    real independent_variance = diagonal_variance - off_diagonal_covariance;
    real mean_mode_variance =
      independent_variance + K * off_diagonal_covariance;
    vector[K] delta = y - mu;
    real quad = dot_self(delta) / independent_variance
      - off_diagonal_covariance * square(sum(delta))
        / (independent_variance * mean_mode_variance);
    return -0.5 * quad;
  }

  real ri_ar1_centered_sufficient_logpdf(
      int N, int T, vector z_bar, vector mu, matrix A,
      matrix L_ri, matrix L_w1, matrix L_e,
      matrix S11, matrix S_pre, matrix S_post, matrix S_lag,
      matrix S_UU, matrix S_VV, matrix S_UV,
      matrix S_1U, matrix S_1V) {
    int P = rows(A);
    matrix[P, P] I_P = identity_matrix(P);
    matrix[P, P] J = I_P - A;
    matrix[P, P] Qw = chol2inv(L_w1);
    matrix[P, P] Qe = chol2inv(L_e);
    matrix[P, P] J_L_ri = J * L_ri;
    matrix[P, P] H =
      I_P
      + quad_form(Qw, L_ri)
      + (T - 1) * quad_form(Qe, J_L_ri);
    matrix[P, P] L_H;
    matrix[P, P] G;
    matrix[P, P] S_innov;
    matrix[P, P] S_ss;
    matrix[P, P] S_1s;
    matrix[P, P] C;
    matrix[P, P] S_h;
    real scatter_quad;
    vector[P * T] d = z_bar - mu;
    vector[P] d_first = segment(d, 1, P);
    vector[P] d_prev = d_first;
    vector[P] innov_sum = rep_vector(0.0, P);
    real mean_base = quad_form(Qw, d_first);
    vector[P] h_bar;
    real mean_quad;
    real half_log_det;

    H = 0.5 * (H + H');
    L_H = cholesky_decompose(H);
    G = mdivide_left_tri_low(L_H, L_ri');

    // Only symmetric-weighted traces of these matrices are used below.
    // Each transpose-paired cross term can therefore be evaluated once.
    S_innov = S_post - 2 * A * S_lag + A * S_pre * A';
    S_ss = S_UU - 2 * S_UV * A' + A * S_VV * A';
    S_1s = S_1U - S_1V * A';
    C = J' * Qe;
    S_h = Qw * S11 * Qw
      + 2 * Qw * S_1s * C'
      + C * S_ss * C';

    scatter_quad =
      trace_dot(Qw, S11)
      + trace_dot(Qe, S_innov)
      - trace_quad_form(S_h, G');

    for (t in 2:T) {
      vector[P] d_curr = segment(d, P * (t - 1) + 1, P);
      vector[P] innov = d_curr - A * d_prev;
      mean_base += quad_form(Qe, innov);
      innov_sum += innov;
      d_prev = d_curr;
    }

    h_bar = Qw * d_first + C * innov_sum;
    mean_quad = mean_base - dot_self(G * h_bar);
    half_log_det =
      sum(log(diagonal(L_w1)))
      + (T - 1) * sum(log(diagonal(L_e)))
      + sum(log(diagonal(L_H)));

    return -N * half_log_det
      - 0.5 * (scatter_quad + N * mean_quad);
  }
}

data {
  int<lower=1> N;
  int<lower=2> P;
  int<lower=2> T;
  int<lower=0, upper=1> prior_only;
  vector[P * T] z_sum;
  matrix[P * T, P * T] z_crossprod;
  vector[P * (P - 1)] slab_mean;
  cov_matrix[P * (P - 1)] slab_cov;
}

transformed data {
  int K = P * T;
  vector[K] z_bar = z_sum / N;
  matrix[P, P] S11 = rep_matrix(0.0, P, P);
  matrix[P, P] S_pre = rep_matrix(0.0, P, P);
  matrix[P, P] S_post = rep_matrix(0.0, P, P);
  matrix[P, P] S_lag = rep_matrix(0.0, P, P);
  matrix[P, P] S_UU = rep_matrix(0.0, P, P);
  matrix[P, P] S_VV = rep_matrix(0.0, P, P);
  matrix[P, P] S_UV = rep_matrix(0.0, P, P);
  matrix[P, P] S_1U = rep_matrix(0.0, P, P);
  matrix[P, P] S_1V = rep_matrix(0.0, P, P);
  real slab_diagonal_variance = slab_cov[1, 1];
  real slab_off_diagonal_covariance = slab_cov[1, 2];

  {
    real slab_tolerance = 1e-10 * fmax(1.0, slab_diagonal_variance);
    matrix[K, K] S = z_crossprod - N * (z_bar * z_bar');
    for (i in 1:(P * (P - 1))) {
      for (j in 1:(P * (P - 1))) {
        real expected = i == j
          ? slab_diagonal_variance : slab_off_diagonal_covariance;
        if (abs(slab_cov[i, j] - expected) > slab_tolerance) {
          reject("The optimized slab kernel requires an exchangeable covariance");
        }
      }
    }
    S = 0.5 * (S + S');
    S11 = block(S, 1, 1, P, P);

    for (t in 1:(T - 1)) {
      int i_pre = P * (t - 1) + 1;
      int i_post = P * t + 1;
      S_pre += block(S, i_pre, i_pre, P, P);
      S_post += block(S, i_post, i_post, P, P);
      S_lag += block(S, i_pre, i_post, P, P);
    }

    for (t in 2:T) {
      int i_t = P * (t - 1) + 1;
      S_1U += block(S, 1, i_t, P, P);
      for (s in 2:T) {
        int i_s = P * (s - 1) + 1;
        S_UU += block(S, i_t, i_s, P, P);
      }
    }

    for (t in 1:(T - 1)) {
      int i_t = P * (t - 1) + 1;
      S_1V += block(S, 1, i_t, P, P);
      for (s in 1:(T - 1)) {
        int i_s = P * (s - 1) + 1;
        S_VV += block(S, i_t, i_s, P, P);
      }
    }

    for (t in 2:T) {
      int i_t = P * (t - 1) + 1;
      for (s in 1:(T - 1)) {
        int i_s = P * (s - 1) + 1;
        S_UV += block(S, i_t, i_s, P, P);
      }
    }
  }
}

parameters {
  vector[P * T] mu;
  matrix[P, P] A;
  vector<lower=0>[P] sd_ri;
  cholesky_factor_corr[P] Lcorr_ri;
  vector<lower=0>[P] sd_w1;
  cholesky_factor_corr[P] Lcorr_w1;
  vector<lower=0>[P] sd_e;
  cholesky_factor_corr[P] Lcorr_e;
}

model {
  vector[P * (P - 1)] cross_lagged;
  int cross_lagged_index = 1;

  for (column in 1:P) {
    for (row in 1:P) {
      if (row != column) {
        cross_lagged[cross_lagged_index] = A[row, column];
        cross_lagged_index += 1;
      }
    }
  }

  mu ~ normal(0, 1.5);
  for (j in 1:P) {
    A[j, j] ~ normal(0, 0.35);
  }
  target += exchangeable_normal_kernel(
    cross_lagged, slab_mean, slab_diagonal_variance,
    slab_off_diagonal_covariance);
  sd_ri ~ normal(0, 1);
  sd_w1 ~ normal(0, 1);
  sd_e ~ normal(0, 1);
  Lcorr_ri ~ lkj_corr_cholesky(2);
  Lcorr_w1 ~ lkj_corr_cholesky(2);
  Lcorr_e ~ lkj_corr_cholesky(2);

  if (prior_only == 0) {
    matrix[P, P] L_ri = diag_pre_multiply(sd_ri, Lcorr_ri);
    matrix[P, P] L_w1 = diag_pre_multiply(sd_w1, Lcorr_w1);
    matrix[P, P] L_e = diag_pre_multiply(sd_e, Lcorr_e);
    target += ri_ar1_centered_sufficient_logpdf(
      N, T, z_bar, mu, A, L_ri, L_w1, L_e,
      S11, S_pre, S_post, S_lag,
      S_UU, S_VV, S_UV, S_1U, S_1V);
  }
}

generated quantities {
  matrix[P, P] Corr_ri =
    multiply_lower_tri_self_transpose(Lcorr_ri);
  matrix[P, P] Corr_w1 =
    multiply_lower_tri_self_transpose(Lcorr_w1);
  matrix[P, P] Corr_e =
    multiply_lower_tri_self_transpose(Lcorr_e);
  real ar_x = A[1, 1];
  real y_to_x = A[1, 2];
  real x_to_y = A[2, 1];
  real ar_y = A[2, 2];
}
