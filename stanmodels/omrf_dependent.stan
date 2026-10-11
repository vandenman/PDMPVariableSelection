functions {
  int pdmp_get_subsample_size();
  int pdmp_get_subsample_index(int n);

  real omrf_row_lp(
      array[] int x,
      row_vector x_real,
      matrix thresholds,
      matrix interactions,
      array[] int seen,
      int P) {
    real result = 0.0;
    row_vector[P] field = x_real * interactions;

    for (i in 1:P) {
      if (x[i] > 0) {
        result += thresholds[i, x[i]];
      }
    }
    if (P > 1) {
      for (i in 1:(P - 1)) {
        for (j in (i + 1):P) {
          result += 2 * interactions[i, j] * x[i] * x[j];
        }
      }
    }
    for (i in 1:P) {
      real max_eta = 0.0;
      real sum_exp;
      if (seen[i] > 1) {
        for (u in 1:(seen[i] - 1)) {
          max_eta = fmax(max_eta, thresholds[i, u] + u * field[i]);
        }
      }
      sum_exp = exp(-max_eta);
      if (seen[i] > 1) {
        for (u in 1:(seen[i] - 1)) {
          sum_exp += exp(thresholds[i, u] + u * field[i] - max_eta);
        }
      }
      result -= max_eta + log(sum_exp);
    }
    return result;
  }
}

data {
  int<lower=1> N;
  int<lower=2> P;
  int<lower=2> K;
  int<lower=0, upper=1> prior_only;
  array[N, P] int<lower=0> X;
  array[P] int<lower=1> seen;
  real<lower=0> prior_threshold_alpha;
  real<lower=0> prior_threshold_beta;
  real<lower=0> prior_node_scale_sd;
  real<lower=0> slab_base_scale;
}

transformed data {
  int E = (P * (P - 1)) %/% 2;
  int noThresholds = sum(seen) - P;
  matrix[N, P] X_real;
  array[N, P] int X_cat;
  array[P] int threshold_start;
  matrix[1, K] category_scores;

  for (n in 1:N) {
    for (i in 1:P) {
      X_real[n, i] = X[n, i];
      X_cat[n, i] = X[n, i] + 1;
    }
  }

  {
    int idx = 1;
    for (i in 1:P) {
      threshold_start[i] = idx;
      idx += seen[i] - 1;
    }
  }

  for (k in 1:K) {
    category_scores[1, k] = k - 1;
  }
}

parameters {
  vector[noThresholds] thresholds_0;
  vector[E] interactions_0;
  vector<lower=0>[P] node_scale;
}

transformed parameters {
  matrix[P, K] thresholds;
  matrix[P, P] interactions = rep_matrix(0.0, P, P);

  {
    int idx = 1;
    for (i in 1:P) {
      if (seen[i] > 1) {
        for (j in 1:(seen[i] - 1)) {
          thresholds[i, j] = thresholds_0[idx];
          idx += 1;
        }
      }
      for (j in seen[i]:K) {
        thresholds[i, j] = negative_infinity();
      }
    }
  }

  {
    int idx = 1;
    if (P > 1) {
      for (i in 1:(P - 1)) {
        for (j in (i + 1):P) {
          interactions[i, j] = interactions_0[idx];
          interactions[j, i] = interactions_0[idx];
          idx += 1;
        }
      }
    }
  }
}

model {
  target += sum(prior_threshold_alpha * thresholds_0
                - (prior_threshold_alpha + prior_threshold_beta) * log1p_exp(thresholds_0));
  node_scale ~ lognormal(0, prior_node_scale_sd);

  {
    int edge = 1;
    for (i in 1:(P - 1)) {
      for (j in (i + 1):P) {
        interactions_0[edge] ~ normal(
          0,
          slab_base_scale * sqrt(node_scale[i] * node_scale[j])
        );
        edge += 1;
      }
    }
  }

  if (prior_only == 0) {
    int m_sub = pdmp_get_subsample_size();
    if (m_sub > 0) {
      for (i in 1:m_sub) {
        int n = pdmp_get_subsample_index(i);
        target += omrf_row_lp(
          X[n], X_real[n], thresholds, interactions, seen, P);
      }
    } else {
      matrix[N, P] field = X_real * interactions;
      for (i in 1:P) {
        if (seen[i] > 1) {
          vector[K] alpha_i = rep_vector(negative_infinity(), K);
          matrix[N, 1] field_i;
          alpha_i[1] = 0;
          for (u in 1:(seen[i] - 1)) {
            alpha_i[u + 1] = thresholds_0[threshold_start[i] + u - 1];
          }
          for (n in 1:N) {
            field_i[n, 1] = field[n, i];
          }
          target += categorical_logit_glm_lpmf(
            X_cat[, i] | field_i, alpha_i, category_scores);
        }
      }
    }
  }
}
