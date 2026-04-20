functions {
  real beta_interval_log(real lo, real hi, real beta_a, real beta_b) {
    if (hi <= lo) {
      return negative_infinity();
    }

    if (lo <= 0 && hi >= 1) {
      return 0;
    }

    if (lo <= 0) {
      return beta_lcdf(hi | beta_a, beta_b);
    }

    if (hi >= 1) {
      return beta_lccdf(lo | beta_a, beta_b);
    }

    return log_diff_exp(
      beta_lcdf(hi | beta_a, beta_b),
      beta_lcdf(lo | beta_a, beta_b)
    );
  }

  real partial_score_sum(
      array[] int score_cell_count_slice,
      int start,
      int end,
      array[] int score_cell_arm,
      array[] int score_cell_index,
      array[] int eta_group_by_arm,
      int H,
      int K,
      array[] vector w,
      array[] vector m_comp,
      array[] vector phi_comp,
      array[] vector eta,
      matrix bin_lower,
      matrix bin_upper,
      array[,] int bin_valid) {
    real lp = 0;

    for (slice_i in 1:size(score_cell_count_slice)) {
      int c = start + slice_i - 1;
      vector[H] component_lp;
      int r = score_cell_arm[c];
      int j = score_cell_index[c];
      int eta_group = eta_group_by_arm[r];

      for (h in 1:H) {
        vector[K] grid_lp;
        real beta_a = m_comp[r][h] * phi_comp[r][h];
        real beta_b = (1 - m_comp[r][h]) * phi_comp[r][h];

        for (k in 1:K) {
          if (bin_valid[k, j] == 1) {
            grid_lp[k] =
              log(eta[eta_group][k]) +
              beta_interval_log(bin_lower[k, j], bin_upper[k, j], beta_a, beta_b);
          } else {
            grid_lp[k] = negative_infinity();
          }
        }

        component_lp[h] = log(w[r][h]) + log_sum_exp(grid_lp);
      }

      lp += score_cell_count_slice[slice_i] * log_sum_exp(component_lp);
    }

    return lp;
  }
}

data {
  int<lower=2> H;
  int<lower=1> J;
  vector[J] score_value;
  int<lower=1> K;
  matrix[K, J] bin_lower;
  matrix[K, J] bin_upper;
  array[K, J] int<lower=0, upper=1> bin_valid;
  int<lower=1, upper=2> eta_groups;
  array[2] int<lower=1, upper=2> eta_group_by_arm;

  int<lower=1> N_score_cells;
  array[N_score_cells] int<lower=1, upper=2> score_cell_arm;
  array[N_score_cells] int<lower=1, upper=J> score_cell_index;
  array[N_score_cells] int<lower=1> score_cell_count;

  array[2] int<lower=0> n_arm;
  array[2] int<lower=0> n_obs_arm;
  int<lower=1> grainsize;

  real atom;
  real<lower=0> rho_prior_alpha;
  real<lower=0> rho_prior_beta;
  real<lower=0> alpha_prior_shape;
  real<lower=0> alpha_prior_rate;
  real<lower=0> m_prior_alpha;
  real<lower=0> m_prior_beta;
  real phi_prior_meanlog;
  real<lower=0> phi_prior_sdlog;
  real<lower=0> gap_prior_alpha;
  real<lower=0> gap_prior_beta;
  int<lower=0, upper=1> use_fixed_v_prior;
  real<lower=0> v_prior_alpha;
  real<lower=0> v_prior_beta;
  vector<lower=0>[K] eta_prior;
}

transformed data {
  real mean_eps = 1e-5;

  if (H != 2) {
    reject("The non-centered-stick IST-3 speed model currently requires H = 2.");
  }
}

parameters {
  array[2] vector<lower=0, upper=1>[H - 1] u_stick;
  array[2] ordered[H] m_logit;
  array[2] vector<lower=log(0.25), upper=log(100)>[H] log_phi_comp;
  vector<lower=0>[2] alpha;
  vector<lower=0, upper=1>[2] rho;
  array[eta_groups] simplex[K] eta;
}

transformed parameters {
  array[2] vector<lower=0, upper=1>[H - 1] v;
  array[2] simplex[H] w;
  array[2] vector<lower=0, upper=1>[H] m_comp;
  array[2] vector<lower=0>[H] phi_comp;
  vector<lower=0, upper=1>[2] pi;

  for (r in 1:2) {
    real remaining = 1;

    for (h in 1:(H - 1)) {
      v[r][h] = 1 - exp(log(u_stick[r][h]) / alpha[r]);
      w[r][h] = v[r][h] * remaining;
      remaining = remaining * (1 - v[r][h]);
    }

    w[r][H] = remaining;
    for (h in 1:H) {
      m_comp[r][h] = mean_eps + (1 - 2 * mean_eps) * inv_logit(m_logit[r][h]);
    }
    phi_comp[r] = exp(log_phi_comp[r]);
    pi[r] = 1 - rho[r];
  }
}

model {
  for (r in 1:2) {
    rho[r] ~ beta(rho_prior_alpha, rho_prior_beta);
    alpha[r] ~ gamma(alpha_prior_shape, alpha_prior_rate);
    for (h in 1:H) {
      real m_unit = inv_logit(m_logit[r][h]);
      target += beta_lpdf(m_comp[r][h] | m_prior_alpha, m_prior_beta) +
        log1m(2 * mean_eps) + log(m_unit) + log1m(m_unit);
    }
    log_phi_comp[r] ~ normal(phi_prior_meanlog, phi_prior_sdlog);
    n_obs_arm[r] ~ binomial(n_arm[r], pi[r]);
  }

  for (g in 1:eta_groups) {
    eta[g] ~ dirichlet(eta_prior);
  }

  target += reduce_sum(
    partial_score_sum,
    score_cell_count,
    grainsize,
    score_cell_arm,
    score_cell_index,
    eta_group_by_arm,
    H,
    K,
    w,
    m_comp,
    phi_comp,
    eta,
    bin_lower,
    bin_upper,
    bin_valid
  );
}
