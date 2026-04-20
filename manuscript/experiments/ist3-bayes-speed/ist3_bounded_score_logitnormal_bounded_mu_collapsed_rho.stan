functions {
  real logitnormal_interval_log(real lo, real hi, real mu, real sigma) {
    if (hi <= lo) {
      return negative_infinity();
    }

    if (lo <= 0 && hi >= 1) {
      return 0;
    }

    if (lo <= 0) {
      return normal_lcdf(logit(hi) | mu, sigma);
    }

    if (hi >= 1) {
      return normal_lccdf(logit(lo) | mu, sigma);
    }

    return log_diff_exp(
      normal_lcdf(logit(hi) | mu, sigma),
      normal_lcdf(logit(lo) | mu, sigma)
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
      array[] vector mu_logit_comp,
      array[] vector sigma_comp,
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
        real mu = mu_logit_comp[r][h];
        real sigma = sigma_comp[r][h];

        for (k in 1:K) {
          if (bin_valid[k, j] == 1) {
            grid_lp[k] =
              log(eta[eta_group][k]) +
              logitnormal_interval_log(bin_lower[k, j], bin_upper[k, j], mu, sigma);
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
  real mu_min = -3;
  real mu_max = 3;
  real mu_range = mu_max - mu_min;

  if (H != 2) {
    reject("The external collapsed-rho bounded-mu logit-normal IST-3 speed model requires H = 2.");
  }
}

parameters {
  array[2] vector<lower=0, upper=1>[H - 1] v;
  vector<lower=0, upper=1>[2] mu_low_unit;
  vector<lower=1e-3, upper=1 - 1e-3>[2] mu_gap_prop;
  array[2] vector<lower=log(0.05), upper=log(10)>[H] log_sigma_comp;
  vector<lower=0>[2] alpha;
  array[eta_groups] simplex[K] eta;
}

transformed parameters {
  array[2] simplex[H] w;
  array[2] vector[H] mu_logit_comp;
  array[2] vector<lower=0>[H] sigma_comp;

  for (r in 1:2) {
    real remaining = 1;
    real mu_low = mu_min + mu_range * mu_low_unit[r];

    for (h in 1:(H - 1)) {
      w[r][h] = v[r][h] * remaining;
      remaining = remaining * (1 - v[r][h]);
    }

    w[r][H] = remaining;
    mu_logit_comp[r][1] = mu_low;
    mu_logit_comp[r][2] = mu_low + (mu_max - mu_low) * mu_gap_prop[r];
    sigma_comp[r] = exp(log_sigma_comp[r]);
  }
}

model {
  for (r in 1:2) {
    alpha[r] ~ gamma(alpha_prior_shape, alpha_prior_rate);
    for (h in 1:(H - 1)) {
      v[r][h] ~ beta(1, alpha[r]);
    }

    target += normal_lpdf(mu_logit_comp[r][1] | 0, 2) + log(mu_range);
    target += normal_lpdf(mu_logit_comp[r][2] | 0, 2) +
      log(mu_max - mu_logit_comp[r][1]);
    log_sigma_comp[r] ~ normal(log(1.2), 0.5);
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
    mu_logit_comp,
    sigma_comp,
    eta,
    bin_lower,
    bin_upper,
    bin_valid
  );
}

generated quantities {
  vector<lower=0, upper=1>[2] rho;
  vector<lower=0, upper=1>[2] pi;

  for (r in 1:2) {
    rho[r] = beta_rng(
      rho_prior_alpha + n_arm[r] - n_obs_arm[r],
      rho_prior_beta + n_obs_arm[r]
    );
    pi[r] = 1 - rho[r];
  }
}
