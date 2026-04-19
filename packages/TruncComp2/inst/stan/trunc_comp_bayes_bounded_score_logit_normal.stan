functions {
  real logitnormal_interval_log(real lower, real upper, real mu, real sigma) {
    if (upper <= lower) {
      return negative_infinity();
    }

    if (lower <= 0 && upper >= 1) {
      return 0;
    }

    if (lower <= 0) {
      return normal_lcdf(logit(upper) | mu, sigma);
    }

    if (upper >= 1) {
      return normal_lccdf(logit(lower) | mu, sigma);
    }

    return log_diff_exp(
      normal_lcdf(logit(upper) | mu, sigma),
      normal_lcdf(logit(lower) | mu, sigma)
    );
  }
}

data {
  int<lower=2> H;
  int<lower=1> J;
  int<lower=0> n_arm[2];
  int<lower=0> n_obs_arm[2];
  int<lower=1> N_score_cells;
  int<lower=1, upper=2> score_cell_arm[N_score_cells];
  int<lower=1, upper=J> score_cell_index[N_score_cells];
  int<lower=1> score_cell_count[N_score_cells];
  vector[J] score_value;
  int<lower=1> K;
  matrix[K, J] bin_lower;
  matrix[K, J] bin_upper;
  int<lower=0, upper=1> bin_valid[K, J];
  int<lower=1, upper=2> eta_groups;
  int<lower=1, upper=2> eta_group_by_arm[2];
  real atom;
  real<lower=0> rho_prior_alpha;
  real<lower=0> rho_prior_beta;
  real<lower=0> alpha_prior_shape;
  real<lower=0> alpha_prior_rate;
  real mu_logit_prior_mean;
  real<lower=0> mu_logit_prior_sd;
  real sigma_logit_prior_meanlog;
  real<lower=0> sigma_logit_prior_sdlog;
  vector<lower=0>[K] eta_prior;
}

parameters {
  vector<lower=0, upper=1>[H - 1] v[2];
  ordered[H] mu_logit_comp[2];
  vector[H] log_sigma_comp[2];
  vector<lower=0>[2] alpha;
  vector<lower=0, upper=1>[2] rho;
  simplex[K] eta[eta_groups];
}

transformed parameters {
  simplex[H] w[2];
  vector<lower=0>[H] sigma_logit_comp[2];
  vector<lower=0, upper=1>[2] pi;

  for (r in 1:2) {
    real remaining;
    remaining = 1;

    for (h in 1:(H - 1)) {
      w[r][h] = v[r][h] * remaining;
      remaining = remaining * (1 - v[r][h]);
    }

    w[r][H] = remaining;
    sigma_logit_comp[r] = exp(log_sigma_comp[r]);
    pi[r] = 1 - rho[r];
  }
}

model {
  for (r in 1:2) {
    rho[r] ~ beta(rho_prior_alpha, rho_prior_beta);
    alpha[r] ~ gamma(alpha_prior_shape, alpha_prior_rate);
    for (h in 1:(H - 1)) {
      v[r][h] ~ beta(1, alpha[r]);
    }
    mu_logit_comp[r] ~ normal(mu_logit_prior_mean, mu_logit_prior_sd);
    log_sigma_comp[r] ~ normal(sigma_logit_prior_meanlog, sigma_logit_prior_sdlog);
  }

  for (g in 1:eta_groups) {
    eta[g] ~ dirichlet(eta_prior);
  }

  for (r in 1:2) {
    n_obs_arm[r] ~ binomial(n_arm[r], pi[r]);
  }

  for (c in 1:N_score_cells) {
    vector[H] component_lp;
    int r;
    int j;
    int eta_group;
    r = score_cell_arm[c];
    j = score_cell_index[c];
    eta_group = eta_group_by_arm[r];

    for (h in 1:H) {
      vector[K] grid_lp;

      for (k in 1:K) {
        if (bin_valid[k, j] == 1) {
          grid_lp[k] =
            log(eta[eta_group][k]) +
            logitnormal_interval_log(
              bin_lower[k, j],
              bin_upper[k, j],
              mu_logit_comp[r][h],
              sigma_logit_comp[r][h]
            );
        } else {
          grid_lp[k] = negative_infinity();
        }
      }

      component_lp[h] = log(w[r][h]) + log_sum_exp(grid_lp);
    }

    target += score_cell_count[c] * log_sum_exp(component_lp);
  }
}

generated quantities {
  real rho_0;
  real rho_1;
  real pi_0;
  real pi_1;
  matrix[2, J] survivor_score_pmf;
  real mu_0_c;
  real mu_1_c;
  real delta_atom;
  real mu_delta;
  real alpha_delta;
  real delta;

  rho_0 = rho[1];
  rho_1 = rho[2];
  pi_0 = pi[1];
  pi_1 = pi[2];

  for (r in 1:2) {
    int eta_group;
    eta_group = eta_group_by_arm[r];

    for (j in 1:J) {
      vector[H] component_lp;

      for (h in 1:H) {
        vector[K] grid_lp;

        for (k in 1:K) {
          if (bin_valid[k, j] == 1) {
            grid_lp[k] =
              log(eta[eta_group][k]) +
              logitnormal_interval_log(
                bin_lower[k, j],
                bin_upper[k, j],
                mu_logit_comp[r][h],
                sigma_logit_comp[r][h]
              );
          } else {
            grid_lp[k] = negative_infinity();
          }
        }

        component_lp[h] = log(w[r][h]) + log_sum_exp(grid_lp);
      }

      survivor_score_pmf[r, j] = exp(log_sum_exp(component_lp));
    }
  }

  mu_0_c = 0;
  mu_1_c = 0;
  for (j in 1:J) {
    mu_0_c = mu_0_c + score_value[j] * survivor_score_pmf[1, j];
    mu_1_c = mu_1_c + score_value[j] * survivor_score_pmf[2, j];
  }

  delta_atom = rho_1 - rho_0;
  mu_delta = mu_1_c - mu_0_c;
  alpha_delta = (pi_1 / (1 - pi_1)) / (pi_0 / (1 - pi_0));
  delta =
    (atom * rho_1 + pi_1 * mu_1_c) -
    (atom * rho_0 + pi_0 * mu_0_c);
}
