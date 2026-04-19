functions {
  real logitnormal_unit_lpdf(real x, real mu, real sigma) {
    if (x <= 0 || x >= 1) {
      return negative_infinity();
    }

    return normal_lpdf(logit(x) | mu, sigma) - log(x) - log1m(x);
  }

  real logitnormal_unit_mean(vector nodes, vector weights, real mu, real sigma) {
    real out;
    out = 0;

    for (q in 1:num_elements(nodes)) {
      out = out + weights[q] * inv_logit(mu + sqrt(2.0) * sigma * nodes[q]);
    }

    return out;
  }
}

data {
  int<lower=1> N;
  int<lower=0, upper=1> A[N];
  int<lower=1, upper=2> arm[N];
  int<lower=2> H;
  int<lower=1> N_obs;
  vector<lower=0, upper=1>[N_obs] x_obs;
  int<lower=1, upper=2> arm_obs[N_obs];
  real score_min;
  real<lower=score_min> score_max;
  real atom;
  real<lower=0> rho_prior_alpha;
  real<lower=0> rho_prior_beta;
  real<lower=0> alpha_prior_shape;
  real<lower=0> alpha_prior_rate;
  real mu_logit_prior_mean;
  real<lower=0> mu_logit_prior_sd;
  real sigma_logit_prior_meanlog;
  real<lower=0> sigma_logit_prior_sdlog;
  int<lower=1> mean_quad_n;
  vector[mean_quad_n] mean_quad_node;
  vector<lower=0>[mean_quad_n] mean_quad_weight;
}

parameters {
  vector<lower=0, upper=1>[H - 1] v[2];
  ordered[H] mu_logit_comp[2];
  vector[H] log_sigma_comp[2];
  vector<lower=0>[2] alpha;
  vector<lower=0, upper=1>[2] rho;
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

  for (n in 1:N) {
    A[n] ~ bernoulli(pi[arm[n]]);
  }

  for (n in 1:N_obs) {
    vector[H] component_lp;

    for (h in 1:H) {
      component_lp[h] =
        log(w[arm_obs[n]][h]) +
        logitnormal_unit_lpdf(
          x_obs[n] |
          mu_logit_comp[arm_obs[n]][h],
          sigma_logit_comp[arm_obs[n]][h]
        ) -
        log(score_max - score_min);
    }

    target += log_sum_exp(component_lp);
  }
}

generated quantities {
  real rho_0;
  real rho_1;
  real pi_0;
  real pi_1;
  real mu_0_c_unit;
  real mu_1_c_unit;
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
  mu_0_c_unit = 0;
  mu_1_c_unit = 0;

  for (h in 1:H) {
    mu_0_c_unit =
      mu_0_c_unit +
      w[1][h] *
      logitnormal_unit_mean(
        mean_quad_node,
        mean_quad_weight,
        mu_logit_comp[1][h],
        sigma_logit_comp[1][h]
      );
    mu_1_c_unit =
      mu_1_c_unit +
      w[2][h] *
      logitnormal_unit_mean(
        mean_quad_node,
        mean_quad_weight,
        mu_logit_comp[2][h],
        sigma_logit_comp[2][h]
      );
  }

  mu_0_c = score_min + (score_max - score_min) * mu_0_c_unit;
  mu_1_c = score_min + (score_max - score_min) * mu_1_c_unit;
  delta_atom = rho_1 - rho_0;
  mu_delta = mu_1_c - mu_0_c;
  alpha_delta = (pi_1 / (1 - pi_1)) / (pi_0 / (1 - pi_0));
  delta =
    (atom * rho_1 + pi_1 * mu_1_c) -
    (atom * rho_0 + pi_0 * mu_0_c);
}
