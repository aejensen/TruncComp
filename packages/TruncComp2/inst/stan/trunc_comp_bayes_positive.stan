data {
  int<lower=1> N;
  int<lower=0, upper=1> A[N];
  int<lower=1, upper=2> arm[N];
  int<lower=2> H;
  int<lower=1> N_obs;
  vector<lower=0>[N_obs] y_obs_scaled;
  int<lower=1, upper=2> arm_obs[N_obs];
  real<lower=0> y_scale;
  real atom;
  real<lower=0> rho_prior_alpha;
  real<lower=0> rho_prior_beta;
  real<lower=0> alpha_prior_shape;
  real<lower=0> alpha_prior_rate;
  real mean_prior_meanlog;
  real<lower=0> mean_prior_sdlog;
  real shape_prior_meanlog;
  real<lower=0> shape_prior_sdlog;
}

parameters {
  vector<lower=0, upper=1>[H - 1] v[2];
  vector<lower=0>[H] mean_comp[2];
  vector<lower=0>[H] shape_comp[2];
  vector<lower=0>[2] alpha;
  vector<lower=0, upper=1>[2] rho;
}

transformed parameters {
  simplex[H] w[2];
  vector<lower=0, upper=1>[2] pi;

  for (r in 1:2) {
    real remaining;
    remaining = 1;

    for (h in 1:(H - 1)) {
      w[r][h] = v[r][h] * remaining;
      remaining = remaining * (1 - v[r][h]);
    }

    w[r][H] = remaining;
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
    mean_comp[r] ~ lognormal(mean_prior_meanlog, mean_prior_sdlog);
    shape_comp[r] ~ lognormal(shape_prior_meanlog, shape_prior_sdlog);
  }

  for (n in 1:N) {
    A[n] ~ bernoulli(pi[arm[n]]);
  }

  for (n in 1:N_obs) {
    vector[H] component_lp;

    for (h in 1:H) {
      component_lp[h] =
        log(w[arm_obs[n]][h]) +
        gamma_lpdf(
          y_obs_scaled[n] |
          shape_comp[arm_obs[n]][h],
          shape_comp[arm_obs[n]][h] / mean_comp[arm_obs[n]][h]
        );
    }

    target += log_sum_exp(component_lp);
  }
}

generated quantities {
  real rho_0;
  real rho_1;
  real pi_0;
  real pi_1;
  real mu_0_c_scaled;
  real mu_1_c_scaled;
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
  mu_0_c_scaled = dot_product(w[1], mean_comp[1]);
  mu_1_c_scaled = dot_product(w[2], mean_comp[2]);
  mu_0_c = y_scale * mu_0_c_scaled;
  mu_1_c = y_scale * mu_1_c_scaled;
  delta_atom = rho_1 - rho_0;
  mu_delta = mu_1_c - mu_0_c;
  alpha_delta = (pi_1 / (1 - pi_1)) / (pi_0 / (1 - pi_0));
  delta =
    (atom * rho_1 + pi_1 * mu_1_c) -
    (atom * rho_0 + pi_0 * mu_0_c);
}
