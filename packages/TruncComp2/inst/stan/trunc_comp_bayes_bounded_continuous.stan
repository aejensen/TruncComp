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
  real<lower=0> m_prior_alpha;
  real<lower=0> m_prior_beta;
  real phi_prior_meanlog;
  real<lower=0> phi_prior_sdlog;
}

parameters {
  vector<lower=0, upper=1>[H - 1] v[2];
  vector<lower=1e-6, upper=1 - 1e-6>[H] m_comp[2];
  vector<lower=1e-6>[H] phi_comp[2];
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
    m_comp[r] ~ beta(m_prior_alpha, m_prior_beta);
    phi_comp[r] ~ lognormal(phi_prior_meanlog, phi_prior_sdlog);
  }

  for (n in 1:N) {
    A[n] ~ bernoulli(pi[arm[n]]);
  }

  for (n in 1:N_obs) {
    vector[H] component_lp;

    for (h in 1:H) {
      real beta_a;
      real beta_b;
      beta_a = m_comp[arm_obs[n]][h] * phi_comp[arm_obs[n]][h];
      beta_b = (1 - m_comp[arm_obs[n]][h]) * phi_comp[arm_obs[n]][h];
      component_lp[h] =
        log(w[arm_obs[n]][h]) +
        beta_lpdf(x_obs[n] | beta_a, beta_b);
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
  mu_0_c_unit = dot_product(w[1], m_comp[1]);
  mu_1_c_unit = dot_product(w[2], m_comp[2]);
  mu_0_c = score_min + (score_max - score_min) * mu_0_c_unit;
  mu_1_c = score_min + (score_max - score_min) * mu_1_c_unit;
  delta_atom = rho_1 - rho_0;
  mu_delta = mu_1_c - mu_0_c;
  alpha_delta = (pi_1 / (1 - pi_1)) / (pi_0 / (1 - pi_0));
  delta =
    (atom * rho_1 + pi_1 * mu_1_c) -
    (atom * rho_0 + pi_0 * mu_0_c);
}
