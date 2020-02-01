// -*- coding: utf-8 -*-

data {
  int<lower=1> N;
  int<lower=2> K;
  // Kb2 = K/2 is here to avoid the (irrelevant) warning about integer
  // division
  int<lower=1> Kb2;
  real criteria_scale;
  int<lower=1> X_delta_ncol;
  row_vector[X_delta_ncol] X_delta[N];
  int<lower=1> X_gamma_ncol;
  row_vector[X_gamma_ncol] X_gamma[N];
  vector[N] stim_sign;
  int<lower=0> counts[N, K];
  int<lower=1> G_%; //common
  int<lower=1,upper=G_%> group_%[N]; //common
  int<lower=1> Z_delta_ncol_%; //delta
  int<lower=1> Z_gamma_ncol_%; //gamma
  row_vector[Z_delta_ncol_%] Z_delta_%[N]; //delta
  row_vector[Z_gamma_ncol_%] Z_gamma_%[N]; //gamma
  // Priors
  vector[X_delta_ncol] delta_fixed_mu;
  vector<lower=0>[X_delta_ncol] delta_fixed_sd;
  matrix[K - 1, X_gamma_ncol] gamma_fixed_mu;
  row_vector<lower=0>[X_gamma_ncol] gamma_fixed_sd[K - 1];
  real<lower=1> lkj_delta_nu_%; //delta
  real<lower=1> lkj_gamma_nu_%; //gamma
  vector<lower=0>[Z_delta_ncol_%] delta_sd_scale_%; //delta
  vector<lower=0>[(K - 1) * Z_gamma_ncol_%] gamma_sd_scale_%; //gamma
}

parameters {
  vector[X_delta_ncol] delta_fixed;
  matrix[K - 1, X_gamma_ncol] gamma_fixed;
  cholesky_factor_corr[Z_delta_ncol_%] L_corr_delta_%; //delta
  vector<lower=0>[Z_delta_ncol_%] delta_sd_%; //delta
  vector[Z_delta_ncol_%] delta_z_%[G_%]; //delta
  cholesky_factor_corr[(K - 1) * Z_gamma_ncol_%] L_corr_gamma_%; //gamma
  vector<lower=0>[(K - 1) * Z_gamma_ncol_%] gamma_sd_%; //gamma
  vector[(K - 1) * Z_gamma_ncol_%] gamma_z_%[G_%]; //gamma
}

transformed parameters {
  vector[Z_delta_ncol_%] delta_random_%[G_%]; //delta
  matrix[K - 1, Z_gamma_ncol_%] gamma_random_%[G_%]; //gamma
  matrix[Z_delta_ncol_%, Z_delta_ncol_%] Corr_delta_%; //delta
  matrix[(K - 1) * Z_gamma_ncol_%, (K - 1) * Z_gamma_ncol_%] Corr_gamma_%; //gamma
  vector[K-1] gamma;
  vector[K-1] criteria;
  real delta;
  real dprim;
  real distr_shift;
  vector[K] multinomial_p[N];
  Corr_delta_% = L_corr_delta_% * L_corr_delta_%'; //delta
  Corr_gamma_% = L_corr_gamma_% * L_corr_gamma_%'; //gamma
  for(g in 1:G_%){ delta_random_%[g] = diag_pre_multiply(delta_sd_%, L_corr_delta_%) * delta_z_%[g]; } //delta
  for(g in 1:G_%){ gamma_random_%[g] = to_matrix(diag_pre_multiply(gamma_sd_%, L_corr_gamma_%) * gamma_z_%[g], K - 1, Z_gamma_ncol_%); } //gamma
  for(n in 1:N){
    delta = X_delta[n] * delta_fixed
                   + Z_delta_%[n] * delta_random_%[group_%[n]] //delta
      ;
    gamma = gamma_fixed * X_gamma[n]'
      + gamma_random_%[group_%[n]] * Z_gamma_%[n]' //gamma
      ;
    dprim = exp(delta); //link-delta
    criteria = criteria_scale * inv_Phi(head(cumulative_sum(softmax(append_row(gamma, 0))), K - 1)); //link-gamma
    distr_shift = -0.5 * (stim_sign[n] * dprim);
    multinomial_p[n, 1] = Phi(criteria[1] + distr_shift);
    for(k in 2:(K - 1))
      multinomial_p[n, k] = Phi(criteria[k] + distr_shift) - Phi(criteria[k - 1] + distr_shift);
    multinomial_p[n, K] = Phi(-(criteria[K - 1] + distr_shift));
  }
}

model {
  for(i in 1:X_delta_ncol)
    delta_fixed[i] ~ normal(delta_fixed_mu[i], delta_fixed_sd[i]); //0, 2
  for(i in 1:(K - 1))
    for(j in 1:X_gamma_ncol)
      gamma_fixed[i, j] ~ normal(gamma_fixed_mu[i, j], gamma_fixed_sd[i, j]); //0, 5
  L_corr_delta_% ~ lkj_corr_cholesky(lkj_delta_nu_%); //delta
  L_corr_gamma_% ~ lkj_corr_cholesky(lkj_gamma_nu_%); //gamma
  for(i in 1:Z_delta_ncol_%){ delta_sd_%[i] ~ cauchy(0, delta_sd_scale_%[i]); } //delta
  for(i in 1:((K - 1) * Z_gamma_ncol_%)){ gamma_sd_%[i] ~ cauchy(0, gamma_sd_scale_%[i]); } //gamma
  for(g in 1:G_%){ delta_z_%[g] ~ normal(0, 1); } //delta
  for(g in 1:G_%){ gamma_z_%[g] ~ normal(0, 1); } //gamma
  for(n in 1:N)
    counts[n] ~ multinomial(multinomial_p[n]);
}

generated quantities{
  int<lower=0> counts_new[N, K];
  for(n in 1:N)
    counts_new[n] = multinomial_rng(multinomial_p[n], sum(counts[n]));
}
