// -*- coding: utf-8 -*-

// Zakładamy, że meta-d powoduje {roz/z}sunięcie rozkładów po każdej
// stronie głównego kryterium. meta-d będzie w ogólnym przypadku
// skorelowany z d' więc te dwa d'my modelujemy jak kryteria.

data {
  int<lower=1> N;
  int<lower=2> K;
  real criteria_scale;
  int<lower=1> X_delta_ncol;
  row_vector[X_delta_ncol] X_delta[N];
  int<lower=1> X_gamma_ncol;
  row_vector[X_gamma_ncol] X_gamma[N];
  vector[N] stim_sign; // = -1, 1 zależnie od bodźca
  int<lower=0> counts[N, K];
  int<lower=1> G_%; //common
  int<lower=1,upper=G_%> group_%[N]; //common
  int<lower=1> Z_delta_ncol_%; //delta
  int<lower=1> Z_gamma_ncol_%; //gamma
  row_vector[Z_delta_ncol_%] Z_delta_%[N]; //delta
  row_vector[Z_gamma_ncol_%] Z_gamma_%[N]; //gamma
  // Priors
  matrix[2, X_delta_ncol] delta_fixed_mu;
  row_vector<lower=0>[X_delta_ncol] delta_fixed_sd[2];
  matrix[K - 1, X_gamma_ncol] gamma_fixed_mu;
  row_vector<lower=0>[X_gamma_ncol] gamma_fixed_sd[K - 1];
  real<lower=1> lkj_delta_nu_%; //delta
  real<lower=1> lkj_gamma_nu_%; //gamma
  vector<lower=0>[2 * Z_delta_ncol_%] delta_sd_scale_%; //delta
  vector<lower=0>[(K - 1) * Z_gamma_ncol_%] gamma_sd_scale_%; //gamma
}

parameters {
  matrix[2, X_delta_ncol] delta_fixed;
  matrix[K - 1, X_gamma_ncol] gamma_fixed;
  cholesky_factor_corr[2 * Z_delta_ncol_%] L_corr_delta_%; //delta
  vector<lower=0>[2 * Z_delta_ncol_%] delta_sd_%; //delta
  vector[2 * Z_delta_ncol_%] delta_z_%[G_%]; //delta
  cholesky_factor_corr[(K - 1) * Z_gamma_ncol_%] L_corr_gamma_%; //gamma
  vector<lower=0>[(K - 1) * Z_gamma_ncol_%] gamma_sd_%; //gamma
  vector[(K - 1) * Z_gamma_ncol_%] gamma_z_%[G_%]; //gamma
}

transformed parameters {
  vector[X_delta_ncol] dprim_fixed;
  vector[X_delta_ncol] mratio_fixed;
  matrix[2, Z_delta_ncol_%] delta_random_%[G_%]; //delta
  matrix[K - 1, Z_gamma_ncol_%] gamma_random_%[G_%]; //gamma
  matrix[2 * Z_delta_ncol_%, 2 * Z_delta_ncol_%] Corr_delta_%; //delta
  matrix[(K - 1) * Z_gamma_ncol_%, (K - 1) * Z_gamma_ncol_%] Corr_gamma_%; //gamma
  vector[K-1] criteria[N];
  vector[2] dprim[N];
  vector[N] shift;
  vector[K] multinomial_p[N];
  vector[2] normalization[N];
  for(i in 1:X_delta_ncol){
    dprim_fixed[i] = exp(delta_fixed[1, i]);
    mratio_fixed[i] = exp(delta_fixed[2, i]) / dprim_fixed[i];
  }
  Corr_delta_% = L_corr_delta_% * L_corr_delta_%'; //delta
  Corr_gamma_% = L_corr_gamma_% * L_corr_gamma_%'; //gamma
  for(g in 1:G_%){ delta_random_%[g] = to_matrix(diag_pre_multiply(delta_sd_%, L_corr_delta_%) * delta_z_%[g], 2, Z_delta_ncol_%); } //delta
  for(g in 1:G_%){ gamma_random_%[g] = to_matrix(diag_pre_multiply(gamma_sd_%, L_corr_gamma_%) * gamma_z_%[g], K - 1, Z_gamma_ncol_%); } //gamma
  for(n in 1:N){
    dprim[n] = exp(delta_fixed * X_delta[n]'
                   + delta_random_%[group_%[n]] * Z_delta_%[n]'  //delta
                   );
    criteria[n] = criteria_scale * inv_Phi(head(cumulative_sum(softmax(append_row(gamma_fixed * X_gamma[n]'
                                                                   + gamma_random_%[group_%[n]] * Z_gamma_%[n]' //gamma
                                                                   , 0))),
                                 K - 1));
  }
  shift = 0.5 * stim_sign; // dla stim = 1 mamy shift = -0.5
  for(n in 1:N){
    normalization[n, 1] = Phi(criteria[n, K/2] - shift[n] * dprim[n, 1]) / Phi(criteria[n, K/2] - shift[n] * dprim[n, 2]);
    normalization[n, 2] = Phi(-(criteria[n, K/2] - shift[n] * dprim[n, 1])) / Phi(-(criteria[n, K/2] - shift[n] * dprim[n, 2]));
    multinomial_p[n, 1] = Phi(criteria[n, 1] - shift[n] * dprim[n, 2]) * normalization[n, 1];
    for(k in 2:(K - 1))
      if(k < (K / 2 + 1)){
        multinomial_p[n, k] = (Phi(criteria[n, k] - shift[n] * dprim[n, 2]) - Phi(criteria[n, k - 1] - shift[n] * dprim[n, 2])) * normalization[n, 1];
      }else{
        multinomial_p[n, k] = (Phi(criteria[n, k] - shift[n] * dprim[n, 2]) - Phi(criteria[n, k - 1] - shift[n] * dprim[n, 2])) * normalization[n, 2];
      }
    multinomial_p[n, K] = Phi(-(criteria[n, K - 1] - shift[n] * dprim[n, 2])) * normalization[n, 2];
  }
}

model {
  for(i in 1:2)
    for(j in 1:X_delta_ncol)
      delta_fixed[i, j] ~ normal(delta_fixed_mu[i, j], delta_fixed_sd[i, j]); //0, 2
  for(i in 1:(K - 1))
    for(j in 1:X_gamma_ncol)
      gamma_fixed[i, j] ~ normal(gamma_fixed_mu[i, j], gamma_fixed_sd[i, j]); //0, 5
  L_corr_delta_% ~ lkj_corr_cholesky(lkj_delta_nu_%); //delta
  L_corr_gamma_% ~ lkj_corr_cholesky(lkj_gamma_nu_%); //gamma
  for(i in 1:(2 * Z_delta_ncol_%)){ delta_sd_%[i] ~ cauchy(0, delta_sd_scale_%[i]); } //delta
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
