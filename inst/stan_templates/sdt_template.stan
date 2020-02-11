// -*- coding: utf-8 -*-

// Every kind of parameter is represented using the type vector, and
// so the effects are represented by a matrix, which may have more
// than one row (if the vector of parameters has length > 1).

data {
  int<lower=0,upper=1> PRINT;
  int<lower=1> N;
  int<lower=2> K;
  // Kb2 = K/2 is here to avoid the (irrelevant) warning about integer
  // division
  int<lower=1> Kb2;
  // for the parsimonious link function
  vector[K-1] unbiased;
  real criteria_scale;
  vector[N] stim_sign;
  int<lower=0> counts[N, K];
  // this will be replaced with delta_size, gamma_size and possibly theta_size
  int<lower=1> PAR_size;
  int<lower=1> dprim_size;
  int<lower=1> X_PAR_ncol;
  row_vector[X_PAR_ncol] X_PAR[N];
  matrix[PAR_size, X_PAR_ncol] PAR_is_fixed;
  matrix[PAR_size, X_PAR_ncol] PAR_fixed_value;
  int<lower=1> group_%_size; //common
  int<lower=1,upper=group_%_size> group_%[N]; //common
  int<lower=1> Z_PAR_ncol_%; //PAR
  row_vector[Z_PAR_ncol_%] Z_PAR_%[N]; //PAR
  // Priors
  matrix[PAR_size, X_PAR_ncol] PAR_fixed_mu;
  row_vector<lower=0>[X_PAR_ncol] PAR_fixed_sd[PAR_size];
  real<lower=1> lkj_PAR_nu_%; //PAR
  vector<lower=0>[PAR_size * Z_PAR_ncol_%] PAR_sd_scale_%; //PAR
}

parameters {
  matrix[PAR_size, X_PAR_ncol] PAR_fixed;
  cholesky_factor_corr[PAR_size * Z_PAR_ncol_%] L_corr_PAR_%; //PAR
  vector<lower=0>[PAR_size * Z_PAR_ncol_%] PAR_sd_%; //PAR
  // Random effects vectors are converted to PAR_size x Z_PAR_ncol
  // matrix in column major order, so for example
  // gamma_1,...,gamma_K-1 for the first column (effect) of the
  // Z_gamma matrix, then gamma_1,...,gamma_K-1 for the second column
  // (effect), etc.
  vector[PAR_size * Z_PAR_ncol_%] PAR_z_%[group_%_size]; //PAR
}

transformed parameters {
  matrix[PAR_size, X_PAR_ncol] PAR_fixed_;
  matrix[PAR_size, Z_PAR_ncol_%] PAR_random_%[group_%_size]; //PAR
  matrix[PAR_size * Z_PAR_ncol_%, PAR_size * Z_PAR_ncol_%] Corr_PAR_%; //PAR
  vector[PAR_size] PAR;
  vector[K - 1] criteria;
  vector[K + 1] multinomial_cum;
  vector[K] multinomial_p[N];
  // used only in the sdt model family
  vector[dprim_size] dprim;
  real shift;
  // used only in the uvsdt model
  real sd_ratio;
  // used only in the metad model
  vector[2] normalization;
  // fixing fixed effects if requested
  for(i in 1:PAR_size)for(j in 1:X_PAR_ncol)if(PAR_is_fixed[i, j] == 1){ PAR_fixed_[i, j] = PAR_fixed_value[i, j]; }else{ PAR_fixed_[i, j] = PAR_fixed[i, j]; }
  Corr_PAR_% = L_corr_PAR_% * L_corr_PAR_%'; //PAR
  for(g in 1:group_%_size)PAR_random_%[g] = to_matrix(diag_pre_multiply(PAR_sd_%, L_corr_PAR_%) * PAR_z_%[g], PAR_size, Z_PAR_ncol_%); //PAR
  if(PRINT == 1){
    print("PRIORS: ");
    print("PAR_fixed_mu "); for(i in 1:PAR_size)print(PAR_fixed_mu[i,]);
    print("PAR_fixed_sd"); for(i in 1:PAR_size)print(PAR_fixed_sd[i,]);
    print("PAR_sd_scale_% = ", PAR_sd_scale_%); //PAR
    print("INITIAL VALUES: ");
    print("PAR_fixed"); for(i in 1:PAR_size)print(PAR_fixed[i,]);
    for(g in 1:group_%_size)print("PAR_z_%[", g, "] = ", PAR_z_%[g]); //PAR
    print("PAR_sd_% = ", PAR_sd_%); //PAR
  }
  for(n in 1:N){
    PAR = PAR_fixed_ * X_PAR[n]';
    PAR = PAR + PAR_random_%[group_%[n]] * Z_PAR_%[n]';  //PAR

    //link-gamma

    //likelihood
  }
}

model {
  for(i in 1:PAR_size)for(j in 1:X_PAR_ncol)PAR_fixed[i, j] ~ normal(PAR_fixed_mu[i, j], PAR_fixed_sd[i, j]);
  L_corr_PAR_% ~ lkj_corr_cholesky(lkj_PAR_nu_%); //PAR
  for(i in 1:(PAR_size * Z_PAR_ncol_%)){ PAR_sd_%[i] ~ cauchy(0, PAR_sd_scale_%[i]); } //PAR
  for(g in 1:group_%_size){ PAR_z_%[g] ~ normal(0, 1); } //PAR
  for(n in 1:N)
    counts[n] ~ multinomial(multinomial_p[n]);
}

generated quantities{
  int<lower=0> counts_new[N, K];
  for(n in 1:N)
    counts_new[n] = multinomial_rng(multinomial_p[n], sum(counts[n]));
}
