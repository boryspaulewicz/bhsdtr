// meta-d' likelihood
shift = 0.5 * stim_sign[n]; // when stim == 1 shift == -0.5
normalization[1] = Phi(criteria[Kb2] - shift * dprim[1]) / Phi(criteria[Kb2] - shift * dprim[2]);
normalization[2] = Phi(-(criteria[Kb2] - shift * dprim[1])) / Phi(-(criteria[Kb2] - shift * dprim[2]));
multinomial_p[n, 1] = Phi(criteria[1] - shift * dprim[2]) * normalization[1];
for(k in 2:(K - 1))
  if(k < (Kb2 + 1)){
    multinomial_p[n, k] = (Phi(criteria[k] - shift * dprim[2]) - Phi(criteria[k - 1] - shift * dprim[2])) * normalization[1];
  }else{
    multinomial_p[n, k] = (Phi(criteria[k] - shift * dprim[2]) - Phi(criteria[k - 1] - shift * dprim[2])) * normalization[2];
  }
multinomial_p[n, K] = Phi(-(criteria[K - 1] - shift * dprim[2])) * normalization[2];
// meta-d' likelihood end
