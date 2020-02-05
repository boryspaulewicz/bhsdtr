// UVSDT likelihood
sd_ratio = exp(theta[1]); //link-theta
shift = -0.5 * (stim_sign[n] * dprim[1]);
if(stim_sign[n] > 0){
  multinomial_p[n, 1] = Phi((criteria[1] + shift) / sd_ratio);
  for(k in 2:gamma_size)
    multinomial_p[n, k] = Phi((criteria[k] + shift) / sd_ratio) - Phi((criteria[k - 1] + shift) / sd_ratio);
  multinomial_p[n, K] = Phi(-(criteria[gamma_size] + shift) / sd_ratio);
 }else{
  multinomial_p[n, 1] = Phi(criteria[1] + shift);
  for(k in 2:(gamma_size))
    multinomial_p[n, k] = Phi(criteria[k] + shift) - Phi(criteria[k - 1] + shift);
  multinomial_p[n, K] = Phi(-(criteria[gamma_size] + shift));
 }
// UVSDT likelihood end
