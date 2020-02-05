// SDT likelihood
shift = -0.5 * (stim_sign[n] * dprim[1]);
multinomial_p[n, 1] = Phi(criteria[1] + shift);
for(k in 2:gamma_size)
  multinomial_p[n, k] = Phi(criteria[k] + shift) - Phi(criteria[k - 1] + shift);
multinomial_p[n, K] = Phi(-(criteria[gamma_size] + shift));
// SDT likelihood end
