// 'parsimonious' gamma link function
for(k in 1:(K - 1)){
  criteria[k] = gamma[1] + exp(gamma[2]) * unbiased[k];
 }
