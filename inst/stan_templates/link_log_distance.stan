// 'log_distance' gamma link function
criteria[Kb2] = gamma[Kb2]; // main criterion i.e. bias
if(K > 2){
  for(k in 1:(Kb2 - 1)){
    criteria[Kb2 + k] = criteria[Kb2 + k - 1] + exp(gamma[Kb2 + k]);
    criteria[Kb2 - k] = criteria[Kb2 - k + 1] - exp(gamma[Kb2 - k]);
  }
 }