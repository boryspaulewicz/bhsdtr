// gamma link twoparameter
criteria[Kb2] = gamma[1]; // main criterion i.e. bias
if(K > 2){
  for(k in 1:(Kb2 - 1)){
    criteria[Kb2 + k] = criteria[Kb2 + k - 1] + exp(gamma[2]);
    criteria[Kb2 - k] = criteria[Kb2 - k + 1] - exp(gamma[2]);
  }
 }
// gamma link end
