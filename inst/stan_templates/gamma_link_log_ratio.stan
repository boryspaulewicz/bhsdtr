// gamma link function log_ratio
criteria[Kb2] = gamma[Kb2]; // main criterion, i.e., bias
if(K > (Kb2 + 1)){
  // spread
  criteria[Kb2 + 1] = criteria[Kb2] + exp(gamma[Kb2 + 1]);
  if(Kb2 > 1){
    // symmetry
    criteria[Kb2 - 1] = criteria[Kb2] - exp(gamma[Kb2 - 1]) * (criteria[Kb2 + 1] - criteria[Kb2]);
    if(K > 4){
      for(k in 1:(Kb2 - 2)){
        // upper consistency
        criteria[Kb2 + k + 1] = criteria[Kb2 + k] + exp(gamma[Kb2 + k + 1]) * (criteria[Kb2 + 1] - criteria[Kb2]);
        // lower consistency
        criteria[Kb2 - k - 1] = criteria[Kb2 - k] - exp(gamma[Kb2 - k - 1]) * (criteria[Kb2] - criteria[Kb2 - 1]);
      }
    }
  }
 }
// gamma_link end
