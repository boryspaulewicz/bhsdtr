## -*- coding: utf-8 -*-

#' Calculates delta = log(d') given accuracy assuming no bias and p(stim = 1) = 0.5.
#' @param acc a vector of average accuracy scores. All the values must
#'     be greater than 0 and lower than 1.
#' @return delta = log(d')
#' @export
acc_to_delta = function(acc){
    if(any((acc <= 0) || (acc >= 1)))
        stop('Some accuracy scores were to extreme')
    log(2 * stats::qnorm(acc))
}
## ok
