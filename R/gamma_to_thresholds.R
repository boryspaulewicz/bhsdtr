## -*- coding: utf-8 -*-

#' Transforms gamma posterior samples or contrasts of gamma posterior samples to thresholds
#'
#' @export
gamma_to_thresholds = function(fit, ..., s = 2, K = NULL){
    contrasts = unlist(list(...))
    gamma = list()
    samples = extract(fit)$gamma_fixed
    for(i in 1:dim(samples)[3])
        gamma[[sprintf('m%d', i)]] = samples[,, i]
    if(is.null(contrasts))
        contrasts = paste('m', 1:length(gamma), sep = '')
    res = list()
    for(i in 1:length(contrasts))
        res[[contrasts[i]]] = with(gamma, eval(rlang::parse_expr(contrasts[i])))
    for(i in 1:length(res))
    res[[i]] = gamma_to_crit(res[[i]], NULL, extract_link(fit, 'gamma'), s, K)
    res
}
