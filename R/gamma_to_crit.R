## -*- coding: utf-8 -*-

#' Transforms gamma posterior samples to midpoint centered criteria.
#'
#' \code{gamma_to_crit} transforms gamma posterior samples to midpoint
#' centered criteria. This will work correctly only if the column of
#' the gamma fixed effects parameter matrix specified by
#' \code{beta_index} represents values of gamma for a given condition,
#' not the regression slope or the difference in gamma values between
#' conditions (see also 'Details'). WARNING: The same gamma link
#' function as the one used when specifying the model has to be used.
#'
#' In general, in an SDT model there is only one sensitivity (there
#' are two or even three in a meta-d' model) but there can be more
#' than one criterion. When there is more than one criterion the
#' regression coefficients that represent the dependence of criteria
#' on additional variables form a matrix, not a vector. Row j of this
#' fixed effects matrix contains regression coefficients for the j-th
#' criterion. Let's say we have data from a simple experiment with
#' experimental and control groups and we want to test if criteria
#' depend on group membership. If we specify the fixed effects model
#' matrix by calling \code{make_stan_data} with \code{fixed =
#' list(gamma = ~ condition, ...)} then the first column of the
#' \code{gamma_fixed} parameter matrix will contain the intercepts for
#' every criterion (in the gamma space), which is the value of gamma
#' in the condition chosen as base level, and the second column will
#' contain the difference in gamma values between experimental and
#' control groups. In this case \code{gamma_to_crit(samples, 1)} will
#' work, but \code{gamma_to_crit(samples, 2)} will not work. If, on
#' the other hand, we wanted to work with the actual criteria, not the
#' gamma values from which criteria are derived, then \code{fixed =
#' list(gamma = ~ -1 + condition, ...)} would allow us to apply the
#' \code{gamma_to_crit} function to both columns of the
#' \code{gamma_fixed} matrix because \code{~ -1 + condition} results
#' in a separate intercepts and slopes parametrization of the gamma
#' fixed effects model matrix. This way we could easily estimate the
#' criteria for both conditions.
#'
#' @param samples a data frame of posterior samples from stanfit
#'     object (e.g., samples = as.data.frame(stanfit)).
#' @param beta_index a fixed effect index (a scalar) = an index of a
#'     \code{gamma_fixed} matrix column. If only the relevant samples
#'     are provided (an appropriate gamma matrix) using beta_index =
#'     NULL will convert the provided samples. See also 'Details' for
#'     which columns of the \code{gamma_fixed} matrix can and which
#'     cannot be converted to criteria using this function.
#' @param gamma_link is either 'softmax' (described in the paper),
#'     'log_distance', 'log_ratio', or 'parsimonious' (See the Readme
#'     file in the github repository)
#' @param s a criteria scaling factor (a scalar), which is relevant
#'     only for the softmax link function. WARNING: must be equal to
#'     the value used when fitting the model.
#' @return a data frame containing posterior samples of midpoint
#'     centered criteria.
#' @examples
#' ## First we apply the inverse of what happens inside the gamma_to_crit function
#' criteria = c(-2, -1, 0, 1, 2)
#' cumprobs = pnorm(criteria)
#' cumprobs = c(pnorm(criteria, sd = 2), 1)
#' areas = c(cumprobs[1], cumprobs[-1] - cumprobs[-length(cumprobs)])
#' gamma = log(areas / areas[length(areas)])
#' ## than we apply the simplified version of gamma_to_crit (link softmax)
#' g_to_c = function(x, s = 2){
#'   x = exp(x)
#'   s * qnorm(cumsum(x/sum(x))[-length(x)])
#' }
#' ## so that we can see that all is well...
#' rbind(g_to_c(gamma), criteria)
#' @export
gamma_to_crit = function(samples, beta_index = 1, gamma_link = 'softmax', s = 2, K = NULL){
    check_link(gamma_link)
    link = gamma_link[1]
    if(is.null(beta_index)){
        nms = 1:ncol(samples)
        beta_index = 1
    }else{
        if (length(beta_index) != 1) 
            warning("Using only the first element of the beta_index vector")
        nms = grep(sprintf("gamma_fixed\\[[0-9]+,%d]", beta_index[1]), 
                   names(samples))
    }
    if (length(nms) == 0) 
        stop(sprintf("Could not find gamma_fixed[.,%d] samples", 
                     beta_index[1]))
    if(link %in% c('parsimonious', 'twoparameter')){
        if(is.null(K))
            stop("K has to be specified for the parsimonious or the twoparameter link functions")
    }else{
        K = length(nms) + 1
    }
    ## as.matrix because if samples is a vector nrow is NULL
    samples = as.matrix(samples[, nms])
    criteria = matrix(nrow = nrow(samples), ncol = K - 1)
    if(link == 'log_ratio'){
        criteria[, K / 2] = samples[, K / 2];
        if(K > 2){
            ## spread
            criteria[, K / 2 + 1] = criteria[, K / 2] + exp(samples[, K / 2 + 1]);
            ## symmetry
            criteria[, K / 2 - 1] = criteria[, K / 2] - exp(samples[, K / 2 - 1]) * (criteria[, K / 2 + 1] - criteria[, K / 2]);
            if(K > 4){
                for(k in 1:(K / 2 - 2)){
                    ## upper consistency
                    criteria[, K / 2 + k + 1] = criteria[, K / 2 + k] + exp(samples[, K / 2 + k + 1]) * (criteria[, K / 2 + 1] - criteria[, K / 2]);
                    ## lower consistency
                    criteria[, K / 2 - k - 1] = criteria[, K / 2 - k] - exp(samples[, K / 2 - k - 1]) * (criteria[, K / 2] - criteria[, K / 2 - 1]);
                }
            }
        }
    }
    if(link == 'log_distance'){
        criteria[, K / 2] = samples[, K / 2];
        if(K > 2){
            for(k in 1:(K / 2 - 1)){
                criteria[, K / 2 + k] = criteria[, K / 2 + k - 1] + exp(samples[, K / 2 + k]);
                criteria[, K / 2 - k] = criteria[, K / 2 - k + 1] - exp(samples[, K / 2 - k]);
            }
        }
    }
    if(link == 'parsimonious'){
        unb = unbiased(K)
        for(k in 1:(K - 1))
            criteria[, k] = samples[, 1] + exp(samples[, 2]) * unb[k]
    }
    if(link == 'twoparameter'){
        criteria[, K / 2] = samples[, 1];
        if(K > 2){
            for(k in 1:(K / 2 - 1)){
                criteria[, K / 2 + k] = criteria[, K / 2 + k - 1] + exp(samples[, 2]);
                criteria[, K / 2 - k] = criteria[, K / 2 - k + 1] - exp(samples[, 2]);
            }
        }
    }
    if(link == 'softmax'){
        criteria = t(apply(exp(cbind(samples, 0)), 1,
                           function(x) s * stats::qnorm(cumsum(x/sum(x))[-length(x)])))
    }
    if(link == 'identity'){
        warning("You want me to convert between identical gamma and criteria vectors. Are you shure you should be sciencing?")
        samples
    }
    colnames(criteria) = gsub("gamma", "criteria", colnames(criteria))
    criteria
}

check_link = function(gamma_link, delta_link = NULL){
    gamma_links = c('softmax', 'log_ratio', 'log_distance', 'parsimonious', 'twoparameter', 'identity')
    if(!(gamma_link %in% gamma_links))
        stop(paste("gamma_link must be one of the following:", gamma_links))
    if(!is.null(delta_link)){
        delta_links = c('log', 'log_ratio', 'identity')
        if(!(delta_link %in% delta_links))
            stop(paste("delta_link must be one of the following:", delta_links))
    }
}

check_model = function(model){
    models = c('sdt', 'uvsdt', 'ordinal', 'uvordinal', 'metad')
    if(!(model %in% models))
        stop(paste("Model must be one of the following:", models))
}

unbiased = function(K){
    log((1:K / K) / (1 - 1:K / K))[1:(K - 1)]
}
