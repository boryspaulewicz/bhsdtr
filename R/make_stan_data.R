## -*- coding: utf-8 -*-

#' Creates the list of data structures required by the stan function.
#'
#' \code{make_stan_data} creates the list of data structures required by the
#' \code{stan} function when fitting the model generated by
#' \code{make_stan_model}.
#'
#' \code{fixed} must be a list of model formulae. It may also contain prior
#' parameter values, if non-default priors on fixed effects are required. This
#' list is composed of the following elements:
#'
#' \describe{
#'
#' \item{delta}{is a model formula that defines the delta fixed effects model
#' matrix, e.g., \code{delta = ~ condition}}
#'
#' \item{gamma}{is a model formula that defines the gamma fixed effects model
#' matrix.}
#'
#' \item{delta_mu}{is na optional vector specifying means of independent normal
#' priors on delta fixed effects. This must be of the same length as the number
#' of delta fixed effects or of length 1, in which case the same value will be
#' used for every delta fixed effect. The default value is
#' \code{acc_to_delta(.75)}}
#'
#' \item{delta_sd}{is na optional vector specifying standard deviations of
#' independent normal priors on delta fixed effects. This must be of the same
#' length as the number of delta fixed effects or of length 1, in which case the
#' same value will be used for every delta fixed effect. The default value is
#' \code{.5 * (acc_to_delta(.99) - acc_to_delta(.51))}}
#'
#' \item{gamma_mu}{is na optional vector specifying means of
#' independent normal priors on gamma fixed effects. The default value
#' is 0. This must be of the same length as the number of gamma fixed
#' effects or of length 1, in which case the same value will be used
#' for every gamma fixed effect. Note that when there is more than one
#' criterion the number of gamma fixed effects is equal to the number
#' I of columns of gamma fixed effects model matrix times the number
#' of criteria K-1.  Internally, the \code{gamma_mu} vector is stored
#' as a K-1xI matrix in column major order: element (K-1)*j+i of the
#' \code{gamma_mu} vector (where 0 < j and 0 < i < K-1) represents
#' fixed effect corresponding to the i-th column of the gamma fixed
#' effects model matrix and the j-th criterion. In other words, the
#' indices for the fixed effects change faster than the indices for
#' the criteria.}
#'
#' \item{gamma_sd}{is na optional vector specifying standard
#' deviations of independent normal priors on gamma fixed
#' effects. This must be of the same length as the number of gamma
#' fixed effects or of length 1, in which case the same value will be
#' used for every gamma fixed effect. See \code{gamma_mu} above for
#' details on how the elements of the \code{gamma_sd} vector
#' correspond to elements of the gamma fixed effects parameter
#' matrix. The default value depends on the chosen gamma link
#' function. For the 'softmax' function it is log(100), which means
#' that a priori areas under the gamma-to-criterion mapping
#' distribution curve delineated by the criteria are expected to vary
#' by a factor of 100 or less but varying by a factor of
#' \code{exp(2*log(100)) = 10000} is highly unlikely. For the
#' 'log_distance' and the 'log_ratio' gamma link functions the default
#' value is 2.}
#'
#' }
#'
#' \code{random} is on optional list of lists of model formulae and
#' optional prior parameter values, if non-default priors on random
#' effects are required. Each list specifies delta and gamma random
#' effects of one grouping factor. Note that the same grouping factor
#' (e.g., subject id) can be reused in different random effects
#' specification lists to force the corresponding random effects to be
#' uncorrelated, which will lower the number of free parameters but
#' may result in interval estimate bias if the correlations exist and
#' are non- negligible:
#'
#' \describe{
#'
#' \item{group}{is a model formula specifying the random grouping
#' factor, e.g., \code{group = ~ subject} indicates that the random
#' effects specified in this list are associated with the subject
#' grouping factor.}
#'
#' \item{delta}{is a model formula that defines the delta random
#' effects model matrix. This must define a submodel of the delta
#' fixed effects model, e.g., if delta depends on f1 but not on f2 in
#' the model and f1 is a within-subject variable then \code{delta = ~
#' f1} or \code{delta = ~ -1 + f1} is a valid delta random effects
#' specification but e.g., \code{delta = ~ f2} is not.}
#'
#' \item{gamma}{is a model formula that defines the gamma random
#' effects model matrix. This must define a submodel of the gamma
#' fixed effects model (see \code{delta} above for some examples.)}
#'
#' \item{delta_nu}{is an optional vector of lkj prior parameters for
#' delta random effects. This must be of the same length as the number
#' of delta random effects or of length 1, in which case the same
#' value will be used for every delta random effect. The default is 1
#' which corresponds to uniform prior on random effects' correlation
#' matrices. The greater the value of this parameter the more emphasis
#' is put on zer off-diagonal correlations, which represents the a
#' priori assumption that the correlations are low or near zero.}
#'
#' \item{delta_sd_scale}{is an optional vector of half-Cauchy prior parameters
#' for delta random effects. This must be of the same length as the number of
#' delta random effects or of length 1, in which case the same value will be
#' used for every delta random effect. The default value is \code{.5 *
#' (acc_to_delta(.99) - acc_to_delta(.51))}}
#'
#' \item{gamma_nu}{is an optional vector of lkj prior parameters for
#' gamma random effects. This must be of the same length as the number
#' of gamma random effects or of length 1, in which case the same
#' value will be used for every gamma random effect. The default is 1
#' which corresponds to uniform prior on random effects' correlation
#' matrices. The greater the value of this parameter the more emphasis
#' is put on zer off-diagonal correlations. See \code{gamma_mu} above
#' for details on how the elements of the \code{gamma} vector
#' correspond to the elements of the gamma parameter matrix.}
#'
#' \item{gamma_sd_scale}{is an optional vector of half-Cauchy prior parameters
#' for gamma random effects. This must be of the same length as the number of
#' gamma random effects or of length 1, in which case the same value will be
#' used for every gamma random effect. The default value is \code{log(100)}. See
#' \code{gamma_mu} above for details on how the elements of the \code{gamma}
#' vector correspond to elements of the gamma parameter matrix.}
#'
#' }
#'
#' @param adata an aggregated data object created by
#'     \code{aggregate_responses} function.
#' @param fixed a list specifying gamma and delta fixed effects and
#'     priors. See also 'Details'.
#' @param random an optional list specifying gamma and delta random
#'     effects and priors. The default is list(), which corresponds to
#'     a non-hierarchical SDT model. See also 'Details'.
#' @param criteria_scale a scaling factor corresponding to mapping
#'     distribution's standard deviation, applies only to the softmax
#'     gamma link function. The default is 2. See also 'Details'.
#' @param gamma_link either 'softmax' (described in the paper),
#'     'log_distance' or 'log_ratio' (See the Readme file in the
#'     github repository)
#' @param model can be either 'sdt' (the default), 'uvsdt', 'metad', 'ordinal',
#' or 'uvordinal'
#' @return a list with response and stimulus data, model matrices,
#'     prior parameter values, and other data required by the stan
#'     model generated using \code{make_stan_model}.
#' @examples
#' data(gabor)
#' gabor$r = combined_response(gabor$stim, gabor$rating, gabor$acc)
#' adata = aggregate_responses(gabor, 'stim', 'r', c('duration', 'id', 'order'))
#' fixed = list(delta = ~ -1 + duration:order, gamma = ~ order)
#' random = list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1))
#' sdata = make_stan_data(adata, fixed, random)
#' sdata
#' @export
make_stan_data = function(adata, fixed, random = list(), criteria_scale = 2, gamma_link = 'softmax', model = 'sdt', delta_link = 'log'){
    sf = sprintf
    check_fixed(fixed)
    check_random(random)
    check_adata(adata)
    check_link(gamma_link)
    check_model(model)
    if((model %in% c('ordinal', 'uvordinal')) & !(gamma_link %in% c('log_distance', 'twoparameter', 'parsimonious')))
        stop('This link function is not implemented for ordinal models')
    default_prior = list(mu = list(delta = acc_to_delta(.75), theta = 0, gamma = 0, eta = 0),
                         sd = list(delta = .5 * (acc_to_delta(.99) - acc_to_delta(.51)), theta = log(2), eta = 5),
                         scale = list(delta = .5 * (acc_to_delta(.99) - acc_to_delta(.51)), theta = log(2), eta = 5))
    if(gamma_link != 'softmax'){
        default_prior$scale$gamma = default_prior$sd$gamma = 2
    }else{
        default_prior$scale$gamma = default_prior$sd$gamma = log(100)
    }
    K = ncol(adata$counts)
    data = list(PRINT = 0,
                N = nrow(adata$counts),
                K = K,
                Kb2 = round(K / 2),
                fixed = fixed,
                random = random,
                gamma_link = gamma_link,
                model = model,
                criteria_scale = criteria_scale,
                unbiased = fix_stan_dim(unbiased(K)),
                eta_size = 1,
                dprim_size = 1,
                delta_size = 1,
                theta_size = 1,
                gamma_size = K - 1,
                criteria_scale = criteria_scale,
                counts = adata$counts)
    if(model %in% c('sdt', 'uvsdt', 'metad')){
        data$delta_size = data$dprim_size = c('sdt' = 1, 'uvsdt' = 1, 'metad' = 2)[model]
        ## in SDT models stim_sign = -1, 1
        data$stim_sign = 2 * as.numeric(as.factor(as.character(adata$stimulus))) - 3
    }else{
        data$stim_sign = rep(0, nrow(adata$counts))
    }
    if(gamma_link %in% c('parsimonious', 'twoparameter'))
        data$gamma_size = 2
    par_types = c('gamma')
    if(model %in% c('ordinal', 'uvordinal')){
        par_types = c('eta', par_types)
    }else{
        par_types = c('delta', par_types)
    }
    if(model %in% c('uvsdt', 'uvordinal'))
        par_types = c(par_types, 'theta')
    for(par_type in par_types){
        ## Adding fixed effects model matrices and their ncols, e.g., X_delta, X_delta_ncol
        data[[sf('X_%s', par_type)]] = remove.zero.cols(stats::model.matrix(fixed[[par_type]], adata$data))
        data[[sf('X_%s_ncol', par_type)]] = ncol(data[[sf('X_%s', par_type)]])
        ## Fixed effects' priors, e.g., delta_fixed_mu / sd
        for(par in c('mu', 'sd'))
            data[[sf('%s_fixed_%s', par_type, par)]] = parse_prior(fixed[[sf('%s_%s', par_type, par)]],
                                                                   default_prior[[par]][[par_type]],
                                                                   ## matrix dimensions
                                                                   c(data[[sf('%s_size', par_type)]], ncol(data[[sf('X_%s', par_type)]])),
                                                                   sf('%s_%s', par_type, par))
        data[[sf('%s_is_fixed', par_type)]] = matrix(0, nrow = data[[sf('%s_size', par_type)]], ncol = ncol(data[[sf('X_%s', par_type)]]))
        data[[sf('%s_fixed_value', par_type)]] = matrix(0, nrow = data[[sf('%s_size', par_type)]], ncol = ncol(data[[sf('X_%s', par_type)]]))
        ## In ordinal models the effect associated with the intercept is fixed at 0 for identifiability
        if((model %in% c('ordinal', 'uvordinal')) & par_type == 'eta'){
            ## We are fixing one element of the fixed effects vector for identifiability
            warning('The first element of the eta fixed effects vector is fixed at 0 for identifiability')
            data$eta_is_fixed[,1] = 1
            data$eta_fixed_value[,1] = 0
        }
    }
    if((gamma_link == 'identity') & (is.null(fixed$gamma_mu))){
        for(i in 1:(data$X_gamma_ncol))
            data$gamma_fixed_mu[, i] = seq(-2, 2, length.out = data$gamma_size)
    }
    ## Random effects
    if(length(random) > 0){ 
        ## converting the formulae (e.g., ~ id) to numeric indices (1, 2, 3, ..., max(id))        
        for(l in 1:length(random)){
            ##! this is not pretty ! Just the group indicator, make sure it is a column vector
            group.mm = stats::model.matrix(random[[l]]$group, adata$data)
            if(ncol(group.mm) == 2){
                ## Probably a numeric variable
                ## warning(sprintf('Grouping variable %s appears to be numeric, group indices in the model may be different than in the data',
                ##                 as.character(random[[l]]$group)))
                random[[l]]$group = as.numeric(as.factor(as.character(group.mm[,2])))
            }else{
                ## Probably a factor
                group.mm = remove.zero.cols(group.mm)
                random[[l]]$group = (group.mm %*% 0:(ncol(group.mm) - 1))[,1] + 1
            }
        }
    }
    ## ... now fill the stan data structure
    if(length(random) > 0){
        for(l in 1:length(random)){
            data[[sf('group_%d_size', l)]] = max(random[[l]]$group)
            data[[sf('group_%d', l)]] = random[[l]]$group
            for(par_type in par_types){
                if(!is.null(random[[l]][[par_type]])){
                    mm = remove.zero.cols(stats::model.matrix(random[[l]][[par_type]], adata$data))
                    data[[sf('Z_%s_ncol_%d', par_type, l)]] = ncol(mm)
                    data[[sf('Z_%s_%d', par_type, l)]] = mm
                    ## just the name of the variable in the random effects specification list
                    v = sf('%s_sd_scale', par_type)
                    data[[sf('%s_sd_scale_%d', par_type, l)]] = parse_prior(random[[l]][[v]],
                                                                            default_prior$scale[[par_type]],
                                                                            ## matrix dimensions
                                                                            c(data[[sf('%s_size', par_type)]], data[[sf('Z_%s_ncol_%d', par_type, l)]]),
                                                                            sf('%s_sd_scale', par_type))
                    v = sf('%s_nu', par_type)
                    data[[sf('lkj_%s_nu_%d', par_type, l)]] = if(is.null(random[[l]][[v]])){ 1 }else{ random[[l]][[v]][1] }
                    ##     v = sf('%s_sd_scale', par_type)
                    ##     ## total prior size
                    ##     s = data[[sf('%s_size', par_type)]] * ncol(mm)
                    ##     if(is.null(random[[l]][[v]])){
                    ##         res = rep(default_prior$scale[[par_type]], s)
                    ##     }else{
                    ##         if(length(random[[l]][[v]]) == 1){
                    ##             res = rep(random[[l]][[v]][1], s)
                    ##         }else{
                    ##             if(length(random[[l]][[v]]) != s)
                    ##                 stop(sf("Vector of prior parameters for %s, grouping factor %d, is not of length 1 or %d",
                    ##                         par_type, l, s))
                    ##             res = random[[l]][[v]]
                    ##         }
                    ##     }
                    ##     data[[sf('%s_sd_scale_%d', par_type, l)]] = fix_stan_dim(res)
                }
            }
        }
    }
    data
}

is.formula = function(x)class(x) == 'formula'

check_fixed = function(fixed){
    for(par in c('eta', 'delta', 'gamma', 'theta'))
        if(!is.null(fixed[[par]]))
           if(!is.formula(fixed[[par]]))
               stop(sprintf('fixed$%s must be of type formula', par))
}

check_random = function(random){
    if(length(random) > 0){
        not_ll = "random must be a list of lists, e.g., list(list(group = ~ id, delta = ~ 1))"
        if(!is.list(random))
            stop(ll)
        for(i in 1:length(random)){
            if(!is.list(random[[i]]))
                stop(ll)
            if(is.null(random[[i]]$group))
                stop("Each list in the random list must contain a 'group' element")
            for(par in c('group', 'eta', 'delta', 'gamma', 'theta'))
                if(!is.null(random[[i]][[par]]))
                    if(!is.formula(random[[i]][[par]]))
                        stop(sprintf('random[[%d]]$%s is not of type formula', i, par))
        }
    }
}

check_adata = function(adata){
    if(!all(c('data', 'counts') %in% names(adata)))
        stop('something is wrong with the aggregated data object')
}

fix_stan_dim = function(x)if(length(x) == 1){ array(x, dim = 1) }else{ x }
fix_stan_dim_m = function(x)if(length(x) == 1){ matrix(x) }else{ x }

remove.zero.cols = function(m)as.matrix(m[,apply(m, 2, function(x)!all(x == 0))])
remove.one.cols = function(m)as.matrix(m[,apply(m, 2, function(x)!all(x == 1))])

parse_prior = function(value = NULL, default, dims, name){
    if(is.null(value)){
        if(length(dims) == 1){
            value = rep(default, dims)
        }else{
            value = matrix(default, nrow = dims[1], ncol = dims[2])
        }
    }else{
        if(is.matrix(value)){
            if(!all(c(nrow(value), ncol(value)) == dims))
                stop(sprintf('Incorrect dimensions of the %s matrix', name))
        }else{
            if(length(value) == 1){
                if(length(dims) == 1){
                    value = rep(value, dims)
                }else{
                    value = matrix(value, nrow = dims[1], ncol = dims[2])
                }
            }else{
                if(length(value) != prod(dims))
                    stop(sprintf("Prior specification %s must contain 1 or %d elements", name, dims))
            }
        }
    }
    fix_stan_dim_m(value)
}
