## -*- coding: utf-8 -*-

#' Creates the SDT stan model code.
#'
#' \code{make_stan_model} function creates the stan model code
#' defining the SDT model with additional regression structure and
#' random effects if these are specified using the \code{random}
#' argument. See \code{\link{make_stan_data}} for details on how to
#' specify a random effects structure. The model is described in detail
#' in the paper (TODO reference)
#'
#' @param random an optional list specifying random effects
#'     structure. See \code{\link{make_stan_data}} for details.
#' @param gamma_link can be either 'softmax' (described in the paper),
#'     'log_distance', or 'log_ratio' (See the Readme file in the
#'     github repository)
#' @param model can be either 'sdt' (the default), 'uvsdt', or 'metad'
#' @return a string containing the full model definition in the stan
#'     modelling language.
#' @examples
#' data(gabor)
#' model = make_stan_model(list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1)))
#' cat(model)
#' @export
make_stan_model = function(random = NULL, gamma_link = 'softmax', model = 'sdt'){
    check_link(gamma_link)
    check_model(model)
    if(model %in% c('ordinal', 'uvordinal')){
        par_types = c('eta', 'gamma')
    }else{
        par_types = c('delta', 'gamma')
    }
    if(model %in% c('uvsdt', 'uvordinal'))par_types = c(par_types, 'theta')
    f = file(paste(path.package('bhsdtr'), '/stan_templates/sdt_template.stan', sep = ''))
    ## We go line by line
    first_pass = NULL
    for(l in readLines(f)){
        if(rmatch('PAR', l)){
            for(par_type in par_types)
                first_pass[length(first_pass) + 1] = gsub('PAR', par_type, l)
        }else{
            if(rmatch('//likelihood', l)){
                first_pass = c(first_pass,
                               readLines(sprintf('%s/stan_templates/likelihood_%s.stan',
                                                 paste(path.package('bhsdtr')), model)))
            }else{
                first_pass = c(first_pass, l)
            }
        }
    }
    model = ''
    for(part in first_pass){
        ## If this is part of the random effects' specification ...
        if(rmatch(sprintf('//(%s)', paste(c('common', par_types), collapse = '|')), part)){
            ## ... and there are some random effects in the model ...
            if(length(random) > 0)
                for(l in 1:length(random)){
                    ## ... then add the line read from the template with % replaced by the grouping factor number ...
                    if(rmatch('//common', part))
                        model[length(model)+1] = gsub('%', l, part)
                    ## ... and do the same with parts of delta/gamma
                    ## random effects' specification if delta/gamma is
                    ## associated with random effects of the given
                    ## grouping factor ...
                    for(par in par_types)
                        if(!is.null(random[[l]][[par]]) & rmatch(sprintf('//%s', par), part))
                            model[length(model)+1] = gsub('%', l, part)
                }
        }else if(rmatch('//link-gamma', part)){
            ## Replace the gamma link function specification with the
            ## chosen link function
            model = c(model,
                      readLines((f2 = file(paste(path.package('bhsdtr'), sprintf('/stan_templates/link_%s.stan', gamma_link), sep = '')))))
            close(f2)
        }else{
            ## This line is not about the random effects' specification
            model[length(model)+1] = part
        }
    }
    close(f)
    paste(model, collapse = '\n')
}

rmatch = function (pattern, vector){
  res = TRUE
  for (i in 1:length(vector)) {
    if (length(grep(pattern, vector[i])) > 0) {
      res[i] = TRUE
    }
    else {
      res[i] = FALSE
    }
  }
  res
}
