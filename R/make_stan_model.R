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
#' @param gamma_link is either 'softmax' (described in the paper), 'log_distance' or 'log_ratio'
#' (See the Readme file in the github repository)
#' @param metad if TRUE meta-d' model code (only with the softmax gamma link function) is created,
#' default is FALSE.
#' @return a string containing the full model definition in the stan
#'     modelling language.
#' @examples
#' data(gabor)
#' model = make_stan_model(list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1)))
#' cat(model)
#' @export
make_stan_model = function(random = NULL, gamma_link = 'softmax', metad = FALSE){
    if(!(gamma_link %in% c('softmax', 'log_ratio', 'log_distance')))
        stop("The gamma_link function must be one of the following: 'softmax', 'log_ratio', 'log_distance'")
    model = ''
    if(!metad){
        f = file(paste(path.package('bhsdtr'), '/stan_templates/sdt_template.stan', sep = ''))
    }else{
        f = file(paste(path.package('bhsdtr'), '/stan_templates/metad_template.stan', sep = ''))
    }
    ## We go line by line
    for(part in readLines(f)){
        ## If this is part of the random effects' specification ...
        if(rmatch('//(common|delta|gamma)', part)){
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
                    for(par in c('delta', 'gamma'))
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
