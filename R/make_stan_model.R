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
#' @return a string containing the full model definition in the stan
#'     modelling language.
#' @examples
#' data(gabor)
#' model = make_stan_model(list(list(group = ~ id, delta = ~ -1 + duration, gamma = ~ 1)))
#' cat(model)
#' @export
make_stan_model = function(random = NULL){
    model = ''
    f = file(paste(path.package('bhsdtr'), '/stan_templates/sdt_template.stan', sep = ''))
    for(part in readLines(f)){
        ## Jeżeli to jest fragment dotyczący efektów losowych ...
        if(rmatch('//(common|delta|gamma)', part)){
            ## ... i w ogóle modelujemy efekty losowe ...
            if(length(random) > 0)
                for(l in 1:length(random)){
                    ## ... to indeksuj część wspólną ...
                    if(rmatch('//common', part))
                        model[length(model)+1] = gsub('%', l, part)
                    ## ... i części specyficzne dla parametrów delta i
                    ## gamma, o ile mają być pod wpływem czynników
                    ## losowych
                    for(par in c('delta', 'gamma'))
                        if(!is.null(random[[l]][[par]]) & rmatch(sprintf('//%s', par), part))
                            model[length(model)+1] = gsub('%', l, part)
                }
        }else{
            ## To nie jest fragment dotyczący efektów losowych, a więc
            ## kopiujemy bez zmian
            model[length(model)+1] = part
        }
    }
    close(f)
    paste(model, collapse = '\n')
}
## ok

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
