## -*- coding: utf-8 -*-

#' Creates the aggregated data object needed by other functions.
#'
#' \code{aggregate_responses} aggregates combined or binary responses and stores
#' them in a list together with additional variables and the stimulus variable.
#' The resulting object is required by \code{plot_sdt_fit} and
#' \code{make_stan_data} functions.
#'
#' @param data a data frame containing stimulus and response variables. If SDT
#'   parameters are regressed on additional variables then these variables have
#'   to be present as well.
#' @param stimulus a name of the stimulus variable in the provided data frame.
#'   The stimulus variable can be of any type but it must contain only two kinds
#'   of values.
#' @param response a name of the response variable in the provided data frame.
#'   This can be a binary classification response or a combined response (see
#'   \code{\link{combined_response}}).
#' @param variables an optional vector of names of additional variables (such as a
#'   variable encoding group membership) that must by preserved in the resulting
#'   aggregated data object.
#' @return a list with three components: \describe{ \item{data}{The data frame
#'   containing all the variables listed in the variables argument excluding the
#'   stimulus and response variables.} \item{stimulus}{The stimulus variable
#'   representing the stimulus class as either 1 or 2.}  \item{counts}{The
#'   response counts matrix with K columns where K is the number of possible
#'   (combined or binary) responses.}}
#' @examples
#' data(gabor)
#' gabor$resp = combined_response(gabor$stim, gabor$rating, gabor$acc)
#' ## See how the lower (higher) combined response values are more common
#' ## for the 1 (2) stimulus class
#' aggregate_responses(gabor, 'stim', 'resp')
#' @export
aggregate_responses = function(data, stimulus = NULL, response, variables = NULL, K = NULL){
    if(length(stimulus) > 1)
        stop('stimulus vector must be of length 1 or NULL')
    if(length(response) != 1)
        stop('response vector must be of length 1')
    ## stimulus classes encoded as 1 or 2
    if(!is.null(stimulus))
        data[[stimulus]] = as.numeric(as.factor(as.character(data[[stimulus]])))
    if(is.null(K))
        K = max(data[[response]], na.rm = T)
    res = plyr::ddply(data, unique(c(variables, stimulus)),
                      function(df)table(c(df[[response]], 1:K)) - 1)
    counts = res[, c((ncol(res)-K+1):ncol(res))]
    ## drop = F is important for the single column data frame case
    adata = list(data = res[, setdiff(variables, response), drop = F],
               counts = counts)
    if(!is.null(stimulus))
        adata$stimulus = res[[stimulus]]
    adata
}
