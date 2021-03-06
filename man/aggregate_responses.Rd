% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_responses.R
\name{aggregate_responses}
\alias{aggregate_responses}
\title{Creates the aggregated data object needed by other functions.}
\usage{
aggregate_responses(
  data,
  stimulus = NULL,
  response,
  variables = NULL,
  K = NULL
)
}
\arguments{
\item{data}{a data frame containing stimulus and response variables. If SDT
parameters are regressed on additional variables then these variables have
to be present as well.}

\item{stimulus}{a name of the stimulus variable in the provided data frame.
The stimulus variable can be of any type but it must contain only two kinds
of values.}

\item{response}{a name of the response variable in the provided data frame.
This can be a binary classification response or a combined response (see
\code{\link{combined_response}}).}

\item{variables}{an optional vector of names of additional variables (such as a
variable encoding group membership) that must by preserved in the resulting
aggregated data object.}
}
\value{
a list with three components: \describe{ \item{data}{The data frame
  containing all the variables listed in the variables argument excluding the
  stimulus and response variables.} \item{stimulus}{The stimulus variable
  representing the stimulus class as either 1 or 2.}  \item{counts}{The
  response counts matrix with K columns where K is the number of possible
  (combined or binary) responses.}}
}
\description{
\code{aggregate_responses} aggregates combined or binary responses and stores
them in a list together with additional variables and the stimulus variable.
The resulting object is required by \code{plot_sdt_fit} and
\code{make_stan_data} functions.
}
\examples{
data(gabor)
gabor$resp = combined_response(gabor$stim, gabor$rating, gabor$acc)
## See how the lower (higher) combined response values are more common
## for the 1 (2) stimulus class
aggregate_responses(gabor, 'stim', 'resp')
}
