#' Perceptual discrimination of short-duration gabor patches with Perceptual Awareness Scale ratings.
#'
#' A small study (n=47) dataset containing Perceptual Awareness Scale
#' ratings of discrimination responses to gabor patches presented for
#' 32 or 64ms.
#'
#' @format A data frame with 4509 rows and 10 variables:
#' \describe{
#'   \item{timestamp}{the time when response was recorded}
#'   \item{duration}{stimulus duration, in ms}
#'   \item{trial}{trial number}
#'   \item{acc}{accuracy, 1 = correct, 0 = incorrect}
#'   \item{id}{unique participant identification number}
#'   \item{order}{between-subject manipulation, DR = rating after decision, RD = decision after rating}
#'   \item{age}{participant's age at the time of the study, in years}
#'   \item{gender}{participant's gender, M = male, F = female}
#'   \item{rating}{Perceptual Awareness Scale rating, 0 = 'did not see nothing', 3 = 'saw it very clearly'}
#'   \item{stim}{stimulus class, 0 = gabor tilted to the left, 1 = gabor tilted to the right}
#' }
"gabor"
