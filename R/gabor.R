## -*- coding: utf-8 -*-

#' Perceptual discrimination of short-duration gabor patches with Perceptual
#' Awareness Scale ratings.
#'
#' The Perceptual Awareness Scale ratings ('no experience', ..., 'absolutely
#' clear image') of discrimination responses to gabor patches presented for 32ms
#' or 64ms before or after the ratings were provided.
#'
#' On each trial participants had to classify the briefly presented gabor patch
#' as either oriented to the left or to the right using the arrow keys. The
#' participants were also asked to rate the stimuli on a 4 point Perceptual
#' Awareness Scale presented at the bottom of the screen. The gabor patch was
#' immediately followed by a mask. The PAS ratings ranged from "no experience"
#' to "absolutely clear image" and were provided either before (RD condition) or
#' after (DR condition) pressing the arrow keys. On each trial the gabor patch
#' was equally likely to be presented for 32ms or 64ms. Decision-rating order
#' was a between-subject variable and gabor patch presentation duration was a
#' within-subject variable. There were 47 participants and 48 trials per
#' condition.
#'
#' @format A data frame with 4509 rows and 10 variables: \describe{
#'   \item{timestamp}{Time when response was recorded} \item{duration}{Stimulus
#'   duration (ms)} \item{trial}{Trial number} \item{acc}{Accuracy: 1 = correct,
#'   0 = incorrect} \item{id}{Unique participant identification number}
#'   \item{order}{Between-subject order manipulation: DR = rating after
#'   decision, RD = decision after rating} \item{age}{Participant's age at the
#'   time of the study (years)} \item{gender}{Participant's gender: M = male, F
#'   = female} \item{rating}{Perceptual Awareness Scale rating: 0 = 'no
#'   experience', 3 = 'absolutely clear image'} \item{stim}{Stimulus class: 0 =
#'   gabor tilted to the left, 1 = gabor tilted to the right.} }
"gabor"
## ok
