# As all the parameters must be non-negative and less than 1, we
# parameterize the model in terms of the logit transformed parameters.

#' logit transformation
#'
#' Convert to p to logit p
#' @param p parameter to transform
#' @export
logitT <- function(p){
  log(p/(1-p))
}

#' back transform logit
#'
#' Convert logit eta to eta
#' @param eta parameter to backtransform
#' @export
ilogitT	<- function(eta){
  1/(1+exp(-eta))
}
