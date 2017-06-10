################################################################################
#' Function to calculate the D metric in a case/control setting
#'
#' @author Emily C Zabor \email{zabore@@mskcc.org}
#'
#' @description \code{ehCalcD} takes a data set and information about class
#' membership and calculates the D metric.
#'
#' @param data a data frame with the covariates
#' @param cls the class variable with values 0 through k, where 0 is for
#' control subjects and 1:k are labels for the subtypes
#' @param k is the number of classes
#' @param formula is the model formula for the polytomous logistic regression
#' model, with "resp" on the lhs and the covariates of interest, located in
#' \code{data} on the rhs
#'
#' @return returns a 3-digit numeric value
#'
#' @references Begg CB, Zabor EC, Bernstein JL, Bernstein L,
#' Press MF, Seshan VE. A conceptual and methodological framework
#' for investigating etiologic heterogeneity. Stat Med 2013; 32(29):5039-52.
#'
#' @export
#'
################################################################################

ehCalcD <- function(data, cls, k, formula) {
  ncontrol <- nrow(data[cls == 0, ])
  ncase <- nrow(data[cls != 0, ])
  fprob <- matrix(NA, ncontrol, k)

  for(i in 1:k) {
    data$resp <- 1 * (cls == i)
    data$csubset <- (cls == i | cls == 0) # labels all controls and all cases of class k as TRUE, all others FALSE
    fit0 <- glm(formula, family = binomial, data, subset = csubset) # Fit model for class k
    fprob[, i] <- predict(fit0, newdata = data[cls == 0, ], type = 'response') # pred prob for each class for controls
  }

  mus <- colSums(fprob) / ncontrol # mus are the average predicted probs for controls for each class
  pis <- mus / sum(mus)

  D <- 0 # initialize D
  for(i in 1:(k - 1)) {
    for(j in (i + 1):k) {
      D <- D + (1 / ncontrol) * pis[i] * pis[j] * (sum(fprob[, i]^2) / mus[i]^2 + sum(fprob[, j]^2) / mus[j]^2 -
                                                     2 * sum(fprob[, i] * fprob[, j]) / (mus[i] * mus[j]))
    }
  }
  return(D = round(D, 3))
}
