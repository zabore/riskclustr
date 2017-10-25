################################################################################
#' Calculate the D metric to measure the extent of etiologic heterogeneity in
#' a case-control study
#'
#' @param formula an mFormula() model formula for a polytmous logistic
#' regression model to be fit with mlogit() using the appropriate variable
#' names from the data as interest e.g. formula = mFormula(class ~ 1 | x1 + x2)
#' for a model with subtype variable class and two individual-specific
#' predictors x1 and x2
#' @param cls the names of the subtype variable in the data, should be
#' supplied in quotes, e.g. cls = "class"
#' @param M the number of subtypes. This could should not include controls, but
#' only the number of subtypes among case subjects.
#' @param data the name of the dataframe that contains the relevant variables
#'
#' @export
#'
################################################################################

ehCalcD <- function(data, cls, M, formula) {

  library(mlogit)

  ncontrol <- nrow(data[data[, cls] == 0, ])
  ncase <- nrow(data[data[, cls] != 0, ])
  fprob <- matrix(NA, ncontrol, M)

  # transform the data for use in mlogit
  data2 <- mlogit.data(data, choice = cls, shape = "wide")

  # fit the polytomous logistic regression model
  mod <- mlogit(formula = formula, data = data2)

  # predicted risk for each class for controls only
  fprob[, 1:M] <- fitted(mod, outcome = FALSE)[, 2:(M + 1)][
    fitted(mod, outcome = FALSE)[, 1] == fitted(mod), ]

  # mus are the average predicted probs for controls for each class
  mus <- colMeans(fprob)
  pis <- mus / sum(mus)

  # calculate D
  D <- 0
  for(i in 1:(M - 1)) {
    for(j in (i + 1):M) {
      D <- D + (1 / ncontrol) * pis[i] * pis[j] *
        (sum(fprob[, i]^2) / mus[i]^2 + sum(fprob[, j]^2) /
           mus[j]^2 - 2 * sum(fprob[, i] * fprob[, j]) / (mus[i] * mus[j]))
    }
  }

  return(D = D)
}
