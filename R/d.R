#' Estimate the incremental explained risk variation
#'
#' @description \code{d} estimates the incremental explained risk variation
#' across a set of pre-specified disease subtypes in a case-control study
#'
#' @param formula an \code{mlogit::mFormula()} model formula for a polytmous
#' logistic regression model to be fit with \code{mlogit()} using the
#' appropriate variable names from the data of interest
#' @param var the name of the subtype variable in the data, where 0 indicates
#' control subjects, should be supplied in quotes, e.g. var = "subtype"
#' @param M the number of subtypes. This should not include controls, but
#' only the number of subtypes among case subjects. For M>=2.
#' @param data the name of the dataframe that contains the relevant variables
#'
#' @examples
#'
#' # specify the model formula for three risk factors, x1, x2 and x3
#' mform <- mlogit::mFormula(subtype ~ 1 | x1 + x2 + x3)
#'
#' # calculate D
#' d(mform, "subtype", 3, subtype_data)
#'
#' @references
#' Begg, C. B., Zabor, E. C., Bernstein, J. L., Bernstein, L., Press, M. F., &
#' Seshan, V. E. (2013). A conceptual and methodological framework for
#' investigating etiologic heterogeneity. Stat Med, 32(29), 5039-5052.
#' doi: 10.1002/sim.5902
#'
#' @export
#'

d <- function(formula, var, M, data) {
  if (any(class(formula) == "mFormula") == FALSE) {
    stop("The formula argument must be of class mFormula. Please correctly specify the model formula and try again.")
  }

  ncontrol <- nrow(data[data[, var] == 0, ])
  ncase <- nrow(data[data[, var] != 0, ])

  # transform the data for use in mlogit
  data2 <- mlogit::mlogit.data(data, choice = var, shape = "wide")

  # fit the polytomous logistic regression model
  mod <- mlogit::mlogit(formula = formula, data = data2)

  # predicted risk for each class for controls only
  fprob <- mod$probabilities[which(data[, var] == 0), -1]

  # mus are the average predicted probs for controls for each class
  mus <- colMeans(fprob)
  pis <- mus / sum(mus)

  # calculate D
  D <- 0
  for (i in 1:(M - 1)) {
    for (j in (i + 1):M) {
      D <- D + (1 / ncontrol) * pis[i] * pis[j] *
        (sum(fprob[, i]^2) / mus[i]^2 + sum(fprob[, j]^2) /
          mus[j]^2 - 2 * sum(fprob[, i] * fprob[, j]) / (mus[i] * mus[j]))
    }
  }

  return(unname(D))
}
