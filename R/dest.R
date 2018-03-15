################################################################################
#' Estimate the incremental explained risk variation
#'
#' @description \code{dest} estimates the incremental explained risk variation
#' across a set of pre-specified disease subtypes in a case-control study
#'
#' @param formula an mFormula() model formula for a polytmous logistic
#' regression model to be fit with mlogit() using the appropriate variable
#' names from the data of interest, see Examples
#' @param cls the name of the subtype variable in the data, where 0 indicates
#' control subjects, should be supplied in quotes, e.g. cls = "class"
#' @param M the number of subtypes. This could should not include controls, but
#' only the number of subtypes among case subjects
#' @param data the name of the dataframe that contains the relevant variables
#'
#' @examples
#'
#' # generate data - 4 classes and 2 risk factors
#' class <- rep(c(0, 1, 2, 3, 4), times = c(1000, 250, 250, 250, 250))
#' x <- matrix(rnorm(2000 * 2), 2000, 2) +
#' model.matrix(~factor(class))[, -1] %*%
#' t(matrix(c(1.5, 0, 0.75, 0.25, 0.25, 0.75, 0, 1.5), ncol = 4))
#' df <- data.frame(class = class, x1 = x[, 1], x2 = x[, 2])
#'
#' # specify the model formula
#' library(mlogit)
#' mform <- mFormula(class ~ 1 | x1 + x2)
#'
#' dest(mform, "class", 4, df)
#'
#' @references
#' Begg, C. B., Zabor, E. C., Bernstein, J. L., Bernstein, L., Press, M. F., &
#' Seshan, V. E. (2013). A conceptual and methodological framework for
#' investigating etiologic heterogeneity. Stat Med, 32(29), 5039-5052.
#' doi: 10.1002/sim.5902
#'
#' @export
#'
################################################################################

dest <- function(formula, cls, M, data) {

  if(any(class(formula) == "mFormula") == FALSE) {
    stop("The formula argument must be of class mFormula. Please correctly specify the model formula and try again.")
  }

  suppressWarnings(suppressMessages(library(mlogit)))

  ncontrol <- nrow(data[data[, cls] == 0, ])
  ncase <- nrow(data[data[, cls] != 0, ])

  # transform the data for use in mlogit
  data2 <- mlogit.data(data, choice = cls, shape = "wide")

  # fit the polytomous logistic regression model
  mod <- mlogit(formula = formula, data = data2)

  # predicted risk for each class for controls only
  fprob <- mod$probabilities[data[, cls] == 0, -1]

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
