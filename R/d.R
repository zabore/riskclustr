#' Estimate the incremental explained risk variation in a case-control study
#'
#' @description \code{d} estimates the incremental explained risk variation
#' across a set of pre-specified disease subtypes in a case-control study.
#' This function takes a model formula and a wide dataset and does the needed
#' transformation on the dataset to get the correct format, fits the polytomous
#' logistic regression model using \code{\link[mlogit]{mlogit}},
#' and calculates D based on the resulting risk predictions.
#'
#' @param label the name of the subtype variable in the data. This should be a
#' numeric variable with values 0 through M, where 0 indicates control subjects.
#' Must be supplied in quotes, e.g. \code{label = "subtype"}.
#' quotes.
#' @param M is the number of subtypes. For M>=2.
#' @param factors a list of the names of the binary or continuous risk factors.
#' For binary risk factors the lowest level will be used as the reference level.
#' e.g. \code{factors = list("age", "sex", "race")}.
#' @param data the name of the dataframe that contains the relevant variables.
#'
#' @examples
#'
#' d(
#'  label = "subtype",
#'  M = 4,
#'  factors = list("x1", "x2", "x2"),
#'  data = subtype_data)
#'
#' @references
#' Begg, C. B., Zabor, E. C., Bernstein, J. L., Bernstein, L., Press, M. F., &
#' Seshan, V. E. (2013). A conceptual and methodological framework for
#' investigating etiologic heterogeneity. Stat Med, 32(29), 5039-5052.
#' doi: 10.1002/sim.5902
#'
#' @export
#'

d <- function(label, M, factors, data) {

  # Check if the label variable is numeric/integer
  if (!class(data[[label]]) %in% c("numeric", "integer")) {
    stop("The argument to label must be numeric or integer. Arguments of type character and factor are not supported, please see the documentation.")
  }

  # Check if the label variable starts with 0
  if (min(data[[label]] != 0)) {
    stop("The argument to label should start with 0. 0 indicates control subjects and cases should be labeled 1 through M, the total number of subtypes.")
  }

  # Check if M is a numeric variable >=2
  if (class(M) != "numeric" | M < 2) {
    stop("The argument to M, the total number of subtypes, must be a numeric value >=2.")
  }

  # Check if M is equal to the number of non-zero levels of label
  if (length(levels(factor(data[[label]][data[[label]] != 0]))) != M) {
    stop("M is not equal to the number of non-zero levels in the variable supplied to label. Please make sure M reflects the number of subtypes in the data.")
  }

  # write the formula
  mform <- mlogit::mFormula(
    stats::as.formula(paste0(label, " ~ 1 |", paste(factors, collapse = " + "))))

  ncontrol <- nrow(data[data[, label] == 0, ])
  ncase <- nrow(data[data[, label] != 0, ])

  # transform the data for use in mlogit
  data2 <- mlogit::mlogit.data(data, choice = label, shape = "wide")

  # fit the polytomous logistic regression model
  mod <- mlogit::mlogit(formula = mform, data = data2)

  # predicted risk for each label for controls only
  fprob <- mod$probabilities[which(data[, label] == 0), -1]

  # mus are the average predicted probs for controls for each label
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
