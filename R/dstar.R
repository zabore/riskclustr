################################################################################
#' Estimate the incremental explained risk variation in case-only study
#'
#' @description \code{dstar} estimates the incremental explained risk variation
#' across a set of pre-specified disease subtypes in a case-only study.
#' This function takes a model formula and a wide dataset and does the needed
#' transformation on the dataset to get the correct format, fits the polytomous
#' logistic regression model, and calculates D based on the resulting risk
#' predictions.
#'
#' @param formula an mFormula() model formula for a polytmous logistic
#' regression model to be fit with mlogit() using the appropriate variable
#' names from the data of interest.
#' @param label the name of the subtype variable in the data with values 1
#' through M. Should be supplied in quotes, e.g. label = "subtype".
#' The highest frequency subtype will be used as the reference subtype in
#' model
#' @param M the number of subtypes among case subjects. For M>=2.
#' @param data the name of the dataframe that contains the relevant variables
#'
#' @examples
#'
#' # specify the model formula for three risk factors, x1, x2 and x3
#' mform <- mlogit::mFormula(subtype ~ 1 | x1 + x2 + x3)
#'
#' # calculate D*
#' # Exclude controls from data as this is a case-only calculation
#' dstar(mform, "subtype", 3, subtype_data[subtype_data$subtype > 0, ])
#'
#' @references
#' Begg, C. B., Seshan, V. E., Zabor, E. C., Furberg, H., Arora, A.,
#' Shen, R., . . . Hsieh, J. J. (2014). Genomic investigation of etiologic
#' heterogeneity: methodologic challenges. BMC Med Res Methodol, 14, 138.
#' doi: 10.1186/1471-2288-14-138
#'
#' @export
#'
################################################################################

dstar <- function(formula, label, M, data) {
  if (any(class(formula) == "mFormula") == FALSE) {
    stop("The formula argument must be of class mFormula. Please correctly specify the model formula and try again.")
  }

  # Check if the label variable is numeric/integer
  if (!class(data[[label]]) %in% c("numeric", "integer")) {
    stop("The argument to label must be numeric or integer. Arguments of type character and factor are not supported, please see the documentation.")
  }

  # Check if the label variable starts with 1
  if (min(data[[label]] != 1)) {
    stop("The argument to label should start with 1. Cases should be labeled with subtypes 1 through M, the total number of subtypes.")
  }

  # Check if M is a numeric variable >=2
  if (class(M) != "numeric" | M < 2) {
    stop("The argument to M, the total number of subtypes, must be a numeric value >=2.")
  }

  # Check if M is equal to the number of levels of label
  if (length(levels(factor(data[[label]]))) != M) {
    stop("M is not equal to the number of levels in the variable supplied to label. Please make sure M reflects the number of subtypes in the data.")
  }

  # sample size (all cases)
  ncase <- nrow(data)

  # make a label variable for use in the model
  # needs to have a reference level of 0
  # use the highest frequency label as reference for stability
  tlabel <- table(data[label])
  data$mlabel <- (data[label] != (which(tlabel == max(tlabel))[1])) * data[label]
  data$origlabel <- data[label]
  data[label] <- data$mlabel

  # transform the data for use in mlogit
  data2 <- mlogit::mlogit.data(data, choice = label, shape = "wide")

  # fit the polytomous logistic regression model
  mod <- mlogit::mlogit(formula = formula, data = data2)

  # predicted risk for each label for controls only
  fprob <- mod$probabilities

  # mus are the average predicted probs for controls for each label
  mus <- colMeans(fprob)
  pis <- mus / sum(mus)

  # calculate D*
  Dstar <- 0
  for (i in 1:(M - 1)) {
    for (j in (i + 1):M) {
      Dstar <- Dstar + (1 / ncase) * pis[i] * pis[j] *
        (sum(fprob[, i]^2) / mus[i]^2 + sum(fprob[, j]^2) /
          mus[j]^2 - 2 * sum(fprob[, i] * fprob[, j]) / (mus[i] * mus[j]))
    }
  }

  return(unname(Dstar))
}
