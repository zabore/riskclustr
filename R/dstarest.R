################################################################################
#' Estimate the incremental explained risk variation in case-only study
#'
#' @description \code{dstarest} estimates the incremental explained risk variation
#' across a set of pre-specified disease subtypes in a case-only study
#'
#' @param formula an mFormula() model formula for a polytmous logistic
#' regression model to be fit with mlogit() using the appropriate variable
#' names from the data of interest, see Examples
#' @param cls the name of the subtype variable in the data with values 1
#' through M. Should be supplied in quotes, e.g. cls = "class"
#' @param M the number of subtypes among case subjects. For M>=2.
#' @param data the name of the dataframe that contains the relevant variables
#'
#' @examples
#'
#' # generate data - 4 classes and 2 risk factors
#' set.seed(20180521)
#' class <- rep(c(0, 1, 2, 3, 4), times = c(1000, 250, 250, 250, 250))
#' x <- matrix(rnorm(2000 * 2), 2000, 2) +
#' model.matrix(~factor(class))[, -1] %*%
#' t(matrix(c(1.5, 0, 0.75, 0.25, 0.25, 0.75, 0, 1.5), ncol = 4))
#' df <- data.frame(class = class, x1 = x[, 1], x2 = x[, 2])
#'
#' # remove the controls - this is a case only analysis example
#' df <- df[df$class != 0, ]
#'
#' library(mlogit)
#' mform <- mFormula(class ~ 1 | x1 + x2)
#'
#' dstarest(mform, "class", 4, df)
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

dstarest <- function(formula, cls, M, data) {

  if(any(class(formula) == "mFormula") == FALSE) {
    stop("The formula argument must be of class mFormula. Please correctly specify the model formula and try again.")
  }

  suppressWarnings(suppressMessages(library(mlogit)))

  # sample size (all cases)
  ncase <- nrow(data)

  # make a class variable for use in the model
  # needs to have a reference level of 0
  # use the highest frequency class as reference for stability
  tcls <- table(data[cls])
  data$mcls <- (data[cls] != (which(tcls == max(tcls))[1])) * data[cls]
  data$origclass <- data[cls]
  data[cls] <- data$mcls

  # transform the data for use in mlogit
  data2 <- mlogit.data(data, choice = cls, shape = "wide")

  # fit the polytomous logistic regression model
  mod <- mlogit(formula = formula, data = data2)

  # predicted risk for each class for controls only
  fprob <- mod$probabilities

  # mus are the average predicted probs for controls for each class
  mus <- colMeans(fprob)
  pis <- mus / sum(mus)

  # calculate D*
  Dstar <- 0
  for(i in 1:(M - 1)) {
    for(j in (i + 1):M) {
      Dstar <- Dstar + (1 / ncase) * pis[i] * pis[j] *
        (sum(fprob[, i]^2) / mus[i]^2 + sum(fprob[, j]^2) /
           mus[j]^2 - 2 * sum(fprob[, i] * fprob[, j]) / (mus[i] * mus[j]))
    }
  }

  return(Dstar = Dstar)
}
