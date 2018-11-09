################################################################################
#' Estimate the overall risk heterogeneity
#'
#' @description \code{ksq} estimates the overall risk heterogeneity in a
#' two-class system (cases versus controls, class 2 versus class 1)
#'
#' @param formula a \code{formula()} object for use in \code{glm()}
#' @param cls the name of the subtype variable in the data, where 0 indicates
#' control subjects, or the reference class level), should be supplied in
#' quotes, e.g. cls = "class". Must be 0/1.
#' @param data the name of the dataframe that contains the relevant variables
#'
#' @examples
#'
#' class <- rep(c(0, 1), times = c(1000, 500))
#' x <- matrix(rnorm(1500 * 2), 1500, 2) +
#' model.matrix(~factor(class))[, -1] * 1.5
#' df <- data.frame(class = class, x1 = x[, 1], x2 = x[, 2])
#' ksq(formula = formula(class ~ x1 + x2), cls = "class", data = df)
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

ksq <- function(formula, cls, data) {

  # number of controls
  ncontrol <- nrow(data[data[, cls] == 0, ])

  # fit the polytomous logistic regression model
  mod <- glm(formula = formula, family = binomial, data = data)

  # predicted risk for controls only
  lp <- mod$linear.predictors[data[, cls] == 0]
  prob <- 1 / (1 + exp(-lp))

  # Want to calculate K^2, the risk predictability of the disease as a single entity
  Ksq <- ((sum(prob^2) / ncontrol) / ((sum(prob) / ncontrol)^2)) - 1

  return(Ksq = Ksq)
}
