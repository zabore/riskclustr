#' Test for etiologic heterogeneity of risk factors according to disease
#' subtypes in a case-control study
#'
#' @author Emily C Zabor \email{zabore@@mskcc.org}
#'
#' @description \code{eh_test_subtype} takes a vector of class labels for
#' pre-specified subtypes and a list of risk factors and returns results
#' related to the
#' question of whether each risk factor differs across levels of the disease
#' subtypes
#'
#' Input is a dataframe that contains the risk factors of interest and a
#' variable containing numeric class labels that is 0 for control subjects.
#' Risk factors can be either binary or continuous. For categorical risk
#' factors, a reference level should be selected and then indicator variables
#' for each remaining level of the risk factor should be created.
#' Categorical risk factors entered as is will be treated as ordinal.
#' Class labels for the cases can be specified as a vector. The multinomial
#' logistic regression model is fit using \code{\link[mlogit]{mlogit}}.
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
#' @param digits the number of digits to round the odds ratios and associated
#' confidence intervals, and the estimates and associated standard errors.
#' Defaults to 2.
#'
#' @return Returns a list.
#'
#' \code{beta} is a matrix containing the raw estimates from the
#' polytomous logistic regression model fit with \code{\link[mlogit]{mlogit}}
#' with a row for each risk factor and a column for each disease subtype.
#'
#' \code{beta_se} is a matrix containing the raw standard errors from the
#' polytomous logistic regression model fit with \code{\link[mlogit]{mlogit}}
#' with a row for each risk factor and a column for each disease subtype.
#'
#' \code{eh_pval} is a vector of unformatted p-values for testing whether each
#' risk factor differs across the levels of the disease subtype.
#'
#' \code{or_ci_p} is a dataframe with the odds ratio (95% CI) for each risk
#' factor/subtype combination, as well as a column of formatted etiologic
#' heterogeneity p-values.
#'
#' \code{beta_se_p} is a dataframe with the estimates (SE) for
#' each risk factor/subtype combination, as well as a column of formatted
#' etiologic heterogeneity p-values.
#'
#' @examples
#'
#' eh_test_subtype(
#'     label = "subtype",
#'     M = 4,
#'     factors = list("x1", "x2", "x3"),
#'     data = subtype_data,
#'     digits = 2)
#'
#' @export
#'

eh_test_subtype <- function(label, M, factors, data, digits = 2) {

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

  # Check if there are any colons in a variable name and stop if so
  if (any(grep("^[^:]+:", factors)) == TRUE) {
    stop("Risk factor names cannot include colons. Please rename the offending risk factor and try again.")
  }

  p <- length(factors) # number of covariates

  # write the formula
  mform <- mlogit::mFormula(
    stats::as.formula(paste0(label, " ~ 1 |", paste(factors, collapse = " + ")))
  )

  # transform the data for use in mlogit
  data2 <- mlogit::mlogit.data(data, choice = label, shape = "wide")

  # fit the polytomous logistic regression model
  fit <- mlogit::mlogit(formula = mform, data = data2)

  coefnames <- unique(sapply(
    strsplit(rownames(summary(fit)$CoefTable), ":"),
    "[[", 2
  ))[-1]
  beta_plr <- matrix(summary(fit)$CoefTable[, 1], ncol = M, byrow = T)[-1, ]
  beta_se <- matrix(summary(fit)$CoefTable[, 2], ncol = M, byrow = T)[-1, ]
  colnames(beta_plr) <- colnames(beta_se) <- levels(as.factor(data[, label]))[-1]
  rownames(beta_plr) <- rownames(beta_se) <- coefnames

  # Calculate the ORs and 95% CIs
  or <- round(exp(beta_plr), digits)
  lci <- round(exp(beta_plr - stats::qnorm(0.975) * beta_se), digits)
  uci <- round(exp(beta_plr + stats::qnorm(0.975) * beta_se), digits)

  ### Calculate the etiologic heterogeneity p-values
  # initiate null vectors
  beta_se_p <- or_ci_p <- pval <- NULL

  # V is a list where each element is the vcov matrix for a diff RF
  vcov_plr <- stats::vcov(fit)
  V <- lapply(coefnames, function(x) {
    vcov_plr[
      which(sapply(strsplit(rownames(vcov_plr), ":"), "[[", 2) == x),
      which(sapply(strsplit(rownames(vcov_plr), ":"), "[[", 2) == x)
    ]
  })

  # Lmat is the contrast matrix to get the etiologic heterogeneity pvalue
  Lmat <- matrix(0, nrow = (M - 1), ncol = M)
  Lmat[row(Lmat) == col(Lmat)] <- 1
  Lmat[row(Lmat) - col(Lmat) == -1] <- -1

  pval <- sapply(1:p, function(i) {
    aod::wald.test(
      b = beta_plr[i, ],
      Sigma = V[[i]], L = Lmat
    )$result$chi2["P"]
  })
  pval <- as.data.frame(pval)
  rownames(pval) <- coefnames
  colnames(pval) <- "p_het"

  # Store results on both beta and OR scale
  beta_se_p <- data.frame(matrix(paste0(
    round(beta_plr, digits), " (",
    round(beta_se, digits), ")"
  ), ncol = M),
  round(pval, 3),
  stringsAsFactors = FALSE
  )

  or_ci_p <- data.frame(matrix(paste0(or, " (", lci, "-", uci, ")"), ncol = M),
    round(pval, 3),
    stringsAsFactors = FALSE
  )

  # Format the resulting dataframes
  rownames(or_ci_p) <- rownames(beta_se_p) <- coefnames
  colnames(or_ci_p) <- colnames(beta_se_p) <-
    c(levels(as.factor(data[, label]))[-1], "p_het")
  or_ci_p$p_het[or_ci_p$p_het == "0"] <- "<.001"
  beta_se_p$p_het[beta_se_p$p_het == "0"] <- "<.001"

  # Returns
  return(list(
    beta = beta_plr,
    beta_se = beta_se,
    eh_pval = pval,
    or_ci_p = or_ci_p,
    beta_se_p = beta_se_p
  ))
}
