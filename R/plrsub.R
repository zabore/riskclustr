################################################################################
#' Conduct an analysis of etiologic heterogeneity for pre-defined subtypes
#' using polytomous logistic regression
#'
#' @author Emily C Zabor \email{zabore@@mskcc.org}
#'
#' @description \code{plrsub} takes a vector of class labels for pre-specified
#' subtypes and a list of risk factors and returns results related to the
#' question of whether each risk factor differs across levels of the disease
#' subtypes
#'
#' Input is a dataframe that contains the risk factors of interest and a
#' variable containing numeric class labels that is 0 for control subjects.
#' Risk factors can be either binary or continuous. For categorical risk
#' factors, a reference level should be selected and then indicator variables
#' for each remaining level of the risk factor should be created.
#' Categorical risk factors entered as is will be treated as ordinal.
#' Class labels for the cases can be specified as a vector.
#'
#' @param cls the name of the variable in the data that contains numeric class
#' labels. This should be 0 for all controls. Argument must be supplied in
#' quotes.
#' @param m is the number of subtypes
#' @param rf a list of the names of the binary or continuous risk factors.
#' For binary risk factors the lowest level will be used as the reference level.
#' @param data the name of the dataframe that contains the tumor markers and risk
#' factors.
#'
#' @return Returns a list.
#'
#' \code{beta} is a matrix containing the estimated beta parameters with a
#' column for each risk factor and a row for each disease subtype.
#'
#' \code{beta_se} contains the associated standard errors.
#'
#' \code{eh_pval} is a vector of p-values for testing whether each risk factor
#' differs across levels of the disease subtype.
#'
#' \code{or_ci_p} is a dataframe with a odds ratio (95% CI) for each risk
#' factor/subtype combination, as well as a column of etiologic heterogeneity
#' p-values.
#'
#' \code{beta_se_p} is a dataframe with the estimated beta parameters (SE) for
#' each risk factor/subtype combination, as well as a column of
#' etiologic heterogeneity p-values.
#'
#' @examples
#'
#' # generate the data - 4 classes and 2 risk factors
#' class <- rep(c(0, 1, 2, 3, 4), times = c(1000, 250, 250, 250, 250))
#' x <- matrix(rnorm(2000 * 2), 2000, 2) +
#' model.matrix(~factor(class))[, -1] %*%
#' t(matrix(c(1.5, 0, 0.75, 0.25, 0.25, 0.75, 0, 1.5), ncol = 4))
#' df <- data.frame(class = class, x1 = x[, 1], x2 = x[, 2])
#'
#' plrsub("class", 4, c("x1", "x2"), df)
#'
#' @export
#'
################################################################################

plrsub <- function(cls, m, rf, data) {

  library(nnet)
  library(aod)

  # Check if cls is a numeric vector
  if(is.numeric(data[, cls]) == FALSE) {
    stop("Class labels must be supplied as a numeric vector")
  }

  # Check if there are any colons in a variable and stop if so
  if(any(grep("^[^:]+:", rf)) == TRUE) {
    stop("Risk factor names cannot include colons. Please rename the offending risk factor and try again.")
  }

  p <- length(rf) # number of covariates

  ### Fit the polytomous logistic regression model
  fit <- multinom(as.formula(paste0(cls, " ~ ", paste(rf, collapse = " + "))),
                  data = data, trace = FALSE)

  beta_plr <- matrix(coef(fit)[, -1], nrow = m)
  beta_se <- matrix(summary(fit)$standard.errors[, -1], nrow = m)
  rownames(beta_plr) <- rownames(beta_se) <- levels(as.factor(data[, cls]))[-1]
  colnames(beta_plr) <- colnames(beta_se) <- fit$coefnames[-1]

  # Calculate the ORs and 95% CIs
  or <- round(exp(beta_plr), 2)
  lci <- round(exp(beta_plr - qnorm(0.975) *
                     summary(fit)$standard.errors[, -1]), 2)
  uci <- round(exp(beta_plr + qnorm(0.975) *
                     summary(fit)$standard.errors[, -1]), 2)

  ### Calculate the etiologic heterogeneity p-values
  # initiate null vectors
  beta_se_p <- or_ci_p <- pval <- NULL

  # V is a list where each element is the vcov matrix for a diff RF
  vcov_plr <- vcov(fit)
  V <- lapply(fit$coefnames[-1], function(x) {
    vcov_plr[which(sapply(strsplit(rownames(vcov_plr), ":"), "[[", 2) == x),
             which(sapply(strsplit(rownames(vcov_plr), ":"), "[[", 2) == x)]})

  # Lmat is the contrast matrix to get the etiologic heterogeneity pvalue
  Lmat <- matrix(0, nrow = (m - 1), ncol = m)
  Lmat[row(Lmat) == col(Lmat)] <- 1
  Lmat[row(Lmat) - col(Lmat) == -1] <- -1

  pval <- sapply(1:p, function(i) {
    wald.test(b = beta_plr[, i], Sigma = V[[i]], L = Lmat)$result$chi2["P"]})

  # Store results on both beta and OR scale
  beta_se_p <- data.frame(t(matrix(paste0(round(beta_plr, 2), " (",
                                          round(beta_se, 2), ")"), nrow = m)),
                          round(pval, 3), stringsAsFactors = FALSE)

  or_ci_p <- data.frame(t(matrix(paste0(or, " (", lci, " - ", uci, ")"),
                                 nrow = m)),
                        round(pval, 3), stringsAsFactors = FALSE)

  # Format the resulting dataframes
  rownames(or_ci_p) <- rownames(beta_se_p) <- fit$coefnames[-1]
  colnames(or_ci_p) <- colnames(beta_se_p) <- c(levels(as.factor(data[, cls]))[-1],
                                                "p_het")
  or_ci_p$p_het[or_ci_p$p_het == "0"] <- "<.001"
  beta_se_p$p_het[beta_se_p$p_het == "0"] <- "<.001"

  # Returns
  return(list(beta = beta_plr,
              beta_se = beta_se,
              eh_pval = pval,
              or_ci_p = or_ci_p,
              beta_se_p = beta_se_p))
}
