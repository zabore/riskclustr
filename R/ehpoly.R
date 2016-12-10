################################################################################
#' Conduct an analysis of etiologic heterogeneity using polytomous logistic
#' regression
#'
#' \code{ehpoly} takes a multinomial logistic regression model fit and returns
#' results related to the question of whether a risk factor differs across
#' levels of a disease subtype, and whether a risk factor differs across levels
#' of each individual tumor marker or which the disease subtypes are comprised.
#'
#' @author Emily C Zabor \email{zabore@@mskcc.org}
#'
#' @param fit a model fit from \code{multinom}
#' @param z design matrix for correspondence between disease subtypes and
#' tumor markers. Number of rows should equal number of subtypes, number of
#' columns should equal number of tumor markers. -1 indicates tumor marker
#' positive within a subtype, 1 indicates tumor marker negative within a
#' subtype.
#'
#' @return Returns a list. \code{beta} is the estimated parameters from
#' the fit. \code{beta_se} is estimates of the standard
#' errors of the parameters from the fit. \code{eh_pval}
#' is a vector of the etiologic heterogeneity p-values. \code{or_ci_p}
#' is a dataframe with a odds ratio (95% CI) for each risk factor/subtype
#' combination, as well as a column of etiologic heterogeneity p-values.
#' \code{beta_p} is a dataframe with the estimated parameters from the
#'  fit as well as a column of etiologic heterogeneity p-values.
#' \code{gamma_p} is a dataframe with estimates of the tumor marker effects and
#' their associated p-values.
#'
#' @export
#'
################################################################################

ehpoly <- function(fit, z) {

  library(nnet)

  if(class(fit)[1] != "multinom") {
    stop("Argument to fit must be of class multinom")
  }

  p <- length(fit$coefnames) - 1 # number of covariates
  m <- length(summary(fit)$lev) - 1 # number of subtypes
  beta_plr <- coef(fit)[, -1]
  vcov_plr <- vcov(fit)

  # Calculate the ORs and 95% CIs
  or <- round(exp(beta_plr), 2)
  lci <- round(exp(beta_plr - qnorm(0.975) *
                     summary(fit)$standard.errors[, -1]), 2)
  uci <- round(exp(beta_plr + qnorm(0.975) *
                     summary(fit)$standard.errors[, -1]), 2)

  # Calculate the etiologic heterogeneity p-values
  # Create two tables: 1) ORs, 95% CIs, p-values; 2) Betas, p-values
  or_ci_p <- beta_p <- matrix(NA, nrow = p, ncol = (m + 1))
  pval <- NA
  for(i in 1:p) {
    if(p == 1) {
      or_ci_p[i, 1:m] <- paste0(or, " (", lci, " - ", uci, ")")
      beta_p[i, 1:m] <- round(beta_plr, 2)
    } else {
      or_ci_p[i, 1:m] <- t(paste0(or[, i], " (", lci[, i], " - ",
                                  uci[, i], ")"))
      beta_p[i, 1:m] <- round(t(beta_plr[, i]), 2)
    }

    B <- beta_plr[, i]
    V <- vcov_plr[grep(fit$coefnames[i + 1], rownames(vcov_plr), fixed = T),
                   grep(fit$coefnames[i + 1], rownames(vcov_plr), fixed = T)]
    Lmat <- matrix(0, nrow = (m - 1), ncol = m)
    for(j in 1:nrow(Lmat)) {
      Lmat[j, j:(j+1)] <- c(1, -1)
    }
    pval[i] <- wald.test(b = B, Sigma = V, L = Lmat)$result$chi2["P"]
    or_ci_p[i, (m + 1)] <- beta_p[i, (m + 1)] <- round(pval[i], 3)
  }

  # Format the resulting dataframes
  or_ci_p <- as.data.frame(or_ci_p, stringsAsFactors = F)
  beta_p <- as.data.frame(beta_p, stringsAsFactors = F)
  rownames(or_ci_p) <- rownames(beta_p) <- fit$coefnames[-1]
  colnames(or_ci_p) <- colnames(beta_p) <- c(fit$lev[-1], "p_het")
  or_ci_p$p_het[or_ci_p$p_het == "0"] <- "<.001"
  beta_p$p_het[beta_p$p_het == 0.000] <- "<.001"

  # Calculate the gammas and their associated p-values
  V <- vcov_plr[-seq(1, nrow(vcov_plr), ncol(beta_plr) + 1),
                -seq(1, nrow(vcov_plr), ncol(beta_plr) + 1)]

  gamma_plr <- gamma_pval_plr <- matrix(NA, ncol = ncol(beta_plr), nrow = 2)
  gamma_p <- NULL
  for(i in 1:ncol(z)) {
    for(j in 1:ncol(beta_plr)) {

      gamma_plr[i, j] <- sum(beta_plr[, j] * z[, i]) / (m / 2)

      gamma_pval_plr[i, j] <- wald.test(b = beta_plr[, j],
                                        Sigma = V[seq(j, nrow(V), ncol(beta_plr)),
                                                  seq(j, nrow(V), ncol(beta_plr))],
                                        L = t(z[, i]))$result$chi2["P"]
    }
    gamma_p <- cbind(gamma_p, round(gamma_plr[i, ], 2),
                     round(gamma_pval_plr[i, ], 3))
  }

  # Format the results
  gamma_p <- as.data.frame(gamma_p, stringsAsFactors = F)
  gamma_p[, seq(2, ncol(gamma_p), 2)][gamma_p[, seq(2, ncol(gamma_p), 2)]
                                      == 0.000] <- "<.001"
  rownames(gamma_p) <- colnames(beta_plr)
  tm <- seq(1:ncol(z))
  colnames(gamma_p) <- as.vector(sapply(tm, function(x)
    {c(paste0("est", tm[x]), paste0("pval", tm[x]))}))

  # Returns
    return(list(beta = coef(fit)[, -1],
              beta_se = summary(fit)$standard.errors[, -1],
              eh_pval = pval,
              or_ci_p = or_ci_p,
              beta_p = beta_p,
              gamma_p = gamma_p))
}
