#' Post-hoc test to obtain overall p-value for a factor variable used in a
#' \code{eh_test_subtype} fit.
#'
#' @author Emily C Zabor \email{zabore@@mskcc.org}
#'
#' @description \code{posthoc_factor_test} takes a \code{eh_test_subtype} fit
#' and returns an overall p-value for a specified factor variable.
#'
#' @param fit the resulting \code{eh_test_subtype} fit.
#' @param factor is the name of the factor variable of interest, supplied
#' in quotes, e.g. \code{factor = "race"}. Only supports a single factor.
#' @param nlevels is the number of levels the factor variable in \code{factor}
#' has.
#'
#' @return Returns a list.
#'
#' \code{pval} is a formatted p-value.
#'
#' \code{pval_raw} is the raw, unformatted p-value.
#'
#' @export
#'

posthoc_factor_test <- function(fit, factor, nlevels) {
  # Determine the number of subtypes based on the provided fitted object
  M <- ncol(fit$beta)

  # Extract the variance-covariance elements for the factor variable
  # Need to grab all variance-covariance elements associated with the single var
  # Uses name matching - unique names required!!
  V <- fit$var_covar[
    which(stringr::str_detect(rownames(fit$var_covar), factor)),
    which(stringr::str_detect(rownames(fit$var_covar), factor))
    ]

  # Extract beta coefficients for the factor variable
  # Uses name matching - unique names required!!
  B <- fit$beta[
    which(stringr::str_detect(rownames(fit$beta), factor)),
    ]

  # Convert B into a vector, row-wise
  Bvec <- as.vector(t(B))

  # Lmat is the contrast matrix to get the etiologic heterogeneity pvalue
  Lmat <- matrix(0, nrow = (M - 1), ncol = M)
  Lmat[row(Lmat) == col(Lmat)] <- 1
  Lmat[row(Lmat) - col(Lmat) == -1] <- -1

  # Now we need to create a block-diagonal matrix based on Lmat
  Lblock <- Matrix::bdiag(rep(list(Lmat), nlevels - 1))
  Lblock <- as.matrix(Lblock)

  pval_raw <- aod::wald.test(
    b = Bvec,
    Sigma = V,
    L = Lblock
  )$result$chi2["P"]

  pval <- ifelse(pval_raw < .001, "<.001", round(pval_raw, 3))

  return(list(
    pval_raw = pval_raw,
    pval = pval
  ))
}




