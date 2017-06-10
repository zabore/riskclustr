################################################################################
#' Function to calculate the test statistic and proportion of variation
#' explained using MDMR
#'
#' @author Emily C Zabor \email{zabore@@mskcc.org}
#'
#' @description \code{fstat} Computes the test statistic and proportion of
#' variation explained using multivariate distance matrix regression. For binary
#' tumor marker data. Uses asymmetric binary distance metric.
#'
#' @param Y an NxK matrix of tumor marker data
#' @param X an NxP matrix of risk factor data
#'
#' @return returns a list containing f, the test statistic, and R, the
#' proportion of variation explained
#'
#' @references Zapala MA, Schork NJ. Multivariate regression analysis of
#' distance matrices for testing associations between gene expression patterns
#' and related variables. PNAS 2006; 103(51):19430-35.
#'
#' @export
#'
################################################################################

fstat_bin <- function(Y, X) {
  N <- nrow(Y)
  K <- ncol(Y)
  P <- ncol(X)

  # asymmetric binary distance matrix
  D <- as.matrix(dist(Y, upper = TRUE, diag = TRUE, method = "binary"))
  A <- (-1/2) * D^2

  # Gower's G
  I <- diag(N)
  one <- rep(1, N)
  G <- (I - (1/N) * one %*% t(one)) %*% A %*% (I - (1/N) * one %*% t(one))

  # Hat matrix
  H <- X %*% solve(t(X) %*% X) %*% t(X)

  # Compute the F-statistic
  f <- (sum(diag(H %*% G %*% H)) / (P - 1)) / (sum(diag((I - H) %*% G %*% (I - H))) / (N - P))

  # Calculate the proportion of variance explained
  R <- sum(diag(H %*% G %*% H)) / sum(diag(G))

  return(list(f = f, R = R))
}
