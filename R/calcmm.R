################################################################################
#' Calculate the minimum misclassification rate in a cross-tabulation
#'
#' @author Emily C Zabor \email{zabore@@mskcc.org}
#'
#' @param tab a cross-tabulation table of the simulation-based class solution
#' that resulted in maximum D to the true class solution
#'
#' @description This function provides a way of reconciling the arbitrary class
#' labels that result from k-means clustering to obtain the misclassification
#' rate of the k-means clustering result to the true class solution
#'
#' @export
#'
################################################################################

calcmm <- function(tab) {
  library(gtools)
  m <- ncol(tab)
  permmat <- permutations(m, m, 1:m)
  misclass <- NULL
  for(i in 1:nrow(permmat)) {
    t2 <- tab[order(permmat[i, ]), ]
    misclass[i] <- 1 - sum(diag(t2)) / sum(t2)
  }
  return(misclass[which.min(misclass)])
}
