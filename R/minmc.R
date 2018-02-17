################################################################################
#' Calculate the minimum misclassification rate in a cross-tabulation
#'
#' @description \code{minmc} provides a way of reconciling the arbitrary class
#' labels that result from k-means clustering to obtain the minimum
#' misclassification rate of a k-means clustering result to a true class
#' solution, after aligning the class labels
#'
#' @author Emily C Zabor \email{zabore@@mskcc.org}
#'
#' @param tab a cross-tabulation table of a class solution from k-means
#' clustering to the true class solution
#'
#' @examples
#'
#' # Example where there is perfect alignment of class results, but the class
#' # labels do not correspond and need to be aligned
#' trueclass <- rep(c(1, 2, 3, 4), times = 20)
#' kclass <- rep(c(3, 2, 1, 4), times = 20)
#' table(kclass, trueclass)
#' minmc(table(kclass, trueclass))
#'
#' @export
#'
################################################################################

minmc <- function(tab) {
  library(gtools)

  if(class(tab) != "table") {
    stop("Argument tab should be the result of table()")
  }

  m <- ncol(tab)
  permmat <- permutations(m, m, 1:m)
  misclass <- NULL
  t2 <- NULL
  for(i in 1:nrow(permmat)) {
    t2[[i]] <- tab[order(permmat[i, ]), ]
    misclass[i] <- 1 - sum(diag(t2[[i]])) / sum(t2[[i]])
  }
  return(list(minmisclass = misclass[which.min(misclass)],
              mincrosstab = t2[[which.min(misclass)]]))
}
