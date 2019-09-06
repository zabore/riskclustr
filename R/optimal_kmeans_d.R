#' Obtain optimal D solution based on k-means clustering of disease marker
#' data in a case-control study
#'
#' @description \code{optimal_kmeans_d} applies k-means clustering using the
#' \code{\link[stats]{kmeans}} function with many random starts. The D value is
#' then calculated for the cluster solution at each random start using the
#' \code{\link{d}} function, and the cluster solution that maximizes D is returned,
#' along with the corresponding value of D. In this way the optimally
#' etiologically heterogeneous subtype solution can be identified from possibly
#' high-dimensional disease marker data.
#'
#' @param markers a vector of the names of the disease markers. These markers
#' should be of a type that is suitable for use with
#' \code{\link[stats]{kmeans}} clustering. All markers will be missing
#' for control subjects. e.g. \code{markers = c("marker1", "marker2")}
#' @param M is the number of clusters to identify using
#' \code{\link[stats]{kmeans}} clustering. For M>=2.
#' @param factors a list of the names of the binary or continuous risk factors.
#' For binary risk factors the lowest level will be used as the reference level.
#' e.g. \code{factors = list("age", "sex", "race")}
#' @param case denotes the variable that contains each subject's status as a
#' case or control. This value should be 1 for cases and 0 for controls.
#' Argument must be supplied in quotes, e.g. \code{case = "status"}.
#' @param data the name of the dataframe that contains the relevant variables.
#' @param nstart the number of random starts to use with
#' \code{\link[stats]{kmeans}} clustering. Defaults to 100.
#' @param seed an integer argument passed to \code{\link{set.seed}}.
#' Default is NULL. Recommended to set in order to obtain reproducible results.
#'
#' @return Returns a list
#'
#' \code{optimal_d} The D value for the optimal D solution
#'
#' \code{optimal_d_data} The original data frame supplied through the
#' \code{data} argument, with a column called \code{optimal_d_label}
#' added for the optimal D subtype label.
#' This has the subtype assignment for cases, and is 0 for all controls.
#'
#' @examples
#' \donttest{
#' # Cluster 30 disease markers to identify the optimally
#' # etiologically heterogeneous 3-subtype solution
#' res <- optimal_kmeans_d(
#'   markers = c(paste0("y", seq(1:30))),
#'   M = 3,
#'   factors = list("x1", "x2", "x3"),
#'   case = "case",
#'   data = subtype_data,
#'   nstart = 100,
#'   seed = 81110224
#' )
#'
#' # Look at the value of D for the optimal D solution
#' res[["optimal_d"]]
#'
#' # Look at a table of the optimal D solution
#' table(res[["optimal_d_data"]]$optimal_d_label)
#' }
#'
#' @references
#' Begg, C. B., Zabor, E. C., Bernstein, J. L., Bernstein, L., Press, M. F., &
#' Seshan, V. E. (2013). A conceptual and methodological framework for
#' investigating etiologic heterogeneity. Stat Med, 32(29), 5039-5052.
#'
#' @export
#'

optimal_kmeans_d <- function(markers, M, factors, case, data,
                             nstart = 100, seed = NULL) {

  # Check if M is a numeric variable >=2
  if (class(M) != "numeric" | M < 2) {
    stop("The argument to M, the total number of subtypes, must be a numeric value >=2.")
  }

  # Check if levels of case-control indicator are all 0 or 1
  if (any(!(as.factor(data[, case]) %in% c("0", "1"))) == TRUE) {
    stop("All elements of case must have values 0 or 1")
  }

  # set the seed, if supplied
  # otherwise print a message
  if (!is.null(seed)) {
    set.seed(seed)
  } else {
    message("When no seed is set the results will not be reproducible.")
  }

  # Add a rowname variable to the dataframe for later merging
  data[["rowname"]] <- rownames(data)

  # Create a matrix of the disease marker data among cases
  y <- data[data[[case]] == 1, markers]

  # Run k-means clustering
  # This returns a list of cluster labels
  kres <- lapply(
    1:nstart,
    function(i) {
      stats::kmeans(
        x = y,
        centers = M,
        algorithm = "MacQueen",
        iter.max = 30
      )$cluster
    }
  )

  # Dalculate D for each solution
  d_kres <- lapply(
    1:nstart,
    function(j) {
      temp <- dplyr::right_join(
        dplyr::data_frame(
          rowname = rownames(data[data[[case]] == 1, ]),
          cls = kres[[j]]
        ),
        data,
        by = "rowname"
      )
      temp$cls[is.na(temp$cls)] <- 0
      d(
        label = "cls",
        M = M,
        factors = factors,
        data = temp
      )
    }
  )

  # Convert d_kres to vector
  d_kres <- unlist(d_kres)

  # Index of optimal D solution
  max_res <- which(d_kres == max(d_kres, na.rm = TRUE))[1]

  # Add optimal class solution into the dataframe
  optimal_d_data <- dplyr::right_join(
    dplyr::data_frame(
      rowname = rownames(data[data[[case]] == 1, ]),
      optimal_d_label = kres[[max_res]]
    ),
    data,
    by = "rowname"
  )
  optimal_d_data$optimal_d_label[is.na(optimal_d_data$optimal_d_label)] <- 0

  # Return results
  return(list(
    optimal_d = d_kres[max_res],
    optimal_d_data = optimal_d_data
  ))
}
