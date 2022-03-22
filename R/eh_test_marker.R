#' Test for etiologic heterogeneity of risk factors according to individual
#' disease markers in a case-control study
#'
#' @author Emily C Zabor \email{zabore@@mskcc.org}
#'
#' @description \code{eh_test_marker} takes a list of individual disease
#' markers,
#' a list of risk factors, a variable name denoting case versus control status,
#' and a dataframe, and returns results related to the question of
#' whether each risk factor differs across levels of the disease subtypes and
#' the question of whether each risk factor differs across levels of each
#' individual disease marker of which the disease subtypes are comprised.
#' Input is a dataframe that contains the individual disease markers, the risk
#' factors of interest, and an indicator of case or control status.
#' The disease markers must be binary and must have levels
#' 0 or 1 for cases. The disease markers should be left missing for control
#' subjects. For categorical disease markers, a reference level should be
#' selected
#' and then indicator variables for each remaining level of the disease marker
#' should be created. Risk factors can be either binary or continuous. For
#' categorical risk factors, a reference level should be selected and then
#' indicator variables for each remaining level of the risk factor should be
#' created.
#'
#' @param markers a list of the names of the binary disease markers.
#' Each must have levels 0 or 1 for case subjects. This value will be missing
#' for all control subjects. e.g. \code{markers = list("marker1", "marker2")}
#' @param factors a list of the names of the binary or continuous risk factors.
#' For binary risk factors the lowest level will be used as the reference level.
#' e.g. \code{factors = list("age", "sex", "race")}
#' @param case denotes the variable that contains each subject's status as a
#' case or control. This value should be 1 for cases and 0 for controls.
#' Argument must be supplied in quotes, e.g. \code{case = "status"}.
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
#' \code{gamma} is a matrix containing the estimated disease marker parameters,
#' obtained as linear combinations of the \code{\link{beta}} estimates,
#' with a row for each risk factor and a column for each disease marker.
#'
#' \code{gamma_se} is a matrix containing the estimated disease marker
#' standard errors, obtained based on a transformation of the \code{\link{beta}}
#' standard errors, with a row for each risk factor and a column for each
#' disease marker.
#'
#' \code{gamma_p} is a matrix of p-values for testing whether each risk factor
#' differs across levels of each disease marker, with a row for each risk
#' factor and a column for each disease marker.
#'
#' \code{or_ci_p} is a dataframe with the odds ratio (95\% CI) for each risk
#' factor/subtype combination, as well as a column of formatted etiologic
#' heterogeneity p-values.
#'
#' \code{beta_se_p} is a dataframe with the estimates (SE) for
#' each risk factor/subtype combination, as well as a column of formatted
#' etiologic heterogeneity p-values.
#'
#' \code{gamma_se_p} is a dataframe with disease marker estimates (SE) and
#' their associated p-values.
#'
#' @examples
#'
#' # Run for two binary tumor markers, which will combine to form four subtypes
#' eh_test_marker(
#'   markers = list("marker1", "marker2"),
#'   factors = list("x1", "x2", "x3"),
#'   case = "case",
#'   data = subtype_data,
#'   digits = 2
#' )
#' @export
#'

eh_test_marker <- function(markers, factors, case, data, digits = 2) {

  # Check if levels of each element in markers have only values 0 or 1
  if (any(sapply(
    lapply(markers, function(x) {
      levels(as.factor(data[[x]]))
    }),
    function(y) {
      any(!(y %in% c("0", "1")))
    }
  ) == TRUE)) {
    stop("All non-missing elements of tm must have values 0 or 1")
  }

  # Check if levels of case-control indicator are all 0 or 1
  if (any(!(as.factor(data[[case]]) %in% c("0", "1"))) == TRUE) {
    stop("All elements of case must have values 0 or 1")
  }

  # Check if there are any colons in a variable name and stop if so
  if (any(grep("^[^:]+:", factors)) == TRUE) {
    stop("Risk factor names cannot include colons. Please rename the offending risk factor and try again.")
  }

  k <- length(markers) # number of tumor markers
  p <- length(factors) # number of covariates
  st <- expand.grid(rep(list(c(0, 1)), k)) # subtype by tm matrix
  colnames(st) <- markers # columns of matrix have names of tumor markers
  st$sub_name <- apply(st, 1, function(x) {
    paste(x, collapse = "/")
  })
  st$sub <- rownames(st)
  m <- nrow(st) # number of subtypes

  # Need to create a subtype variable in the data based on matrix st
  data <- merge(data, st, by = unlist(markers), all = T)
  data$sub[data[[case]] == 0] <- 0 # value of subtype is 0 for controls

  # Order the subtype names to correspond to the numbered order
  data$sub_name <- factor(data$sub_name, levels = unique(st$sub_name))

  # write the formula
  mform <- mlogit::mFormula(
    stats::as.formula(paste0("sub ~ 1 |", paste(factors, collapse = " + ")))
  )

  # transform the data for use in mlogit
  data2 <- mlogit::mlogit.data(data, choice = "sub", shape = "wide")

  # fit the polytomous logistic regression model
  fit <- mlogit::mlogit(formula = mform, data = data2)

  coefnames <- unique(sapply(
    strsplit(rownames(summary(fit)$CoefTable), ":"),
    "[[", 1
  ))[-1]
  beta_plr <- matrix(summary(fit)$CoefTable[, 1], ncol = m, byrow = T)[-1, , drop = FALSE]
  beta_se <- matrix(summary(fit)$CoefTable[, 2], ncol = m, byrow = T)[-1, , drop = FALSE]
  colnames(beta_plr) <- colnames(beta_se) <- levels(as.factor(data[["sub"]]))[-1]
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
      which(sapply(strsplit(rownames(vcov_plr), ":"), "[[", 1) == x),
      which(sapply(strsplit(rownames(vcov_plr), ":"), "[[", 1) == x)
    ]
  })

  # Lmat is the contrast matrix to get the etiologic heterogeneity pvalue
  Lmat <- matrix(0, nrow = (m - 1), ncol = m)
  Lmat[row(Lmat) == col(Lmat)] <- 1
  Lmat[row(Lmat) - col(Lmat) == -1] <- -1

  pval <- sapply(1:p, function(i) {
    chisq1 <- t(Lmat %*% beta_plr[i, ]) %*%
      solve(Lmat %*% V[[i]] %*% t(Lmat)) %*%
      (Lmat %*% beta_plr[i, ])
    1 - stats::pchisq(chisq1, df = nrow(Lmat))
  })
  names(pval) <- coefnames

  # Store results on both beta and OR scale
  beta_se_p <- data.frame(t(matrix(paste0(
    round(beta_plr, digits), " (",
    round(beta_se, digits), ")"
  ), nrow = m)),
  round(pval, 3),
  stringsAsFactors = FALSE
  )

  or_ci_p <- data.frame(t(matrix(paste0(or, " (", lci, "-", uci, ")"),
    nrow = m
  )),
  round(pval, 3),
  stringsAsFactors = FALSE
  )

  # Format the resulting dataframes
  rownames(or_ci_p) <- rownames(beta_se_p) <- coefnames
  colnames(or_ci_p) <- colnames(beta_se_p) <-
    c(levels(as.factor(data[["sub"]]))[-1], "p_het")
  or_ci_p$p_het[or_ci_p$p_het == "0"] <- "<.001"
  beta_se_p$p_het[beta_se_p$p_het == "0"] <- "<.001"

  ### Calculate the gammas and their SEs and associated p-values
  # Need to alter the st matrix to remove subtype and replace 0 with -1
  gamma_plr <- gamma_plr_se <- gamma_plr_pval <- gamma_se_p <- NULL

  d <- as.matrix(st[, 1:k])
  d[d == 0] <- -1

  # Get the parameter ests as linear combo of betas
  gamma_plr <- matrix(beta_plr %*% d / (m / 2), nrow = p)

  # Get the standard error ests
  gamma_plr_se <- t(sapply(V, function(x) sqrt((m / 2)^(-2) *
      diag(t(d) %*% x %*% d))))

  # Calculate the p-values for each tumor marker
  gamma_plr_pval <- matrix(sapply(1:k, function(j) {
    sapply(1:p, function(i) {
      chisq1 <- t(t(d[, j]) %*% beta_plr[i, ]) %*%
        solve(t(d[, j]) %*% V[[i]] %*% t(t(d[, j]))) %*%
        (t(d[, j]) %*% beta_plr[i, ])
      1 - stats::pchisq(chisq1, df = nrow(t(d[, j])))
    })
  }), nrow = p)

  # Put results together into data frame
  gamma_se_p <- data.frame(
    matrix(
      unlist(
        lapply(1:k, function(k) {
          cbind(
            matrix(paste0(
              round(gamma_plr, 2),
              " (",
              round(gamma_plr_se, 2),
              ")"
            ),
            nrow = p
            )[, k],
            round(gamma_plr_pval[, k], 3)
          )
        })
      ),
      nrow = p
    ),
    stringsAsFactors = FALSE
  )

  # Format the results
  rownames(gamma_plr) <-
    rownames(gamma_plr_se) <-
    rownames(gamma_plr_pval) <-
    rownames(gamma_se_p) <-
    coefnames

  colnames(gamma_plr) <-
    colnames(gamma_plr_se) <-
    colnames(gamma_plr_pval) <-
    unlist(markers)

  # format p-values < 0.001
  gamma_se_p[, seq(2, ncol(gamma_se_p), 2)][gamma_se_p[, seq(2, ncol(gamma_se_p), 2)]
  == "0"] <- "<.001"

  colnames(gamma_se_p) <- as.vector(
    sapply(1:k, function(x) {
      c(
        paste(markers[[x]], "est"),
        paste(markers[[x]], "pval")
      )
    })
  )

  # Returns
  return(list(
    beta = beta_plr,
    beta_se = beta_se,
    eh_pval = pval,
    gamma = gamma_plr,
    gamma_se = gamma_plr_se,
    gamma_pval = gamma_plr_pval,
    or_ci_p = or_ci_p,
    beta_se_p = beta_se_p,
    gamma_se_p = gamma_se_p
  ))
}
