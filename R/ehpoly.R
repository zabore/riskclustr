################################################################################
#' Conduct an analysis of etiologic heterogeneity using polytomous logistic
#' regression
#'
#' @author Emily C Zabor \email{zabore@@mskcc.org}
#'
#' @description \code{ehpoly} takes a list of individual tumor markers and a
#' list of risk factors and returns results related to the question of whether
#' each risk factor differs across levels of the disease subtypes and the
#' question of whether each risk factor differs across levels of each individual
#' tumor marker of which the disease subtypes are comprised.
#'
#' Input is a dataframe that contains the individual tumor markers, the risk
#' factors of interest, and an indicator of case or control status.
#' The tumor markers must be binary and must have levels
#' 0 or 1 for cases. The tumor markers should be left missing for control
#' subjects. For categorical tumor markers, a reference level should be selected
#' and then indicator variables for each remaining level of the tumor marker
#' should be created. For continuous tumor markers, categories should be formed
#' and then indicator variables can be constructed as in the case of categorical
#' tumor markers. Risk factors can be either binary or continuous. For
#' categorical risk factors, a reference level should be selected and then
#' indicator variables for each remaining level of the risk factor should be
#' created. Categorical risk factors entered as is will be treated as ordinal.
#'
#' @param tm a list of the names of the binary tumor markers.
#' Each must have levels 0 or 1 for case subjects. This value will be missing
#' for all control subjects.
#' @param rf a list of the names of the binary or continuous risk factors.
#' For binary risk factors the lowest level will be used as the reference level.
#' @param case denotes the variable that contains each subject's status as a
#' case or control. This value should be 1 for cases and 0 for controls.
#' Argument must be supplied in quotes.
#' @param df the name of the dataframe that contains the tumor markers and risk
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
#' \code{gamma} is a matrix containing the estimated gamma parameters,
#' obtained as linear combinations of the beta parameters with a column for each
#' risk factor and a row for each tumor marker.
#'
#' \code{gamma_se} contains the associated standard errors.
#'
#' \code{gamma_p} is a matrix of p-values for testing whether each risk factor
#' differs across levels of each tumor marker, with a column for each risk
#' factor and a row for each tumor marker.
#'
#' \code{or_ci_p} is a dataframe with a odds ratio (95% CI) for each risk
#' factor/subtype combination, as well as a column of etiologic heterogeneity
#' p-values.
#'
#' \code{beta_se_p} is a dataframe with the estimated beta parameters (SE) for
#' each risk factor/subtype combination, as well as a column of
#' etiologic heterogeneity p-values.
#'
#' \code{gamma_se_p} is a dataframe with estimates of the gamma tumor marker
#' effects (SE) and their associated p-values.
#'
#'
#' @export
#'
################################################################################

ehpoly <- function(tm, rf, case, df) {

  library(nnet)
  library(aod)

  # Check if levels of each element in tm have only values 0 or 1
  if(any(sapply(lapply(tm, function(x) {levels(as.factor(df[, x]))}),
                function(y) {any(!(y %in% c("0", "1")))
  }) == TRUE)) {
    stop("All non-missing elements of tm must have values 0 or 1")
  }
  
  # Check if levels of case-control indicator are all 0 or 1
  if(any(!(as.factor(df[, "case"]) %in% c("0", "1"))) == TRUE) {
    stop("All elements of case must have values 0 or 1")
  }

  k <- length(tm) # number of tumor markers
  p <- length(rf) # number of covariates
  st <- expand.grid(rep(list(c(0, 1)), k)) #subtype by tm matrix
  colnames(st) <- tm # columns of matrix have names of tumor markers
  st$sub_name <- apply(st, 1, function(x) {paste(x, collapse = "/")})
  st$sub <- rownames(st)
  m <- nrow(st) # number of subtypes

  # Need to create a subtype variable in the data based on matrix st
  df <- merge(df, st, by = tm, all = T)
  df$sub[df[, case] == 0] <- 0 # value of subtype is 0 for controls

  # Order the subtype names to correspond to the numbered order
  df$sub_name <- factor(df$sub_name, levels = unique(st$sub_name))

  ### Now fit the polytomous logistic regression model
  fit <- multinom(as.formula(paste0("sub ~ ", paste(rf, collapse = " + "))),
                  data = df, trace = FALSE)

  beta_plr <- matrix(coef(fit)[, -1], nrow = m)
  beta_se <- matrix(summary(fit)$standard.errors[, -1], nrow = m)
  rownames(beta_plr) <- rownames(beta_se) <- levels(df$sub_name)
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
    vcov_plr[grep(x, rownames(vcov_plr), fixed = T),
             grep(x, rownames(vcov_plr), fixed = T)]})
  
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
  colnames(or_ci_p) <- colnames(beta_se_p) <- c(levels(df$sub_name), "p_het")
  or_ci_p$p_het[or_ci_p$p_het == "0"] <- "<.001"
  beta_se_p$p_het[beta_se_p$p_het == "0"] <- "<.001"
  
  ### Calculate the gammas and their SEs and associated p-values
  # Need to alter the st matrix to remove subtype and replace 0 with -1
  gamma_plr <- gamma_plr_se <- gamma_plr_pval <- gamma_se_p <- NULL
  
  d <- as.matrix(st[, 1:k])
  d[d == 0] <- -1
  
  # Get the parameter ests as linear combo of betas
  gamma_plr <- matrix(t(beta_plr) %*% d / (m / 2), nrow = p)
  
  # Get the standard error ests
  gamma_plr_se <- t(sapply(V, function(x) sqrt((m / 2)^(-2) * 
                                                 diag(t(d) %*% x %*% d))))
  
  # Calculate the p-values for each tumor marker
  gamma_plr_pval <- matrix(sapply(1:k, function(j) {sapply(1:p, function(i) {
    wald.test(b = beta_plr[, i], Sigma = V[[i]], 
              L = t(d[, j]))$result$chi2["P"]})}), nrow = p)
  
  # Put results together into data frame
  gamma_se_p <- data.frame(matrix(unlist(lapply(1:k, function(k) {
    cbind(matrix(paste0(round(gamma_plr, 2), " (", round(gamma_plr_se, 2), ")"), 
                 nrow = p)[, k], round(gamma_plr_pval[, k], 3))})), nrow = p),
    stringsAsFactors = FALSE)

  # Format the results
  rownames(gamma_plr) <- rownames(gamma_plr_se) <-
    rownames(gamma_plr_pval) <- rownames(gamma_se_p) <- fit$coefnames[-1]
  colnames(gamma_plr) <- colnames(gamma_plr_se) <-
    colnames(gamma_plr_pval) <- tm
  gamma_se_p[, seq(2, ncol(gamma_se_p), 2)][gamma_se_p[, seq(2, ncol(gamma_se_p), 2)]
                                            == "0"] <- "<.001"
  colnames(gamma_se_p) <- as.vector(sapply(1:k, function(x){
    c(paste(tm[x], "est"), paste(tm[x], "pval"))}))
  
  # Returns
  return(list(beta = beta_plr,
              beta_se = beta_se,
              eh_pval = pval,
              gamma = gamma_plr,
              gamma_se = gamma_plr_se,
              gamma_pval = gamma_plr_pval,
              or_ci_p = or_ci_p,
              beta_se_p = beta_se_p,
              gamma_se_p = gamma_se_p))
}
