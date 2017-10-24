################################################################################
#' Estimate the true population D
#'
#' @author Emily C Zabor \email{zabore@@mskcc.org}
#'
#' @description Adapted from code written by Venkat Seshan and made
#' generalizable to any number of subtypes, any number of risk factors, and
#' any risk factor means for the cases
#'
#' @param N population sample size
#' @param M number of disease subtypes
#' @param pi vector of control and subtype prevalences
#' i.e. c(pi0, pi1, pi2, pi3) for controls and 3 subtypes
#' @param P number of risk factors
#' @param mu_m P x M matrix of risk factor means for the cases, by default all
#' risk factors have mean 0 for control subjects
#'
#' @export
#'
################################################################################

trueDsim <- function(N, M, pi, P, mu_m) {
  # subtype indicator matrix
  cls <- rep(1:M, N * pi[2:(M + 1)])
  t <- matrix(0, nrow = length(cls), ncol = length(unique(cls)))
  t[cbind(seq_along(cls), cls)] <- 1
  t <- rbind(matrix(rep(0, N * pi[1] * M), N * pi[1], M), t)

  # risk factor matrix
  x <- matrix(rnorm(N * P), N, P) + t %*% t(mu_m)

  # likelihood contribution
  er <- exp(x %*% mu_m)

  # disease risks
  r0 <- as.vector(cbind(1, er) %*% pi)
  r <- pi[2:(M + 1)] * er / r0
  rtot <- rowSums(r)

  # coefficient of variations
  K2 <- var(rtot) / mean(rtot)^2
  Km2 <- apply(r, 2, function(x) var(x) / mean(x)^2)

  # D measure
  D <- sum(tp[2:(M + 1)] * Km2 / (1 - pi[1])) - K2

  return(D)
}



# ### New version
# set.seed(20171024)
# N <- 1e3
# M <- 3
# pi <- c(.4, .2, .2, .2)
# P <- 2
# mu_m <- matrix(c(1, 0, sqrt(0.5), sqrt(0.5), 0, 1), 2, 3)
#
# # subtype indicator matrix
# cls <- rep(1:M, N * pi[2:(M + 1)])
# t <- matrix(0, nrow = length(cls), ncol = length(unique(cls)))
# t[cbind(seq_along(cls), cls)] <- 1
# t <- rbind(matrix(rep(0, N * pi[1] * M), N * pi[1], M), t)
#
# # risk factor matrix
# x <- matrix(rnorm(N * P), N, P) + t %*% t(mu_m)
#
# # likelihood contribution
# er <- exp(x %*% mu_m)
#
# # disease risks
# r0 <- as.vector(cbind(1, er) %*% pi)
# r <- pi[2:(M + 1)] * er / r0
# rtot <- rowSums(r)
#
# # coefficient of variations
# K2 <- var(rtot) / mean(rtot)^2
# Km2 <- apply(r, 2, function(x) var(x) / mean(x)^2)
#
# # D measure
# D <- sum(tp[2:(M + 1)] * Km2 / (1 - pi[1])) - K2
# D2 <- D
#
#
# ### Old version
# set.seed(20171024)
# nn <- 1e3
# tp <- c(.4, .2, .2, .2)
# mu <- 1
#
# # subtype indicator
# t1 <- rep(c(0, 1, 0, 0), nn*tp)
# t2 <- rep(c(0, 0, 1, 0), nn*tp)
# t3 <- rep(c(0, 0, 0, 1), nn*tp)
#
# # 2 risk factors
# x1 <- rnorm(nn) + mu*(t1 + sqrt(0.5)*t2 + 0*t3)
# x2 <- rnorm(nn) + mu*(0*t1 + sqrt(0.5)*t2 + t3)
#
# # likelihood contribution
# er1 <- exp(mu*x1)
# er2 <- exp(mu*(sqrt(0.5)*(x1+x2)))
# er3 <- exp(mu*x2)
#
# # disease risks
# r0 <- cbind(1, er1, er2, er3) %*% tp
# r1 <- tp[2]*er1/r0
# r2 <- tp[3]*er2/r0
# r3 <- tp[4]*er3/r0
# r0 <- tp[1]/r0
# # overall risk
# r <- r1+r2+r3
#
# # coefficient of variations
# K2 <- var(r)/mean(r)^2
#
# Ka2 <- var(r1)/mean(r1)^2
# Kb2 <- var(r2)/mean(r2)^2
# Kc2 <- var(r3)/mean(r3)^2
#
# # D measure
# D <- sum(tp[2:4]*c(Ka2, Kb2, Kc2)/(1-tp[1])) - K2
# D1 <- D
#
#
# D1
# D2
