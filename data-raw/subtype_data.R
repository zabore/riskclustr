# Generate data for use testing the riskclustr package
# These data have 2 continuous risk factors (x1, x2) and one binary variable (x3)
# There is a true known subtype label and also 30 continuous tumor markers

set.seed(20181109)

M <- 3
N <- 2000
pi <- c(0.4, 0.2, 0.2, 0.2)
P1 <- 2
mu_m <- matrix(c(1.5, 0.75, 0, 0, 0.75, 1.5), ncol = 3, byrow = T)
K_A <- 15
K_C <- 15
a <- 2.1
b <- 1.8
lambda_Am <- matrix(rep(c(a, 0, a, 0, a), times = c(5, 15, 5, 15, 5)), ncol = 3)
lambda_Bj <- matrix(rep(c(b, 0, b, 0, b), times = c(5, 15, 5, 15, 5)), ncol = 3)

# generate etiologically distinct subtypes indicator matrix
tAcls <- rep(1:M, N * pi[2:(M + 1)])
tA <- matrix(0, nrow = length(tAcls), ncol = length(unique(tAcls)))
tA[cbind(seq_along(tAcls), tAcls)] <- 1
tA <- rbind(matrix(rep(0, N * pi[1] * M), N * pi[1], M), tA)

# generate continuous risk factor data
x1 <- matrix(rnorm(N * P1), N, P1) + tA %*% t(mu_m)

# generate binary risk factor data
x2 <- c(
  rbinom(N * pi[1], 1, 0.1),
  rbinom(N * pi[2], 1, 0.2),
  rbinom(N * pi[3], 1, 0.4),
  rbinom(N * pi[4], 1, 0.6)
)

x <- cbind(x1, x2)
colnames(x) <- paste0("x", seq(1:ncol(x)))

# matrix of EH-related tumor marker data (N x K_A)
tmA <- matrix(rnorm(N * K_A), N, K_A) + tA %*% t(lambda_Am)


# matrix of additional noise tumor marker data (N x K_C)
tmC <- matrix(rnorm(N * K_C), N, K_C)

# matrix of all tumor marker data for CASES ONLY with numeric colnames
y1 <- tail(cbind(tmA, tmC), -N * pi[1])

# then add NA for all controls
y <- rbind(
  matrix(rep(NA, ncol(y1) * N * pi[1]), N * pi[1], ncol(y1)),
  y1
)
colnames(y) <- paste0("y", seq(1:ncol(y)))


# combine all data and save
subtype_data <- as.data.frame(
  cbind(
    subtype = c(rep(0, N * pi[1]), tAcls),
    x,
    y
  )
)

subtype_data <- dplyr::mutate(subtype_data,
  subtype_name = dplyr::case_when(
    subtype == 0 ~ "Control",
    subtype == 1 ~ "Subtype A",
    subtype == 2 ~ "Subtype B",
    subtype == 3 ~ "Subtype C"
  )
)

subtype_data <- dplyr::select(
  subtype_data,
  subtype, subtype_name, dplyr::everything()
)

usethis::use_data(subtype_data, overwrite = TRUE)
