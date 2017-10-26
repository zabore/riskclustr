source("C:/Users/zabore/Google Drive/Columbia/Dissertation/code/ethet/R_testing/ehCalcD_simData.R")
source("C:/Users/zabore/Google Drive/Columbia/Dissertation/code/ethet/R/ehCalcD.R")


set.seed(20171024)
xx <- simData(2000)

ehCalcD(data = xx, cls = "cls", M = 3, formula = mFormula(cls ~ 1 | x1 + x2))


suppressMessages(library(mlogit))
data <- xx
cls <- "cls"
M <- 3
formula <- mFormula(cls ~ 1 | x1 + x2)


ncontrol <- nrow(data[data[, cls] == 0, ])
ncase <- nrow(data[data[, cls] != 0, ])
# fprob <- matrix(NA, ncontrol, M)

# transform the data for use in mlogit
data2 <- mlogit.data(data, choice = cls, shape = "wide")

# fit the polytomous logistic regression model
mod <- mlogit(formula = formula, data = data2)

# predicted risk for each class for controls only
# fprob[, 1:M] <- fitted(mod, outcome = FALSE)[, 2:(M + 1)][
#     fitted(mod, outcome = FALSE)[, 1] == fitted(mod), ]

fprob <- mod$probabilities[data[, cls] == 0, -1]

# mus are the average predicted probs for controls for each class
mus <- colMeans(fprob)
pis <- mus / sum(mus)

# calculate D
D <- 0
for(i in 1:(M - 1)) {
    for(j in (i + 1):M) {
        D <- D + (1 / ncontrol) * pis[i] * pis[j] *
            (sum(fprob[, i]^2) / mus[i]^2 + sum(fprob[, j]^2) /
                 mus[j]^2 - 2 * sum(fprob[, i] * fprob[, j]) / (mus[i] * mus[j]))
    }
}