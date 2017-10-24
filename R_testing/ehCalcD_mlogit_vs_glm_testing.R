# Comparing D results from glm() and mlogit()
source("C:/Users/zabore/Google Drive/Columbia/Dissertation/code/ethet/R/ehCalcD_simData.R")
library(mlogit)

ehCalcD <- function(data, cls, M, formula, formula2) {
    ncontrol <- nrow(data[data[, cls] == 0, ])
    ncase <- nrow(data[data[, cls] != 0, ])
    fprob <- fprob0 <- matrix(NA, ncontrol, M)

    for(m in 1:M) {
        # make an indicator for whether a case is in the class of interest
        data$resp <- 1 * (data[, cls] == m)
        # labels all controls and all cases of class m as TRUE, all others FALSE
        data$csubset <- (data[, cls] == m | data[, cls] == 0)
        # Fit model for class m
        fit0 <- glm(formula, family = binomial, data, subset = csubset)
        # linear predictor for each class for controls only
        fprob[, m] <- exp(fit0$linear.predictors[fit0$y == 0])
    }

    # predicted risk for each class for controls only
    fprob[, 1:M] <- fprob[, 1:M] / (1 + rowSums(fprob[, 1:M]))

    # transform the data for use in mlogit
    xx1 <- mlogit.data(xx, choice = cls, shape = "wide")

    # fit the polytomous logistic regression model
    mod <- mlogit(formula = formula2, data = xx1)

    # predicted risk for each class for controls only
    fprob0[, 1:M] <- fitted(mod, outcome = FALSE)[, 2:(M + 1)][
        fitted(mod, outcome = FALSE)[, 1] == fitted(mod), ]

    # mus are the average predicted probs for controls for each class
    mus <- colMeans(fprob)
    pis <- mus / sum(mus)

    D <- 0 # initialize D
    for(i in 1:(M - 1)) {
        for(j in (i + 1):M) {
            D <- D + (1 / ncontrol) * pis[i] * pis[j] * (sum(fprob[, i]^2) / mus[i]^2 + sum(fprob[, j]^2) / mus[j]^2 - 2 * sum(fprob[, i] * fprob[, j]) / (mus[i] * mus[j]))
        }
    }

    D0 <- 0 # initialize D
    for(i in 1:(M - 1)) {
        for(j in (i + 1):M) {
            D0 <- D0 + (1 / ncontrol) * pis[i] * pis[j] * (sum(fprob0[, i]^2) / mus[i]^2 + sum(fprob0[, j]^2) / mus[j]^2 - 2 * sum(fprob0[, i] * fprob0[, j]) / (mus[i] * mus[j]))
        }
    }

    return(c(D = round(D, 3), D0 = round(D0, 3)))
}


### Running once for testing
set.seed(20171024)
xx <- simData(2000)
colnames(xx) <- c("x1", "x2", "class")
data <- xx
cls <- "class"
M <- 3
formula <- as.formula(resp ~ x1 + x2)
formula2 <- mFormula(class ~ 1 | x1 + x2)

ehCalcD(xx, "class", 3, as.formula(resp ~ x1 + x2),
        mFormula(class ~ 1 | x1 + x2))


### Running 1000 simulated datasets to compare D from glm() and mlogit()
set.seed(20171024)

dd <- matrix(0, 1000, 2)
for(i in 1:1000) {
    xx <- simData(2000)
    dd[i,] <- ehCalcD(xx, xx$cls, 3, as.formula(resp ~ x1 + x2),
                      mFormula(cls ~ 1 | x1 + x2))
}



