ehCalcD <- function(data, cls, k, formula) {
    ncontrol <- nrow(data[cls == 0, ])
    ncase <- nrow(data[cls != 0, ])
    fprob <- fprob0 <- matrix(NA, ncontrol, k)

    for(i in 1:k) {
        data$resp <- 1 * (cls == i)
        data$csubset <- (cls == i | cls == 0) # labels all controls and all cases of class k as TRUE, all others FALSE
        fit0 <- glm(formula, family = binomial, data, subset = csubset) # Fit model for class k
        fprob0[,i] <- exp(fit0$linear.predictors[fit0$y==0])
        fprob[, i] <- predict(fit0, newdata = data[cls == 0, ], type = 'response') # pred prob for each class for controls
    }

    mus <- colSums(fprob) / ncontrol # mus are the average predicted probs for controls for each class
    pis <- mus / sum(mus)

    D <- 0 # initialize D
    for(i in 1:(k - 1)) {
        for(j in (i + 1):k) {
            D <- D + (1 / ncontrol) * pis[i] * pis[j] * (sum(fprob[, i]^2) / mus[i]^2 + sum(fprob[, j]^2) / mus[j]^2 - 2 * sum(fprob[, i] * fprob[, j]) / (mus[i] * mus[j]))
        }
    }

    # D calculated by VES
    fprob0[,1:k] <- fprob0[,1:k]/(1+rowSums(fprob0[,1:k]))
    mus0 <- colSums(fprob0) / ncontrol # mus are the average predicted probs for controls for each class
    pis0 <- mus0 / sum(mus0)

    D0 <- 0 # initialize D
    for(i in 1:(k - 1)) {
        for(j in (i + 1):k) {
            D0 <- D0 + (1 / ncontrol) * pis0[i] * pis0[j] * (sum(fprob0[, i]^2) / mus0[i]^2 + sum(fprob0[, j]^2) / mus0[j]^2 - 2 * sum(fprob0[, i] * fprob0[, j]) / (mus0[i] * mus0[j]))
        }
    }
    
    return(c(D = round(D, 3), D0 = round(D0, 3)))
}

simData <- function(n, mu=1) {
    nn <- round(n*c(0.4, 0.2, 0.2, 0.2))
    # subtype indicator
    t1 <- rep(c(0, 1, 0, 0), nn)
    t2 <- rep(c(0, 0, 1, 0), nn)
    t3 <- rep(c(0, 0, 0, 1), nn)

    # 2 risk factors
    x1 <- rnorm(n) + mu*(t1 + sqrt(0.5)*t2 + 0*t3)
    x2 <- rnorm(n) + mu*(0*t1 + sqrt(0.5)*t2 + t3)
    cls <- rep(0:3, nn)
    data.frame(x1=x1, x2=x2, cls=cls)
}

set.seed(1234);
dd <- matrix(0, 1000, 2)
for(i in 1:1000) {
    xx <- simData(2000)
    dd[i,] <- ehCalcD(xx, xx$cls, 3, as.formula(resp ~ x1 + x2))
}

dd1 <- matrix(0, 1000, 2)
for(i in 1:1000) {
    xx <- simData(5000)
    dd1[i,] <- ehCalcD(xx, xx$cls, 3, as.formula(resp ~ x1 + x2))
}
