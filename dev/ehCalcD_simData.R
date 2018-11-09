# This function generates simulated data with 2 continuous risk factors
# and 3 subtypes
# This data is for use testing the ehCalcD() function
# It will NOT help produce good warnings, etc, as the data are perfectly clean

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

