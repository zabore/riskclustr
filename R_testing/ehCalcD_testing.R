source("C:/Users/zabore/Google Drive/Columbia/Dissertation/code/ethet/R/ehCalcD_simData.R")
source("C:/Users/zabore/Google Drive/Columbia/Dissertation/code/ethet/R/ehCalcD.R")

set.seed(20171024)
xx <- simData(2000)

ehCalcD(data = xx, cls = "cls", M = 3, formula = mFormula(cls ~ 1 | x1 + x2))
