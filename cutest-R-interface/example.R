

source("cutest.R")
library(dfsaneacc)

cutest_init("BOOTH")

n <- cutest_getn()
x0 <- cutest_getx0()
ret <- dfsaneacc(x=x0, evalr=cutest_evalr, nhlim=6,
                 epsf=1.0e-6*sqrt(n),iprint=0)
ret

cutest_end()
