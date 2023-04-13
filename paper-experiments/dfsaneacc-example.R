installed <- require("dfsaneacc")

if (!installed) {
    install.packages("dfsaneacc")
    library("dfsaneacc")
}

n <- 3
x0 <- rep(1 / n^2, n)

expfun2 <- function(x) {
    n <- length(x)
    f <- rep(NA, n)
    f[1] <- exp(x[1]) - 1.0
    f[2:n] <- (2:n) / 10.0 * (exp(x[2:n]) + x[1:n - 1] - 1)
    f
}

ret <- dfsaneacc(x = x0, evalr = expfun2, nhlim = 6, epsf = 1.0e-6 * sqrt(n),
                 iprint = 0)
ret
