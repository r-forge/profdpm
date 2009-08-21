library(profdpm)

x  <- matrix(c(rep(1,200),runif(200, -10, 10)),200,2)
y1 <- x[1:100,2] + rnorm(100, 0, 1)
y2 <- -x[101:200,2] + rnorm(100, 0, 1)

y  <- c(y1, y2)
g  <- c(rep(1:2, 50), rep(3:4, 50))
p  <- list(alpha=5, s0=1, m0=c(0,0), a0=1, b0=0)

fit <- profLinear(y, x, g, p, iter=1000)
