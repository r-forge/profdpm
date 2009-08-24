\name{profLinear}
\alias{profLinear}
\title{Dirichlet Process Mixture of Linear Models using the Profile Method}
\description{
  This function performs inference in a Dirichlet process mixture of linear models using the profile method.
}
\usage{
  profLinear(y, x, group, parm, iter = 1000, crit = 0.001, verb=FALSE)
}
\arguments{
  \item{y}{a numeric vector observations}
  \item{x}{a numeric matrix of covariates}
  \item{group}{an optional vector of integers indicating grouping among observations}
  \item{parm}{a list containing the named elements \sQuote{alpha}, \sQuote{a0}, \sQuote{b0}, \sQuote{m0}, and \sQuote{s0} corresponding to the prior parameters of the normal-gamma Dirichlet process mixture.}
  \item{iter}{an integer value specifying the maximum number of iterations for the optimization algorithm}
  \item{crit}{an numeric value specifying the minimum proportional gain in the log posterior relative to  the cumulative gain over all previous iterations. This value constitutes a 'stopping criterion' for the optimization algorithm.}
  \item{verb}{a logical value indicating whether the routine should be verbose in printing}
}
\value{
  An instance of the class profLinear
}
\details{
This function fits a Dirichlet process mixture of linear models using the profile method.

The prior parameters of the normal-gamma Dirichlet process mixture should all be scalars except \sQuote{m0} which should be a vector of length equal to the number of columns in x.
}

\examples{
library(profdpm)
x0 <- matrix(rnorm(100, 0, 3), 50, 2)
x1 <- matrix(rnorm(100, 0, 3), 50, 2)
y0 <- x0[,1] + rnorm(50, 0, 1)
y1 <- 10 - x1[,1] + rnorm(50, 1, 1)
x  <- cbind(rep(1,100), rbind(x0, x1))
y  <- c(y0, y1)
g  <- c( rep(seq(0,4),10), rep(seq(5,9),10) )
p  <- list(alpha=0.001,a0=0.001,b0=0.001,m0=c(0,0,0),s0=1)
fit<- profLinear(y, x, g, p, iter=100)
}
\author{Matt Shotwell}
