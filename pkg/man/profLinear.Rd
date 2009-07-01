\name{profLinear}
\alias{profLinear}
\title{Dirichlet Process Mixture of Linear Models using the Profile Method}
\description{
  This function performs inference in a Dirichlet process mixture of linear models using the profile method.
}
\usage{
  profLinear(form, data, group, parm, iter = 1000, crit = 0.001)
}
\arguments{
  \item{form}{a formula object}
  \item{data}{a dataframe}
  \item{group}{an optional vector of integers indicating grouping among observations}
  \item{parm}{a list containing the named elements \sQuote{alpha}, \sQuote{a0}, \sQuote{b0}, \sQuote{m0}, and \sQuote{s0} corresponding to the prior parameters of the normal-gamma dirichlet process mixture (See Details)}
  \item{iter}{an integer value specifying the maximum number of iterations for the optimization algorithm}
  \item{crit}{an numeric value specifying the minimum proportional gain in the log posterior relative to  the cumulative gain over all previous iterations. This value constitutes a 'stopping criterion' for the optimization algorithm.}
}
\value{
  An instance of the class profLinear
}
\details{
This function will fit a Dirichlet process mixture of linear models using the profile method. The function does nothing at the moment and simply returns the an empty instance of the profLinear class.
}

\examples{
library(profdpm)
data(crbr08)
fit <- profLinear(a~b, crbr08, list(c=1,d=1))
}
\author{Matt Shotwell}
