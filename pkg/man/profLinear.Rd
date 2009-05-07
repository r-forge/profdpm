\name{profLinear}
\alias{profLinear}
\title{Dirichlet Process Mixture of Linear Models using the Profile Method}
\description{
  This function performs inference in a Dirichlet process mixture of linear models using the profile method.
}
\usage{
  profDensity(form, data, parm, iter = 1000, crit = 0.001)
}
\arguments{
  \item{form}{a formula object}
  \item{data}{a dataframe}
  \item{parm}{a list containing the named elements \sQuote{alpha}, \sQuote{a0}, \sQuote{b0}, \sQuote{m0}, and \sQuote{s0} corresponding to the prior parameters of the normal-gamma dirichlet process mixture (See Details)}
  \item{iter}{an integer value specifying the maximum number of iterations for the optimization algorithm}
  \item{crit}{an numeric value specifying the minimum proportional gain in the log posterior relative to  the cumulative gain over all previous iterations. This value constitutes a 'stopping criterion' for the optimization algorithm.}
}
\value{
  the formula object passed into the function
}
\details{
This function will fit a Dirichlet process mixture of linear models using the profile method. The function does nothing at the moment and simply returns the formula object that is initially passed in as the \sQuote{form} argument.

#The Dirichlet process mixture (\acronym{DPM}) is often used to model an error distribution nonparametrically, for example, in a linear model. The DPM groups sample observations with similar values in to clusters and treats these clusters independently. Within clusters, sample observations are modeled parametrically, for example, the normal distribution. The parameters of the within-cluster model are given a \sQuote{base prior} distribution. The base prior distribution conjugate to the within-cluster likelihood is often selected for simplicity. The function \code{profileDensity} attempts to find the maximum posterior (\acronym{MAP}) estimate of the indicator variables that determine cluster membership. In this way, \dQuote{profile} inference may be conducted on the within-cluster models. \cr \cr

#In \code{profDensity}, the within-cluster sample observations are assumed to be normally distributed with mean \eqn{\mu} and precision \eqn{\lambda}. The parameters \eqn{\mu} and \eqn{\lambda} are assumed to take the conjugate normal-gamma form with prior degrees of freedom \sQuote{a0}, sum of squared errors \sQuote{b0}, mean \sQuote{m0}, and sample size \sQuote{s0}. Once cluster membership is estimated, the within-cluster posterior distribution for \eqn{\mu} and \eqn{\lambda} also has the normal-gamma form. The vectors \sQuote{a}, \sQuote{b}, \sQuote{m}, and \sQuote{s} returned by \code{profileDensity} are posterior degrees of freedom, sum of squared errors, mean, and sample size for each cluster. \cr \cr

#Nonparametric density estimation proceeds by computing the within-cluster predictive densities at the posterior mode. The final density estimate is given by mixing the within-cluster densities in proportion to the corresponding posterior sample sizes. This step is not actually performed by \code{profDensity}. However, the S3 \code{plot} method for the \code{profDensity} class computes and plots this estimate against the histogram of the original sample observations.
}

\examples{
library(profdpm)
data(crbr08)
fit <- profLinear(a~b, crbr08, list(c=1,d=1))
}
\author{Matt Shotwell}
