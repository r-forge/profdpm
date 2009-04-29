\name{profDensity}
\alias{profDensity}
\title{Nonparametric Density Estimation Using the Profile Method for the Dirichlet Process Mixture}
\description{
  This function performs nonparametric density estimation using the profile method for the normal-gamma Dirichlet process mixture model.
}
\usage{
  profDensity(data, 
                 parm=list(alpha=0.001,a0=0.001,b0=0.001,m0=0,s0=1), 
                 iter = 1000, crit = 0.001)
}
\arguments{
  \item{data}{a numeric vector}
  \item{parm}{a list containing the named elements \sQuote{alpha}, \sQuote{a0}, \sQuote{b0}, \sQuote{m0}, and \sQuote{s0} corresponding to the prior parameters of the normal-gamma dirichlet process mixture (See Details)}
  \item{iter}{an integer value specifying the maximum number of iterations for the optimization algorithm}
  \item{crit}{an numeric value specifying the minimum proportional gain in the log posterior relative to  the cumulative gain over all previous iterations. This value constitutes a 'stopping criterion' for the optimization algorithm.}
}
\value{
  A list of class type 'profDensity' containing the following
  \item{data}{a copy of the data vector}
  \item{parm}{a copy of the parm list}
  \item{index}{estimated cluster membership indices}
  \item{a}{ posterior degrees of freedom \sQuote{a} for each cluster}
  \item{b}{ posterior sum of squared errors \sQuote{b} for each cluster}
  \item{m}{ posterior mean \sQuote{m} for each cluster}
  \item{s}{ posterior sample size \sQuote{s} for each cluster}
}
\details{
The Dirichlet process mixture (\acronym{DPM}) is often used to model an error distribution nonparametrically, for example, in a linear model. The DPM groups sample observations with similar values in to clusters and treats these clusters independently. Within clusters, sample observations are modeled parametrically, for example, the normal distribution. The parameters of the within-cluster model are given a \sQuote{base prior} distribution. The base prior distribution conjugate to the within-cluster likelihood is often selected for simplicity. The function \code{profileDensity} attempts to find the maximum posterior (\acronym{MAP}) estimate of the indicator variables that determine cluster membership. In this way, \dQuote{profile} inference may be conducted on the within-cluster models. \cr \cr

In \code{profDensity}, the within-cluster sample observations are assumed to be normally distributed with mean \eqn{\mu} and precision \eqn{\lambda}. The parameters \eqn{\mu} and \eqn{\lambda} are assumed to take the conjugate normal-gamma form with prior degrees of freedom \sQuote{a0}, sum of squared errors \sQuote{b0}, mean \sQuote{m0}, and sample size \sQuote{s0}. Once cluster membership is estimated, the within-cluster posterior distribution for \eqn{\mu} and \eqn{\lambda} also has the normal-gamma form. The vectors \sQuote{a}, \sQuote{b}, \sQuote{m}, and \sQuote{s} returned by \code{profileDensity} are posterior degrees of freedom, sum of squared errors, mean, and sample size for each cluster. \cr \cr

Nonparametric density estimation proceeds by computing the within-cluster predictive densities at the posterior mode. The final density estimate is given by mixing the within-cluster densities in proportion to the corresponding posterior sample sizes. This step is not actually performed by \code{profDensity}. However, the S3 \code{plot} method for the \code{profDensity} class computes and plots this estimate against the histogram of the original sample observations.
}
\seealso{
plot.profDensity
}
\examples{
library(profdpm)
data(crbr08)
fit <- profDensity(crbr08[["chip"]])
plot(fit, breaks=40)
}
\author{Matt Shotwell}
