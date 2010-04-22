\name{profBinary}
\alias{profBinary}
\title{Clustering with Dirichlet Process Mixture of Binary Models using the Profile Method}

\description{
This function finds the most probable cluster partition in a Dirichlet process mixture of binary models using the profile method.
}

\usage{
profLinear(y, clust, param, method="stochastic",
           maxiter=1000, crit=1e-5, verbose=FALSE)
}

\arguments{
\item{y}{numeric matrix, where each row corresponds to a binary vector observation.}
\item{clust}{optional vector of factors (or coercible to factors) indicating initial clustering among observations.}
\item{param}{optional list containing the any of the named elements \sQuote{alpha}, \sQuote{a0}, and \sQuote{b0} corresponding to the prior parameters of the beta-binary Dirichlet process mixture. The prior parameters of the beta-binary Dirichlet process mixture should all be scalars.}
\item{method}{character string indicating the optimization method to be used. Meaningful values for this string are \sQuote{stochastic} (default), \sQuote{gibbs}, \sQuote{agglomerative}, and \sQuote{none}. 

The \sQuote{stochastic} method is still under investigation and may be unsuitable for some datasets. Optimization involves computing a moving average of the relative change in the marginal posterior distribution over the possible clusters after each iteration. The optimization stopping criterion is the minumim value this quantity can take before stopping the optimization cycle. If the optimization cycle reaches the maximum allowable iterations before meeting the stopping criterion, a warning is issued. Even if the optimization criterion is met, the resulting clusters may not represent the global maximum under the marginal posterior distribution over the possible clusters.

The \sQuote{gibbs} method implements the Polya urn Gibbs sampler. This method draws samples from the posterior distribution over the cluster partition in a sequential Gibbs fashion. The sample value with greatest posterior mass is returned. See MacEachern(1994) for details.

The \sQuote{agglomerative} method initially places each observation into seperate clusters. At each iteration, two of the remaining clusters are merged, where the merged clusters are chosen such that the resulting increase in the posterior mass function is maximized. This is repeated until only one cluster remains. The MAP estimate is the cluster partition, among those considered, which maximizes the posterior mass function over the possible cluster partitions. See Ward (1963) for additional details.  

The \sQuote{none} method is typically used in conjunction with the \sQuote{clust} option to specify an initial cluster partition. If \sQuote{none} is specified without \sQuote{clust}, a simple algorithm is used to initialize the cluster partition. Otherwise, the cluster partition is initialized using the \sQuote{clust} argument. The posterior statistics are then computed for initialized clusters.
  }
\item{maxiter}{integer value specifying the maximum number of iterations for the optimization algorithm.}
\item{crit}{numeric scalar constituting a stopping criterion for the \sQuote{stochastic} optimization method.}
\item{verbose}{logical value indicating whether the routine should be verbose in printing.}
}

\value{
An instance of the class profLinear containing the following objects
\item{y}{the numeric matrix of observations, where rows with missing observations (NA) are removed}
\item{param}{the list of prior parameters}
\item{clust}{a numeric vector of integers indicating cluster membership for each non-missing observation}
\item{a}{a list of numeric vectors containing the posterior vector \eqn{ a }{ a } for each cluster}
\item{b}{a list of numeric vectors containing the posterior vector \eqn{ b }{ b } for each cluster}
\item{logp}{the log value of the marginal posterior mass function for the cluster partition evaluated at \sQuote{clust}}
}

\details{
This function fits a Dirichlet process mixture of binary models (DPMBM) using the profile method. This method will group the observations into clusters. The clusters are determined by maximizing the marginal posterior distribution over the space of possible clusters. Each cluster has an associated binary model. The binary model assigns Bernoulli probabilities independently to each binary valued outcome, corresponding to the columns of \eqn{y}{y}. Each cluster has a unique set of outcome probabilities. 

Missing observations (NA) are removed automatically and a warning is issued. The return value contains the reduced observation vector and covariate matrix.

}

\seealso{
\code{\link{pci}}
}

\examples{
library(profdpm)
}
\author{Matt Shotwell}
\references{
Ward, J. H. (1963) Heirarchical Grouping to Optimize an Objective Function. Journal of the American Statistical Association 58:236-244 

MacEachern, S. N. (1994) Estimating Normal Means with Conjugate Style Dirichlet Process Prior. Communications in Statistics B 23:727-741 
}

