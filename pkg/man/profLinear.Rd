\name{profLinear}
\alias{profLinear}
\title{Clustering with Dirichlet Process Mixture of Linear Models using the Profile Method}

\description{
This function finds the most probable cluster partition in a Dirichlet process mixture of linear models using the profile method.
}

\usage{
profLinear(y, x, group, clust, param, method="stochastic",
           maxiter=1000, crit=1e-5, verbose=FALSE)
}

\arguments{
\item{y}{numeric vector of observations.}
\item{x}{numeric matrix of covariates.}
\item{group}{optional vector of factors (or coercible to factors) indicating grouping among observations. Observations that are grouped will be clustered together. This is useful if several values form a single longitudinal observation.}
\item{clust}{optional vector of factors (or coercible to factors) indicating initial clustering among observations. Grouped observations (see \sQuote{group}) must have the same \sQuote{clust} value.}
\item{param}{optional list containing the any of the named elements \sQuote{alpha}, \sQuote{a0}, \sQuote{b0}, \sQuote{m0}, and \sQuote{s0} corresponding to the prior parameters of the normal-gamma Dirichlet process mixture. The prior parameters of the normal-gamma Dirichlet process mixture should all be scalars except \sQuote{m0} which should be a vector of length equal to the number of columns in \sQuote{x}.}
\item{method}{character string indicating the optimization method to be used. Meaningful values for this string are \sQuote{stochastic} (default), \sQuote{agglomerative}, and \sQuote{none}. 

The \sQuote{stochastic} method is still under investigation and may be unsuitable for some datasets. Optimization involves computing a moving average of the relative change in the marginal posterior distribution over the possible clusters after each iteration. The optimization stopping criterion is the minumim value this quantity can take before stopping the optimization cycle. If the optimization cycle reaches the maximum allowable iterations before meeting the stopping criterion, a warning is issued. Even if the optimization criterion is met, the resulting clusters may not represent the global maximum under the marginal posterior distribution over the possible clusters.

The \sQuote{agglomerative} method initially places each observation into seperate clusters. At each iteration, two of the remaining clusters are merged, where the merged clusters are chosen such that the resulting increase in the posterior mass function is maximized. This is repeated until only one cluster remains. The MAP estimate is the cluster partition, among those considered, which maximizes the posterior mass function over the possible cluster partitions. See Ward (1963) for additional details.  

The \sQuote{none} method is typically used in conjunction with the \sQuote{clust} option to specify an initial cluster partition. If \sQuote{none} is specified without \sQuote{clust}, a simple algorithm is used to initialize the cluster partition. Otherwise, the cluster partition is initialized using the \sQuote{clust} argument. The posterior statistics are then computed for initialized clusters.
  }
\item{maxiter}{integer value specifying the maximum number of iterations for the optimization algorithm.}
\item{crit}{numeric scalar constituting a stopping criterion for the \sQuote{stochastic} optimization method.}
\item{verbose}{logical value indicating whether the routine should be verbose in printing.}
}

\value{
An instance of the class profLinear containing the following objects
\item{y}{the numeric vector of observations, where missing observations (NA) are removed}
\item{x}{the numeric matrix of covariates, where missing covariates (NA) are removed}
\item{group}{the grouping vector, where missing group values (NA) are removed}
\item{param}{the list of prior parameters}
\item{clust}{a numeric vector of integers indicating cluster membership for each non-missing observation}
\item{a}{a numeric vector containing the posterior parameter \eqn{ a_k }{ a_k } for each cluster}
\item{b}{a numeric vector containing the posterior parameter \eqn{ b_k }{ b_k } for each cluster}
\item{m}{a list of numeric vectors containing the posterior parameter \eqn{ m_k }{ m_k } for each cluster}
\item{s}{a list of numeric matrices containing the posterior parameter \eqn{ s_k }{ s_k } for each cluster}
\item{logp}{the log value of the marginal posterior mass function for the cluster partition evaluated at \sQuote{clust}}
}

\details{
This function fits a Dirichlet process mixture of linear models (DPMLM) using the profile method. This method will group the observations into clusters. The clusters are determined by maximizing the marginal posterior distribution over the space of possible clusters. Each cluster has an associated linear model. Notationally, the linear model for cluster \eqn{ k }{ k } has the form \deqn{ y_k = \gamma_k\prime x_k + \epsilon_k }{ y_k = \g_k' x_k + e_k }, where \eqn{ y_k }{ y_k } and \eqn{ x_k }{ x_k } are the observation vector and covariate matrix for a particular cluster, \eqn{ \epsilon_k }{ e_k } has a multivariate normal distribution with mean zero and precision matrix \eqn{ \tau_kI }{ t_k I }, and \eqn{ \gamma_k }{ g_k } is the vector of linear coefficients, about which inference is made. In the DPMLM, conditional on the clustering, \eqn{ \gamma_k }{ g_k } and \eqn{ \tau_k }{ t_k } have a joint normal-gamma posterior distribution of the form \deqn{ p(m_k, \tau_k| y, x) = n( \gamma_k | m_k, \tau_k s_k ) g( \tau_k | a_k, b_k ) }{ p(m_k, t_k| y, x) = n( g_k | m_k, t_k s_k ) g( t_k | a_k, b_k ) }, where \eqn{ n(\cdot) }{ n() } is the multivariate normal density function with mean vector \eqn{ m_k }{ m_k } and precision matrix \eqn{ \tau_k s_k }{ t_k s_k } and \eqn{ g(\cdot) }{ g() } is the gamma density function with shape and scale parameters \eqn{ a_k }{ a_k } and \eqn{ b_k }{ b_k }. In addition to the cluster indicators, the posterior quantities \eqn{ s_k }{ s_k }, \eqn{ m_k }{ m_k }, \eqn{ a_k }{ a_k }, and \eqn{ b_k }{ b_k } are provided for each cluster in the return value.


Missing observations (NA) are removed automatically and a warning is issued. The return value contains the reduced observation vector and covariate matrix.

}

\seealso{
\code{\link{profRand}}
}

\examples{
library(profdpm)

# set up some data
# linear model 0
x0  <- rnorm(50, 0, 3)
y0  <- x0 + rnorm(50, 0, 1) 

# linear model 1
x1  <- rnorm(50, 0, 3)
y1  <- 10 - x1 + rnorm(50, 0, 1)

# add a column of ones to the covariate matrix (intercept)
x   <- cbind(rep(1,100), c(x0, x1))
y   <- c(y0, y1)

# indicate grouping within each linear model
g   <- c( rep(seq(0,4),10), rep(seq(5,9),10) )

# fit the DPMLM
fit <- profLinear(y, x, g)

# plot the resulting fit(s)
plot(x[,2], y)
for(i in 1:length(fit$m)) {
  abline( a=fit$m[[i]][1], b=fit$m[[i]][2] )
}
}
\author{Matt Shotwell}
\references{
Ward, J. H. (1963), \sQuote{Heirarchical Grouping to Optimize an Objective Function,} Journal of the American Statistical Association, 58, 236-244.
}

