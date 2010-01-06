\name{profLinear}
\alias{profLinear}
\title{Dirichlet Process Mixture of Linear Models using the Profile Method}

\description{
  This function performs inference in a Dirichlet process mixture of linear models using the profile method.
}

\usage{
  profLinear(y, x, group, clust, param, method="Shotwell", stop=FALSE,
             maxiter=1000, crit=1e-5, prior="Dirichlet", verbose=FALSE)
}

\arguments{
  \item{y}{a numeric vector of observations.}
  \item{x}{a numeric matrix of covariates.}
  \item{group}{an optional vector of factors (or coercible to factors) indicating grouping among observations. Observations that are grouped will be clustered together. This is useful if several observations form a single longitudinal observation.}
  \item{clust}{an optional vector factors indicating initial cluster membership.}
  \item{param}{an optional list containing the any of the named elements \sQuote{alpha}, \sQuote{a0}, \sQuote{b0}, \sQuote{m0}, and \sQuote{s0} corresponding to the prior parameters of the normal-gamma Dirichlet process mixture. The prior parameters of the normal-gamma Dirichlet process mixture should all be scalars except \sQuote{m0} which should be a vector of length equal to the number of columns in x.}
  \item{method}{a character string indicating the optimization method to be used. Meaningful values for this string are \sQuote{Shotwell} (default), \sQuote{agglomerative}, and \sQuote{none}.}
  \item{stop}{either the logical FALSE (default), or a positive integer specifying the iteration at which the optimization algorithm should stop.}
  \item{maxiter}{an integer value specifying the maximum number of iterations for the optimization algorithm.}
  \item{crit}{a numeric scalar constituting a 'stopping criterion' for the optimization algorithm.}
  \item{prior}{a character string indicating the type of prior to be used for the clustering indicators. Meaningful values for this string are \sQuote{Dirichlet} (default) and \sQuote{cluster}. }
  \item{verbose}{a logical value indicating whether the routine should be verbose in printing.}
}

\value{
  An instance of the class profLinear containing the following objects
  \item{y}{the numeric vector of observations, where missing observations (NA) are removed}
  \item{x}{the numeric matrix of covariates, where missing covariates (NA) are removed}
  \item{group}{the grouping vector, where missing group values (NA) are removed}
  \item{parm}{the list of prior parameters}
  \item{clust}{a numeric vector of integers indicating cluster membership for each non-missing observation}
  \item{a}{a numeric vector containing the posterior parameter \eqn{ a_k }{ a_k } for each cluster}
  \item{b}{a numeric vector containing the posterior parameter \eqn{ b_k }{ b_k } for each cluster}
  \item{m}{a list of numeric vectors containing the posterior parameter \eqn{ m_k }{ m_k } for each cluster}
  \item{s}{a list of numeric matrices containing the posterior parameter \eqn{ s_k }{ s_k } for each cluster}
  \item{logp}{the log value of the marginal posterior mass function for the cluster partition evaluated at \sQuote{clust}}
}

\details{
This function fits a Dirichlet process mixture of linear models (DPMLM) using the profile method. This method will group the observations into clusters. The clusters are determined by maximizing the marginal posterior distribution over the space of possible clusters. Each cluster has an associated linear model. Notationally, the linear model for cluster \eqn{ k }{ k } has the form \deqn{ y_k = \gamma_k\prime x_k + \epsilon_k }{ y_k = \g_k' x_k + e_k }, where \eqn{ y_k }{ y_k } and \eqn{ x_k }{ x_k } are the observation vector and covariate matrix for a particular cluster, \eqn{ \epsilon_k }{ e_k } has a multivariate normal distribution with mean zero and precision matrix \eqn{ \tau_kI }{ t_k I }, and \eqn{ \gamma_k }{ g_k } is the vector of linear coefficients, about which inference is made. In the DPMLM, conditional on the clustering, \eqn{ \gamma_k }{ g_k } and \eqn{ \tau_k }{ t_k } have a joint normal-gamma posterior distribution of the form \deqn{ p(m_k, \tau_k| y, x) = n( \gamma_k | m_k, \tau_k s_k ) g( \tau_k | a_k, b_k ) }{ p(m_k, t_k| y, x) = n( g_k | m_k, t_k s_k ) g( t_k | a_k, b_k ) }, where \eqn{ n(\cdot) }{ n() } is the multivariate normal density function with mean vector \eqn{ m_k }{ m_k } and precision matrix \eqn{ \tau_k s_k }{ t_k s_k } and \eqn{ g(\cdot) }{ g() } is the gamma density function with shape and scale parameters \eqn{ a_k }{ a_k } and \eqn{ b_k }{ b_k }. In addition to the cluster indicators, the posterior quantities \eqn{ s_k }{ s_k }, \eqn{ m_k }{ m_k }, \eqn{ a_k }{ a_k }, and \eqn{ b_k }{ b_k } are provided for each cluster in the return value.

The augmented Dirichlet process yields a natural prior distribution for the cluster indicators. By default, this \sQuote{Dirichlet} prior distribution is used. However, this prior tends to force observations into a small number of clusters with non-homogeneous sizes. In contrast, the \sQuote{cluster} prior distribution allows for more clusters and promotes homogeneity is cluster size.

The optimization strategy for this model is still under investigation and may be unsuitable for some datasets. Optimization involves computing a moving average of the relative change in the marginal posterior distribution over the possible clusters after each iteration. The optimization stopping criterion is the minumim value this quantity can take before stopping the optimization cycle. If the optimization cycle reaches the maximum allowable iterations before meeting the stopping criterion, a warning is issued. Even if the optimization criterion is met, the resulting clusters may not represent the global maximum under the marginal posterior distribution over the possible clusters. 

This method will perform clustering without any specific data grouping. However, specifying the appropriate grouping structure may greatly reduce the computational burden of the algorithm by reducing the complexity of the marginal posterior over the possible clusters. 

Observations with missing data (NA) are removed automatically and a warning is issued. The return value contains the reduced observation vector and covariate matrix.

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
