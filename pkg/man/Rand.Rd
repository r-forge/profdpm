\name{Rand}
\alias{Rand}
\title{The Rand Statistic}

\description{
  This function computes the Rand statistic (1971), a metric that quantifies the similarity between two vectors that may be coerced into factors.
}

\usage{
  Rand(x1, x2)
}

\arguments{
  \item{x1}{a vector}
  \item{x2}{a vector}
}

\value{
  A real number between 0.0 and 1.0, where 0.0 indicates no similarity and 1.0 indicates identical clustering.
}

\details{
  This function computes the similarity between two vectors representing the cluster partition of \eqn{ n }{ n } items. The two vectors should be of the same length. If observation \eqn{ k }{ k } is clustered together with observation \eqn{ j }{ j } in both \sQuote{x1} and \sQuote{x2}, then this constitutes an equivalent clustering. Hence, there are \eqn{ n(n-1)/2 }{ n(n-1)/2 } possible equivalent clusterings, where \eqn{ n }{ n } is the length of \sQuote{x1} and \sQuote{x2}. The metric is computed as the number of equivalent clusterings divided by the number possible equivalent clusterings. See Rand (1971) for additional details.
}

\references{
Rand, W. (1971) Objective Criteria for the Evaluation of Clustering Methods. Journal of the American Statistical Association 66:846-850
}

\author{Matt Shotwell}
