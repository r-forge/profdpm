\name{plot.profDensity}
\alias{plot.profDensity}
\title{An S3 'plot' Method for the 'profDensity' class}
\description{
  Generate plots associated with the \code{profDensity} class
}
\usage{
  \method{plot}{profDensity}(x, index, ...)
}
\arguments{
  \item{x}{a instance of the 'profDensity' class.}
  \item{index}{index of the cluster to be plotted. If missing, all clusters are plotted.}
  \item{...}{additional arguments passed to hist().}
}
\value{
  None
}
\details{
  This function plots the histogram of the original sample observations and overlays the nonparametric density estimate using the output of \code{profDensity}.
}
\seealso{
profDensity
}
\author{Matt Shotwell}
