# This R script tests a variety of erroneous profLinear invocations 
# using a simple dataset containing three groups of observations, 
# each exhibiting a (different) linear relationship with a covariate.

# test an expression that should generate a warning
testWarning <- function( expr, test ) {
  res <- tryCatch( expr, warning=function(w){} )
  if( !is.null(res) ) { stop( paste("test failed:", test) ) }
}
  
library(profdpm)
library(lattice)

x   <- matrix( c( rep( 1,600 ), rep( 1:10, 60 ) ), 600, 2 )
y   <- 2*x[ 1:200, 2 ] + rnorm( 200, 0, 1 )
y   <- c( y, -2*x[ 201:400, 2 ] + rnorm( 200, 0, 1 ) )
y   <- c( y, rnorm( 200, 0, 1 ) )
g   <- rep( 1:60, rep( 10, 60 ) )

# maxiter is only meaningful for positive integers
testWarning( profLinear( y, x, g, maxiter=-1 ), "negative maxiter")
testWarning( profLinear( y, x, g, maxiter=0.1 ), "non-integer maxiter")
testWarning( profLinear( y, x, g, maxiter=FALSE ), "logical maxiter (FALSE)")
testWarning( profLinear( y, x, g, maxiter=c(1,2) ), "length(maxiter) = 2")

# crit is only meaningful for positive reals
testWarning( profLinear( y, x, g, crit=-0.1 ), "negative crit")
testWarning( profLinear( y, x, g, crit=FALSE ), "logical crit (FALSE)")
testWarning( profLinear( y, x, g, crit=c(0.1,0.2) ), "length(crit) = 2")
