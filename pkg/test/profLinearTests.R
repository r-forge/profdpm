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
