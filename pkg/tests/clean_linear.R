# This R script tests the profLinear function on a dataset containing
# two groups of observation, each exhibiting a (different) linear relationship
# with a covariate.

library(profdpm)

x   <- matrix( c( rep( 1,100 ), rnorm( 200, 0, 3 ) ), 100, 3 )
y   <- x[ 1:50, 2 ] + rnorm( 50, 0, 1 )
y   <- c( y, 10 - x[ 51:100, 2 ] + rnorm( 50, 0, 1 ) )
g   <- c( rep( seq( 0,4 ), 10 ), rep( seq( 5,9 ), 10 ) )
p   <- list( alpha=0.001, a0=0.001, b0=0.001, m0=c( 0,0,0 ), s0=1 )
fit <- profLinear( y, x, g, p, iter=1000 )

plot( fit$x[ , 2 ], fit$y )
for( cl in 1:length( fit$m ) ) {
  abline( a=fit$m[[ cl ]][ 1 ], b=fit$m[[ cl ]][ 2 ] )
}
