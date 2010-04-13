#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/*
   The Rand function computes the Rand statistic
   for two vectors of integers. The computing time of this 
   function increases with the square of the length of the 
   two vectors.

   Rand, W. (1971). Objective Criteria for the Evaluation
   of Clustering Methods. JASA 66:846-850.
 */

SEXP Rand( SEXP c1, SEXP c2 ) {
  unsigned int i, j, n, con = 0;
  unsigned int *v1 = (unsigned int *) INTEGER( c1 ), *v2 = (unsigned int *) INTEGER( c2 );
  SEXP ret;
  n = (unsigned int) LENGTH( c1 );
  for( i = 0; i < LENGTH( c1 ) - 1; i++ ) {
    for( j = i + 1; j < LENGTH( c1 ); j++ ) {
      //alternative 0 - times for n=10k - mean: 0.63908s sd: 0.00755s
      //if( !( (v1[ i ] ^ v1[ j ] ? 1 : 0) ^ (v2[ i ] ^ v2[ j ] ? 1 : 0) ) ) { con++; }
      //alternative 1 - times for n=10k - mean: 0.54451s sd: 0.00565s
      //if( v1[ i ] == v1[ j ] && v2[ i ] == v2[ j ] ) { con++; continue; }
      //if( v1[ i ] != v1[ j ] && v2[ i ] != v2[ j ] ) { con++; }
      //alternative 2 - times for n=10k - mean: 0.54287s sd: 0.00765s 
      if( v1[ i ] == v1[ j ] ) {
        if( v2[ i ] == v2[ j ] ) { con++; }
      } else {
        if( v2[ i ] != v2[ j ] ) { con++; }
      }
    }
  }
  PROTECT( ret = allocVector( REALSXP, 1 ) );
  REAL( ret )[ 0 ] = ( (double) con ) / ( ( n * ( n - 1 ) ) / 2 );
  UNPROTECT( 1 );
  return ret;
}
        

