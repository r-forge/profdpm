#include "util.h"

/* stats */

/* 
   The lfactorial function computes the log (base e)
   factorial for an integer through a series expansion
   approximation given by Jolly (1961). This approximation 
   was found to significantly improve performance of this
   computation with little loss of accuracy.

   Jolley, L.B.W. (1961) Summation of Series. pp. 28
   ISBN 0-486-60023-8 
*/

double lfactorial(unsigned int x) {
  if( x == 0 ) { return 0; }
  return ( LN_SQRT_2PI + (0.5 + x)*log(x) - x );
}

/* R */

SEXP getListElementByName(SEXP list, const char * name) {
  SEXP elem = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  unsigned int i;
  for (i = 0; i < LENGTH(list); i++) {
    if(strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      elem = VECTOR_ELT(list, i);
      break;
    }
  }
  return(elem);
}
