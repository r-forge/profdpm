#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP sandbox() {
  SEXP mat;
  int i, j;
  PROTECT(mat = allocMatrix(INTSXP, 5, 5));
  for(i = 0; i < 5; i++) {
    for(j = 0; j < 5; j++) {
      INTEGER(mat)[i*5 + j] = i*5 + j;
    }
  }
  UNPROTECT(1); 
  return mat;
}
