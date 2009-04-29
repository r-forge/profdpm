#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "pdpm.h"

pdpm_t * pdpm_R_alloc(unsigned int size);
SEXP getListElementByName(SEXP list, const char * name);

SEXP profDensity(SEXP data, SEXP parm, SEXP iter, SEXP crit) {
  SEXP retval, elem, names, class, index;
  pdpm_t * obj;
  unsigned int i;
  double s, m, a, b; 
 
  PROTECT(retval = allocVector(VECSXP, 7));
  PROTECT(names = allocVector(STRSXP, 7));
  PROTECT(class = allocVector(STRSXP, 1));
  PROTECT(index = allocVector(INTSXP, LENGTH(data)));
  SET_STRING_ELT(names, 0, mkChar("data"));
  SET_STRING_ELT(names, 1, mkChar("parm"));
  SET_STRING_ELT(names, 2, mkChar("index"));
  SET_STRING_ELT(names, 3, mkChar("a"));
  SET_STRING_ELT(names, 4, mkChar("b"));
  SET_STRING_ELT(names, 5, mkChar("m"));
  SET_STRING_ELT(names, 6, mkChar("s"));
  SET_STRING_ELT(class, 0, mkChar("profDensity"));
  setAttrib(retval, R_NamesSymbol, names);
  setAttrib(retval, R_ClassSymbol, class);
  SET_VECTOR_ELT(retval, 0, data);
  SET_VECTOR_ELT(retval, 1, parm);
  SET_VECTOR_ELT(retval, 2, index);
 
  obj        = pdpm_R_alloc((unsigned int) LENGTH(data));
  obj->data  = (double *) REAL(VECTOR_ELT(retval, 0));
  obj->index = (unsigned int *) INTEGER(VECTOR_ELT(retval, 2));
  elem       = getListElementByName(parm, "alpha");
  obj->alpha = elem == R_NilValue ? DEFAULT_ALPHA : *(REAL(elem));
  elem       = getListElementByName(parm, "a0");
  obj->a0    = elem == R_NilValue ? DEFAULT_A0 : *(REAL(elem));
  elem       = getListElementByName(parm, "b0");
  obj->b0    = elem == R_NilValue ? DEFAULT_B0 : *(REAL(elem));
  elem       = getListElementByName(parm, "m0");
  obj->m0    = elem == R_NilValue ? DEFAULT_M0 : *(REAL(elem));
  elem       = getListElementByName(parm, "s0");
  obj->s0    = elem == R_NilValue ? DEFAULT_S0 : *(REAL(elem));

  pdpm_chunk(obj, *(INTEGER(iter)), *(REAL(crit))); /* optimization */

  pdpm_sort(obj);
  for(i = 0; i < obj->size; i++) {
    obj->index[i]++;
  }
  SET_VECTOR_ELT(retval, 3, allocVector(REALSXP, obj->usize));
  SET_VECTOR_ELT(retval, 4, allocVector(REALSXP, obj->usize));
  SET_VECTOR_ELT(retval, 5, allocVector(REALSXP, obj->usize));
  SET_VECTOR_ELT(retval, 6, allocVector(REALSXP, obj->usize));
  for(i = 0; i < obj->usize; i++) {
    s = obj->s0 + obj->unique[i];
    m = ( obj->s0 * obj->m0 + obj->sum[i] ) / s;
    a = obj->a0 + obj->unique[i];
    b = obj->b0 + obj->sumsq[i] + obj->s0 * pow(obj->m0,2) - s * pow(m,2);   
    REAL(VECTOR_ELT(retval, 3))[i] = a;
    REAL(VECTOR_ELT(retval, 4))[i] = b;
    REAL(VECTOR_ELT(retval, 5))[i] = m;
    REAL(VECTOR_ELT(retval, 6))[i] = s;
  }

  UNPROTECT(4);
  return(retval);
}

pdpm_t * pdpm_R_alloc(unsigned int size) {
  pdpm_t * obj;
  obj         = (pdpm_t *) R_alloc(1, sizeof(pdpm_t));
  obj->datasq = (double *) R_alloc(size, sizeof(double));
  obj->sum    = (double *) R_alloc(size, sizeof(double));
  obj->sumsq  = (double *) R_alloc(size, sizeof(double));
  obj->unique = (unsigned int *) R_alloc(size, sizeof(unsigned int));
  obj->indmap = (unsigned int *) R_alloc(size, sizeof(unsigned int));
  obj->size   = size;
  return(obj);
}

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
