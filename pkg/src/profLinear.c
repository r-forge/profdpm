#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "util.h"
#include "pdpmlm.h"

pdpmlm_t * pdpmlm_R_alloc(unsigned int size);

//SEXP profDensity(SEXP data, SEXP parm, SEXP iter, SEXP crit) {

SEXP profLinear(SEXP form, SEXP data, SEXP group, SEXP parm, SEXP iter, SEXP crit) {
 
  SEXP retval, elem, names, class, index;
  pdpmlm_t * obj;
  unsigned int i;
  double s, m, a, b; 
 
  PROTECT(retval = allocVector(VECSXP, 8));
  PROTECT(names = allocVector(STRSXP, 8));
  PROTECT(class = allocVector(STRSXP, 1));
  PROTECT(index = allocVector(INTSXP, LENGTH(data)));
  SET_STRING_ELT(names, 0, mkChar("data"));
  SET_STRING_ELT(names, 1, mkChar("parm"));
  SET_STRING_ELT(names, 2, mkChar("group"));
  SET_STRING_ELT(names, 3, mkChar("index"));
  SET_STRING_ELT(names, 4, mkChar("a"));
  SET_STRING_ELT(names, 5, mkChar("b"));
  SET_STRING_ELT(names, 6, mkChar("m"));
  SET_STRING_ELT(names, 7, mkChar("s"));
  SET_STRING_ELT(class, 0, mkChar("profLinear"));
  setAttrib(retval, R_NamesSymbol, names);
  setAttrib(retval, R_ClassSymbol, class);
  SET_VECTOR_ELT(retval, 0, data);
  SET_VECTOR_ELT(retval, 1, parm);
  SET_VECTOR_ELT(retval, 2, group);
  SET_VECTOR_ELT(retval, 3, index);
 
  obj        = pdpmlm_R_alloc((unsigned int) LENGTH(data));
  obj->data  = (double *) REAL(VECTOR_ELT(retval, 0));
  obj->group = (unsigned int *) INTEGER(VECTOR_ELT(retval, 2));
  obj->index = (unsigned int *) INTEGER(VECTOR_ELT(retval, 3));
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

  //pdpmlm_chunk(obj, *(INTEGER(iter)), *(REAL(crit))); /* optimization */

  //pdpmlm_sort(obj);
  //for(i = 0; i < obj->size; i++) {
  //  obj->index[i]++;
  //}
  SET_VECTOR_ELT(retval, 4, allocVector(REALSXP, obj->usize));
  SET_VECTOR_ELT(retval, 5, allocVector(REALSXP, obj->usize));
  SET_VECTOR_ELT(retval, 6, allocVector(REALSXP, obj->usize));
  SET_VECTOR_ELT(retval, 7, allocVector(REALSXP, obj->usize));
  for(i = 0; i < obj->usize; i++) {
    s = obj->s0 + obj->unique[i];
    m = ( obj->s0 * obj->m0 + obj->sum[i] ) / s;
    a = obj->a0 + obj->unique[i];
    b = obj->b0 + obj->sumsq[i] + obj->s0 * pow(obj->m0,2) - s * pow(m,2);   
    REAL(VECTOR_ELT(retval, 4))[i] = a;
    REAL(VECTOR_ELT(retval, 5))[i] = b;
    REAL(VECTOR_ELT(retval, 6))[i] = m;
    REAL(VECTOR_ELT(retval, 7))[i] = s;
  }

  UNPROTECT(4);
  return(retval);
}

pdpmlm_t * pdpmlm_R_alloc(unsigned int size) {
  pdpmlm_t * obj;
  obj         = (pdpmlm_t *) R_alloc(1, sizeof(pdpmlm_t));
  obj->datasq = (double *) R_alloc(size, sizeof(double));
  obj->sum    = (double *) R_alloc(size, sizeof(double));
  obj->sumsq  = (double *) R_alloc(size, sizeof(double));
  obj->unique = (unsigned int *) R_alloc(size, sizeof(unsigned int));
  obj->indmap = (unsigned int *) R_alloc(size, sizeof(unsigned int));
  obj->size   = size;
  return(obj);
}
