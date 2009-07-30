#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "util.h"
#include "pdpmlm.h"

void debug( char * msg ) { Rprintf("debug: %s\n", msg); }
void memerror() { error("failed to allocate memory"); }

SEXP profLinear(SEXP y, SEXP x, SEXP group, SEXP parm, SEXP iter, SEXP crit) {
  SEXP retval, elem, names, class, clust, dim;
  pdpmlm_t * obj;
  unsigned int i, j, k, cls, onei=1;
  double s, m, a, b; 
  double *xp, *yp, oned=1.0, zerod=0.0;

  //0. setup the return value 
  PROTECT(retval = allocVector(VECSXP, 9));
  PROTECT(names = allocVector(STRSXP, 9));
  PROTECT(class = allocVector(STRSXP, 1));
  PROTECT(clust = allocVector(INTSXP, LENGTH(y)));
  SET_STRING_ELT(names, 0, mkChar("y"));
  SET_STRING_ELT(names, 1, mkChar("x"));
  SET_STRING_ELT(names, 2, mkChar("group"));
  SET_STRING_ELT(names, 3, mkChar("parm"));
  SET_STRING_ELT(names, 4, mkChar("clust"));
  SET_STRING_ELT(names, 5, mkChar("a"));
  SET_STRING_ELT(names, 6, mkChar("b"));
  SET_STRING_ELT(names, 7, mkChar("m"));
  SET_STRING_ELT(names, 8, mkChar("s"));
  SET_STRING_ELT(class, 0, mkChar("profLinear"));
  setAttrib(retval, R_NamesSymbol, names);
  setAttrib(retval, R_ClassSymbol, class);
  SET_VECTOR_ELT(retval, 0, y);
  SET_VECTOR_ELT(retval, 1, x);
  SET_VECTOR_ELT(retval, 2, group);
  SET_VECTOR_ELT(retval, 3, parm);
  SET_VECTOR_ELT(retval, 4, clust);

  //1. Allocate memory for obj and pgr
  obj = (pdpmlm_t *) pdpmlm_alloc( 1, sizeof(pdpmlm_t) );
  if( obj == NULL ) { memerror(); }
  obj->pgr = (unsigned int *) pdpmlm_alloc( obj->p, sizeof(unsigned int) );
  if( obj->pgr == NULL ) { memerror(); }
  debug("step 1 ok");

  //2. Make assignments, check priors
  obj->y     = REAL(y);
  obj->x     = REAL(x);
  obj->vgr   = INTEGER(group);
  dim        = getAttrib(x, R_DimSymbol); 
  obj->p     = INTEGER(dim)[ 1 ];
  obj->q     = INTEGER(dim)[ 0 ];
  elem       = getListElementByName(parm, "alpha");
  if( elem == R_NilValue ) {
    warning( "list item \"alpha\" missing from parm, using default value" );
    obj->alp = DEFAULT_ALP;
  } else { obj->alp = *(REAL(elem)); }
  elem       = getListElementByName(parm, "s0");
  if( elem == R_NilValue ) {
    warning( "list item \"s0\" missing from parm, using default value" );
    obj->s0 = DEFAULT_S0;
  } else { obj->s0 = *(REAL(elem)); }
  elem       = getListElementByName(parm, "m0");
  if( elem == R_NilValue ) {
    warning( "list item \"m0\" missing from parm, using default values" );
    obj->m0 = (double *) pdpmlm_alloc( obj->q, sizeof(double) );
    for( i = 0; i < obj->q; i++ ) { obj->m0[i] = DEFAULT_M0; }
  } else if ( LENGTH(elem) < obj->q ) {
    warning( "list item \"m0\" should be of length ncol(x), using default values" );
    obj->m0 = (double *) pdpmlm_alloc( obj->q, sizeof(double) );
    for( i = 0; i < obj->q; i++ ) { obj->m0[i] = DEFAULT_M0; }
  } else { obj->m0 = REAL(elem); }
  elem       = getListElementByName(parm, "a0");
  if( elem == R_NilValue ) { 
    warning( "list item \"a0\" missing, using default value" );
    obj->a0 = DEFAULT_A0;
  } else { obj->a0 = *(REAL(elem)); }
  elem       = getListElementByName(parm, "b0");
  if( elem == R_NilValue ) {
    warning( "list item \"b0\" missing, using default value" );
    obj->b0 = DEFAULT_B0;
  } else { obj->b0 = *(REAL(elem)); }
  debug("step 2 ok");

  //3. Compute pgr, ngr
  obj->ngr = 0;
  for( i = 0; i < obj->p; i++ ) { obj->pgr[ i ] = 0; }
  for( i = 0; i < obj->p; i++ ) { obj->pgr[ obj->vgr[ i ] ]++; }
  for( i = 0; i < obj->p; i++ ) { if( obj->pgr[ i ] > 0 ) { obj->ngr++; } }
  debug("step 3 ok");
 
  //4. Allocate and zero memory vcl, pcl, ncl
  obj->ncl = 0;
  obj->vcl = (unsigned int *) pdpmlm_alloc( obj->ngr, sizeof(unsigned int) );
  obj->pcl = (unsigned int *) pdpmlm_alloc( obj->ngr, sizeof(unsigned int) );
  if( obj->pcl == NULL ) { memerror(); }
  for( i = 0; i < obj->ngr; i++ ) { 
    obj->vcl[ i ] = BAD_CLS;
    obj->pcl[ i ] = 0; 
  }
  debug("step 4 ok");

  //5. Allocate and zero memory for xxgr xygr, and yygr
  obj->xxgr = (double **) pdpmlm_alloc( obj->ngr, sizeof(double *) );
  if( obj->xxgr == NULL ) { memerror(); }
  obj->xygr = (double **) pdpmlm_alloc( obj->ngr, sizeof(double *) );
  if( obj->xygr == NULL ) { memerror(); }
  obj->yygr = (double *) pdpmlm_alloc( obj->ngr, sizeof(double) );
  if( obj->yygr == NULL ) { memerror(); }
  for( i = 0; i < obj->ngr; i++ ) {
    obj->xxgr[ i ] = (double *) pdpmlm_alloc( obj->q * obj->q, sizeof(double) );
    if( obj->xxgr[ i ] == NULL ) { memerror(); }
    obj->xygr[i] = (double *) pdpmlm_alloc( obj->q, sizeof(double) );
    if( obj->xygr[ i ] == NULL ) { memerror(); }
    obj->yygr[i] = 0.0;
    for( j = 0; j < obj->q; j++ ) {
      obj->xygr[ i ][ j ] = 0.0;
      for( k = 0; k < obj->q; k++ ) {
        obj->xxgr[ i ][ k + j*obj->q ] = 0.0;
      }
    }
  }
  debug("step 5 ok");
  
  //6. Compute xxgr, xygr, yygr
  for( i = 0; i < obj->p; i++ ) {
       xp = obj->x + i*obj->q;
       yp = obj->y + i;
       F77_CALL(dgemm)("N","T",
                      &obj->q,&obj->q,&onei,&oned,
                      xp,&obj->q,xp,&obj->q,
                      &oned,obj->xxgr[ obj->vgr[ i ] ],&obj->q); 
       F77_CALL(daxpy)(&obj->q, yp, xp, &onei, obj->xygr[ obj->vgr[ i ] ], &onei); 
       obj->yygr[ obj->vgr[ i ] ] += (*yp) * (*yp);
  }
  debug("step 6 ok");

  //7. allocate and zero xxcl, xycl, yycl
  obj->xxcl = (double **) pdpmlm_alloc( obj->ngr, sizeof(double *) );
  if( obj->xxcl == NULL ) { memerror(); }
  obj->xycl = (double **) pdpmlm_alloc( obj->ngr, sizeof(double *) );
  if( obj->xycl == NULL ) { memerror(); }
  obj->yycl = (double *) pdpmlm_alloc( obj->ngr, sizeof(double) );
  if( obj->yycl == NULL ) { memerror(); }
  for( i = 0; i < obj->ngr; i++ ) {
    obj->xxcl[ i ] = NULL;
    obj->xycl[ i ] = NULL;
    obj->yycl[ i ] = 0.0;
  }
  debug("step 7 ok");

  //8. allocate s, m, and buf
  obj->s = (double *) pdpmlm_alloc( obj->q * obj->q, sizeof(double) );
  if( obj->s == NULL ) { memerror(); }
  obj->m = (double *) pdpmlm_alloc( obj->q, sizeof(double) );
  if( obj->m == NULL ) { memerror(); }
  obj->buf = (double *) pdpmlm_alloc( obj->q, sizeof(double) );
  if( obj->buf == NULL ) { memerror(); }
  debug("step 8 ok");

  pdpmlm_add( obj, 0, 0 );
  pdpmlm_add( obj, 1, 1 );
  pdpmlm_add( obj, 2, 2 );
  pdpmlm_add( obj, 3, 3 );
  pdpmlm_Rdump( obj );

  SET_VECTOR_ELT(retval, 5, allocVector(REALSXP, obj->ncl)); //a
  SET_VECTOR_ELT(retval, 6, allocVector(REALSXP, obj->ncl)); //b
  SET_VECTOR_ELT(retval, 7, allocVector(VECSXP, obj->ncl)); //m
  SET_VECTOR_ELT(retval, 8, allocVector(VECSXP, obj->ncl)); //s
  for( i = 0; i < obj->p; i++ ) { 
    INTEGER(VECTOR_ELT(retval, 4))[i] = obj->vcl[ obj->vgr[ i ] ]; 
  }
  cls = 0;
  for( i = 0; i < obj->ncl; i++) {
    while( obj->pcl[ cls ] == 0 ) { cls++; }
    SET_VECTOR_ELT(VECTOR_ELT(retval, 8), i, allocVector(REALSXP, obj->q*obj->q));
    SET_VECTOR_ELT(VECTOR_ELT(retval, 7), i, allocVector(REALSXP, obj->q));
    dim = allocVector(INTSXP, 2);
    INTEGER(dim)[0] = obj->q;
    INTEGER(dim)[1] = obj->q;
    setAttrib(VECTOR_ELT(VECTOR_ELT(retval, 8), i), R_DimSymbol, dim);
    pdpmlm_parm( obj, cls,
        REAL(VECTOR_ELT(VECTOR_ELT(retval, 8), i)),
        REAL(VECTOR_ELT(VECTOR_ELT(retval, 7), i)),
        REAL(VECTOR_ELT(retval, 5))+i,
        REAL(VECTOR_ELT(retval, 6))+i
    );
    cls++;
  }

  UNPROTECT(4);
  return(retval);
}
