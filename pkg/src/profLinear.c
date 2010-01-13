#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "util.h"
#include "pdpmlm.h"


SEXP profLinear(SEXP y, SEXP x, SEXP group, SEXP clust, SEXP param, SEXP method, SEXP maxiter, SEXP crit, SEXP prior, SEXP verbose) {
  SEXP retval, elem, names, class, dim;
  pdpmlm_t * obj;
  int i, j, k, cls, onei=1; 
  double *xp, *yp, oned=1.0;

  //0. setup the return value 
  PROTECT(retval = allocVector(VECSXP, 10));
  PROTECT(names = allocVector(STRSXP, 10));
  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(names, 0, mkChar("y"));
  SET_STRING_ELT(names, 1, mkChar("x"));
  SET_STRING_ELT(names, 2, mkChar("group"));
  SET_STRING_ELT(names, 3, mkChar("param"));
  SET_STRING_ELT(names, 4, mkChar("clust"));
  SET_STRING_ELT(names, 5, mkChar("a"));
  SET_STRING_ELT(names, 6, mkChar("b"));
  SET_STRING_ELT(names, 7, mkChar("m"));
  SET_STRING_ELT(names, 8, mkChar("s"));
  SET_STRING_ELT(names, 9, mkChar("logp"));
  SET_STRING_ELT(class, 0, mkChar("profLinear"));
  setAttrib(retval, R_NamesSymbol, names);
  setAttrib(retval, R_ClassSymbol, class);
  SET_VECTOR_ELT(retval, 0, y);
  SET_VECTOR_ELT(retval, 1, x);
  SET_VECTOR_ELT(retval, 2, group);
  SET_VECTOR_ELT(retval, 3, param);
  SET_VECTOR_ELT(retval, 4, allocVector(INTSXP, LENGTH(y)));

  //1. Allocate obj, make assignments, check priors
  obj = (pdpmlm_t *) R_alloc( 1, sizeof(pdpmlm_t) );
  obj->mem = sizeof(pdpmlm_t);

  //1.1 Set flags
  obj->flags   = 0;
  if( LOGICAL(verbose)[0] )    { obj->flags |= FLAG_VERBOSE; }
  if( INTEGER(prior)[0] == 1 ) { obj->flags |= FLAG_PRICLUS; }

  //1.2 Set pointers to data
  obj->y     = REAL(y);
  obj->x     = REAL(x);
  obj->vgr   = (unsigned int *) INTEGER(group);
  dim        = getAttrib(x, R_DimSymbol); 
  obj->p     = INTEGER(dim)[ 1 ];
  obj->q     = INTEGER(dim)[ 0 ];

  //1.3 Check values in prior list
  elem       = getListElementByName(param, "alpha");
  if( elem == R_NilValue ) {
    warning( "list item \"alpha\" missing from param, using default value" );
    obj->alp = DEFAULT_ALP;
  } else { obj->alp = REAL(elem)[0]; }
  elem       = getListElementByName(param, "s0");
  if( elem == R_NilValue ) {
    warning( "list item \"s0\" missing from param, using default value" );
    obj->s0 = DEFAULT_S0;
  } else { obj->s0 = REAL(elem)[0]; }
  elem       = getListElementByName(param, "m0");
  if( elem == R_NilValue ) {
    warning( "list item \"m0\" missing from param, using default values" );
    obj->m0 = (double *) pdpmlm_alloc( obj, obj->q, sizeof(double) ); 
    for( i = 0; i < obj->q; i++ ) { obj->m0[i] = DEFAULT_M0; }
  } else if ( LENGTH(elem) < obj->q ) {
    warning( "list item \"m0\" should be of length ncol(x), using default values" );
    obj->m0 = (double *) pdpmlm_alloc( obj, obj->q, sizeof(double) );
    for( i = 0; i < obj->q; i++ ) { obj->m0[i] = DEFAULT_M0; }
  } else { obj->m0 = REAL(elem); }
  elem       = getListElementByName(param, "a0");
  if( elem == R_NilValue ) { 
    warning( "list item \"a0\" missing, using default value" );
    obj->a0 = DEFAULT_A0;
  } else { obj->a0 = REAL(elem)[0]; }
  elem       = getListElementByName(param, "b0");
  if( elem == R_NilValue ) {
    warning( "list item \"b0\" missing, using default value" );
    obj->b0 = DEFAULT_B0;
  } else { obj->b0 = REAL(elem)[0]; }
 
  //2. Allocate memory for pgr
  obj->pgr = (unsigned int *) pdpmlm_alloc( obj, obj->p, sizeof(unsigned int) );

  //3. Compute pgr, ngr
  obj->ngr = 0;
  for( i = 0; i < obj->p; i++ ) { obj->pgr[ i ] = 0; }
  for( i = 0; i < obj->p; i++ ) { obj->pgr[ obj->vgr[ i ] ]++; }
  for( i = 0; i < obj->p; i++ ) { if( obj->pgr[ i ] > 0 ) { obj->ngr++; } }
 
  //4. Allocate and zero memory vcl, pcl, ncl
  obj->ncl = 0;
  obj->vcl = (unsigned int *) pdpmlm_alloc( obj, obj->ngr, sizeof(unsigned int) );
  obj->pcl = (unsigned int *) pdpmlm_alloc( obj, obj->ngr, sizeof(unsigned int) );
  for( i = 0; i < obj->ngr; i++ ) { 
    obj->vcl[ i ] = BAD_CLS;
    obj->pcl[ i ] = 0; 
  }

  //5. Allocate and zero memory for xxgr xygr, and yygr
  obj->xxgr = (double **) pdpmlm_alloc( obj, obj->ngr, sizeof(double *) );
  obj->xygr = (double **) pdpmlm_alloc( obj, obj->ngr, sizeof(double *) );
  obj->yygr = (double *)  pdpmlm_alloc( obj, obj->ngr, sizeof(double) );
  for( i = 0; i < obj->ngr; i++ ) {
    obj->xxgr[ i ] = (double *) pdpmlm_alloc( obj, obj->q * obj->q, sizeof(double) );
    obj->xygr[ i ] = (double *) pdpmlm_alloc( obj, obj->q, sizeof(double) );
    obj->yygr[i] = 0.0;
    for( j = 0; j < obj->q; j++ ) {
      obj->xygr[ i ][ j ] = 0.0;
      for( k = 0; k < obj->q; k++ ) {
        obj->xxgr[ i ][ k + j*obj->q ] = 0.0;
      }
    }
  }
  
  //6. Compute xxgr, xygr, yygr
  for( i = 0; i < obj->p; i++ ) {
       xp = obj->x + i*obj->q;
       yp = obj->y + i;
       F77_CALL(dgemm)("N", "T",
                      (int *) &obj->q, (int *) &obj->q, &onei, &oned,
                      xp, (int *) &obj->q, xp, (int *) &obj->q,
                      &oned, obj->xxgr[ obj->vgr[ i ] ], (int *) &obj->q); 
       F77_CALL(daxpy)((int *) &obj->q, yp, xp, &onei, obj->xygr[ obj->vgr[ i ] ], &onei); 
       obj->yygr[ obj->vgr[ i ] ] += (*yp) * (*yp);
  }
  
  //7. allocate and zero xxcl, xycl, yycl
  obj->xxcl = (double **) pdpmlm_alloc( obj, obj->ngr, sizeof(double *) );
  obj->xycl = (double **) pdpmlm_alloc( obj, obj->ngr, sizeof(double *) );
  obj->yycl = (double *)  pdpmlm_alloc( obj, obj->ngr, sizeof(double) );
  for( i = 0; i < obj->ngr; i++ ) {
    obj->xxcl[ i ] = NULL;
    obj->xycl[ i ] = NULL;
    obj->yycl[ i ] = 0.0;
  }

  //8. allocate s, m, fbuf, and pbuf
  obj->s = (double *) pdpmlm_alloc( obj, obj->q * obj->q, sizeof(double) );
  obj->m = (double *) pdpmlm_alloc( obj, obj->q, sizeof(double) );
  obj->fbuf = (double *) pdpmlm_alloc( obj, obj->q, sizeof(double) );
  obj->pbuf = (unsigned int *) pdpmlm_alloc( obj, obj->ngr, sizeof(unsigned int) );

  //9. distribute clusters initially and perform optimization
  if( isInteger(clust) ) {
    i = 0;
    for( j = 0; j < obj->ngr; j++ ) {
      pdpmlm_add( obj, j, INTEGER(clust)[i] );
      i += obj->pgr[ j ];
    }
  } 

  if( INTEGER(method)[0] == METHOD_NONE ) {
    if( isLogical(clust) ) { pdpmlm_divy( obj ); }
  }
  else if( INTEGER(method)[0] == METHOD_STOCH ) {
    if( isLogical(clust) ) { pdpmlm_divy( obj ); }
    pdpmlm_stoch( obj, INTEGER(maxiter)[0], REAL(crit)[0] );
  }
  else if( INTEGER(method)[0] == METHOD_AGGLO ) {
    if( isLogical(clust) ) { for( i = 0; i < obj->ngr; i++ ) { pdpmlm_add( obj, i, i ); } }
    pdpmlm_agglo( obj, INTEGER(maxiter)[0] );
  }

  if( obj->flags & FLAG_VERBOSE ) {
    pdpmlm_printf( "allocated memory: %fMb\n", obj->mem/1000000.0 );
  }

  //10. complete the return value
  SET_VECTOR_ELT(retval, 5, allocVector(REALSXP, obj->ncl)); //a
  SET_VECTOR_ELT(retval, 6, allocVector(REALSXP, obj->ncl)); //b
  SET_VECTOR_ELT(retval, 7, allocVector(VECSXP, obj->ncl)); //m
  SET_VECTOR_ELT(retval, 8, allocVector(VECSXP, obj->ncl)); //s
  SET_VECTOR_ELT(retval, 9, allocVector(REALSXP, 1)); //logp
  REAL(VECTOR_ELT(retval, 9))[0] = pdpmlm_logp( obj );

  for( i = 0; i < obj->ngr; i++ ) { obj->pbuf[ i ] = BAD_CLS; }
  cls = 1;
  for( i = 0; i < obj->p; i++ ) { 
    if( obj->pbuf[ obj->vcl[ obj->vgr[ i ] ] ] == BAD_CLS ) { 
      obj->pbuf[ obj->vcl[ obj->vgr[ i ] ] ] = cls++;
    }
    INTEGER(VECTOR_ELT(retval, 4))[i] = obj->pbuf[ obj->vcl[ obj->vgr[ i ] ] ];
  }
  cls = 0;
  for( i = 0; i < obj->ncl; i++) {
    while( obj->pcl[ cls ] == 0 ) { cls++; }
    SET_VECTOR_ELT(VECTOR_ELT(retval, 8), obj->pbuf[ cls ]-1, allocVector(REALSXP, obj->q*obj->q));
    SET_VECTOR_ELT(VECTOR_ELT(retval, 7), obj->pbuf[ cls ]-1, allocVector(REALSXP, obj->q));
    PROTECT(dim = allocVector(INTSXP, 2));
    INTEGER(dim)[0] = obj->q;
    INTEGER(dim)[1] = obj->q;
    setAttrib(VECTOR_ELT(VECTOR_ELT(retval, 8), obj->pbuf[ cls ]-1), R_DimSymbol, dim);

    pdpmlm_parm( obj, cls,
        REAL(VECTOR_ELT(VECTOR_ELT(retval, 8), obj->pbuf[ cls ]-1)),
        REAL(VECTOR_ELT(VECTOR_ELT(retval, 7), obj->pbuf[ cls ]-1)),
        REAL(VECTOR_ELT(retval, 5))+obj->pbuf[ cls ]-1,
        REAL(VECTOR_ELT(retval, 6))+obj->pbuf[ cls ]-1
    );
    cls++;
  }

  UNPROTECT(3+obj->ncl);
  return retval;
}
