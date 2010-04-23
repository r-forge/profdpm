#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "util.h"
#include "pdpmb.h"

SEXP profBinary(SEXP y, SEXP clust, SEXP param, SEXP method,\
  SEXP maxiter, SEXP crit, SEXP verbose) {
  SEXP retval, elem, names, class, dim;
  pdpmb_t * obj;
  int i, j, k, cls, onei=1; 
  double *xp, *yp, oned=1.0;
  unsigned int *buf;

  //0. setup the return value 
  PROTECT(retval = allocVector(VECSXP, 6));
  PROTECT(names = allocVector(STRSXP, 6));
  PROTECT(class = allocVector(STRSXP, 1));
  SET_STRING_ELT(names, 0, mkChar("y"));
  SET_STRING_ELT(names, 1, mkChar("param"));
  SET_STRING_ELT(names, 2, mkChar("clust"));
  SET_STRING_ELT(names, 3, mkChar("a"));
  SET_STRING_ELT(names, 4, mkChar("b"));
  SET_STRING_ELT(names, 5, mkChar("logp"));
  SET_STRING_ELT(class, 0, mkChar("profBinary"));
  setAttrib(retval, R_NamesSymbol, names);
  setAttrib(retval, R_ClassSymbol, class);
  SET_VECTOR_ELT(retval, 0, y);
  SET_VECTOR_ELT(retval, 1, param);
  dim = getAttrib(y, R_DimSymbol); 
  SET_VECTOR_ELT(retval, 2, allocVector(INTSXP, INTEGER(dim)[ 0 ]));

  //1. Allocate obj, make assignments, check priors
  obj = (pdpmb_t *) R_alloc( 1, sizeof(pdpmb_t) );
  obj->mem = sizeof(pdpmb_t);

  //1.1 Set flags
  obj->flags   = 0;
  if( LOGICAL(verbose)[0] )    { obj->flags |= FLAG_VERBOSE; }

  //1.2 Set pointers to data
  obj->y     = INTEGER(y);
  obj->ngr   = INTEGER(dim)[ 0 ];
  obj->q     = INTEGER(dim)[ 1 ];

  //1.3 Check values in param list
  elem       = getListElementByName(param, "lambda");
  if( elem == R_NilValue ) { obj->flags != FLAG_DIRICHL; obj->lam = 0; }
  else if( REAL(elem)[0] < 0 || REAL(elem)[0] > 1 ) {
    warning( "list item \"lambda\" must be between zero and one, using default value" );
    obj->lam = DEFAULT_LAM;
  } else { obj->lam = REAL(elem)[0]; }
  elem       = getListElementByName(param, "alpha");
  if( elem == R_NilValue ) {
    obj->alp = DEFAULT_ALP;
  } else if( REAL(elem)[0] <= 0 ) {
    warning( "list item \"alpha\" must be positive, using default value" );
    obj->alp = DEFAULT_ALP;
  } else { obj->alp = REAL(elem)[0]; }
  elem       = getListElementByName(param, "a0");
  if( elem == R_NilValue ) { 
    obj->a0 = DEFAULT_A0;
  } else if( REAL(elem)[0] <= 0 ) {
    warning( "list item \"a0\" must be positive, using default value" );
    obj->a0 = DEFAULT_A0;
  } else { obj->a0 = REAL(elem)[0]; }
  elem       = getListElementByName(param, "b0");
  if( elem == R_NilValue ) {
    obj->b0 = DEFAULT_B0;
  } else if( REAL(elem)[0] < 0 ) {
    warning( "list item \"b0\" must be nonnegative, using default value" );
    obj->b0 = DEFAULT_B0;
  } else { obj->b0 = REAL(elem)[0]; }

  //2. (Placeholder) 

  //3. (Placeholder)

  //4. Allocate and zero memory vcl, gcl, and ncl
  obj->ncl = 0;
  obj->vcl = (unsigned int *) pdpmb_alloc( obj, obj->ngr, sizeof(unsigned int) );
  obj->gcl = (unsigned int *) pdpmb_zalloc( obj, obj->ngr, sizeof(unsigned int) );
  for( i = 0; i < obj->ngr; i++ ) { obj->vcl[ i ] = BAD_VCL; }

  //5. (Placeholder)
  
  //6. (Placeholder)
  
  //7. allocate and zero gqcl
  obj->gqcl = (unsigned int *) pdpmb_zalloc( obj, obj->ngr * obj->q, sizeof(double) );

  //8. (Placeholder)
  obj->pbuf = (unsigned int *) pdpmb_alloc( obj, obj->ngr, sizeof(unsigned int) ); 

  //9. distribute clusters initially and perform optimization
  if( isInteger(clust) ) {
    for( j = 0; j < obj->ngr; j++ ) {
      pdpmb_add( obj, j, INTEGER(clust)[j] );
    }
  } 
  
  //***********************************************************
  if( INTEGER(method)[0] == METHOD_NONE ) {
    if( isLogical(clust) ) { pdpmb_divy( obj ); }
    obj->logp = pdpmb_logp( obj );
  }
  else if( INTEGER(method)[0] == METHOD_STOCH ) {
    if( isLogical(clust) ) { pdpmb_divy( obj ); }
    GetRNGstate();
    pdpmb_stoch( obj, INTEGER(maxiter)[0], REAL(crit)[0] );
    PutRNGstate();
  }
  else if( INTEGER(method)[0] == METHOD_GIBBS ) {
    if( isLogical(clust) ) { pdpmb_divy( obj ); }
    GetRNGstate();
    pdpmb_gibbs( obj, INTEGER(maxiter)[0], REAL(crit)[0] );
    PutRNGstate();
  }
  else if( INTEGER(method)[0] == METHOD_AGGLO ) {
    if( isLogical(clust) ) { for( i = 0; i < obj->ngr; i++ ) { pdpmb_add( obj, i, i ); } }
    pdpmb_agglo( obj, INTEGER(maxiter)[0] );
  }
  //***********************************************************

  if( !(obj->flags & FLAG_OPTCRIT) ) { warning("optimization criterion not met"); }
  if( obj->flags & FLAG_VERBOSE ) {
    pdpmb_printf( "allocated memory: %fMb\n", obj->mem/1000000.0 );
  }

  //10. complete the return value
  SET_VECTOR_ELT(retval, 3, allocVector(VECSXP, obj->ncl)); //a
  SET_VECTOR_ELT(retval, 4, allocVector(VECSXP, obj->ncl)); //b
  SET_VECTOR_ELT(retval, 5, allocVector(REALSXP, 1)); //logp
  REAL(VECTOR_ELT(retval, 5))[0] = obj->logp;

  for( i = 0; i < obj->ngr; i++ ) { obj->pbuf[ i ] = BAD_VCL; }
  cls = 1;
  for( i = 0; i < obj->ngr; i++ ) { 
    if( obj->pbuf[ obj->vcl[ i ] ] == BAD_VCL ) { 
      obj->pbuf[ obj->vcl[ i ] ] = cls++;
    }
    INTEGER(VECTOR_ELT(retval, 2))[i] = obj->pbuf[ obj->vcl[ i ] ];
  }
  cls = 0;
  for( i = 0; i < obj->ncl; i++ ) {
    while( obj->gcl[ cls ] == 0 ) { cls++; }
    //FIXME obj->pbuf[ cls ] instead of i
    SET_VECTOR_ELT(VECTOR_ELT(retval, 3), i, allocVector(REALSXP, obj->q));
    SET_VECTOR_ELT(VECTOR_ELT(retval, 4), i, allocVector(REALSXP, obj->q));
    for( j = 0; j < obj->q; j++ ) {
      REAL(VECTOR_ELT(VECTOR_ELT(retval, 3), i))[j] = obj->a0 +\
        obj->gqcl[ cls * obj->q + j ];
      REAL(VECTOR_ELT(VECTOR_ELT(retval, 4), i))[j] = obj->b0 +\
        obj->gcl[ cls ] - obj->gqcl[ cls * obj->q + j ];
    }
  }
  UNPROTECT(3);
  return retval;
}
