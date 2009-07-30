#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
SEXP sandbox(SEXP x) {
  SEXP out; 
  int i;
 
  PROTECT(out = allocMatrix(REALSXP, 2, 2));
  for(i = 0; i < 4; i++) { REAL(out)[i] = i; }
  UNPROTECT(1);
  return(out);
}
/*
SEXP sandbox(SEXP x) {
  double *dx = REAL(x);
  int i, lx;
  int  m, n, k, incx, incy;
  double alpha;
  
  lx = LENGTH(x);

  m      = lx;   // rows of op( a ) and out
  n      = lx;   // cols of op( b ) and out
  k      = 1;    // cols of op( a ) and rows of op( b )
  alpha  = 1;
  incx   = 1;
  incy   = 1;
  F77_CALL(daxpy)(&n, &alpha, dx, &incx, dx, &incy); // obj->xx += x[p,]'y[p]

  return x;
}
*/
/*
  double *dx = REAL(x), *dout;
  int i, lx;
  SEXP out;
  char transa, transb;
  int  m, n, k, lda, ldb, ldout;
  double alpha, beta;
  
  lx = LENGTH(x);

  PROTECT(out = allocMatrix(REALSXP, lx, lx));
  dout   = REAL(out);
  transa = 'N';
  transb = 'T';  // 'N' op( a ) = a, 'T' op( a ) = a'
  m      = lx;   // rows of op( a ) and out
  n      = lx;   // cols of op( b ) and out
  k      = 1;    // cols of op( a ) and rows of op( b )
  alpha  = 1;
  lda    = lx;   // first dimension of a 
  ldb    = lx;   // first dimension of b 
  beta   = 0;    
  ldout  = lx;   // first dimension of c
  F77_CALL(dgemm)(&transa,&transb,&m,&n,&k,&alpha,dx,&lda,dx,&ldb,&beta,dout,&ldout);

  // out = alpha*op( y )*op( x ) + beta*out
  UNPROTECT(1); 
  return out;
*/

/*
SEXP sandbox(SEXP x) {
  double *dx = REAL(x), *dout;
  int i, lx;
  SEXP out;
  char transa, transb;
  int  m, n, k, lda, ldb, ldout;
  double alpha, beta;
  
  lx = LENGTH(x);

  PROTECT(out = allocMatrix(REALSXP, lx, lx));
  dout   = REAL(out);
  transa = 'N';
  transb = 'T';  // 'N' op( a ) = a, 'T' op( a ) = a'
  m      = lx;   // rows of op( a ) and out
  n      = lx;   // cols of op( b ) and out
  k      = 1;    // cols of op( a ) and rows of op( b )
  alpha  = 1;
  lda    = lx;   // first dimension of a 
  ldb    = lx;   // first dimension of b 
  beta   = 0;    
  ldout  = lx;   // first dimension of c
  F77_CALL(dgemm)(&transa,&transb,&m,&n,&k,&alpha,dx,&lda,dx,&ldb,&beta,dout,&ldout);

  // out = alpha*op( y )*op( x ) + beta*out
  UNPROTECT(1); 
  return out;
}
*/
/*
SEXP sandbox(SEXP y, SEXP x) {
  double *dy = REAL(y), *dx = REAL(x), *dout;
  int i, ry, cy, rx, cx;
  SEXP dim, out;
  char transy, transx;
  int  m, n, k, ldy, ldx, ldout;
  double alpha, beta;

  dim = getAttrib(y, R_DimSymbol);
  ry  = INTEGER(dim)[0];
  cy  = INTEGER(dim)[1];
  dim = getAttrib(x, R_DimSymbol);
  rx  = INTEGER(dim)[0];
  cx  = INTEGER(dim)[1]; 

  PROTECT(out = allocMatrix(REALSXP, ry, cx));
  dout   = REAL(out);
  transy = 'N';
  transx = 'N';  // 'N' op( y ) = y, 'T' op( y ) = y'
  m      = ry;   // rows of op( y ) and out
  n      = cx;   // cols of op( x ) and out
  k      = cy;   // cols of op( y ) and rows of op( x )
  alpha  = 1;
  ldy    = ry;   // first dimension of y (ry if 'N', cy if 'T')
  ldx    = rx;   // first dimension of x (rx if 'N', cx if 'T')
  beta   = 0;    
  ldout  = ry;   // first dimension of out (ry if 'N', cy if 'T')
  F77_CALL(dgemm)(&transy,&transx,&m,&n,&k,&alpha,dy,&ldy,dx,&ldx,&beta,dout,&ldout);

  // out = alpha*op( y )*op( x ) + beta*out
  UNPROTECT(1); 
  return out;
}
*/
