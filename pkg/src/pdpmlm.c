#include "pdpmlm.h"


void memerror() { error("failed to allocate memory"); }

void pdpmlm_divy( pdpmlm_t * obj, unsigned int ncl ) {
  unsigned int i;
  for( i = 0; i < obj->ngr; i++ ) { 
    pdpmlm_add( obj, i, i % ncl );
  }
}

void pdpmlm_add( pdpmlm_t * obj, unsigned int grp, unsigned int cls ) {
  unsigned int i, j;
 
  if( grp >= obj->ngr ) { pdpmlm_printf("grp: %u\n", grp); error( "pdpmlm_add: invalid argument" ); }
  if( cls >= obj->ngr ) { pdpmlm_printf("cls: %u\n", cls); error( "pdpmlm_add: invalid argument" ); }
 
  // 1. set vcl, recompute pcl, and possibly ncl
  obj->vcl[ grp ] = cls;
  if( obj->pcl[ cls ] == 0 ) { obj->ncl++; }
  obj->pcl[ cls ] += obj->pgr[ grp ];

  // 2. allocate memory for xxcl and xycl if necessary, zero xxcl, xycl
  if( obj->xxcl[ cls ] == NULL ) {
    obj->xxcl[ cls ] = (double *) pdpmlm_alloc( obj->q * obj->q, sizeof(double) );
    if( obj->xxcl[ cls ] == NULL ) { error("memory allocation failed"); }
    obj->xycl[ cls ] = (double *) pdpmlm_alloc( obj->q, sizeof(double) );
    if( obj->xycl[ cls ] == NULL ) { error("memory allocation failed"); }
    obj->yycl[ cls ] = 0.0;
    for( i = 0; i < obj->q; i++ ) {
      obj->xycl[ cls ][ i ] = 0.0;
      for( j = 0; j < obj->q; j++ ) {
        obj->xxcl[ cls ][ j + i*obj->q ] = 0.0;
      }
    }
  }

  // 3. recompute xxcl, xycl yycl for cls using xxgr, xygr, and yygr from grp
  obj->yycl[ cls ] += obj->yygr[ grp ];
  for( i = 0; i < obj->q; i++ ) {
    obj->xycl[ cls ][ i ] += obj->xygr[ grp ][ i ];
    for( j = 0; j < obj->q; j++ ) {
      obj->xxcl[ cls ][ j + i * obj->q ] += obj->xxgr[ grp ][ j + i * obj->q ];
    }
  }
}

void pdpmlm_sub( pdpmlm_t * obj, unsigned grp, unsigned int cls ) {
  unsigned int i, j;

  if( grp >= obj->ngr || cls >= obj->ngr ) { error( "pdpmlm_sub: invalid argument" ); } 

  // 1. set vcl, recompute pcl, and possibly ncl
  obj->vcl[ grp ] = BAD_CLS; // comment this out after debug
  obj->pcl[ cls ] -= obj->pgr[ grp ];
  if( obj->pcl[ cls ] == 0 ) { obj->ncl--; }

  // 2. recompute xxcl, xycl yycl for cls using xxgr, xygr, and yygr from grp
  obj->yycl[ cls ] -= obj->yygr[ grp ];
  for( i = 0; i < obj->q; i++ ) {
    obj->xycl[ cls ][ i ] -= obj->xygr[ grp ][ i ];
    for( j = 0; j < obj->q; j++ ) {
      obj->xxcl[ cls ][ j + i * obj->q ] -= obj->xxgr[ grp ][ j + i * obj->q ];
    }
  }
}

void pdpmlm_move( pdpmlm_t * obj, unsigned int grp, unsigned int cls ) {
  unsigned int old = obj->vcl[ grp ];
  if( old == cls ) { return; }
  pdpmlm_sub( obj, grp, old );
  pdpmlm_add( obj, grp, cls );
}

void pdpmlm_parm( pdpmlm_t * obj, unsigned int cls, double * s, double * m, double * a, double * b ) {
  unsigned int i, j, d, ione=1, inf, info;
  double *x, *y, done=1.0, dzero=0.0;

  if( obj->pcl[ cls ] == 0 ) { return; }

  if( cls >= obj->ngr ) { error( "pdpmlm_parm: invalid argument" ); }

  // 1. load s with s0I + x'x, load m with s0m0 + x'y
  d = 0;
  for( i = 0; i < obj->q; i++ ) {
    m[ i ] = obj->s0 * obj->m0[ i ] + obj->xycl[ cls ][ i ];
    for( j = 0; j < obj->q; j++ ) {
      s[ i * obj->q + j ] = obj->xxcl[ cls ][ i * obj->q + j ];
    }
    s[ i * obj->q + (d++) ] += obj->s0;
  }
     
  // 2. m = (s)^(-1) * m
  // dgesv overwrites the matrix passed in. Hence, we must load s again
  // when the call is finished. obj->buf holds some temporary data.
  //  FIXME use dpbsv instead (for sym,pd matrices)
  F77_CALL(dgesv)(&obj->q, &ione, s, &obj->q, (unsigned int *) obj->buf, m, &obj->q, &info);
  if( info > 0 ) { error("dgesv: system is singular"); }
  if( info < 0 ) { error("dgesv: invalid argument"); }

  // 3. reload s
  d = 0;
  for( i = 0; i < obj->q; i++ ) {
    for( j = 0; j < obj->q; j++ ) {
      s[ i * obj->q + j ] = obj->xxcl[ cls ][ i * obj->q + j ];
    }
    s[ i * obj->q + (d++) ] += obj->s0;
  }

  // 4. b = y'y + s0*m0'm0 + m'sm
  *b = obj->yycl[ cls ];
  *b += obj->s0*F77_CALL(ddot)(&obj->q, obj->m0, &ione, obj->m0, &ione);  // b += s0*m0'm0
  F77_CALL(dgemv)( "N", &obj->q, &obj->q, &done, s, &obj->q, m, &ione, &dzero, obj->buf, &ione );  // obj->buf = s*m
  *b += F77_CALL(ddot)( &obj->q, m, &ione, obj->buf, &ione ); // b += m'obj->buf

  // 5. a = a0 + nk;
  *a = obj->a0 + obj->pcl[ cls ];
}


double pdpmlm_logp( pdpmlm_t * obj ) {
  unsigned int i, cls = 0;
  double logp = obj->ncl * log( obj->alp ) - lfactorial( obj->ncl );
  for( i = 0; i < obj->ncl; i++ ) {
    while( obj->pcl[ cls ] == 0 ) { cls++; }
    pdpmlm_parm( obj, cls, obj->s, obj->m, &obj->a, &obj->b );
    logp += lgamma( obj->a / 2 ) - ( obj->a / 2 ) * log( obj->b / 2 );
  //  logp += lfactorial( obj->pcl[ cls ] - 1 );
    cls++;
  }
  return logp;
}

unsigned int pdpmlm_free( pdpmlm_t * obj ) {
  unsigned int cls = 0;
  while( cls < obj->ngr && obj->pcl[ cls ] > 0 ) { cls++; }
  if( cls == obj->ngr ) { cls = BAD_CLS; }
  return cls;
}

void pdpmlm_best( pdpmlm_t * obj, unsigned int grp ) {
  unsigned int i, test_cls, best_cls;
  double test_logp, best_logp;

  if( grp >= obj->ngr ) { error( "pdpmlm_best: invalid argument" ); }

  best_logp = pdpmlm_logp( obj );
  best_cls = obj->vcl[ grp ];

  if( obj->pcl[ best_cls ] > obj->pgr[ grp ] ) {
    test_cls = pdpmlm_free( obj );
    if( test_cls == BAD_CLS ) { error("pdpmlm_best: test_cls should not == BAD_CLS"); }
    pdpmlm_move( obj, grp, test_cls );
    test_logp = pdpmlm_logp( obj );
    if( test_logp > best_logp ) { 
      best_logp = test_logp;
      best_cls  = test_cls;
    }
  }

  test_cls = 0;
  for( i = 0; i < obj->ncl; i++ ) {
    while( obj->pcl[ test_cls ] == 0 ) { test_cls++; }
    pdpmlm_move( obj, grp, test_cls );
    test_logp = pdpmlm_logp( obj );
    if( test_logp > best_logp ) { 
      best_logp = test_logp;
      best_cls  = test_cls;
    }
    test_cls++;
  }

  if( obj->vcl[ grp ] != best_cls ) { pdpmlm_move( obj, grp, best_cls ); }
}
    

void pdpmlm_chunk( pdpmlm_t * obj, unsigned int itermax) {
  unsigned int i, *vcl_old, *grps, ngrps, cls, iter = 0;
  double logp_old, logp;
  // 0. select number of groups to be shuffled (%50 of total)
  ngrps = obj->ngr / 2;
  ngrps = ngrps == 0 ? 1 : ngrps;

  // 1. allocate memory for vcl_old, grps
  vcl_old = (unsigned int *) pdpmlm_alloc( ngrps, sizeof( unsigned int ) );
  if( vcl_old == NULL ) { memerror(); }
  grps    = (unsigned int *) pdpmlm_alloc( ngrps, sizeof( unsigned int ) );
  if( grps == NULL ) { memerror(); }

  while( iter++ < itermax ) {
 
    // 2. randomly select ngrps groups to shuffle, save indicators
    GetRNGstate();
    for( i = 0; i < ngrps; i++ ) {
      grps[ i ] = (unsigned int) floor( obj->ngr * runif( 0.0, 1.0 ) );
      vcl_old[ i ] = obj->vcl[ grps[ i ] ];
    }
  
    // 3. compute old logp, move groups to same/new cluster
    logp_old = pdpmlm_logp( obj ); 
    //cls = pdpmlm_free( obj );
    //if( cls == BAD_CLS ) { cls = obj->vcl[ grps[ 0 ] ]; }
    cls = (unsigned int) floor( obj->ngr * runif( 0.0, 1.0 ) );
    for( i = 0; i < ngrps; i++ ) { pdpmlm_move( obj, grps[ i ], cls ); }
    PutRNGstate();
         
    // 4. move each group to best cluster
    for( i = 0; i < ngrps; i++ ) { pdpmlm_best( obj, grps[ i ] ); }

    // 5. compute logp, keep new clustering if better, else revert to old
    logp = pdpmlm_logp( obj );
    if( logp < logp_old ) { 
      for( i = 0; i < ngrps; i++ ) { pdpmlm_move( obj, grps[ i ], vcl_old[ i ] ); }
    }
    //pdpmlm_printf("iter: %u, ncl: %u, logp: %f\n", iter, obj->ncl, pdpmlm_logp( obj ) );
  }
}


void pdpmlm_Rdump( pdpmlm_t * obj ) {
  unsigned int i, j, cls;
  // 1. print variables  
//  pdpmlm_printf("vgr\n");
//  printIntegerVector( obj->vgr, obj->p, 1 );
//  pdpmlm_printf("pgr\n");
//  printIntegerVector( obj->pgr, obj->ngr, 1 ); // p are allocated
//  pdpmlm_printf("ngr\n");
//  printIntegerVector( &obj->ngr, 1, 1 );
  pdpmlm_printf("vcl\n");
  printIntegerVector( obj->vcl, obj->ngr, 1 );
  pdpmlm_printf("pcl\n");
  printIntegerVector( obj->pcl, obj->ngr, 1 );
  pdpmlm_printf("ncl\n");
  printIntegerVector( &obj->ncl, 1, 1 );
//  pdpmlm_printf("y\n");
//  printRealVector( obj->y, obj->p, 1 );
//  pdpmlm_printf("x\n");
//  printRealVector( obj->x, obj->p*obj->q, 1 );
//  pdpmlm_printf("p\n");
//  printIntegerVector( &obj->p, 1, 1 );
//  pdpmlm_printf("q\n");
//  printIntegerVector( &obj->q, 1, 1 );
  pdpmlm_printf("s\n");
  printRealVector( obj->s, obj->q*obj->q, 1 );
  pdpmlm_printf("m\n");
  printRealVector( obj->m, obj->q, 1 );
  pdpmlm_printf("a\n");
  printRealVector( &obj->a, 1, 1 );
  pdpmlm_printf("b\n");
  printRealVector( &obj->b, 1, 1 );
  
/*
  for(i = 0; i < obj->ngr; i++) {
    pdpmlm_printf( "group %u\n", i );
    pdpmlm_printf( "xxgr\n" );
    printRealVector( obj->xxgr[ i ], obj->q * obj->q, 1 );
    pdpmlm_printf( "xygr\n" );
    printRealVector( obj->xygr[ i ], obj->q, 1 );
    pdpmlm_printf( "yygr\n" );
    printRealVector( &obj->yygr[ i ], 1, 1 );
  }
*/
/*
  cls = 0;
  for( i = 0; i < obj->ncl; i++ ) {
    while( obj->pcl[ cls ] == 0 ) { cls++; }
    pdpmlm_printf( "cluster %u\n", cls );    
    pdpmlm_printf( "xxcl\n" );
    printRealVector( obj->xxcl[ cls ], obj->q * obj->q, 1 );
    pdpmlm_printf( "xycl\n" );
    printRealVector( obj->xycl[ cls ], obj->q, 1 );
    pdpmlm_printf( "yycl\n" );
    printRealVector( &obj->yycl[ cls ], 1, 1 );
    cls++; 
  }
*/
}

