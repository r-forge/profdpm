#include "pdpmlm.h"

void pdpmlm_divy( pdpmlm_t * obj, unsigned int ncl ) {
  unsigned int i, ccl = 0;
  
  // 1. initialize vcl, mcl, compute pcl, ncl
  obj->ncl = ncl;
  for( i = 0; i < obj->ngr; i++ ) { 
    obj->vcl[ i ] = i % ncl; 
    obj->pcl[ i ]++;
  }
    
  // 2. make calls to pdpmlm_add
  for( i = 0; i < obj->ngr; i++ ) {
    pdpmlm_add( obj, i, obj->vcl[ i ] );
  }
}

void pdpmlm_add( pdpmlm_t * obj, unsigned grp, unsigned int cls ) {
  unsigned int i, j;

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

  // 1. set vcl, recompute pcl, and possibly ncl and mcl
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
  F77_CALL(dgemv)( "N", &obj->q, &obj->q, &done, s, &obj->q, m, &obj->q, &dzero, obj->buf, &obj->q );  // obj->buf = s*m
  *b += F77_CALL(ddot)( &obj->q, m, &ione, obj->buf, &ione ); // b += m'obj->buf

  // 5. a = a0 + nk;
  *a = obj->a0 + obj->pcl[ cls ];
}


double pdpmlm_logp( pdpmlm_t * obj ) {
  unsigned int i, cls = 0;
  double logp = obj->alp * log( obj->ncl ) - lfactorial( obj->ncl );
  for( i = 0; i < obj->ncl; i++ ) {
    while( obj->pcl[ cls ] == 0 ) { cls++; }
    pdpmlm_parm( obj, cls, obj->s, obj->m, &obj->a, &obj->b );
    logp += lgamma( obj->a / 2 ) - ( obj->a / 2 ) * log( obj->b / 2 );
    cls++;
  }
  return logp;
}

void pdpmlm_best( pdpmlm_t * obj, unsigned int grp ) {
  unsigned int i, free = 1, tried = 0, test_cls = 0, best_cls;
  double test_logp, best_logp;

  best_logp = pdpmlm_logp( obj );
  best_cls = obj->vcl[ grp ];

  // If the current cluster is composed only of the group in question,
  // then a free cluster need not be tried (tried == 1 indicates a free
  // cluster has bee tried. free add to the total tried clusters )
  if( obj->pcl[ best_cls ] == obj->pgr[ grp ] ) { tried = 1; free = 0; }

  for( i = 0; i < obj->ncl + free; i++ ) {
    while( test_cls == best_cls || ( obj->pcl[ test_cls ] == 0 && tried == 1 ) ) { test_cls++; }
    if( obj->pcl[ test_cls ] == 0 ) { tried == 1; }
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
  unsigned int *vcl_old, *grps, ngrps, iter = 0;
  double lpost_old;

  // 0. select number of groups to be shuffled (%20 of total)
  ngrps = obj->ngr / 5;
  ngrps = ngrps == 0 ? 1 : ngrps;

  // 1. compute old logp, allocate memory for vcl_old, grps
  lpost_old = pdpmlm_lpost( obj ); 
  vcl_old = pdpmlm_alloc( ngrps, sizeof( unsigned int ) );
  grps    = pdpmlm_alloc( ngrps, sizeof( unsigned int ) );

  while( iter++ < itermax ) {
  // 2. randomly select ngrps groups to shuffle, save indicators
  // 3. move groups to same cluster
  // 4. move each group to best cluster
  // 5. compute logp, keep new clustering if better, else revert to old
  }
}


void pdpmlm_Rdump( pdpmlm_t * obj ) {
  unsigned int i, j, cls;
  // 1. print variables  
  pdpmlm_printf("vgr\n");
  printIntegerVector( obj->vgr, obj->p, 1 );
  pdpmlm_printf("pgr\n");
  printIntegerVector( obj->pgr, obj->ngr, 1 ); // p are allocated
  pdpmlm_printf("ngr\n");
  printIntegerVector( &obj->ngr, 1, 1 );
  pdpmlm_printf("vcl\n");
  printIntegerVector( obj->vcl, obj->ngr, 1 );
  pdpmlm_printf("pcl\n");
  printIntegerVector( obj->pcl, obj->ncl, 1 );
  pdpmlm_printf("ncl\n");
  printIntegerVector( &obj->ncl, 1, 1 );
  pdpmlm_printf("y\n");
  printRealVector( obj->y, obj->p, 1 );
  pdpmlm_printf("x\n");
  printRealVector( obj->x, obj->p*obj->q, 1 );
  pdpmlm_printf("p\n");
  printIntegerVector( &obj->p, 1, 1 );
  pdpmlm_printf("q\n");
  printIntegerVector( &obj->q, 1, 1 );
  pdpmlm_printf("s\n");
  printRealVector( obj->s, obj->q*obj->q, 1 );
  pdpmlm_printf("m\n");
  printRealVector( obj->m, obj->q, 1 );
  pdpmlm_printf("a\n");
  printRealVector( &obj->a, 1, 1 );
  pdpmlm_printf("b\n");
  printRealVector( &obj->b, 1, 1 );
  
  for(i = 0; i < obj->ngr; i++) {
    pdpmlm_printf( "group %u\n", i );
    pdpmlm_printf( "xxgr\n" );
    printRealVector( obj->xxgr[ i ], obj->q * obj->q, 1 );
    pdpmlm_printf( "xygr\n" );
    printRealVector( obj->xygr[ i ], obj->q, 1 );
    pdpmlm_printf( "yygr\n" );
    printRealVector( &obj->yygr[ i ], 1, 1 );
  }

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
}
