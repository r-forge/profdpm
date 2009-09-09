#include "pdpmlm.h"

void memerror() { error("failed to allocate memory"); }

void pdpmlm_divy( pdpmlm_t * obj ) {
  unsigned int i, grp = 0, cls = 0, set = 0;
  pdpmlm_add( obj, grp, cls );
  for( grp = 1; grp < obj->ngr; grp++ ) { 
    
    set = 0;
    for( i = 0; i < cls; i++ ) {
      pdpmlm_add( obj, grp, cls );
      pdpmlm_parm( obj, cls, obj->s, obj->m, &obj->a, &obj->b );
      if( ( obj->b / obj->yycl[ cls ] ) < 0.98 ) { set = 1; break; }
      else { pdpmlm_sub( obj, grp, cls ); }
    }
    if( set == 0 ) { pdpmlm_add( obj, grp, ++cls ); }
  }
  if( obj->flags & FLAG_VERBOSE ) {
    pdpmlm_printf("iter: 0, ncl: %u, logp: %f\n", obj->ncl, pdpmlm_logp( obj ) );
  }
}

void pdpmlm_add( pdpmlm_t * obj, unsigned int grp, unsigned int cls ) {
  unsigned int i, j;
 
  if( grp >= obj->ngr ) { pdpmlm_printf("grp: %u\n", grp); error( "pdpmlm_add: invalid grp argument" ); }
  if( cls >= obj->ngr ) { pdpmlm_printf("cls: %u\n", cls); error( "pdpmlm_add: invalid cls argument" ); }
 
  // 1. set vcl, recompute pcl, and possibly ncl
  obj->vcl[ grp ] = cls;
  if( obj->pcl[ cls ] == 0 ) { obj->ncl++; }
  obj->pcl[ cls ] += 1;

  // 2. allocate memory for xxcl and xycl if necessary, zero xxcl, xycl
  if( obj->xxcl[ cls ] == NULL ) {
    obj->xxcl[ cls ] = (double *) pdpmlm_alloc( obj->q * obj->q, sizeof(double) );
    if( obj->xxcl[ cls ] == NULL ) { memerror(); }
    else { obj->mem += obj->q * obj->q * sizeof(double); }
    obj->xycl[ cls ] = (double *) pdpmlm_alloc( obj->q, sizeof(double) );
    if( obj->xycl[ cls ] == NULL ) { memerror(); }
    else { obj->mem += obj->q * sizeof(double); }
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
  obj->pcl[ cls ] -= 1;
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
  if( old != BAD_CLS ) {
    pdpmlm_sub( obj, grp, old );
  }
  pdpmlm_add( obj, grp, cls );
}

double pdpmlm_movep( pdpmlm_t * obj, unsigned int grp, unsigned int cls ) {
  double logp = 0.0;
  unsigned int oldcls = obj->vcl[ grp ];
  unsigned int oldncl = obj->ncl;
  if( oldcls == cls ) { return logp; }
  if( oldcls != BAD_CLS ) {
    logp -= pdpmlm_logpcls( obj, oldcls );
    pdpmlm_sub( obj, grp, oldcls );
    logp += pdpmlm_logpcls( obj, oldcls );
  }
  logp -= pdpmlm_logpcls( obj, cls );
  pdpmlm_add( obj, grp, cls );
  logp += pdpmlm_logpcls( obj, cls );
  if( obj->ncl > oldncl ) { logp += log( obj->alp ) - log( obj->ncl ); }
  else if( oldncl > obj->ncl ) { logp -= log( obj->alp ) - log( oldncl ); }
  return logp;
}

double pdpmlm_logpcls( pdpmlm_t * obj, unsigned int cls ) {
  double logp;
  if( obj->pcl[ cls ] == 0 ) { return 0.0; }
  pdpmlm_parm( obj, cls, obj->s, obj->m, &obj->a, &obj->b );
  logp = lgamma( obj->a / 2 ) - ( obj->a / 2 ) * log( obj->b / 2 );
  // "cluster" form of prior - allows more similar sized groups
  if( obj->flags & FLAG_PRICLUS ) { logp += log( obj->pcl[ cls ] ); }
  // "Dirichlet" form of prior - promotes fewer groups of different size
  else { logp += lfactorial( obj->pcl[ cls ] - 1 ); }
  return logp;
}

double pdpmlm_logp( pdpmlm_t * obj ) {
  unsigned int i, cls = 0;
  double logp;
  logp = obj->ncl * log( obj->alp ) - lfactorial( obj->ncl );
  for( i = 0; i < obj->ncl; i++ ) {
    while( obj->pcl[ cls ] == 0 ) { cls++; }
    logp += pdpmlm_logpcls( obj, cls );
    cls++;
  }
  return logp;
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

  // 4. b = y'y + s0*m0'm0 - m'sm
  *b = obj->yycl[ cls ]; // b = y'y
  *b += obj->s0*F77_CALL(ddot)(&obj->q, obj->m0, &ione, obj->m0, &ione);  // b += s0*m0'm0
  F77_CALL(dgemv)( "N", &obj->q, &obj->q, &done, s, &obj->q, m, &ione, &dzero, obj->buf, &ione );  // obj->buf = s*m
  *b -= F77_CALL(ddot)( &obj->q, m, &ione, obj->buf, &ione ); // b -= m'obj->buf
 

  // 5. a = a0 + nk;
  *a = obj->a0 + obj->pcl[ cls ];
}


unsigned int pdpmlm_free( pdpmlm_t * obj ) {
  unsigned int cls = 0;
  while( cls < obj->ngr && obj->pcl[ cls ] > 0 ) { cls++; }
  if( cls == obj->ngr ) { cls = BAD_CLS; }
  return cls;
}

void pdpmlm_away( pdpmlm_t * obj, unsigned int grp ) {
  unsigned int i, test_cls, old_cls, best_cls;
  double test_delp=0, best_delp=0;

  old_cls = obj->vcl[ grp ];

  if( grp >= obj->ngr ) { error( "pdpmlm_away: invalid argument" ); }

  if( obj->pcl[ best_cls ] > obj->pgr[ grp ] ) {
    best_cls = pdpmlm_free( obj );
    if( test_cls == BAD_CLS ) { error("pdpmlm_best: test_cls should not == BAD_CLS"); }
    best_delp = pdpmlm_movep( obj, grp, best_cls );
  }

  test_cls = 0;
  for( i = 0; i < obj->ncl; i++ ) {
    while( obj->pcl[ test_cls ] == 0 ) { test_cls++; }
    test_delp += pdpmlm_movep( obj, grp, test_cls );
    if( test_delp > best_delp && test_cls != old_cls) { 
      best_delp = test_delp;
      best_cls  = test_cls;
    }
    test_cls++;
  }

  if( obj->vcl[ grp ] != best_cls ) { pdpmlm_move( obj, grp, best_cls ); }
}

void pdpmlm_best( pdpmlm_t * obj, unsigned int grp ) {
  unsigned int i, test_cls, best_cls;
  double test_delp=0, best_delp=0;

  if( grp >= obj->ngr ) { error( "pdpmlm_best: invalid argument" ); }

  best_cls = obj->vcl[ grp ];

  if( obj->pcl[ best_cls ] > 1 ) {
    test_cls = pdpmlm_free( obj );
    if( test_cls == BAD_CLS ) { error("pdpmlm_best: test_cls should not == BAD_CLS"); }
    test_delp += pdpmlm_movep( obj, grp, test_cls );
    if( test_delp > best_delp ) { 
      best_delp = test_delp;
      best_cls  = test_cls;
    }
  }

  test_cls = 0;
  for( i = 0; i < obj->ncl; i++ ) {
    while( obj->pcl[ test_cls ] == 0 ) { test_cls++; }
    test_delp += pdpmlm_movep( obj, grp, test_cls );
    if( test_delp > best_delp ) { 
      best_delp = test_delp;
      best_cls  = test_cls;
    }
    test_cls++;
  }

  if( obj->vcl[ grp ] != best_cls ) { pdpmlm_move( obj, grp, best_cls ); }
}
    
void pdpmlm_chunk( pdpmlm_t * obj, unsigned int itermax, double crit) {
  unsigned int i, *vcl_old, *grps, ngrps, cls, iter = 0;
  double logp_old, logp, pdel = 1, pcum = 0;

  // 0. allocate memory for vcl_old, grps
  vcl_old = (unsigned int *) pdpmlm_alloc( obj->ngr, sizeof( unsigned int ) );
  if( vcl_old == NULL ) { memerror(); }
  grps    = (unsigned int *) pdpmlm_alloc( obj->ngr, sizeof( unsigned int ) );
  if( grps == NULL ) { memerror(); }

  // 1. select the number of groups to shuffle
  ngrps = obj->ngr / 3;
  ngrps = ngrps == 0 ? 1 : ngrps;
  if( obj->flags & FLAG_VERBOSE && (iter % 1) == 0 ) {
    pdpmlm_printf( "updating %u of %u groups per iteration\n", ngrps, obj->ngr );
  }

  GetRNGstate();
  while( iter++ < itermax ) {

    // 2. randomly select ngrps groups to shuffle, save indicators
    for( i = 0; i < ngrps; i++ ) {
      grps[ i ] = (unsigned int) floor( obj->ngr * runif( 0.0, 1.0 ) );
      vcl_old[ i ] = obj->vcl[ grps[ i ] ];
    }
  
    // 3. compute old logp, move groups to same/random cluster
    logp_old = pdpmlm_logp( obj ); 
    //cls = pdpmlm_free( obj );
    //if( cls == BAD_CLS ) { cls = obj->vcl[ grps[ 0 ] ]; }
    cls = (unsigned int) floor( obj->ngr * runif( 0.0, 1.0 ) );
    for( i = 0; i < ngrps; i++ ) { pdpmlm_move( obj, grps[ i ], cls ); }
         
    // 4. move each group to best cluster
    for( i = 0; i < ngrps; i++ ) { pdpmlm_best( obj, grps[ i ] ); }

    // 5. compute logp, keep new clustering if better, else revert to old
    logp = pdpmlm_logp( obj );
    if( logp <= logp_old ) {    
        for( i = 0; i < ngrps; i++ ) { pdpmlm_move( obj, grps[ i ], vcl_old[ i ] ); }
        pdel *= 0.9;
    }
    
    // 6. update the stopping criterion
    else{ 
      pdel = 0.5 * (logp - logp_old) + 0.5 * pdel;
      logp_old = logp;
    }
    pcum += pdel;

    // 7. print summary if requested every 20 iterations
    if( obj->flags & FLAG_VERBOSE && (iter % 1) == 0 ) {
      pdpmlm_printf("iter: %u, ncl: %u, logp: %f, crit: %f\n", iter, obj->ncl, logp_old, pdel / pcum );
    }

    // 8. check stopping criterion, break the optimization loop if met, print if requested
    if( pcum > 0 && ( pdel / pcum ) < crit ) {
      obj->flags |= FLAG_OPTCRIT;
      if( obj->flags & FLAG_VERBOSE ) { 
        pdpmlm_printf( "iter: %u, ncl: %u, logp: %f, crit: %f\nstopping criterion met\n", iter, obj->ncl, logp_old, pdel / pcum ); 
      }
      break;
    }

  }
  if( !(obj->flags & FLAG_OPTCRIT) ) { warning("optimization criterion not met"); }
  PutRNGstate();
}
