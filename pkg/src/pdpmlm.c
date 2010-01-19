#include "pdpmlm.h"

void * pdpmlm_alloc( pdpmlm_t * obj, unsigned int count, unsigned int size ) {
  obj->mem += count * size;
  return R_alloc( count, size );
} 

void pdpmlm_divy( pdpmlm_t * obj ) {
  unsigned int i, grp = 0, cls = 0, set = 0;
  double a, b;
  pdpmlm_add( obj, grp, cls );
  for( grp = 1; grp < obj->ngr; grp++ ) {
    set = 0;
    for( i = 0; i < cls; i++ ) {
      pdpmlm_parm( obj, cls, obj->s, obj->m, &a, &b );
      pdpmlm_add( obj, grp, cls );
      pdpmlm_parm( obj, cls, obj->s, obj->m, &obj->a, &obj->b );
      if( ( (obj->a / obj->b) / (a / b) ) < 1.10 ) { set = 1; break; }
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
   
  if( grp >= obj->ngr || cls >= obj->ngr ) { error( "pdpmlm_add: invalid argument" ); } 
 
  // 1. set vcl, recompute pcl, and possibly ncl
  obj->vcl[ grp ] = cls;
  if( obj->pcl[ cls ] == 0 ) { obj->ncl++; }
  obj->pcl[ cls ] += 1;

  // 2. allocate memory for xxcl and xycl if necessary, zero xxcl, xycl
  if( obj->xxcl[ cls ] == NULL ) {
    obj->xxcl[ cls ] = (double *) pdpmlm_alloc( obj, obj->q * obj->q, sizeof(double) );
    obj->xycl[ cls ] = (double *) pdpmlm_alloc( obj, obj->q, sizeof(double) );
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
  if( obj->flags & FLAG_DIRICHL ) { logp += lfactorial( obj->pcl[ cls ] - 1 ); }
  else { logp += lfactorial( obj->pcl[ cls ] ) - obj->gam * lfactorial( obj->pcl[ cls ] - 1 ); }
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
  int i, j, d, ione=1, info;
  double done=1.0, dzero=0.0;

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
  // when the call is finished. If this could be avoided, would save some time
  // obj->fbuf holds some temporary data.
  // FIXME use dposv or dppsv instead (for sym,pd matrices), may be faster, latter requires changing
  // all matrices to triangular packed storatge, rather than full storage
  // FIXME do not 'error' here (although I have never observed these errors)
  F77_CALL(dgesv)((int *) &obj->q, &ione, s, (int *) &obj->q,
                  (int *) obj->fbuf, m, (int *) &obj->q, &info);
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
  // b = y'y
  *b = obj->yycl[ cls ];   
  // b += s0*m0'm0
  *b += obj->s0*F77_CALL(ddot)( (int *) &obj->q, obj->m0, &ione, obj->m0, &ione);
  // obj->fbuf = s*m 
  F77_CALL(dgemv)( "N", (int *) &obj->q, (int *) &obj->q, &done, s, (int *) &obj->q,
                    m, &ione, &dzero, obj->fbuf, &ione );
  // b -= m'obj->fbuf
  *b -= F77_CALL(ddot)( (int *) &obj->q, m, &ione, obj->fbuf, &ione );
 

  // 5. a = a0 + nk;
  *a = obj->a0 + obj->pcl[ cls ];
}

unsigned int pdpmlm_free( pdpmlm_t * obj ) {
  unsigned int cls = 0;
  while( cls < obj->ngr && obj->pcl[ cls ] > 0 ) { cls++; }
  if( cls == obj->ngr ) { cls = BAD_CLS; }
  return cls;
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

void pdpmlm_stoch( pdpmlm_t * obj, int maxiter, double crit) {
  unsigned int i, *vcl_old, *grps, ngrps, cls, iter = 0, spmercls;
  double logp_old, logp, pdel = 1, pcum = 0;

  // 0. allocate memory for vcl_old, grps
  vcl_old = (unsigned int *) pdpmlm_alloc( obj, obj->ngr, sizeof( unsigned int ) );
  grps    = (unsigned int *) pdpmlm_alloc( obj, obj->ngr, sizeof( unsigned int ) );

  GetRNGstate();
  while( maxiter-- != 0 ) {
  
    // 1. select the number of groups to shuffle
    ngrps = (unsigned int) floor( obj->ngr * runif( 0.0, 1.0 ) );
    ngrps = ngrps == 0 ? 1 : ngrps;
    
    // 2. randomly select ngrps groups to shuffle, save indicators
    for( i = 0; i < ngrps; i++ ) {
      grps[ i ] = (unsigned int) floor( obj->ngr * runif( 0.0, 1.0 ) );
      vcl_old[ i ] = obj->vcl[ grps[ i ] ];
    }
  
    // 3. compute old logp, move groups to random cluster
    logp_old = pdpmlm_logp( obj ); 
    cls = (unsigned int) floor( obj->ngr * runif( 0.0, 1.0 ) );
    for( i = 0; i < ngrps; i++ ) { pdpmlm_move( obj, grps[ i ], cls ); }

    logp = pdpmlm_logp( obj );
    if( logp <= logp_old ) {  

      // 4. move each group to best cluster
      for( i = 0; i < ngrps; i++ ) { pdpmlm_best( obj, grps[ i ] ); }

      // 5. compute logp, keep new clustering if better, else revert to old
      logp = pdpmlm_logp( obj );
      if( logp <= logp_old ) {    
        for( i = 0; i < ngrps; i++ ) { pdpmlm_move( obj, grps[ i ], vcl_old[ i ] ); }
        pdel *= 0.9;
      }

    }
   
    // 6. update the stopping criterion
    else{ 
      pdel = 0.5 * (logp - logp_old) + 0.5 * pdel;
      logp_old = logp;
    }
    pcum += pdel;

    // 7. print summary if requested every 20 iterations
    if( obj->flags & FLAG_VERBOSE && (iter % 1) == 0 ) {
      pdpmlm_printf("iter: %u, ncl: %u, logp: %f, ngrps: %u, crit: %f\n", iter, obj->ncl, logp_old, ngrps, pdel / pcum );
    }

    // 8. check stopping criterion, break the optimization loop if met, print if requested
    if( pcum > 0 && ( pdel / pcum ) < crit ) {
      obj->flags |= FLAG_OPTCRIT;
      if( obj->flags & FLAG_VERBOSE ) { 
        pdpmlm_printf( "iter: %u, ncl: %u, logp: %f, ngrps: %u, crit: %f\nstopping criterion met\n", iter, obj->ncl, logp_old, ngrps, pdel / pcum ); 
      }
      break;
    }

  }
  if( !(obj->flags & FLAG_OPTCRIT) ) { warning("optimization criterion not met"); }
  PutRNGstate();
}

void pdpmlm_merge( pdpmlm_t * obj, unsigned int cls1, unsigned int cls2 ) {
  unsigned int i, grp = 0, size;

  //0. cannot merge an empty group
  if( obj->pcl[ cls1 ] == 0 || obj->pcl[ cls2 ] == 0 ) { return; }
  size = obj->pcl[ cls1 ];

  //1. merge the cluster
  for( i = 0; i < size; i++ ) {
    while( obj->vcl[ grp ] != cls1 ) { grp++; }
    pdpmlm_move( obj, grp, cls2 );
  }
}    


double pdpmlm_mergep( pdpmlm_t * obj, unsigned int cls1, unsigned int cls2 ) {
  unsigned int i, grp = 0, size;
  double del;

  //0. cannot merge an empty group
  if( obj->pcl[ cls1 ] == 0 || obj->pcl[ cls2 ] == 0 ) { return 0.0; }
  size = obj->pcl[ cls1 ];

  //1. account for loss of cls1 in logp 
  //del = log( obj->ncl ) - log( obj->alp );
  //above would be mathematically correct, but not numerically
  //since lfactorial is an approximation 
  del = lfactorial( obj->ncl ) - lfactorial( obj->ncl - 1 ) - log( obj->alp );
  
  //2. merge the cluster
  del -= pdpmlm_logpcls( obj, cls1 );
  del -= pdpmlm_logpcls( obj, cls2 );
  for( i = 0; i < size; i++ ) {
    while( obj->vcl[ grp ] != cls1 ) { grp++; }
    pdpmlm_move( obj, grp, cls2 );
  }
  del += pdpmlm_logpcls( obj, cls1 );
  del += pdpmlm_logpcls( obj, cls2 );

  return del;
}    

double pdpmlm_testmergep( pdpmlm_t * obj, unsigned int cls1, unsigned int cls2 ) {
  unsigned int grp = 0, testgrp, size;
  double del;

  //0. cannot merge an empty group
  if( obj->pcl[ cls1 ] == 0 || obj->pcl[ cls2] == 0 ) { return 0.0; }
  size = obj->pcl[ cls1 ];

  //1. account for loss of cls1 in logp 
  //del = log( obj->ncl ) - log( obj->alp );
  //above would be mathematically correct, but not numerically
  //since lfactorial is an approximation 
  del = lfactorial( obj->ncl ) - lfactorial( obj->ncl - 1 ) - log( obj->alp );
  
  //2. enumerate groups in cls1
  for( testgrp = 0; testgrp < obj->pcl[ cls1 ]; testgrp++ ) {
    while( obj->vcl[ grp ] != cls1 ) { grp++; }
    obj->pbuf[ testgrp ] = grp++;
  }
  
  //3. merge cls1 to cls2
  del -= pdpmlm_logpcls( obj, cls1 );
  del -= pdpmlm_logpcls( obj, cls2 );
  for( testgrp = 0; testgrp < size; testgrp++ ) {
    pdpmlm_move( obj, obj->pbuf[ testgrp ], cls2 );
  }
  del += pdpmlm_logpcls( obj, cls1 );
  del += pdpmlm_logpcls( obj, cls2 );

  //4. unmerge
  for( testgrp = 0; testgrp < size; testgrp++ ) {
    pdpmlm_move( obj, obj->pbuf[ testgrp ], cls1 );
  }
  return del;
}

void pdpmlm_agglo( pdpmlm_t * obj, int maxiter ) {
  //delp[ i, j ] - change in logp by merging cluster j with i, stored as a lower
  //triangular matrix with ( ngr * ( ngr - 1 ) / 2 ) items. The value is accessed
  //via array according to the following.
  //delp[ i, j ] = del[ i * ( ngr - 1 ) - i * ( i - 1 ) / 2 - i + j - 1 ] for 0 <= i < j < ngr 
  double * delp;
  double   delp_temp, delp_best, logp_best = -DBL_MAX;
  unsigned int * vcl_best;
  unsigned int   delp_ind, i, j, icls, jcls, icls_best = BAD_CLS, jcls_best = BAD_CLS;
  unsigned int   icls_last = BAD_CLS, jcls_last = BAD_CLS;
  int calcs = 0, calcs_cent = obj->ngr * ( obj->ngr - 1 ) + 1;
  calcs_cent = calcs_cent > 100 ? calcs_cent / 100 : 1;

  //0. allocate some additional memory
  delp = (double *) pdpmlm_alloc( obj, ( obj->ngr * ( obj->ngr - 1 ) / 2 ), sizeof( double ) );
  vcl_best = (unsigned int *) pdpmlm_alloc( obj, obj->ngr, sizeof( unsigned int ) );

  //1. compute initial logp 
  obj->logp = pdpmlm_logp( obj );  

  // repeat until all clusters are merged into one
  while( obj->ncl > 1 && maxiter-- != 0 ) {

    //2. update and record best delp among all pairs of clusters
    //testmerge all pairs at first loop (i.e. obj->ncl == obj->ngr)
    //for subsequent loops, testmerge only the pairs involving jcls_best
    //other pairs need only be updated with the following:
    //lfactorial( obj->ncl ) - lfactorial( obj->ncl - 1 ) - log( obj->alp );
    delp_temp = 2 * lfactorial( obj->ncl ) - lfactorial( obj->ncl + 1 ) - lfactorial( obj->ncl - 1 );
    delp_best = -DBL_MAX;
    icls = 0;
    for( i = 0; i < obj->ncl - 1; i++ ) {
      while( obj->pcl[ icls ] == 0 ) { icls++; }
      jcls = icls + 1;
      for( j = 0; j < obj->ncl - i - 1; j++ ) {
        while( obj->pcl[ jcls ] == 0 ) { jcls++; }
        delp_ind = icls * ( obj->ngr - 1 ) - icls * ( icls - 1 ) / 2 - icls + jcls - 1;
        if( obj->ncl == obj->ngr || icls == jcls_last || jcls == jcls_last ) { 
          delp[ delp_ind ] = pdpmlm_testmergep( obj, icls, jcls ); 
          calcs++;
          if( (obj->flags & FLAG_VERBOSE) && (calcs % calcs_cent == 0) ) {
            pdpmlm_printf("\rpercent complete: %d%", calcs / calcs_cent); 
          }
        } else { delp[ delp_ind ] += delp_temp; }
        if( delp[ delp_ind ] > delp_best ) {
          delp_best = delp[ delp_ind ];
          icls_best = icls;
          jcls_best = jcls;
        }
        jcls++;
      }
      icls++;
    }
    
    //3. merge the clusters corresponding to the largest delp
    pdpmlm_merge( obj, icls_best, jcls_best );
    icls_last = icls_best;
    jcls_last = jcls_best;
    obj->logp += delp_best;

    //4. if highest logp (so far) is obtained by merge, save the vcl
    if( obj->logp > logp_best ) {
      logp_best = obj->logp;
      for( i = 0; i < obj->ngr; i++ ) { vcl_best[ i ] = obj->vcl[ i ]; }
    }  
  }

  if( obj->flags & FLAG_VERBOSE ) { pdpmlm_printf("\rpercent complete: 100%\n"); }
  //5. cluster the groups according to the highest obtained logp
  obj->logp = logp_best;
  for( i = 0; i < obj->ngr; i++ ) {
    pdpmlm_move( obj, i, vcl_best[ i ] );
  }
}

