#include "pdpmlm.h"

void * pdpmlm_alloc( pdpmlm_t * obj, unsigned int count, unsigned int size ) {
  obj->mem += count * size;
  return R_alloc( count, size );
} 

void * pdpmlm_zalloc( pdpmlm_t * obj, unsigned int count, unsigned int size ) {
  void * data = R_alloc( count, size );
  char * start = (char *) data;
  char * end = start + count * size;
  do { *(start++) = 0; } while( start < end );
  obj->mem += count * size;
  return data;
} 

void pdpmlm_divy( pdpmlm_t * obj ) {
  unsigned int i, grp = 0, cls = 0, set = 0;
  double a, b;
  pdpmlm_add( obj, grp, cls );
  for( grp = 1; grp < obj->ngr; grp++ ) {
    for( cls = 0; cls < obj->ncl; cls++ ) {
      pdpmlm_parm( obj, cls, obj->s, obj->m, &a, &b, &obj->d );
      pdpmlm_add( obj, grp, cls );
      pdpmlm_parm( obj, cls, obj->s, obj->m, &obj->a, &obj->b, &obj->d );
      //Generally precision (a/b) is decreased with addition 
      //of a group to a cluster. We accept the addition when
      //the the precision changes by more than 0.95 fold
      if( ( (obj->a / obj->b) / (a / b) ) > 0.95 ) { break; }
      else { pdpmlm_sub( obj, grp, cls ); }
    }
    if( cls == obj->ncl ) { pdpmlm_add( obj, grp, cls ); } 
  }
  obj->logp = pdpmlm_logp( obj );
  if( obj->flags & FLAG_VERBOSE ) {
    pdpmlm_printf("initialized: logp: %f\n", obj->logp );
  } 
}

void pdpmlm_add( pdpmlm_t * obj, unsigned int grp, unsigned int cls ) {
  unsigned int i, j, index;
   
  if( grp >= obj->ngr || cls >= obj->ngr ) { 
    error( "pdpmlm_add: invalid argument: grp = %u cls = %u", grp, cls ); 
  } 
 
  //set vcl, recompute gcl, pcl, and possibly ncl
  obj->vcl[ grp ] = cls;
  if( obj->gcl[ cls ] == 0 ) { obj->ncl++; }
  obj->pcl[ cls ] += obj->pgr[ grp ];
  obj->gcl[ cls ] += 1;

  //allocate and zero memory for xxcl and xycl if necessary
  if( obj->xxcl[ cls ] == NULL ) {
    obj->xxcl[ cls ] = (double *) pdpmlm_zalloc( obj, ( obj->q * ( obj->q + 1 ) ) / 2, sizeof(double) );
    obj->xycl[ cls ] = (double *) pdpmlm_zalloc( obj, obj->q, sizeof(double) );
    obj->yycl[ cls ] = 0.0;
  }

  //(re)compute xxcl, xycl yycl for cls using xxgr, xygr, and yygr from grp
  //use index here rather then UMAT(i, j) to avoid unnecessary index computation
  index = 0;
  obj->yycl[ cls ] += obj->yygr[ grp ];
  for( i = 0; i < obj->q; i++ ) {
    obj->xycl[ cls ][ i ] += obj->xygr[ grp ][ i ];
    for( j = i; j < obj->q; j++ ) {
      obj->xxcl[ cls ][ index ] += obj->xxgr[ grp ][ index ];
      index++;
    }
  }
}

void pdpmlm_sub( pdpmlm_t * obj, unsigned grp, unsigned int cls ) {
  unsigned int i, j, index;

  if( grp >= obj->ngr || cls >= obj->ngr ) { 
    error( "pdpmlm_sub: invalid argument: grp = %u", grp ); 
  } 

  //set vcl, recompute gcl, and possibly ncl
  obj->vcl[ grp ] = BAD_VCL;
  obj->pcl[ cls ] -= obj->pgr[ grp ];
  obj->gcl[ cls ] -= 1;
  if( obj->gcl[ cls ] == 0 ) { obj->ncl--; }

  //recompute xxcl, xycl yycl for cls using xxgr, xygr, and yygr from grp
  //use index here rather then UMAT(i, j) to avoid unnecessary index computation
  index = 0;
  obj->yycl[ cls ] -= obj->yygr[ grp ];
  for( i = 0; i < obj->q; i++ ) {
    obj->xycl[ cls ][ i ] -= obj->xygr[ grp ][ i ];
    for( j = i; j < obj->q; j++ ) {
      obj->xxcl[ cls ][ index ] -= obj->xxgr[ grp ][ index ];
      index++;
    }
  }
}

void pdpmlm_move( pdpmlm_t * obj, unsigned int grp, unsigned int cls ) {
  unsigned int old = obj->vcl[ grp ];
  if( old == cls ) { return; }
  if( old != BAD_VCL ) {
    pdpmlm_sub( obj, grp, old );
  }
  pdpmlm_add( obj, grp, cls );
}

double pdpmlm_movep( pdpmlm_t * obj, unsigned int grp, unsigned int cls ) {
  double logp = 0.0;
  unsigned int only[2];
  if( obj->vcl[ grp ] == cls ) { return logp; }
  only[0] = obj->vcl[ grp ];
  only[1] = cls;
  logp -= pdpmlm_logponly( obj, only, 2 );
  pdpmlm_move( obj, grp, cls );
  logp += pdpmlm_logponly( obj, only, 2 );
  return logp;
}

double pdpmlm_logpcls( pdpmlm_t * obj, unsigned int cls ) {
  double logp;
  if( obj->gcl[ cls ] == 0 ) { return 0.0; }
  pdpmlm_parm( obj, cls, obj->s, obj->m, &obj->a, &obj->b, &obj->d );
  logp = lgamma( obj->a / 2 ) - ( obj->a / 2 ) * log( obj->b / 2 ) - obj->d;
  if( obj->flags & FLAG_DIRICHL ) { logp += lgamma( obj->gcl[ cls ] ); }
  else { logp += lgamma( obj->gcl[ cls ] + 1 ) - obj->lam * lgamma( obj->gcl[ cls ] ); }
  return logp;
}


double pdpmlm_logp( pdpmlm_t * obj ) {
  unsigned int i, cls = 0;
  double logp;
  for( i = 0; i < obj->ncl; i++ ) {
    while( obj->gcl[ cls ] == 0 ) { cls++; }
    logp += pdpmlm_logpcls( obj, cls );
    cls++;
  }
  return logp;
}

double pdpmlm_logponly( pdpmlm_t * obj, unsigned int * only, unsigned int size ) {
  unsigned int i;
  double logp;
  logp = obj->ncl * log( obj->alp ) - obj->lam * lgamma( obj->ncl + 1 );
  for( i = 0; i < size; i++ ) {
    logp += pdpmlm_logpcls( obj, only[ i ] );
  }
}

void pdpmlm_parm( pdpmlm_t * obj, unsigned int cls, double * s, double * m, double * a, double * b, double * d ) {
  int i, j, dj, ione=1, info, *ipiv;
  double done=1.0, dzero=0.0;

  if( obj->gcl[ cls ] == 0 ) { return; }
  if( cls >= obj->ngr ) { error( "pdpmlm_parm: invalid argument" ); }

  //1. load s with s0I + x'x, load m with s0m0 + x'y
  dj = 0;
  for( i = 0; i < obj->q; i++ ) {
    m[ i ] = obj->s0 * obj->m0[ i ] + obj->xycl[ cls ][ i ];
    for( j = i; j < obj->q; j++ ) {
      s[ dj ] = obj->xxcl[ cls ][ dj ];
      if( j == i ) { s[ dj ] += obj->s0; }
      dj++;
    }
  }
     
  //2. m = (s)^(-1) * m
  //dspsv overwrites s, must reload s aftward.
  ipiv = (int *) obj->fbuf;
  F77_CALL(dspsv)("U", (int *) &obj->q, &ione, s, ipiv, m, (int *) &obj->q, &info);
  if( info > 0 ) { warning("dppsv: system is singular"); }
  if( info < 0 ) { error("dppsv: invalid argument"); }

  //2.5 d = 0.5*log|det(s)| (see dspsv/dsptrf documentation)
  *d = 0.0;
  for( i = 0; i < obj->q; i++ ) {
    if( ipiv[ i ] > 0 ) {
      *d += log( ABS( s[ UMAT(i, i) ] ) );
    } else if( i > 0 && ipiv[ i ] < 0 && ipiv[ i-1 ] == ipiv[ i ] ) {
      *d += log( ABS(
        s[ UMAT(i-1,i-1) ] * s[ UMAT(i,i) ] -\
        s[ UMAT(i-1,i) ] * s[ UMAT(i-1,i) ] )\
      );
    }
  }
  *d *= 0.5;

  //3. reload s
  dj = 0;
  for( i = 0; i < obj->q; i++ ) {
    for( j = i; j < obj->q; j++ ) {
      s[ dj ] = obj->xxcl[ cls ][ dj ];
      if( j == i ) { s[ dj ] += obj->s0; }
      dj++;
    }
  }

  //4. b = b0 + y'y + s0*m0'm0 - m'sm
  //b = b0 + y'y
  *b = obj->b0 + obj->yycl[ cls ];   
  //b += s0*m0'm0
  *b += obj->s0*F77_CALL(ddot)( (int *) &obj->q, obj->m0, &ione, obj->m0, &ione);
  //obj->fbuf = s*m
  if( info == 0 ) { 
    F77_CALL(dspmv)("U", (int *) &obj->q, &done, s, m, &ione, &dzero, obj->fbuf, &ione);
    //b -= m'obj->fbuf
    *b -= F77_CALL(ddot)( (int *) &obj->q, m, &ione, obj->fbuf, &ione );
  }
 
  //5. a = a0;
  *a = obj->a0 + obj->pcl[ cls ];
}

unsigned int pdpmlm_free( pdpmlm_t * obj ) {
  unsigned int cls = 0;
  while( cls < obj->ngr && obj->gcl[ cls ] > 0 ) { cls++; }
  if( cls == obj->ngr ) { cls = BAD_VCL; }
  return cls;
}

void pdpmlm_best( pdpmlm_t * obj, unsigned int grp ) {
  unsigned int i, test_cls, best_cls;
  double test_delp=0, best_delp=0;

  if( grp >= obj->ngr ) { error( "pdpmlm_best: invalid argument" ); }

  best_cls = obj->vcl[ grp ];

  if( obj->gcl[ best_cls ] > 1 ) {
    test_cls = pdpmlm_free( obj );
    if( test_cls == BAD_VCL ) { error("pdpmlm_best: test_cls should not == BAD_VCL"); }
    test_delp += pdpmlm_movep( obj, grp, test_cls );
    if( test_delp > best_delp ) { 
      best_delp = test_delp;
      best_cls  = test_cls;
    }
  }

  test_cls = 0;
  for( i = 0; i < obj->ncl; i++ ) {
    while( obj->gcl[ test_cls ] == 0 ) { test_cls++; }
    test_delp += pdpmlm_movep( obj, grp, test_cls );
    if( test_delp > best_delp ) { 
      best_delp = test_delp;
      best_cls  = test_cls;
    }
    test_cls++;
  }

  if( obj->vcl[ grp ] != best_cls ) { pdpmlm_move( obj, grp, best_cls ); }
}

static unsigned int rmulti( double * logp, unsigned int n ) {
  unsigned int i, j, ret = 0;
  double pl, ph, u;
  u = pdpmlm_runif(0.0, 0.1);
  pl = 0; 
  for( i = 0; i < n; i++ ) {
    ph = 0;
    for( j = 0; j < n; j++ ) {
      ph += exp( logp[ j ] - logp[ i ] );
    }
    ph = pl + 1/ph;
    if( u > pl && u < ph ) { ret = i; }
    pl = ph;
  }
  return ret;
}

void pdpmlm_gibbs( pdpmlm_t * obj, int maxiter, double crit) {
  unsigned int i, *vcl_best, cls, grp, iter = 0;
  unsigned int cls_old, cls_new, *proposal_cls, proposal_ncl;
  double stopcrit = DBL_MAX, logp_best, *proposal_logp;

  //0. compute initial logp, save initial partition
  obj->logp = pdpmlm_logp( obj );
  logp_best = obj->logp;
  vcl_best = (unsigned int *) pdpmlm_alloc( obj, obj->ngr, sizeof( unsigned int ) );
  for( i = 0; i < obj->ngr; i++ ) { vcl_best[ i ] = obj->vcl[ i ]; }

  proposal_logp = (double*) pdpmlm_alloc( obj, obj->ngr, sizeof( unsigned int ) );
  proposal_cls  = (unsigned int*) pdpmlm_alloc( obj, obj->ngr, sizeof( unsigned int ) );

  while( iter++ < maxiter && stopcrit > crit ) {
  
    for( grp = 0; grp < obj->ngr; grp++ ) {

      //1. save old cluster
      cls_old = obj->vcl[ grp ];

      //2. compute change in logp for possible moves
      cls = 0;
      proposal_ncl = 0;
      for( i = 0; i < obj->ncl; i++ ) {
        while( obj->gcl[ cls ] == 0 ) { cls++; }
        proposal_logp[ i ] = pdpmlm_movep( obj, grp, cls );
        proposal_cls[ i ] = cls;
        proposal_ncl++;
        pdpmlm_move( obj, grp, cls_old );
        cls++;
      }
      if( obj->gcl[ obj->vcl[ grp ] ] > 1 ) { 
        proposal_cls[ proposal_ncl ] = pdpmlm_free( obj );
        proposal_logp[ proposal_ncl ] = pdpmlm_movep( obj, grp, proposal_cls[ proposal_ncl ] );
        proposal_ncl++;
        pdpmlm_move( obj, grp, cls_old );
      }

      //3. draw from the conditional
      cls_new = proposal_cls[ rmulti( proposal_logp, proposal_ncl ) ];
      obj->logp += pdpmlm_movep( obj, grp, cls_new );

    }
   
    //4. save best partition so far
    if( obj->logp > logp_best ) {
      if( stopcrit == DBL_MAX ) { stopcrit = obj->logp - logp_best; }
      else { stopcrit += obj->logp - logp_best; }
      logp_best = obj->logp;
      for( i = 0; i < obj->ngr; i++ ) { vcl_best[ i ] = obj->vcl[ i ]; }
    }
     
    //5. print summary if requested every iteration
    if( obj->flags & FLAG_VERBOSE && (iter % 1) == 0 ) {
      pdpmlm_printf("iter: %u, ncl: %u, logp: %f, logp_best: %f\n", iter, obj->ncl, obj->logp, logp_best );
    }
  
    //6. update stopping criterion
    if(stopcrit != DBL_MAX) { stopcrit *= 0.95; }
  }
   
  //7. restore the best partition
  for( grp = 0; grp < obj->ngr; grp++ ) { 
    pdpmlm_move( obj, grp, vcl_best[ grp ] );  
  }
  obj->logp = pdpmlm_logp( obj );

  if( !(obj->flags & FLAG_OPTCRIT) ) { warning("optimization criterion not met"); }
}


void pdpmlm_stoch( pdpmlm_t * obj, int maxiter, double crit) {
  unsigned int i, *vcl_old, *grps, ngrps, cls, iter = 0, spmercls;
  double logp_old, pdel = 1, pcum = 0;

  //0. allocate memory for vcl_old, grps
  vcl_old = (unsigned int *) pdpmlm_alloc( obj, obj->ngr, sizeof( unsigned int ) );
  grps    = (unsigned int *) pdpmlm_alloc( obj, obj->ngr, sizeof( unsigned int ) );

  //0.5 compute initial logp
  obj->logp = pdpmlm_logp( obj );

  while( iter++ < maxiter ) {
  
    //1. select the number of groups to shuffle
    ngrps = (unsigned int) floor( obj->ngr * pdpmlm_runif( 0.0, 1.0 ) );
    ngrps = ngrps == 0 ? 1 : ngrps;
    
    //2. randomly select ngrps groups to shuffle, save indicators
    for( i = 0; i < ngrps; i++ ) {
      grps[ i ] = (unsigned int) floor( obj->ngr * pdpmlm_runif( 0.0, 1.0 ) );
      vcl_old[ i ] = obj->vcl[ grps[ i ] ];
    }
  
    //3. compute old logp, move groups to random cluster 
    logp_old = obj->logp;
    cls = (unsigned int) floor( obj->ngr * pdpmlm_runif( 0.0, 1.0 ) );
    for( i = 0; i < ngrps; i++ ) { pdpmlm_move( obj, grps[ i ], cls ); }

    obj->logp = pdpmlm_logp( obj );
    if( obj->logp <= logp_old ) {  

      //4. move each group to best cluster
      for( i = 0; i < ngrps; i++ ) { pdpmlm_best( obj, grps[ i ] ); }

      //5. compute logp, keep new clustering if better, else revert to old
      obj->logp = pdpmlm_logp( obj );
      if( obj->logp <= logp_old ) {    
        for( i = 0; i < ngrps; i++ ) { pdpmlm_move( obj, grps[ i ], vcl_old[ i ] ); }
        pdel *= 0.9;
        obj->logp = logp_old;
      }

    }
   
    //6. update the stopping criterion
    else{ 
      pdel = 0.5 * (obj->logp - logp_old) + 0.5 * pdel;
      logp_old = obj->logp;
    }
    pcum += pdel;

    //7. print summary if requested every 20 iterations
    if( obj->flags & FLAG_VERBOSE && (iter % 1) == 0 ) {
      pdpmlm_printf("iter: %u, ncl: %u, logp: %f, ngrps: %u, crit: %f\n",\
                    iter, obj->ncl, logp_old, ngrps, pdel / pcum );
    }

    //8. check stopping criterion, break the optimization loop if met, print if requested
    if( pcum > 0 && ( pdel / pcum ) < crit ) {
      obj->flags |= FLAG_OPTCRIT;
      if( obj->flags & FLAG_VERBOSE ) { 
        pdpmlm_printf( "iter: %u, ncl: %u, logp: %f, ngrps: %u, crit: %f\n\
                        stopping criterion met\n",\
                        iter, obj->ncl, logp_old, ngrps, pdel / pcum ); 
      }
      break;
    }

  }
  if( !(obj->flags & FLAG_OPTCRIT) ) { warning("optimization criterion not met"); }
}

void pdpmlm_merge( pdpmlm_t * obj, unsigned int cls1, unsigned int cls2 ) {
  unsigned int i, grp = 0, size;

  //0. cannot merge an empty group
  if( obj->gcl[ cls1 ] == 0 || obj->gcl[ cls2 ] == 0 ) { return; }
  size = obj->gcl[ cls1 ];

  //1. merge the cluster
  for( i = 0; i < size; i++ ) {
    while( obj->vcl[ grp ] != cls1 ) { grp++; }
    pdpmlm_move( obj, grp, cls2 );
  }
}    


double pdpmlm_mergep( pdpmlm_t * obj, unsigned int cls1, unsigned int cls2 ) {
  unsigned int i, grp = 0, size;
  unsigned int only[2];
  double del = 0;
  
  //cannot merge an empty group
  if( obj->gcl[ cls1 ] == 0 || obj->gcl[ cls2 ] == 0 ) { return 0.0; }
  size = obj->gcl[ cls1 ];

  //merge the cluster
  only[0] = cls1;
  only[1] = cls2;
  del -= pdpmlm_logponly( obj, only, 2 );
  for( i = 0; i < size; i++ ) {
    while( obj->vcl[ grp ] != cls1 ) { grp++; }
    pdpmlm_move( obj, grp, cls2 );
  }
  del += pdpmlm_logponly( obj, only, 2 );
  return del;
}    

double pdpmlm_testmergep( pdpmlm_t * obj, unsigned int cls1, unsigned int cls2 ) {
  unsigned int grp = 0, testgrp, size;
  double del = 0;

  //cannot merge an empty group
  if( obj->gcl[ cls1 ] == 0 || obj->gcl[ cls2] == 0 ) { return 0.0; }
  size = obj->gcl[ cls1 ];
  
  //enumerate groups in cls1 (cannot use pbuf elsewhere!!!)
  for( testgrp = 0; testgrp < size; testgrp++ ) {
    while( obj->vcl[ grp ] != cls1 ) { grp++; }
    obj->pbuf[ testgrp ] = grp++;
  }
 
  //merge clusters 
  del = pdpmlm_mergep( obj, cls1, cls2 );

  //unmerge clusters
  for( testgrp = 0; testgrp < size; testgrp++ ) {
    pdpmlm_move( obj, obj->pbuf[ testgrp ], cls1 );
  }
  return del;
}

void pdpmlm_agglo( pdpmlm_t * obj, int maxiter ) {
  //delp[ i, j ] - change in logp by merging cluster j with i, stored as an upper
  //triangular packed storage matrix with ( ngr * ( ngr - 1 ) / 2 ) items. The
  //value is accessed with UMAT(i, j)
  double * delp;
  double   delp_temp, delp_best, logp_best = -DBL_MAX;
  unsigned int * vcl_best;
  unsigned int   delp_ind, i, j, icls, jcls, icls_best = BAD_VCL, jcls_best = BAD_VCL;
  unsigned int   icls_last = BAD_VCL, jcls_last = BAD_VCL;
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
    delp_temp = obj->lam * ( 2 * lgamma( obj->ncl + 1 ) - lgamma( obj->ncl + 2 ) - lgamma( obj->ncl ) );
    delp_best = -DBL_MAX;
    icls = 0;
    for( i = 0; i < obj->ncl - 1; i++ ) {
      while( obj->gcl[ icls ] == 0 ) { icls++; }
      jcls = icls + 1;
      for( j = 0; j < obj->ncl - i - 1; j++ ) {
        while( obj->gcl[ jcls ] == 0 ) { jcls++; }
        delp_ind = UMAT(icls, jcls);
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
