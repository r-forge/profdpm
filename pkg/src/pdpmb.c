#include "pdpmb.h"

void * pdpmb_alloc( pdpmb_t * obj, unsigned int count, unsigned int size ) {
  obj->mem += count * size;
  return R_alloc( count, size );
} 

void * pdpmb_zalloc( pdpmb_t * obj, unsigned int count, unsigned int size ) {
  void * data = R_alloc( count, size );
  char * start = (char *) data;
  char * end = start + count * size;
  do { *(start++) = 0; } while( start < end );
  obj->mem += count * size;
  return data;
} 

void pdpmb_divy( pdpmb_t * obj ) {
  unsigned int i, grp = 0;
  for( grp = 0; grp < obj->ngr; grp++ ) {
    pdpmb_add( obj, grp, grp );
  }
  obj->logp = pdpmb_logp( obj );
  if( obj->flags & FLAG_VERBOSE ) {
    pdpmb_printf("initialized: logp: %f\n", obj->logp );
  } 
}

void pdpmb_add( pdpmb_t * obj, unsigned int grp, unsigned int cls ) {
  unsigned int i, j, index;
  if( grp >= obj->ngr || cls >= obj->ngr ) { 
    error( "pdpmb_add: invalid argument: grp = %u cls = %u", grp, cls ); 
  } 
  //set vcl, recompute gcl, and possibly ncl
  obj->vcl[ grp ] = cls;
  if( obj->gcl[ cls ] == 0 ) { obj->ncl++; }
  obj->gcl[ cls ] += 1;
  //(re)compute gqcl
  for( i = 0; i < obj->q; i++ ) {
    obj->gqcl[ FMAT(cls, i, obj->ngr) ] += obj->y[ FMAT(grp, i, obj->ngr) ];
  }
}

void pdpmb_sub( pdpmb_t * obj, unsigned grp, unsigned int cls ) {
  unsigned int i, j, index;
  if( grp >= obj->ngr || cls >= obj->ngr ) { 
    error( "pdpmb_sub: invalid argument: grp = %u", grp ); 
  } 
  //set vcl, recompute gcl, and possibly ncl
  obj->vcl[ grp ] = BAD_VCL;
  obj->gcl[ cls ] -= 1;
  if( obj->gcl[ cls ] == 0 ) { obj->ncl--; }
  //recompute gqcl
  for( i = 0; i < obj->q; i++ ) {
    obj->gqcl[ FMAT(cls, i, obj->ngr) ] -= obj->y[ FMAT(grp, i, obj->ngr) ];
  }
}

void pdpmb_move( pdpmb_t * obj, unsigned int grp, unsigned int cls ) {
  unsigned int old = obj->vcl[ grp ];
  if( old == cls ) { return; }
  if( old != BAD_VCL ) {
    pdpmb_sub( obj, grp, old );
  }
  pdpmb_add( obj, grp, cls );
}

double pdpmb_movep( pdpmb_t * obj, unsigned int grp, unsigned int cls ) {
  double logp = 0.0;
  unsigned int only[2];
  unsigned int old = obj->vcl[ grp ];
  if( old == cls ) { return logp; }
  only[0] = old;
  only[1] = cls;
  logp -= pdpmb_logponly( obj, only, 2 );
  pdpmb_move( obj, grp, cls );
  logp += pdpmb_logponly( obj, only, 2 );
  return logp;
}

double pdpmb_logpcls( pdpmb_t * obj, unsigned int cls ) {
  unsigned int i;
  double logp = 0.0;
  if( obj->gcl[ cls ] == 0 ) { return logp; }
  //compute posterior mass
  for( i = 0; i < obj->q; i++ ) {
    logp += lgamma( obj->a0 + (double) obj->gqcl[ FMAT(cls, i, obj->ngr) ] ) +\
            lgamma( obj->b0 + (double) obj->gcl[ cls ] -\
            (double) obj->gqcl[ FMAT(cls, i, obj->ngr) ] ) -\
            lgamma( (double) obj->gcl[ cls ] + obj->a0 + obj->b0 );
  }
  if( obj->flags & FLAG_DIRICHL ) { logp += lgamma( obj->gcl[ cls ] ); }
  else { logp += lgamma( obj->gcl[ cls ] + 1 ) - obj->lam * lgamma( obj->gcl[ cls ] ); }
  return logp;
}


double pdpmb_logp( pdpmb_t * obj ) {
  unsigned int i, cls = 0;
  double logp;
  logp = obj->ncl * log( obj->alp ) - obj->lam * lgamma( obj->ncl + 1 );
  for( i = 0; i < obj->ncl; i++ ) {
    while( obj->gcl[ cls ] == 0 ) { cls++; }
    logp += pdpmb_logpcls( obj, cls );
    cls++;
  }
  return logp;
}

double pdpmb_logponly( pdpmb_t * obj, unsigned int * only, unsigned int size ) {
  unsigned int i;
  double logp;
  logp = obj->ncl * log( obj->alp ) - obj->lam * lgamma( obj->ncl + 1 );
  for( i = 0; i < size; i++ ) {
    logp += pdpmb_logpcls( obj, only[ i ] );
  }
  return logp;
}

unsigned int pdpmb_free( pdpmb_t * obj ) {
  unsigned int cls = 0;
  while( cls < obj->ngr && obj->gcl[ cls ] > 0 ) { cls++; }
  if( cls == obj->ngr ) { cls = BAD_VCL; }
  return cls;
}

void pdpmb_best( pdpmb_t * obj, unsigned int grp ) {
  unsigned int i, test_cls, best_cls;
  double test_delp=0, best_delp=0;
  best_cls = obj->vcl[ grp ];

  //try and empty cluster 
  if( obj->gcl[ best_cls ] > 1 ) {
    test_cls = pdpmb_free( obj );
    test_delp += pdpmb_movep( obj, grp, test_cls );
    if( test_delp > best_delp ) { 
      best_delp = test_delp;
      best_cls  = test_cls;
    }
  }
  test_cls = 0;

  //try existing clusters
  for( i = 0; i < obj->ncl; i++ ) {
    while( obj->gcl[ test_cls ] == 0 ) { test_cls++; }
    test_delp += pdpmb_movep( obj, grp, test_cls );
    if( test_delp > best_delp ) { 
      best_delp = test_delp;
      best_cls  = test_cls;
    }
    test_cls++;
  }
  if( obj->vcl[ grp ] != best_cls ) { pdpmb_move( obj, grp, best_cls ); }
}

static unsigned int rmulti( double * logp, unsigned int n ) {
  unsigned int i, j, ret = 0;
  double pl, ph, u;
  u = pdpmb_runif(0.0, 1.0);
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

void pdpmb_gibbs( pdpmb_t * obj, int maxiter, double crit) {
  unsigned int i, *vcl_best, cls, grp, iter = 0, test;
  unsigned int cls_old, cls_new, *proposal_cls, proposal_ncl;
  double stopcrit = 1.0, logp_best, *proposal_logp;

  //compute initial logp, save initial partition
  obj->logp = pdpmb_logp( obj );
  logp_best = obj->logp;
  vcl_best = (unsigned int *) pdpmb_alloc( obj, obj->ngr, sizeof( unsigned int ) );
  for( i = 0; i < obj->ngr; i++ ) { vcl_best[ i ] = obj->vcl[ i ]; }
  proposal_logp = (double*) pdpmb_alloc( obj, obj->ngr, sizeof( double ) );
  proposal_cls  = (unsigned int*) pdpmb_alloc( obj, obj->ngr, sizeof( unsigned int ) );

  while( iter++ < maxiter && stopcrit > crit ) {
    for( grp = 0; grp < obj->ngr; grp++ ) {
      cls_old = obj->vcl[ grp ];
      //compute change in logp for possible moves
      cls = 0;
      proposal_ncl = 0;
      for( i = 0; i < obj->ncl; i++ ) {
        while( obj->gcl[ cls ] == 0 ) { cls++; }
        proposal_logp[ i ] = pdpmb_movep( obj, grp, cls );
        proposal_cls[ i ] = cls;
        proposal_ncl++;
        pdpmb_move( obj, grp, cls_old );
        cls++;
      }
      if( obj->gcl[ obj->vcl[ grp ] ] > 1 ) { 
        proposal_cls[ proposal_ncl ] = pdpmb_free( obj );
        proposal_logp[ proposal_ncl ] = pdpmb_movep( obj, grp, proposal_cls[ proposal_ncl ] );
        proposal_ncl++;
        pdpmb_move( obj, grp, cls_old );
      }
      //draw from the conditional
      cls_new = proposal_cls[ rmulti( proposal_logp, proposal_ncl ) ];
      obj->logp += pdpmb_movep( obj, grp, cls_new );
    }
    //save
    if( obj->logp > logp_best ) {
      stopcrit += obj->logp - logp_best;
      logp_best = obj->logp;
      for( i = 0; i < obj->ngr; i++ ) { vcl_best[ i ] = obj->vcl[ i ]; }
    }
    //print summary if requested every iteration
    if( obj->flags & FLAG_VERBOSE && (iter % 1) == 0 ) {
      pdpmb_printf("iter: %u, ncl: %u, logp: %f, best: %f, crit: %f\n",\
                     iter, obj->ncl, obj->logp, logp_best, stopcrit );
    }
    //update stopping criterion
    if(stopcrit != DBL_MAX) { stopcrit *= 0.95; }
  }
  //restore
  obj->logp = logp_best;
  for( grp = 0; grp < obj->ngr; grp++ ) { 
    pdpmb_move( obj, grp, vcl_best[ grp ] );  
  }
}


void pdpmb_stoch( pdpmb_t * obj, int maxiter, double crit) {
  unsigned int i, *vcl_old, *grps, ngrps, cls, iter = 0, spmercls;
  double logp_old, pdel = 1, pcum = 0;

  //allocate memory for vcl_old, grps
  vcl_old = (unsigned int *) pdpmb_alloc( obj, obj->ngr, sizeof( unsigned int ) );
  grps    = (unsigned int *) pdpmb_alloc( obj, obj->ngr, sizeof( unsigned int ) );

  //compute initial logp
  obj->logp = pdpmb_logp( obj );
  while( iter++ < maxiter ) {

    //select the number of groups to shuffle
    ngrps = (unsigned int) floor( obj->ngr * pdpmb_runif( 0.0, 1.0 ) );
    ngrps = ngrps == 0 ? 1 : ngrps;

    //randomly select ngrps groups to shuffle, save indicators
    for( i = 0; i < ngrps; i++ ) {
      grps[ i ] = (unsigned int) floor( obj->ngr * pdpmb_runif( 0.0, 1.0 ) );
      vcl_old[ i ] = obj->vcl[ grps[ i ] ];
    }

    //compute old logp, move groups to random cluster 
    logp_old = obj->logp;
    cls = (unsigned int) floor( obj->ngr * pdpmb_runif( 0.0, 1.0 ) );
    for( i = 0; i < ngrps; i++ ) { pdpmb_move( obj, grps[ i ], cls ); }
    obj->logp = pdpmb_logp( obj );

    if( obj->logp <= logp_old ) {  
      //move each group to best cluster
      for( i = 0; i < ngrps; i++ ) { pdpmb_best( obj, grps[ i ] ); }
      //compute logp, keep new clustering if better, else revert to old
      obj->logp = pdpmb_logp( obj );
      if( obj->logp <= logp_old ) {    
        for( i = 0; i < ngrps; i++ ) { pdpmb_move( obj, grps[ i ], vcl_old[ i ] ); }
        pdel *= 0.9;
        obj->logp = logp_old;
      }
    }

    //update the stopping criterion
    else{ 
      pdel = 0.5 * (obj->logp - logp_old) + 0.5 * pdel;
      logp_old = obj->logp;
    }
    pcum += pdel;

    //print summary if requested every 20 iterations
    if( obj->flags & FLAG_VERBOSE && (iter % 1) == 0 ) {
      pdpmb_printf("iter: %u, ncl: %u, logp: %f, exp: %u, crit: %f\n",\
                    iter, obj->ncl, logp_old, ngrps, pdel / pcum );
    }

    //check stopping criterion, break the optimization loop if met
    if( pcum > 0 && ( pdel / pcum ) < crit ) {
      obj->flags |= FLAG_OPTCRIT;
      break;
    }

  }
}

void pdpmb_merge( pdpmb_t * obj, unsigned int cls1, unsigned int cls2 ) {
  unsigned int i, grp = 0, size;
  //cannot merge an empty group
  if( obj->gcl[ cls1 ] == 0 || obj->gcl[ cls2 ] == 0 ) { return; }
  size = obj->gcl[ cls1 ];
  //merge the cluster
  for( i = 0; i < size; i++ ) {
    while( obj->vcl[ grp ] != cls1 ) { grp++; }
    pdpmb_move( obj, grp, cls2 );
  }
}    


double pdpmb_mergep( pdpmb_t * obj, unsigned int cls1, unsigned int cls2 ) {
  unsigned int only[2];
  double del = 0.0;
  //cannot merge an empty group
  if( obj->gcl[ cls1 ] == 0 || obj->gcl[ cls2 ] == 0 ) { return 0.0; }
  //merge the cluster
  only[0] = cls2;
  only[1] = cls1;
  del -= pdpmb_logponly( obj, only, 2 );
  pdpmb_merge( obj, cls1, cls2 );
  del += pdpmb_logponly( obj, only, 1 );
  return del;
}    

double pdpmb_testmergep( pdpmb_t * obj, unsigned int cls1, unsigned int cls2 ) {
  unsigned int grp = 0, testgrp, size;
  double del = 0.0;

  //enumerate groups in cls1 (cannot use pbuf elsewhere!!!)
  size = obj->gcl[ cls1 ];
  for( testgrp = 0; testgrp < size; testgrp++ ) {
    while( obj->vcl[ grp ] != cls1 ) { grp++; }
    obj->pbuf[ testgrp ] = grp++;
  }

  //merge clusters 
  del = pdpmb_mergep( obj, cls1, cls2 );

  //unmerge clusters
  for( testgrp = 0; testgrp < size; testgrp++ ) {
    pdpmb_move( obj, obj->pbuf[ testgrp ], cls1 );
  }
  return del;
}

void pdpmb_agglo( pdpmb_t * obj, int maxiter ) {
  double *del, del_best, logp_best = -DBL_MAX;
  unsigned int *vcl_best, i, j, index;
  unsigned int icls, jcls, icls_best = BAD_VCL, jcls_best = BAD_VCL;
  unsigned int icls_last = BAD_VCL, jcls_last = BAD_VCL;
  int calcs = 0, cent;

  //allocate some additional memory
  //del[ i, j ] - upper triangular packed storage
  del = (double *) pdpmb_alloc( obj, ( obj->ngr * ( obj->ngr + 1 ) / 2 ), sizeof( double ) );
  vcl_best = (unsigned int *) pdpmb_alloc( obj, obj->ngr, sizeof( unsigned int ) );

  //compute initial logp 
  obj->logp = pdpmb_logp( obj );
  //compute required calculations (and 1%)
  cent = obj->ngr * ( obj->ngr - 1 ) + 1;
  cent = cent > 100 ? cent / 100 : 1;

  //repeat until all clusters are merged into one
  //while( obj->ncl > 1 && maxiter-- != 0 ) {
  while( obj->ncl > 1 ) {
    //compute best merge
    del_best = -DBL_MAX;
    icls = 0;
    for( i = 0; i < obj->ncl - 1; i++ ) {
      while( obj->gcl[ icls ] == 0 ) { icls++; }
      jcls = icls + 1;
      for( j = 0; j < obj->ncl - i - 1; j++ ) {
        while( obj->gcl[ jcls ] == 0 ) { jcls++; }
        index = UMAT(icls, jcls);
        if( obj->ncl == obj->ngr || icls == jcls_last || jcls == jcls_last ) { 
          del[ index ] = pdpmb_testmergep( obj, icls, jcls ); 
          calcs++;
          if( (obj->flags & FLAG_VERBOSE) && (calcs % cent == 0) ) {
            pdpmb_printf("\rpercent complete: %d%", calcs / cent); 
          }
        }
        if( del[ index ] >= del_best ) {
          del_best = del[ index ];
          icls_best = icls;
          jcls_best = jcls;
        }
        jcls++;
      }
      icls++;
    }

    //merge
    pdpmb_merge( obj, icls_best, jcls_best );
    icls_last = icls_best;
    jcls_last = jcls_best;
    obj->logp += del_best;

    //save
    if( obj->logp > logp_best ) {
      logp_best = obj->logp;
      for( i = 0; i < obj->ngr; i++ ) { vcl_best[ i ] = obj->vcl[ i ]; }
    }  
  }

  //restore
  obj->logp = logp_best;
  obj->flags |= FLAG_OPTCRIT;
  for( i = 0; i < obj->ngr; i++ ) {
    pdpmb_move( obj, i, vcl_best[ i ] );
  }

  if( obj->flags & FLAG_VERBOSE ) { 
    pdpmb_printf("\rpercent complete: 100%\n");
    pdpmb_printf("ncl: %u logp: %f\n", obj->ncl, obj->logp);
  }
}
