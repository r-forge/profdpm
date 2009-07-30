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
  Rprintf("step 1 ok\n");

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
  Rprintf("step 2 ok\n");
   
  // 3. recompute xxcl, xycl yycl for cls using xxgr, xygr, and yygr from grp
  obj->yycl[ cls ] += obj->yygr[ grp ];
  for( i = 0; i < obj->q; i++ ) {
    obj->xycl[ cls ][ i ] += obj->xygr[ grp ][ i ];
    for( j = 0; j < obj->q; j++ ) {
      obj->xxcl[ cls ][ j + i * obj->q ] += obj->xxgr[ grp ][ j + i * obj->q ];
    }
  }
  Rprintf("step 3 ok\n");
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

void pdpmlm_move( pdpmlm_t * obj, unsigned int grp, unsigned int new ) {
  unsigned int old = obj->vcl[ grp ];
  if( old == new ) { return; }
  pdpmlm_sub( obj, grp, old );
  pdpmlm_add( obj, grp, new );
}

double pdpmlm_logp( pdpmlm_t * obj ) {
  // 1. compute the log posterior value from s, m, a, and b
  return 0.0;
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
    pdpmlm_printf( "xygr\n" );
    printRealVector( obj->xycl[ cls ], obj->q, 1 );
    pdpmlm_printf( "yygr\n" );
    printRealVector( &obj->yycl[ cls ], 1, 1 );
  }
}

