#ifndef PDPMLM_H
#define PDPMLM_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "util.h"

// The following should be changed if ported to another interface
// In addition, allocated memory may need to be freed if ported to
// another interface. This is not necessary if allocated using R's
// memory manager.
#define pdpmlm_printf Rprintf

#define DEFAULT_ALP    1.000
#define DEFAULT_A0     0.001
#define DEFAULT_B0     0.001
#define DEFAULT_M0     0.000
#define DEFAULT_S0     1.000

#define BAD_CLS INT_MAX

// bit masks for flags
#define FLAG_VERBOSE  1<<0  // should routine be verbose
#define FLAG_OPTCRIT  1<<1  // has optimization criterion been met
#define FLAG_OPTSTOP  1<<2  // a stopping iteration has been specified
#define FLAG_PRICLUS  1<<3  // use cluster prior instead of Dirichlet prior
#define FLAG_EMPTY_4  1<<4  // not used
#define FLAG_EMPTY_5  1<<5  // not used
#define FLAG_EMPTY_6  1<<6  // not used
#define FLAG_EMPTY_7  1<<7  // not used


typedef struct {

unsigned char   flags;  // some options

double          alp;  // prior alpha parameter
double          s0;   // prior s0 parameter
double        * m0;   // prior m0 parameter
double          a0;   // prior a0 parameter
double          b0;   // prior b0 parameter

// Each of the p entries in y has a corresponding entry in vgr
// that indicates group membership. The values in vgr begin
// at zero and are increasing. The largest value is ngr-1.
// y and x are sorted before being passed to the C code such
// that the values in vgr are in order
unsigned int  * vgr;  // group vector (array of length p)
unsigned int  * pgr;  // number in each group (array of length ngr)
unsigned int    ngr;  // number of groups

// Each of the ngr groups has an entry in vcl indicating
// group membership. The group indicator is used to index
// the values in vcl. Hence, 'vcl[ 0 ]' would yield the 
// cluster indicator for group 0. The values in vcl may 
// range from 0 to ngr-1. These values are not ordered,
// and may not be continous. However, there will always 
// only be ncl distinct values other than BAD_CLS.
unsigned int  * vcl;  // cluster vector (array of length ngr)
unsigned int  * pcl;  // number in each cluster (array of length ngr)
unsigned int    ncl;  // number of clusters

double        * y;    // y vector (array of length p)
double        * x;    // x matrix (array of length p*q)
unsigned int    p;    // nrow(x)
unsigned int    q;    // ncol(x)

double       ** xxgr;   // x'x matrix (array of ngr arrays of length q*q)
double       ** xygr;   // x'y vector (array of ngr arrays of length q)
double        * yygr;   // y'y scalar (array of length ngr)

double       ** xxcl;   // x'x matrix (array of at least ncl arrays of length q*q)
double       ** xycl;   // x'y matrix (array of at least ncl arrays of length q)
double        * yycl;   // y'y scalar (array of at least length ncl)

double        * s;      // storage for an s matrix (array of length q*q)
double        * m;      // storage for an m vector (array of length q)
double          a;      // storage for an a scalar
double          b;      // storage for an b scalar

double          logp;   // log posterior value

double        * fbuf;   // temporary storage for fortran routines (array of length q)
unsigned int  * pbuf;   // temporary storage for pdpmlm routines (array of length ngr)
unsigned int    mem;    // memory usage counter

} pdpmlm_t;

// Allocate memory and count usage in obj->mem
void *       pdpmlm_alloc( pdpmlm_t * obj, unsigned int count, unsigned int size );

// Assign the observations/groups according to a simple algorithm
void         pdpmlm_divy( pdpmlm_t * obj );

// Add an observation/group to a cluster cls
void         pdpmlm_add( pdpmlm_t * obj, unsigned int grp, unsigned int cls );

// Remove an observation/group from cluster cls
void         pdpmlm_sub( pdpmlm_t * obj, unsigned int grp, unsigned int cls );

// Move an observation/group to cluster cls, should only be called after pdpmlm_add
void         pdpmlm_move( pdpmlm_t * obj, unsigned int grp, unsigned int cls );

// Move an observation/group to cluster cls, return the resulting change in logp
double       pdpmlm_movep( pdpmlm_t * obj, unsigned int grp, unsigned int cls );

// Get the index of an empty cluster, or BAD_CLS if none exist
unsigned int pdpmlm_free( pdpmlm_t * obj );

// Compute the change in logp that would result from a merge of cls1 and cls2
double       pdpmlm_testmergep( pdpmlm_t * obj, unsigned int cls1, unsigned int cls2 );
 
// Merge cluster cls1 into cls2 and return the resulting change in logp
double       pdpmlm_mergep( pdpmlm_t * obj, unsigned int cls1, unsigned int cls2 );
 
// Merge cluster cls1 into cls2
void         pdpmlm_merge( pdpmlm_t * obj, unsigned int cls1, unsigned int cls2 );

// Perform agglomerative clustering
void         pdpmlm_agglo( pdpmlm_t * obj, int stop );

// Try to split cluster cls
double       pdpmlm_splitbest( pdpmlm_t * obj, unsigned int cls );

// Try to merge cluster cls
double       pdpmlm_mergebest( pdpmlm_t * obj, unsigned int cls );

// Move an observation/group to the cluster that minimizes the logp
void         pdpmlm_best( pdpmlm_t * obj, unsigned int grp );

// Compute the posterior parameters s, m, a, and b for a given cluster
void         pdpmlm_parm( pdpmlm_t * obj, unsigned int cls, double * s, double * m, double * a, double * b );

// Compute part of the log posterior value for a particular cluster
double       pdpmlm_logpcls( pdpmlm_t * obj, unsigned int cls );

// Compute the log posterior value for the model 
double       pdpmlm_logp( pdpmlm_t * obj );

// Perform chunk style optimization
void         pdpmlm_chunk( pdpmlm_t * obj, unsigned int itermax, double crit );
#endif
