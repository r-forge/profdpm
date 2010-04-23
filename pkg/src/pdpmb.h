#ifndef PDPMB_H
#define PDPMB_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "util.h"

//some R specific macros
#define pdpmb_printf Rprintf
#define pdpmb_runif  runif

//address upper triangular packed storage matrix
//a00, a01, a11, a02, a12, a22, a03, a13, a23, a33, ...
//i <= j
#define UMAT(i, j) (i + j * ( j + 1 ) / 2)

//absolute value
#define ABS(x) (x < 0 ? -x : x)

#define DEFAULT_LAM    0.000
#define DEFAULT_ALP    1.000
#define DEFAULT_A0     1.000
#define DEFAULT_B0     1.000

#define BAD_VCL UINT_MAX

//bit masks for flags
#define FLAG_VERBOSE  1<<0  //should routine be verbose
#define FLAG_OPTCRIT  1<<1  //has optimization criterion been met
#define FLAG_DIRICHL  1<<2  //use cluster Dirichlet prior
#define FLAG_EMPTY_3  1<<3  //not used
#define FLAG_EMPTY_4  1<<4  //not used
#define FLAG_EMPTY_5  1<<5  //not used
#define FLAG_EMPTY_6  1<<6  //not used
#define FLAG_EMPTY_7  1<<7  //not used

//optimization methods
#define METHOD_NONE   0
#define METHOD_STOCH  1
#define METHOD_AGGLO  2
#define METHOD_GIBBS  3

typedef struct {

unsigned char   flags;  //some options

double          lam;  //prior lambda parameter
double          alp;  //prior alpha parameter
double          a0;   //prior a0 parameter
double          b0;   //prior b0 parameter

//ncl - number of clusters. 
//Each group is assigned to exactly one cluster. Each 
//cluster is made up of one or more groups. ncl is the current 
//number of clusters. ncl must be less than or equal to ngr.
unsigned int    ncl;

//vcl - cluster membership vector
//Each of the ngr groups has an entry in vcl indicating
//group membership. The group indicator is used to index
//the values in vcl. Hence, 'vcl[ 0 ]' would yield the 
//cluster indicator for group 0. The values in vcl may 
//range from 0 to ngr-1. These values are not ordered.
//However, there will always be ncl distinct values 
//other than BAD_CLS. (length ngr)
unsigned int  * vcl;

//gcl - number of groups (rows of y) in each 
//cluster, corresponds to nk in model notation
unsigned int  * gcl;

unsigned int  * y;    //y vector (array of length p*q)
unsigned int    ngr;  //nrow(y)
unsigned int    q;    //ncol(y)

//gqcl - number of 1s in the q'th column of the submatrix
//of y corresponding to each cluster, corresponds to nkj
//in model notation (length ngr*q)
unsigned int  * gqcl;

//a and b are temporary storage variables used to hold
//posterior quantities
double        * a;      //storage for a vector (length q)
double        * b;      //storage for a vector (length q)

double          logp;   //log posterior value
unsigned int  * pbuf;   //temporary storage
unsigned int    mem;    //memory usage counter (using sizeof)

} pdpmb_t;

//Allocate memory and count usage in obj->mem
void *       pdpmb_alloc( pdpmb_t * obj, unsigned int count, unsigned int size );

//Allocate and zero memory and count usage in obj->mem
void *       pdpmb_zalloc( pdpmb_t * obj, unsigned int count, unsigned int size );

//Assign the observations/groups according to a simple algorithm
void         pdpmb_divy( pdpmb_t * obj );

//Add an observation/group to a cluster cls
void         pdpmb_add( pdpmb_t * obj, unsigned int grp, unsigned int cls );

//Remove an observation/group from cluster cls
void         pdpmb_sub( pdpmb_t * obj, unsigned int grp, unsigned int cls );

//Move an observation/group to cluster cls, should only be called after pdpmb_add
void         pdpmb_move( pdpmb_t * obj, unsigned int grp, unsigned int cls );

//Move an observation/group to cluster cls, return the resulting change in logp
double       pdpmb_movep( pdpmb_t * obj, unsigned int grp, unsigned int cls );

//Get the index of an empty cluster, or BAD_CLS if none exist
unsigned int pdpmb_free( pdpmb_t * obj );

//Compute the change in logp that would result from a merge of cls1 and cls2
double       pdpmb_testmergep( pdpmb_t * obj, unsigned int cls1, unsigned int cls2 );
 
//Merge cluster cls1 into cls2 and return the resulting change in logp
double       pdpmb_mergep( pdpmb_t * obj, unsigned int cls1, unsigned int cls2 );
 
//Merge cluster cls1 into cls2
void         pdpmb_merge( pdpmb_t * obj, unsigned int cls1, unsigned int cls2 );

//Move an observation/group to the cluster that minimizes the logp
void         pdpmb_best( pdpmb_t * obj, unsigned int grp );

//Compute the posterior parameters s, m, a, and b for a given cluster
void         pdpmb_parm( pdpmb_t * obj, unsigned int cls, double * a, double * b );

//Compute part of the log posterior value for a particular cluster
double       pdpmb_logpcls( pdpmb_t * obj, unsigned int cls );

//Compute the log posterior value for the model 
double       pdpmb_logp( pdpmb_t * obj );

//Compute the log posterior value only for the clusters given by only[0:(size-1)]
double       pdpmb_logponly( pdpmb_t * obj, unsigned int * only, unsigned int size );

//Perform Gibbs sampler optimization
void         pdpmb_gibbs( pdpmb_t * obj, int maxiter, double crit );

//Perform stochastic optimization
void         pdpmb_stoch( pdpmb_t * obj, int maxiter, double crit );

//Perform agglomerative optimization
void         pdpmb_agglo( pdpmb_t * obj, int maxiter );

#endif
