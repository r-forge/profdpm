#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <R_ext/PrtUtil.h>

// The following should be changed if ported to another interface
#define pdpmlm_alloc(count, size)  R_alloc(count, size)
#define pdpmlm_printf Rprintf

#define DEFAULT_ALP    0.001
#define DEFAULT_A0     0.001
#define DEFAULT_B0     0.001
#define DEFAULT_M0     0
#define DEFAULT_S0     1

#define BAD_CLS INT_MAX

typedef struct {

double          alp;  // prior alpha parameter
double          s0;   // prior s0 parameter
double        * m0;   // prior m0 parameter
double          a0;   // prior a0 parameter
double          b0;   // prior b0 parameter

unsigned int  * vgr;  // group vector {0,0,...,ngr-1,ngr-1}
unsigned int  * pgr;  // number in each group (ngr)
unsigned int    ngr;  // number of groups

unsigned int  * vcl;  // cluster vector {0,...,ncl-1} 
unsigned int  * pcl;  // number in each cluster (ngr)
unsigned int    ncl;  // number of clusters

double        * y;    // y vector (px1)
double        * x;    // x matrix (pxq)
unsigned int    p;    // nrow(x)
unsigned int    q;    // ncol(x)

double       ** xxgr;   // x'x matrix (qxq) for each group
double       ** xygr;   // x'y vector (qx1) for each group
double        * yygr;   // y'y scalar for each group

double       ** xxcl;   // x'x matrix (qxq) for each cluster
double       ** xycl;   // x'y vector (qx1) for each cluster
double        * yycl;   // y'y scalar for each cluster

double        * s;      // storage for an s matrix (qxq)
double        * m;      // storage for an m vector (qx1)
double          a;      // storage for an a scalar
double          b;      // storage for an b scalar

double        * buf;    // temporary storage (qx1)

} pdpmlm_t;


// Assign the observations/groups in sequence among ncl groups
void pdpmlm_divy( pdpmlm_t * obj, unsigned int ncl );

// Add an observation/group to a cluster new 
void pdpmlm_add( pdpmlm_t * obj, unsigned grp, unsigned int cls );

// Remove an observation/group from cluster clt 
void pdpmlm_sub( pdpmlm_t * obj, unsigned grp, unsigned int cls );

// Move an observation/group from its current cluster to cluster clt
void pdpmlm_move( pdpmlm_t * obj, unsigned int grp, unsigned int cls );

// Compute the posterior parameters s, m, a, and b for a given cluster
void pdpmlm_parm( pdpmlm_t * obj, unsigned int cls, double * s, double * m, double * a, double * b );

// Compute the log posterior value for the model 
double pdpmlm_logp( pdpmlm_t * obj );

// (Note: This function is R specific) This function prints a 
// representation of the pdpmlm_t to the R terminal
void pdpmlm_Rdump( pdpmlm_t * obj );
