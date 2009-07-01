#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Memory.h>

#include "util.h"

#define DEFAULT_ALPHA  0.001
#define DEFAULT_A0     0.001
#define DEFAULT_B0     0.001
#define DEFAULT_M0     0
#define DEFAULT_S0     1

#define BAD_INDEX      INT_MAX

/**
 * Struct to contain the pdpmlm data and state of the algorithm.
 * This structure is specififc to a Normal-Gamma base distribution
 * in the DP.
 */
typedef struct {
  double          alpha;    /* prior DP precision */
  double          a0;       /* prior rate parameter */
  double          b0;       /* prior sum of squares */
  double          m0;       /* prior mean */
  double          s0;       /* prior number of `observations' */
  double        * data;     /* data storage vector */
  double        * datasq;   /* data squared storage vector */
  double        * sum;      /* sum for each group */
  double        * sumsq;    /* sum of squares for each group */
  unsigned int  * group;    /* data group vector */
  unsigned int  * index;    /* cluster indicator vector, index = 0 .. size */  
  unsigned int  * unique;   /* number of observations for each index */
  unsigned int  * indmap;   /* index map, contains list of non-empty indices */
  unsigned int    usize;    /* number of non-zero items in unique */
  unsigned int    size;     /* size of data and index vectors */
  double          del;      /* change in posterior for index */
} pdpmlm_t;

/**
 * Initialize indicies sequentially to form grps groups
 * and compute sum, sumsq, and unique. If grps == zero,
 * indices are not assigned before computing sum, sumsq
 * and unique (useful when you need to recompute these
 * values for a given index).
 */
void          pdpmlm_init(pdpmlm_t * obj, unsigned int grps);

/**
 * Initialize index using a simple clustering algorithm
 * and compute sum, sumsq, and unique
 */
void          pdpmlm_init_smart(pdpmlm_t * obj);


/**
 * Get an index such that obj->unique[index] == 0
 */
unsigned int  pdpmlm_empty(pdpmlm_t * obj);

/** 
 * Move an observation to another group
 */
void          pdpmlm_item_move(pdpmlm_t * obj, unsigned int item, unsigned int dst);

/** 
 * Make the best move, return the resulting change in
 * log posterior
 */ 
double        pdpmlm_item_move_best(pdpmlm_t * obj, unsigned int item);

/** 
 * Make the best move away from the current item's index,
 * return the resulting change in log posterior
 */
double        pdpmlm_item_move_best_away(pdpmlm_t * obj, unsigned int item);

/**
 * Compute the change in the value of the 
 * log posterior density function for the move
 */
double        pdpmlm_delpost(pdpmlm_t * obj, unsigned int item, unsigned int dst);

/**
 * Sort group indices to occur consecutively starting
 */
void          pdpmlm_sort(pdpmlm_t * obj);

/**
 * Perform cycle style optimization
 * Go through each observation sequentially,
 * moving it to the most favorable cluster.
 */
void          pdpmlm_cycle(pdpmlm_t * obj, unsigned int iter);

/**
 * Perform shuffle style maximization
 */
void          pdpmlm_shuffle(pdpmlm_t * obj, unsigned int iter, double crit);

/**
 * Perform chunk style maximization
 */
void          pdpmlm_chunk(pdpmlm_t * obj, unsigned int iter, double crit);
