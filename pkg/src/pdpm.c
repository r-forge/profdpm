#include "pdpm.h"

static void pdpm_indmap_add(pdpm_t * obj, unsigned int index) {
  obj->indmap[obj->usize++] = index;
}

static void pdpm_indmap_sub(pdpm_t * obj, unsigned int index) {
  unsigned int i, found = INT_MAX;
  for(i = 0; i < obj->usize; i++) {
    if(obj->indmap[i] == index) { found =  i; break; }
  }
  if(found == INT_MAX) { return; }
  obj->usize--;
  for(i = found; i < obj->usize; i++) {
    obj->indmap[i] = obj->indmap[i+1];
  }
}
    
static void pdpm_init_datasq(pdpm_t * obj) {
  unsigned int i;
  for(i = 0; i < obj->size; i++) {
    obj->datasq[i] = obj->data[i] * obj->data[i];
  }
}

static void pdpm_item_sub(pdpm_t * obj, unsigned int item, unsigned int src) {
  if(obj->unique[src] == 1) { pdpm_indmap_sub(obj, src); }
  obj->sum[src]    -= obj->data[item];
  obj->sumsq[src]  -= obj->datasq[item];
  obj->unique[src] -= 1;
  obj->index[item]  = BAD_INDEX;
}

static void pdpm_item_add(pdpm_t * obj, unsigned int item, unsigned int dst) {
  if(obj->unique[dst] == 0) { pdpm_indmap_add(obj, dst); }
  obj->sum[dst]    += obj->data[item];
  obj->sumsq[dst]  += obj->datasq[item];
  obj->unique[dst] += 1;
  obj->index[item]  = dst;
}

void pdpm_init(pdpm_t * obj, unsigned int grps) {
  unsigned int i, dst;
  obj->usize = 0;
  obj->del = 0;
  memset((void *) obj->sum, 0x00, obj->size * sizeof(double));
  memset((void *) obj->sumsq, 0x00, obj->size * sizeof(double));
  memset((void *) obj->unique, 0x00, obj->size * sizeof(unsigned int));
  pdpm_init_datasq(obj);
  grps = grps > obj->size ? obj->size : grps; 
  for(i = 0; i < obj->size; i++) {
    dst = grps == 0 ? obj->index[i] : i % grps;
    pdpm_item_add(obj, i, dst);
  }
}

void pdpm_init_smart(pdpm_t * obj) {
  unsigned int i, j;
  double mse, mse_max = obj->b0 / obj->a0; 
  obj->usize = 0;
  obj->del = 0;
  memset((void *) obj->sum, 0x00, obj->size * sizeof(double));
  memset((void *) obj->sumsq, 0x00, obj->size * sizeof(double));
  memset((void *) obj->unique, 0x00, obj->size * sizeof(unsigned int));
  pdpm_init_datasq(obj);
  for(i = 0; i < obj->size; i++) {
    for(j = 0; j <= obj->usize; j++) {
      mse = obj->sumsq[j] + obj->datasq[i];
      mse -= pow(obj->sum[j] + obj->data[i],2) /  ( obj->unique[j] + 1 ) ;
      mse /= ( obj->unique[j] + 1 );
      if( mse < mse_max ) {
        pdpm_item_add(obj, i, j);
        break;
      }
    }
  }
}

unsigned int pdpm_empty(pdpm_t * obj) {
  unsigned int i;
  if(obj->usize == obj->size) { return(0); }
  for(i = 0; i < obj->size - 1; i++) {
    if( obj->unique[i] == 0 ) { return(i); }
  }
}

void pdpm_item_move(pdpm_t * obj, unsigned int item, unsigned int dst) {
  unsigned int src;
  src = obj->index[item];
  if( dst >= obj->size || src == dst ) { return; } 
  pdpm_item_sub(obj, item, src);
  pdpm_item_add(obj, item, dst);
}

double pdpm_item_move_best(pdpm_t * obj, unsigned int item) {
  unsigned int i;
  double try_del, best_del = 0, best_int = obj->index[item];
  for(i = 0; i < obj->usize; i++) {
    try_del = pdpm_delpost(obj, item, obj->indmap[i]);
    if(try_del > best_del) { best_del = try_del, best_int = obj->indmap[i]; }
  }
  if(best_del > 0) { pdpm_item_move(obj, item, best_int); }
  return(best_del);
}  

double pdpm_item_move_best_away(pdpm_t * obj, unsigned int item) {
  unsigned int i;
  double try_del, best_del = DBL_MIN, best_int;
  if(obj->usize == 1) { return(0.0); }
  for(i = 0; i < obj->usize; i++) {
    if(obj->indmap[i] != obj->index[item]) { 
      try_del = pdpm_delpost(obj, item, obj->indmap[i]);
      if(try_del > best_del) { best_del = try_del, best_int = obj->indmap[i]; }
    }
  }
  pdpm_item_move(obj, item, best_int);
  return(best_del);
}  

double pdpm_delpost(pdpm_t * obj, unsigned int item, unsigned int dst) {
  double new_sum[2], new_sumsq[2], lpost = 0;
  double mul, s, m, a, b;
  unsigned int new_unique[2], ind[2], i;
  ind[0] = obj->index[item];
  ind[1] = dst;
  if( ind[1] >= obj->size || ind[0] == ind[1] ) { return(0.0); }
  if( obj->unique[ind[0]] == 1 ) { lpost -= log(obj->alpha) - log(obj->usize); }
  if( obj->unique[ind[1]] == 0 ) { lpost += log(obj->alpha) - log(obj->usize+1); }
  for(i = 0; i < 2; i++) {
    mul = i == 0 ? -1 : 1;
    new_sum[i]    = obj->sum[ind[i]] + mul * obj->data[item];
    new_sumsq[i]  = obj->sumsq[ind[i]] + mul * obj->datasq[item];
    new_unique[i] = obj->unique[ind[i]] + mul;
    if( i == 0 || obj->unique[ind[i]] > 0 ) {
      s = obj->s0 + obj->unique[ind[i]];
      m = (obj->s0 * obj->m0 + obj->sum[ind[i]]) / s;
      a = obj->a0 + obj->unique[ind[i]];
      b = obj->b0 + obj->sumsq[ind[i]] + obj->s0*pow(obj->m0,2) - s*pow(m,2);
      lpost -= lfactorial(s - 1) + lgammafn(a/2) - a * log(b/2) / 2 - log(s) / 2;
    }
    if( i == 1 || new_unique[i] > 0 ) {
      s = obj->s0 + new_unique[i];
      m = (obj->s0 * obj->m0 + new_sum[i]) / s;
      a = obj->a0 + new_unique[i];
      b = obj->b0 + new_sumsq[i] + obj->s0*pow(obj->m0,2) - s*pow(m,2);
      lpost += lfactorial(s - 1) + lgammafn(a/2) - a * log(b/2) / 2 - log(s) / 2;
    }
  }
  return(lpost);
}
  

void pdpm_sort(pdpm_t * obj) {
  unsigned int i, count = 0;
  memset((void *) obj->unique, 0xff, obj->size * sizeof(unsigned int));
  for(i = 0; i < obj->size; i++) {
    if(obj->unique[obj->index[i]] > obj->size) {
      obj->unique[obj->index[i]] = count++;
    }
    obj->index[i] = obj->unique[ obj->index[i] ];
  } 
  pdpm_init(obj, 0);
}

void pdpm_cycle(pdpm_t * obj, unsigned int iter) {
  unsigned int i, icount;
  double del = 0;
  for(icount = 0; icount < iter; icount++) {
    for(i = 0; i < obj->size; i++) {
    obj->del +=  pdpm_item_move_best(obj, i);
    }
  }
}

void pdpm_shuffle(pdpm_t * obj, unsigned int iter, double crit) {
  unsigned int i, to, icount, nitems, nsmall, cent, *items, *oldind, dec;
  unsigned int ccount = 0;
  double olddel, del1, del2;
  dec = obj->size >= 200 ? obj->size / 20 : 10;

  /* shuffle 50% of obj->size */
  nitems = (obj->size / 50 + 1);
  items = (unsigned int *) R_alloc(obj->size, sizeof(unsigned int));
  oldind = (unsigned int *) R_alloc(obj->size, sizeof(unsigned int));


  /* set print progress frequency */
  cent = iter >= 100 ? iter / 100 : 1;


  for(icount = 0; icount < iter; icount++) {
    del1 = 0;
    to = pdpm_empty(obj);

    /* save old indices */
    GetRNGstate();
    for(i = 0; i < nitems; i++) {
      items[i] = (unsigned int) (obj->size * runif(0.0, 1.0));
      oldind[items[i]] = obj->index[items[i]];
    }
    PutRNGstate();

    /* move items to new group, compute del associated with moves */
    for(i = 0; i < nitems; i++) {
      del1 += pdpm_delpost(obj, items[i], to); 
      pdpm_item_move(obj, items[i], to);
    }

    /* take the best move back, compute del associated with moves */
    for(i = 0; i < nitems; i++) {
      del1 += pdpm_item_move_best(obj, items[i]);
    }

    /* if the overall del is not positive, move items back to 
     * original index. 
     */
    if( del1 <= 0 ) {
      for(i = 0; i < nitems; i++) {
        pdpm_item_move(obj, items[i], oldind[items[i]]);
      }
    }
    else { 
      /* try to take out small groups */
      del2 = 0;
      nsmall = 0;
      for(i = 0; i < obj->size; i++) {
        if(obj->unique[obj->index[i]] < dec) {
          items[nsmall]    = i;
          oldind[nsmall++] = obj->index[i];
          del2 += pdpm_item_move_best_away(obj, i);
        }
      }
      if( del2 < 0 ) {
        for(i = 0; i < nsmall; i++) {
          pdpm_item_move(obj, items[i], oldind[i]);
        }
      }
      else { del1 += del2; }

      /* check stopping criterion */
      olddel = obj->del;
      obj->del += del1;
      if(olddel == 0 || (del1 / olddel) < crit) {
        if(++ccount > 10) { break; }
      }
    }

    /* print iteration information */
    if(icount % cent == 0) { 
      Rprintf("\roptimizing: iter %6d, change %6f groups %3d", icount, obj->del, obj->usize); 
    } 
  }
  Rprintf("\roptimizing: iter %6d, change %6f groups %3d\n", icount,  obj->del, obj->usize);
}

void pdpm_chunk(pdpm_t * obj, unsigned int iter, double crit) {
  unsigned int i, itemp, nitems, *revert, truesize = obj->size; 
  double dtemp, ran;

  pdpm_init_datasq(obj); 
  nitems = truesize > 1000 ? 1000 : truesize;
  revert = (unsigned int *) R_alloc(nitems, sizeof(unsigned int));
  
  GetRNGstate();
  for(i = 0; i < nitems; i++) {
    ran = runif(0.0, 1.0);
    revert[i] = (unsigned int) (obj->size * ran);
    dtemp = obj->data[i];
    obj->data[i] = obj->data[revert[i]];
    obj->data[revert[i]] = dtemp;
    dtemp = obj->datasq[i];
    obj->datasq[i] = obj->datasq[revert[i]];
    obj->datasq[revert[i]] = dtemp;
  } 
  PutRNGstate();

  obj->size = nitems;
  pdpm_init_smart(obj);
  pdpm_shuffle(obj, iter * 10, crit);

  obj->size = truesize;
  itemp = ( (truesize - nitems) / 100 ) > 0 ? (truesize - nitems) / 100 : 1;
  for(i = nitems; i < truesize; i++) {
    pdpm_item_add(obj, i, pdpm_empty(obj));
    pdpm_item_move_best(obj, i);
  } 

  for(i = 0; i < nitems; i++) {
    dtemp = obj->data[revert[i]];
    obj->data[revert[i]] = obj->data[i]; 
    obj->data[i] = dtemp;
    dtemp = obj->datasq[revert[i]];
    obj->datasq[revert[i]] = obj->datasq[i]; 
    obj->datasq[i] = dtemp;
    itemp = obj->index[revert[i]];
    obj->index[revert[i]] = obj->index[i];
    obj->index[i] = itemp;
  }
 
  pdpm_shuffle(obj, iter, crit); 
  
}

  

  
