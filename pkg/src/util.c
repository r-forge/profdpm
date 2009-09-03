#include "util.h"

/* stats */

double lfactorial(unsigned int x) {
  /* Jolley, L.B.W. (1961) Summation of Series. pp. 28
     ISBN 0-486-60023-8  */
  if( x == 0 ) { return 0; }
  return ( LN_SQRT_2PI + (0.5 + x)*log(x) - x );
}

double factorial(unsigned int x) {
  return(exp(lfactorial(x)));
}

double mean(double * vec, unsigned int size) {
  double *it, sum = 0;
  for(it = vec; it < vec + size; it++) { sum += *it; }
  return( sum / size );
}

double sumsq(double * vec, unsigned int size) {
  double *it, avg, dif, sum = 0;
  avg = mean(vec, size);
  for(it = vec; it < vec + size; it++) { 
    dif = *it-avg; 
    sum += dif * dif; 
  }
  return( sum );
}


double stdev(double * vec, unsigned int size) {
  return( sumsq(vec, size) / (size - 1) );
}
  

/* string */

void strip_all(char * str) {
  unsigned int cnt = 0;
  do {
    if( ' ' == *str || '\t' == *str ) { cnt++; }
    else { *(str - cnt) = *str; } 
  } while( '\0' != *(str++) );
}

void strip_front(char * str) {
  int cnt;
  if( ' ' != *str && '\t' != *str ) { return; }
  else {str++; cnt = 1;}
  do { 
    if( ( ' ' == *str || '\t' == *str ) && cnt > 0) { cnt++; }
    else { if(cnt > 0) { cnt = -cnt; } *(str + cnt) = *str; }
  } while ('\0' != *(str++) );
} 

void strip_rear(char * str) {
  unsigned int cnt = 0;
  while( '\0' != *str ) {
    if( ' ' != *str && '\t' != *str ) { cnt = 0; }
    else { cnt++; }
    str++;
  } 
  *(str - cnt) = '\0';
}

void strip_both(char * str) {
  strip_front(str);
  strip_rear(str);
} 
/* R */

SEXP getListElementByName(SEXP list, const char * name) {
  SEXP elem = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  unsigned int i;
  for (i = 0; i < LENGTH(list); i++) {
    if(strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      elem = VECTOR_ELT(list, i);
      break;
    }
  }
  return(elem);
}
