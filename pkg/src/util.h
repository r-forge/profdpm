#include <math.h>

/* log(sqrt(2*pi)) == log(2*pi)/2 */
#define LN_SQRT_2PI 0.918938533204672741780329736406

/* stats */
double lfactorial(unsigned int x);
double factorial(unsigned int x);
double mean(double * vec, unsigned int size);
double sumsq(double * vec, unsigned int size);
double stdev(double * vec, unsigned int size);

/* string */

void strip_all(char * str);
void strip_front(char * str);
void strip_rear(char * str);
void strip_both(char * str); 

