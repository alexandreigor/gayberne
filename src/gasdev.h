#include "ran1.h"

#ifndef GASDEV_H_
#define GASDEV_H_

void getGasDevStatus(int *_iset, float *_gset);
void reinitGasDev(int iset, float gset);

/* Returns a normally distributed deviate
 * with zero mean and unit variance,
 * using ran1(*idum) as the source of
 * uniform deviates.*/
float gasdev(long *idum);

#endif /* GASDEV_H_ */
