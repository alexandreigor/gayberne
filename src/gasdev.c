#include <math.h>
#include <stdio.h>

static int iset = 0;
static float gset;

void getGasDevStatus(int *_iset, float *_gset) {
	*_iset = iset;
	*_gset = gset;
}

void reinitGasDev(int _iset, float _gset) {
	iset = _iset;
	gset = _gset;
}

/* Returns a normally distributed deviate
 * with zero mean and unit variance,
 * using ran1(*idum) as the source of
 * uniform deviates.*/
float gasdev(long *idum) {
	float ran1(long *idum);
	float fac, rsq, v1, v2;
	float ret;

	if (iset == 0) {
		do {
			v1 = 2.0 * ran1(idum) - 1.0;
			v2 = 2.0 * ran1(idum) - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq);
		gset = v1 * fac;
		iset = 1;
		ret = v2 * fac;
	} else {
		iset = 0;
		ret = gset;
	}

//	printf("gasdev(*idum = %li) = %f; // iset=%i, gset=%f\n", *idum, ret, iset, gset);
	return ret;
}
/* (C) Copr. 1986-92 Numerical Recipes Software . */
