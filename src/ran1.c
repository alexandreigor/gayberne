#include <stdlib.h>
#include <stdio.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

static long iy = 0;
static long iv[NTAB];

int getNTab(){
	return NTAB;
}

void getRan1Status(long *_iy, long _iv[NTAB]) {
	*_iy = iy;
	register int i;
	for (i = 0; i < NTAB; i++) {
		printf("%li\t",iv[i]);
		_iv[i] = iv[i];
	}
	printf("\n");
}

void reinitRan1(const long _iy, const long _iv[NTAB]) {
	iy = _iy;
	register int i;
	for (i = 0; i < NTAB; i++) {
		iv[i] = _iv[i];
	}
}

/*
 * "Minimal" random number generator of
 * Park and Miller with Bays-Durham shuffle
 * and added safeguards. Returns a uniform
 * random deviate between 0.0 and 1.0
 * (exclusive of the endpoint values).
 * Call with idum a negative integer to
 * initialize; thereafter, do not alter
 * idum between successive deviates in a
 * sequence. RNMX should approximate the
 * largest floating value that is less than 1.
 */

float ran1(long *idum) {
	int j;
	long k;
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1)
			*idum = 1;
		else
			*idum = -(*idum);
		for (j = NTAB + 7; j >= 0; j--) {
			k = (*idum) / IQ;
			*idum = IA * (*idum - k * IQ) - IR * k;
			if (*idum < 0)
				*idum += IM;
			if (j < NTAB)
				iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ;
	*idum = IA * (*idum - k * IQ) - IR * k;
	if (*idum < 0)
		*idum += IM;
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = *idum;
	if ((temp = AM * iy) > RNMX)
		return RNMX;
	else
		return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

/*
 * "Minimal" random number generator of
 * Park and Miller. Returns a uniform
 * random deviate between 0.0 and 1.0.
 * Set or reset idum to any integer
 * value (except the unlikely value MASK)
 * to initialize the sequence; idum must
 * not be altered between calls for
 * successive deviates in a sequence.
 * */
float ran0(long *idum) {
	long k;
	float ans;
	*idum ^= MASK;
	k = (*idum) / IQ;
	*idum = IA * (*idum - k * IQ) - IR * k;
	if (*idum < 0)
		*idum += IM;
	ans = AM * (*idum);
	*idum ^= MASK;
//	printf("ran0(*idum = %li) = %f;\n", *idum, ans);
	return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK
/* (C) Copr. 1986-92 Numerical Recipes Software . */
