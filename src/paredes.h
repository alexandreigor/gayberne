/*
 * paredes.h
 *
 *  Created on: 31/03/2012
 *      Author: igor
 */

#ifndef PAREDES_H_
#define PAREDES_H_

#include "gb_main.h"
#include "vetor.h"
#include "math_comp.h"

typedef void
(*funcaoInteracaoParede)(particle part, numericType *u, vetor *forca,
		numericType *ta, double sigma0, double epsilon0);

void paredesCorredor(particle part, numericType *u, vetor *forca,
		numericType *ta, double sigma0, double epsilon0);
void paredesTodas(particle part, numericType *u, vetor *forca,
		numericType *ta, double sigma0, double epsilon0);
void paredeNenhuma(particle part, numericType *u, vetor *forca,
		numericType *ta, double sigma0, double epsilon0);
void paredesRugosas(particle part, numericType *u, vetor *forca,
		numericType *ta, double sigma0, double epsilon0);
void paredesCorredor2(particle part, numericType *u, vetor *forca,
		numericType *ta, double sigma0, double epsilon0);

#endif /* PAREDES_H_ */
