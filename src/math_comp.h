/*
 * math_comp.h
 *
 *  Created on: 03/01/2012
 *      Author: igor
 */

#ifndef MATH_COMP_H_
#define MATH_COMP_H_
#define M_2PI						6.28318530717958647692L
#define FATOR_RADIANO_PARA_GRAU		5.72957795130823208767E1L
#define FATOR_GRAU_PARA_RADIANO		1.74532925199432957692E-2L

double uniformDev(long *seed);

double uniformDev2(double ave, double stdDev, long *semente);

// Números aleatórios com distribuição normal
// com média ave e desvio padrão stdDev.
double gasDev(double ave, double stdDev, long *semente);

long double sqr(long double x);

long double cube(long double x);

long double powInt(const long double x, const unsigned int y);

double rad2grau(double rad);

double grau2rad(double grau);

void pontoMaisProximo(const long double x, const long double y,
		const long double x0, const long double y0, const long double x1,
		const long double y1, const char segmento, long double *px,
		long double *py);

long double distanciaDePontoAteReta(const long double x,
		const long double y, const long double x0, const long double y0,
		const long double x1, const long double y1, const char segmento,
		const char retornarPontoMaisProximo, long double *px, long double *py);

#endif /* MATH_COMP_H_ */
