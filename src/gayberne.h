#include "vetor.h"

#ifndef GAYBERNE_H_
#define GAYBERNE_H_

#define RAIZ_6_DE_2		1.122462048309372981433
#define RAIZ_6_DE_2_2	0.5612310241546864907165
#define EQUILIBRIO_UTIL 0.122462048309372981433 // = 2^(1/6) - 1
typedef enum potencial_range {
	POTENCIAL_TOTAL_RANGE = 0,
	POTENCIAL_SIGMA_RANGE = 1,
	POTENCIAL_EQUIL_RANGE = 2
} potencial_range;

numericType distanciaEquilibrioGBV(const vetor r,
		const numericType theta_i, const numericType theta_j,
		const numericType chi, const numericType sigma0);

numericType distanciaContatoGB(const numericType r_x,
		const numericType r_y, const numericType theta_i,
		const numericType theta_j, const numericType chi,
		const numericType sigma0);

numericType distanciaContatoGBV(const vetor r,
		const numericType theta_i, const numericType theta_j,
		const numericType chi, const numericType sigma0);

char overlapGB(const numericType r_x, const numericType r_y,
		const numericType theta_i, const numericType theta_j,
		const numericType chi, const numericType sigma0);

char overlapGBV(const vetor r, const numericType theta_i,
		const numericType theta_j, const numericType chi,
		const numericType sigma0);

void interacaoGayBerne(const numericType r_x, const numericType r_y,
		const numericType theta_i, const numericType theta_j,
		const numericType chi, const numericType chi_,
		const numericType sigma0, const numericType epsilon0, const double nu,
		const double mu, numericType *U, numericType *forca_i_x,
		numericType *forca_i_y, numericType *forca_j_x, numericType *forca_j_y,
		numericType *torque_i, numericType *torque_j);

void interacaoGayBerneV(const vetor rij, const numericType theta_i,
		const numericType theta_j, const numericType chi,
		const numericType chi_, const numericType sigma0,
		const numericType epsilon0, const double nu, const double mu,
		numericType *U, vetor *forca_i, vetor *forca_j, numericType *torque_i,
		numericType *torque_j);

void
interacaoGayBerneV_Rep(const vetor rij, const numericType theta_i,
		const numericType theta_j, const numericType chi,
		const numericType chi_, const numericType sigma0,
		const numericType epsilon0, const double nu, const double mu,
		numericType *U, vetor *forca_i, vetor *forca_j, numericType *torque_i,
		numericType *torque_j, numericType *profPoco,
		const potencial_range apenasRepulsao);

numericType profundidadePocoGBV(const vetor r,
		const numericType theta_i, const numericType theta_j,
		const numericType epsilon0, const numericType chi,
		const numericType chi_, const double nu, const double mu);

numericType profundidadePocoGB(const numericType r_x,
		const numericType r_y, const numericType theta_i,
		const numericType theta_j, const numericType epsilon0,
		const numericType chi, const numericType chi_, const double nu,
		const double mu);

void interacaoNova(const vetor rij, const numericType theta_i,
		const numericType theta_j, const numericType chi,
		const numericType chi_, const numericType sigma0,
		const numericType epsilon0, const double nu, const double mu,
		numericType *U, vetor *forca_i, vetor *forca_j, numericType *torque_i,
		numericType *torque_j);

void interacaoGayBerneV_WallV(numericType y_i,
		const numericType theta_i, const numericType chi,
		const numericType chi_, const numericType sigma0,
		const numericType epsilon0, const double nu, const double mu,
		numericType *U, vetor *forca_i, numericType *torque_i,
		numericType *profPoco, const numericType wTop, const numericType wBot);

void interacaoGayBerneV_WallH(numericType x_i,
		const numericType theta_i, const numericType chi,
		const numericType chi_, const numericType sigma0,
		const numericType epsilon0, const double nu, const double mu,
		numericType *U, vetor *forca_i, numericType *torque_i,
		numericType *profPoco, const numericType wRight,
		const numericType wLeft);
//
//void interacaoNova_2(const vetor rij, const numericType theta_i,
//		const numericType theta_j, const numericType chi,
//		const numericType chi_, const numericType sigma0,
//		const numericType epsilon0, const double nu, const double mu,
//		numericType *U, vetor *forca_i, vetor *forca_j, numericType *torque_i,
//		numericType *torque_j);
#endif /* GAYBERNE_H_ */
