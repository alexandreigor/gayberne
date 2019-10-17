#include "gayberne.h"
#include "vetor.h"
#include "math_comp.h"
#include <math.h>

extern long seedLong;

// Distância em que a força é zero (apenas parte radial)
inline numericType distanciaEquilibrioGBV(const vetor r,
		const numericType theta_i, const numericType theta_j,
		const numericType chi, const numericType sigma0) {
	return distanciaContatoGBV(r, theta_i, theta_j, chi, sigma0) + sigma0
			* EQUILIBRIO_UTIL;
}

// Distância na qual o potencial, é zero.
inline numericType distanciaContatoGB(const numericType r_x,
		const numericType r_y, const numericType theta_i,
		const numericType theta_j, const numericType chi,
		const numericType sigma0) {
	return distanciaContatoGBV(newVetor(r_x, r_y), theta_i, theta_j, chi,
			sigma0);
}

inline numericType distanciaContatoGBV(const vetor r,
		const numericType theta_i, const numericType theta_j,
		const numericType chi, const numericType sigma0) {
	numericType chiHalf = 0.5 * chi;
	numericType gamma = atan2(r.y, r.x);
	numericType cos_theta_i_menos_j = cos(theta_i - theta_j);
	numericType cos_gamma_menos_theta_i = cos(gamma - theta_i);
	numericType cos_gamma_menos_theta_j = cos(gamma - theta_j);
	numericType somaCos2 = sqr(
			cos_gamma_menos_theta_i + cos_gamma_menos_theta_j);
	numericType difCos2 =
			sqr(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j);

	numericType denSomaChi = 1 + chi * cos_theta_i_menos_j;
	numericType denDifChi = 1 - chi * cos_theta_i_menos_j;
	numericType denSomaChi_inv = 1.0 / denSomaChi;
	numericType denDifChi_inv = 1.0 / denDifChi;
	numericType R_chi = 1 - chiHalf * (somaCos2 * denSomaChi_inv + difCos2
			* denDifChi_inv);
	return sigma0 / sqrt(R_chi);
}

inline numericType profundidadePocoGB(const numericType r_x,
		const numericType r_y, const numericType theta_i,
		const numericType theta_j, const numericType epsilon0,
		const numericType chi, const numericType chi_, const double nu,
		const double mu) {
	return profundidadePocoGBV(newVetor(r_x, r_y), theta_i, theta_j, epsilon0,
			chi, chi_, nu, mu);
}

inline numericType profundidadePocoGBV(const vetor r,
		const numericType theta_i, const numericType theta_j,
		const numericType epsilon0, const numericType chi,
		const numericType chi_, const double nu, const double mu) {

	numericType gamma = atan2(r.y, r.x);
	numericType chi_Half = 0.5 * chi_;

	numericType cos_theta_i_menos_j = cos(theta_i - theta_j);
	numericType cos_gamma_menos_theta_i = cos(gamma - theta_i);
	numericType cos_gamma_menos_theta_j = cos(gamma - theta_j);

	numericType somaCos2 = sqr(
			cos_gamma_menos_theta_i + cos_gamma_menos_theta_j);
	numericType difCos2 =
			sqr(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j);

	numericType denSomaChi = 1 + chi_ * cos_theta_i_menos_j;
	numericType denDifChi = 1 - chi_ * cos_theta_i_menos_j;
	numericType denSomaChi_inv = 1.0 / denSomaChi;
	numericType denDifChi_inv = 1.0 / denDifChi;

	numericType R_chi_ = 1 - chi_Half * (somaCos2 * denSomaChi_inv + difCos2
			* denDifChi_inv);

	numericType epsilonLinha = R_chi_;
	numericType denInv = 1.0 / (1.0 - sqr(chi * cos_theta_i_menos_j));
	numericType epsilon2 = sqrt(denInv);

	numericType epsilon2_nu = pow(epsilon2, nu);
	numericType epsilonLinha_mu = pow(epsilonLinha, mu);

	numericType epsilon = epsilon0 * epsilon2_nu * epsilonLinha_mu;
	return -epsilon;
}

inline char overlapGB(const numericType r_x, const numericType r_y,
		const numericType theta_i, const numericType theta_j,
		const numericType chi, const numericType sigma0) {
	return overlapGBV(newVetor(r_x, r_y), theta_i, theta_j, chi, sigma0);
}

inline char overlapGBV(const vetor r, const numericType theta_i,
		const numericType theta_j, const numericType chi,
		const numericType sigma0) {
	numericType r_ = length(r);
	numericType sigma =
			distanciaEquilibrioGBV(r, theta_i, theta_j, chi, sigma0);

	return (r_ <= sigma);
}

inline void interacaoGayBerne(const numericType r_x, const numericType r_y,
		const numericType theta_i, const numericType theta_j,
		const numericType chi, const numericType chi_,
		const numericType sigma0, const numericType epsilon0, const double nu,
		const double mu, numericType *U, numericType *forca_i_x,
		numericType *forca_i_y, numericType *forca_j_x, numericType *forca_j_y,
		numericType *torque_i, numericType *torque_j) {
	vetor rij = newVetor(r_x, r_y);
	vetor forca_i, forca_j;
	interacaoGayBerneV(rij, theta_i, theta_j, chi, chi_, sigma0, epsilon0, nu,
			mu, U, &forca_i, &forca_j, torque_i, torque_j);
	*forca_i_x = forca_i.x;
	*forca_i_y = forca_i.y;
	*forca_j_x = forca_j.x;
	*forca_j_y = forca_j.y;
}

inline void interacaoGayBerneV(const vetor rij, const numericType theta_i,
		const numericType theta_j, const numericType chi,
		const numericType chi_, const numericType sigma0,
		const numericType epsilon0, const double nu, const double mu,
		numericType *U, vetor *forca_i, vetor *forca_j, numericType *torque_i,
		numericType *torque_j) {

	numericType chiHalf = 0.5 * chi;
	numericType chi_Half = 0.5 * chi_;
	numericType dij = length(rij);
	vetor eij = divEscalar(dij, rij);
	numericType gamma = atan2(rij.y, rij.x);

	//numericType f_ = 0.06;
	numericType cos_theta_i_menos_j = cos(theta_i - theta_j);
	numericType cos_gamma_menos_theta_i = cos(gamma - theta_i);
	numericType cos_gamma_menos_theta_j = cos(gamma - theta_j);

	numericType sin_theta_i_menos_j = sin(theta_i - theta_j);
	numericType sin_gamma_menos_theta_i = sin(gamma - theta_i);
	numericType sin_gamma_menos_theta_j = sin(gamma - theta_j);

	numericType somaCos2 = sqr(
			cos_gamma_menos_theta_i + cos_gamma_menos_theta_j);
	numericType difCos2 =
			sqr(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j);

	numericType sigma;
	numericType dSigma_i, dSigma_j, dSigma_g;
	{
		numericType denSomaChi = 1 + chi * cos_theta_i_menos_j;
		numericType denDifChi = 1 - chi * cos_theta_i_menos_j;
		numericType denSomaChi_inv = 1.0 / denSomaChi;
		numericType denDifChi_inv = 1.0 / denDifChi;
		numericType denSomaChi_inv2 = sqr(denSomaChi_inv);
		numericType denDifChi_inv2 = sqr(denDifChi_inv);
		numericType R_chi = 1 - chiHalf * (somaCos2 * denSomaChi_inv + difCos2
				* denDifChi_inv);
		sigma = sigma0 / sqrt(R_chi);
		numericType dRi, dRj, dRg;
		dRi = -chiHalf * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_i * denSomaChi + chi * somaCos2
						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(2
								* (cos_gamma_menos_theta_i
										- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_i * denDifChi - chi
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRj = -chiHalf * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_j * denSomaChi - chi * somaCos2
						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(-2 * (cos_gamma_menos_theta_i
								- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_j * denDifChi + chi
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRg = chi * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
				* denSomaChi_inv + //
				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
						* denDifChi_inv);

		dSigma_i = -0.5 * (sigma / R_chi) * dRi;
		dSigma_j = -0.5 * (sigma / R_chi) * dRj;
		dSigma_g = -0.5 * (sigma / R_chi) * dRg;
	}
	numericType epsilonLinha;
	numericType dEpsilonLinha_i, dEpsilonLinha_j, dEpsilonLinha_g;
	{

		numericType denSomaChi = 1 + chi_ * cos_theta_i_menos_j;
		numericType denDifChi = 1 - chi_ * cos_theta_i_menos_j;
		numericType denSomaChi_inv = 1.0 / denSomaChi;
		numericType denDifChi_inv = 1.0 / denDifChi;
		numericType denSomaChi_inv2 = sqr(denSomaChi_inv);
		numericType denDifChi_inv2 = sqr(denDifChi_inv);

		numericType R_chi_ = 1 - chi_Half * (somaCos2 * denSomaChi_inv
				+ difCos2 * denDifChi_inv);
		epsilonLinha = R_chi_;

		numericType dRi, dRj, dRg;
		dRi = -chi_Half * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_i * denSomaChi + chi_
						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(2
								* (cos_gamma_menos_theta_i
										- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_i * denDifChi - chi_
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRj = -chi_Half * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_j * denSomaChi - chi_
						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(-2 * (cos_gamma_menos_theta_i
								- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_j * denDifChi + chi_
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRg = chi_ * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
				* denSomaChi_inv + //
				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
						* denDifChi_inv);
		dEpsilonLinha_i = dRi;
		dEpsilonLinha_j = dRj;
		dEpsilonLinha_g = dRg;
	}
	numericType epsilon2;
	numericType dEpsilon2_i, dEpsilon2_j, dEpsilon2_g;
	{
		numericType denInv = 1.0 / (1.0 - sqr(chi * cos_theta_i_menos_j));
		epsilon2 = sqrt(denInv);

		dEpsilon2_g = 0.0;
		dEpsilon2_i = -epsilon2 * sqr(chi) * cos_theta_i_menos_j
				* sin_theta_i_menos_j * denInv;
		dEpsilon2_j = -dEpsilon2_i;
	}
	numericType epsilon;
	numericType dEpsilon_i, dEpsilon_j, dEpsilon_g;
	{
		numericType epsilon2_nu_menos_1 = pow(epsilon2, nu - 1);
		numericType epsilonLinha_nu_menos_1 = pow(epsilonLinha, mu - 1);
		numericType prodComum = epsilon0 * epsilon2_nu_menos_1
				* epsilonLinha_nu_menos_1;

		epsilon = prodComum * epsilon2 * epsilonLinha;

		dEpsilon_i = prodComum * (nu * epsilonLinha * dEpsilon2_i + mu
				* epsilon2 * dEpsilonLinha_i);
		dEpsilon_j = prodComum * (nu * epsilonLinha * dEpsilon2_j + mu
				* epsilon2 * dEpsilonLinha_j);
		dEpsilon_g = prodComum * (nu * epsilonLinha * dEpsilon2_g + mu
				* epsilon2 * dEpsilonLinha_g);
	}
	numericType fra, fra6;
	fra = sigma0 / (dij - sigma + sigma0);
	fra6 = sqr(cube(fra));

	numericType dUi, dUj, dUr, dUg, dU_den;

	*U = 4 * epsilon * fra6 * (fra6 - 1.0);

	dU_den = 48.0 * epsilon / sigma0 * fra * fra6 * (fra6 - 0.5);

	dUr = -dU_den;
	dUi = 4 * fra6 * (fra6 - 1.0) * dEpsilon_i + dU_den * dSigma_i;
	dUj = 4 * fra6 * (fra6 - 1.0) * dEpsilon_j + dU_den * dSigma_j;
	dUg = 4 * fra6 * (fra6 - 1.0) * dEpsilon_g + dU_den * dSigma_g;

	numericType f_x = -dUg * eij.y / dij + dUr * eij.x;
	numericType f_y = +dUg * eij.x / dij + dUr * eij.y;

	*forca_i = newVetor(f_x, f_y);
	*forca_j = newVetor(-f_x, -f_y);

	*torque_i = -dUi;
	*torque_j = -dUj;
}

inline void interacaoGayBerneV_Rep(const vetor rij, const numericType theta_i,
		const numericType theta_j, const numericType chi,
		const numericType chi_, const numericType sigma0,
		const numericType epsilon0, const double nu, const double mu,
		numericType *U, vetor *forca_i, vetor *forca_j, numericType *torque_i,
		numericType *torque_j, numericType *profPoco,
		const potencial_range potencialRange) {

	numericType chiHalf = 0.5 * chi;
	numericType chi_Half = 0.5 * chi_;
	numericType dij = length(rij);
	vetor eij = divEscalar(dij, rij);
	numericType gamma = atan2(rij.y, rij.x);

	//numericType f_ = 0.06;
	numericType cos_theta_i_menos_j = cos(theta_i - theta_j);
	numericType cos_gamma_menos_theta_i = cos(gamma - theta_i);
	numericType cos_gamma_menos_theta_j = cos(gamma - theta_j);

	numericType sin_theta_i_menos_j, sin_gamma_menos_theta_i,
			sin_gamma_menos_theta_j; // Inicialização adiada

	numericType somaCos2 = sqr(
			cos_gamma_menos_theta_i + cos_gamma_menos_theta_j);
	numericType difCos2 =
			sqr(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j);

	numericType sigma;
	numericType dSigma_i, dSigma_j, dSigma_g;
	{
		numericType denSomaChi = 1 + chi * cos_theta_i_menos_j;
		numericType denDifChi = 1 - chi * cos_theta_i_menos_j;
		numericType denSomaChi_inv = 1.0 / denSomaChi;
		numericType denDifChi_inv = 1.0 / denDifChi;
		numericType denSomaChi_inv2, denDifChi_inv2; // inicialização adiada;
		numericType R_chi = 1 - chiHalf * (somaCos2 * denSomaChi_inv + difCos2
				* denDifChi_inv);
		sigma = sigma0 / sqrt(R_chi);

		if (potencialRange) {
			numericType distCorte = 0.0;
			if (potencialRange == POTENCIAL_SIGMA_RANGE) { // Corte em sigma
				distCorte = sigma;
			} else { // Corte em r_p
				distCorte = sigma + sigma0 * EQUILIBRIO_UTIL;
			}
			if (dij > distCorte) {
				*U = 0;

				*forca_i = newVetor(0.0, 0.0);
				*forca_j = newVetor(0.0, 0.0);

				*torque_i = 0.0;
				*torque_j = 0.0;

				*profPoco = 0.0;
				return;
			}
		}

		denSomaChi_inv2 = sqr(denSomaChi_inv);
		denDifChi_inv2 = sqr(denDifChi_inv);

		sin_theta_i_menos_j = sin(theta_i - theta_j);
		sin_gamma_menos_theta_i = sin(gamma - theta_i);
		sin_gamma_menos_theta_j = sin(gamma - theta_j);

		numericType dRi, dRj, dRg;
		dRi = -chiHalf * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_i * denSomaChi + chi * somaCos2
						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(2
								* (cos_gamma_menos_theta_i
										- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_i * denDifChi - chi
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRj = -chiHalf * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_j * denSomaChi - chi * somaCos2
						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(-2 * (cos_gamma_menos_theta_i
								- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_j * denDifChi + chi
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRg = chi * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
				* denSomaChi_inv + //
				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
						* denDifChi_inv);

		dSigma_i = -0.5 * (sigma / R_chi) * dRi;
		dSigma_j = -0.5 * (sigma / R_chi) * dRj;
		dSigma_g = -0.5 * (sigma / R_chi) * dRg;
	}
	numericType epsilonLinha;
	numericType dEpsilonLinha_i, dEpsilonLinha_j, dEpsilonLinha_g;
	{

		numericType denSomaChi = 1 + chi_ * cos_theta_i_menos_j;
		numericType denDifChi = 1 - chi_ * cos_theta_i_menos_j;
		numericType denSomaChi_inv = 1.0 / denSomaChi;
		numericType denDifChi_inv = 1.0 / denDifChi;
		numericType denSomaChi_inv2 = sqr(denSomaChi_inv);
		numericType denDifChi_inv2 = sqr(denDifChi_inv);

		numericType R_chi_ = 1 - chi_Half * (somaCos2 * denSomaChi_inv
				+ difCos2 * denDifChi_inv);
		epsilonLinha = R_chi_;

		numericType dRi, dRj, dRg;
		dRi = -chi_Half * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_i * denSomaChi + chi_
						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(2
								* (cos_gamma_menos_theta_i
										- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_i * denDifChi - chi_
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRj = -chi_Half * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_j * denSomaChi - chi_
						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(-2 * (cos_gamma_menos_theta_i
								- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_j * denDifChi + chi_
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRg = chi_ * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
				* denSomaChi_inv + //
				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
						* denDifChi_inv);
		dEpsilonLinha_i = dRi;
		dEpsilonLinha_j = dRj;
		dEpsilonLinha_g = dRg;
	}
	numericType epsilon2;
	numericType dEpsilon2_i, dEpsilon2_j, dEpsilon2_g;
	{
		numericType denInv = 1.0 / (1.0 - sqr(chi * cos_theta_i_menos_j));
		epsilon2 = sqrt(denInv);

		dEpsilon2_g = 0.0;
		dEpsilon2_i = -epsilon2 * sqr(chi) * cos_theta_i_menos_j
				* sin_theta_i_menos_j * denInv;
		dEpsilon2_j = -dEpsilon2_i;
	}
	numericType epsilon;
	numericType dEpsilon_i, dEpsilon_j, dEpsilon_g;
	{
		numericType epsilon2_nu_menos_1 = pow(epsilon2, nu - 1);
		numericType epsilonLinha_nu_menos_1 = pow(epsilonLinha, mu - 1);
		numericType prodComum = epsilon0 * epsilon2_nu_menos_1
				* epsilonLinha_nu_menos_1;

		epsilon = prodComum * epsilon2 * epsilonLinha;

		dEpsilon_i = prodComum * (nu * epsilonLinha * dEpsilon2_i + mu
				* epsilon2 * dEpsilonLinha_i);
		dEpsilon_j = prodComum * (nu * epsilonLinha * dEpsilon2_j + mu
				* epsilon2 * dEpsilonLinha_j);
		dEpsilon_g = prodComum * (nu * epsilonLinha * dEpsilon2_g + mu
				* epsilon2 * dEpsilonLinha_g);
	}
	numericType fra, fra6;
	fra = sigma0 / (dij - sigma + sigma0);
	fra6 = sqr(cube(fra));

	numericType dUi, dUj, dUr, dUg, dU_den;

	*U = 4 * epsilon * fra6 * (fra6 - 1.0);

	dU_den = 48.0 * epsilon / sigma0 * fra * fra6 * (fra6 - 0.5);

	dUr = -dU_den;
	dUi = 4 * fra6 * (fra6 - 1.0) * dEpsilon_i + dU_den * dSigma_i;
	dUj = 4 * fra6 * (fra6 - 1.0) * dEpsilon_j + dU_den * dSigma_j;
	dUg = 4 * fra6 * (fra6 - 1.0) * dEpsilon_g + dU_den * dSigma_g;

	numericType f_x = -dUg * eij.y / dij + dUr * eij.x;
	numericType f_y = +dUg * eij.x / dij + dUr * eij.y;

	*forca_i = newVetor(f_x, f_y);
	*forca_j = newVetor(-f_x, -f_y);

	*torque_i = -dUi;
	*torque_j = -dUj;

	if (1 == 0) {
		//TODO: random
		vetor eij2 = newVetorByDirection(1.0, *forca_j);
		double theta = 2E2;
		double f = theta * uniformDev(&seedLong);

		vetor nov = multEscalar(f, eij2);
		inc(forca_j, nov);
		dec(forca_i, nov);
	}

	if (potencialRange == POTENCIAL_EQUIL_RANGE) {
		*profPoco = -epsilon;
	} else {
		*profPoco = 0.0;
	}
}

inline void interacaoGayBerneV_WallV(numericType y_i,
		const numericType theta_i, const numericType chi,
		const numericType chi_, const numericType sigma0,
		const numericType epsilon0, const double nu, const double mu,
		numericType *U, vetor *forca_i, numericType *torque_i,
		numericType *profPoco, const numericType wTop, const numericType wBot) {

	typedef enum _dir {
		CIMA = 687897, BAIXO = -129837
	} dir;
	dir pos;
	{
		numericType r_y_1 = y_i - wBot;
		numericType r_y_2 = y_i - wTop;
		if (fabs(r_y_1) > fabs(r_y_2)) {
			pos = CIMA;
		} else {
			pos = BAIXO;
		}
	}
	numericType theta_j = -theta_i;

	numericType chiHalf = 0.5 * chi;
	numericType chi_Half = 0.5 * chi_;
	numericType dij;
	vetor rij, eij;
	numericType gamma = M_PI_2;

	//numericType f_ = 0.06;
	numericType cos_theta_i_menos_j = cos(theta_i - theta_j);
	numericType cos_gamma_menos_theta_i = cos(gamma - theta_i);
	numericType cos_gamma_menos_theta_j = cos(gamma - theta_j);

	numericType sin_theta_i_menos_j, sin_gamma_menos_theta_i,
			sin_gamma_menos_theta_j; // Inicialização adiada

	numericType somaCos2 = sqr(
			cos_gamma_menos_theta_i + cos_gamma_menos_theta_j);
	numericType difCos2 =
			sqr(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j);

	numericType sigma;
	numericType dSigma_i, dSigma_j, dSigma_g;

	{

		numericType denSomaChi = 1 + chi * cos_theta_i_menos_j;
		numericType denDifChi = 1 - chi * cos_theta_i_menos_j;
		numericType denSomaChi_inv = 1.0 / denSomaChi;
		numericType denDifChi_inv = 1.0 / denDifChi;
		numericType R_chi = 1 - chiHalf * (somaCos2 * denSomaChi_inv + difCos2
				* denDifChi_inv);
		numericType denSomaChi_inv2, denDifChi_inv2; // inicialização adiada;
		sigma = sigma0 / sqrt(R_chi);

		rij = sub(
				newVetor(
						0,
						(pos == CIMA) ? (wTop + 0.5 * sigma) : (wBot - 0.5
								* sigma)), newVetor(0, y_i));
		dij = length(rij);
		eij = divEscalar(dij, rij);

		numericType distCorte = sigma + sigma0 * EQUILIBRIO_UTIL;
		if (dij > distCorte) {
			*U = 0;

			*forca_i = newVetor(0.0, 0.0);

			*torque_i = 0.0;

			*profPoco = 0.0;
			return;
		}

		denSomaChi_inv2 = sqr(denSomaChi_inv);
		denDifChi_inv2 = sqr(denDifChi_inv);

		sin_theta_i_menos_j = sin(theta_i - theta_j);
		sin_gamma_menos_theta_i = sin(gamma - theta_i);
		sin_gamma_menos_theta_j = sin(gamma - theta_j);

		numericType dRi, dRj, dRg;
		dRi = -chiHalf * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_i * denSomaChi + chi * somaCos2
						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(2
								* (cos_gamma_menos_theta_i
										- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_i * denDifChi - chi
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRj = -chiHalf * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_j * denSomaChi - chi * somaCos2
						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(-2 * (cos_gamma_menos_theta_i
								- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_j * denDifChi + chi
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRg = chi * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
				* denSomaChi_inv + //
				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
						* denDifChi_inv);

		dSigma_i = -0.5 * (sigma / R_chi) * dRi;
		dSigma_j = -0.5 * (sigma / R_chi) * dRj;
		dSigma_g = -0.5 * (sigma / R_chi) * dRg;
	}
	numericType epsilonLinha;
	numericType dEpsilonLinha_i, dEpsilonLinha_j, dEpsilonLinha_g;
	{

		numericType denSomaChi = 1 + chi_ * cos_theta_i_menos_j;
		numericType denDifChi = 1 - chi_ * cos_theta_i_menos_j;
		numericType denSomaChi_inv = 1.0 / denSomaChi;
		numericType denDifChi_inv = 1.0 / denDifChi;
		numericType denSomaChi_inv2 = sqr(denSomaChi_inv);
		numericType denDifChi_inv2 = sqr(denDifChi_inv);

		numericType R_chi_ = 1 - chi_Half * (somaCos2 * denSomaChi_inv
				+ difCos2 * denDifChi_inv);
		epsilonLinha = R_chi_;

		numericType dRi, dRj, dRg;
		dRi = -chi_Half * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_i * denSomaChi + chi_
						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(2
								* (cos_gamma_menos_theta_i
										- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_i * denDifChi - chi_
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRj = -chi_Half * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_j * denSomaChi - chi_
						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(-2 * (cos_gamma_menos_theta_i
								- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_j * denDifChi + chi_
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRg = chi_ * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
				* denSomaChi_inv + //
				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
						* denDifChi_inv);
		dEpsilonLinha_i = dRi;
		dEpsilonLinha_j = dRj;
		dEpsilonLinha_g = dRg;
	}
	numericType epsilon2;
	numericType dEpsilon2_i, dEpsilon2_j, dEpsilon2_g;
	{
		numericType denInv = 1.0 / (1.0 - sqr(chi * cos_theta_i_menos_j));
		epsilon2 = sqrt(denInv);

		dEpsilon2_g = 0.0;
		dEpsilon2_i = -epsilon2 * sqr(chi) * cos_theta_i_menos_j
				* sin_theta_i_menos_j * denInv;
		dEpsilon2_j = -dEpsilon2_i;
	}
	numericType epsilon;
	numericType dEpsilon_i, dEpsilon_j, dEpsilon_g;
	{
		numericType epsilon2_nu_menos_1 = pow(epsilon2, nu - 1);
		numericType epsilonLinha_nu_menos_1 = pow(epsilonLinha, mu - 1);
		numericType prodComum = epsilon0 * epsilon2_nu_menos_1
				* epsilonLinha_nu_menos_1;

		epsilon = prodComum * epsilon2 * epsilonLinha;

		dEpsilon_i = prodComum * (nu * epsilonLinha * dEpsilon2_i + mu
				* epsilon2 * dEpsilonLinha_i);
		dEpsilon_j = prodComum * (nu * epsilonLinha * dEpsilon2_j + mu
				* epsilon2 * dEpsilonLinha_j);
		dEpsilon_g = prodComum * (nu * epsilonLinha * dEpsilon2_g + mu
				* epsilon2 * dEpsilonLinha_g);
	}
	numericType fra, fra6;
	fra = sigma0 / (dij - sigma + sigma0);
	fra6 = sqr(cube(fra));

	numericType dUi, dUj, dUr, dUg, dU_den;

	*U = 4 * epsilon * fra6 * (fra6 - 1.0);

	dU_den = 48.0 * epsilon / sigma0 * fra * fra6 * (fra6 - 0.5);

	dUr = -dU_den;
	dUi = 4 * fra6 * (fra6 - 1.0) * dEpsilon_i + dU_den * dSigma_i;
	dUj = 4 * fra6 * (fra6 - 1.0) * dEpsilon_j + dU_den * dSigma_j;
	dUg = 4 * fra6 * (fra6 - 1.0) * dEpsilon_g + dU_den * dSigma_g;

	numericType f_x = -dUg * eij.y / dij + dUr * eij.x;
	numericType f_y = +dUg * eij.x / dij + dUr * eij.y;

	*forca_i = newVetor(f_x, f_y);
	*torque_i = -dUi;
	*profPoco = -epsilon;
}

inline void interacaoGayBerneV_WallH(numericType x_i,
		const numericType theta_i, const numericType chi,
		const numericType chi_, const numericType sigma0,
		const numericType epsilon0, const double nu, const double mu,
		numericType *U, vetor *forca_i, numericType *torque_i,
		numericType *profPoco, const numericType wRight,
		const numericType wLeft) {

	typedef enum _dir {
		DIREITA = 62897, ESQUERDA = -1439807
	} dir;
	dir pos;
	{
		numericType r_y_1 = x_i - wLeft;
		numericType r_y_2 = x_i - wRight;
		if (fabs(r_y_1) > fabs(r_y_2)) {
			pos = DIREITA;
		} else {
			pos = ESQUERDA;
		}
	}
	numericType theta_j = -theta_i;

	numericType chiHalf = 0.5 * chi;
	numericType chi_Half = 0.5 * chi_;
	numericType dij;
	vetor rij, eij;
	numericType gamma = 0;

	//numericType f_ = 0.06;
	numericType cos_theta_i_menos_j = cos(theta_i - theta_j);
	numericType cos_gamma_menos_theta_i = cos(gamma - theta_i);
	numericType cos_gamma_menos_theta_j = cos(gamma - theta_j);

	numericType sin_theta_i_menos_j, sin_gamma_menos_theta_i,
			sin_gamma_menos_theta_j; // Inicialização adiada

	numericType somaCos2 = sqr(
			cos_gamma_menos_theta_i + cos_gamma_menos_theta_j);
	numericType difCos2 =
			sqr(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j);

	numericType sigma;
	numericType dSigma_i, dSigma_j, dSigma_g;

	{

		numericType denSomaChi = 1 + chi * cos_theta_i_menos_j;
		numericType denDifChi = 1 - chi * cos_theta_i_menos_j;
		numericType denSomaChi_inv = 1.0 / denSomaChi;
		numericType denDifChi_inv = 1.0 / denDifChi;
		numericType R_chi = 1 - chiHalf * (somaCos2 * denSomaChi_inv + difCos2
				* denDifChi_inv);
		numericType denSomaChi_inv2, denDifChi_inv2; // inicialização adiada;
		sigma = sigma0 / sqrt(R_chi);

		rij = sub(
				newVetor(
						(pos == DIREITA) ? (wRight + 0.5 * sigma) : (wLeft
								- 0.5 * sigma), 0), newVetor(x_i, 0));
		dij = length(rij);
		eij = divEscalar(dij, rij);

		numericType distCorte = sigma + sigma0 * EQUILIBRIO_UTIL;
		if (dij > distCorte) {
			*U = 0;

			*forca_i = newVetor(0.0, 0.0);

			*torque_i = 0.0;

			*profPoco = 0.0;
			return;
		}

		denSomaChi_inv2 = sqr(denSomaChi_inv);
		denDifChi_inv2 = sqr(denDifChi_inv);

		sin_theta_i_menos_j = sin(theta_i - theta_j);
		sin_gamma_menos_theta_i = sin(gamma - theta_i);
		sin_gamma_menos_theta_j = sin(gamma - theta_j);

		numericType dRi, dRj, dRg;
		dRi = -chiHalf * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_i * denSomaChi + chi * somaCos2
						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(2
								* (cos_gamma_menos_theta_i
										- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_i * denDifChi - chi
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRj = -chiHalf * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_j * denSomaChi - chi * somaCos2
						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(-2 * (cos_gamma_menos_theta_i
								- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_j * denDifChi + chi
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRg = chi * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
				* denSomaChi_inv + //
				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
						* denDifChi_inv);

		dSigma_i = -0.5 * (sigma / R_chi) * dRi;
		dSigma_j = -0.5 * (sigma / R_chi) * dRj;
		dSigma_g = -0.5 * (sigma / R_chi) * dRg;
	}
	numericType epsilonLinha;
	numericType dEpsilonLinha_i, dEpsilonLinha_j, dEpsilonLinha_g;
	{

		numericType denSomaChi = 1 + chi_ * cos_theta_i_menos_j;
		numericType denDifChi = 1 - chi_ * cos_theta_i_menos_j;
		numericType denSomaChi_inv = 1.0 / denSomaChi;
		numericType denDifChi_inv = 1.0 / denDifChi;
		numericType denSomaChi_inv2 = sqr(denSomaChi_inv);
		numericType denDifChi_inv2 = sqr(denDifChi_inv);

		numericType R_chi_ = 1 - chi_Half * (somaCos2 * denSomaChi_inv
				+ difCos2 * denDifChi_inv);
		epsilonLinha = R_chi_;

		numericType dRi, dRj, dRg;
		dRi = -chi_Half * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_i * denSomaChi + chi_
						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(2
								* (cos_gamma_menos_theta_i
										- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_i * denDifChi - chi_
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRj = -chi_Half * ( //
				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
						* sin_gamma_menos_theta_j * denSomaChi - chi_
						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
						(-2 * (cos_gamma_menos_theta_i
								- cos_gamma_menos_theta_j)
								* sin_gamma_menos_theta_j * denDifChi + chi_
								* difCos2 * sin_theta_i_menos_j)
								* denDifChi_inv2);
		dRg = chi_ * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
				* denSomaChi_inv + //
				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
						* denDifChi_inv);
		dEpsilonLinha_i = dRi;
		dEpsilonLinha_j = dRj;
		dEpsilonLinha_g = dRg;
	}
	numericType epsilon2;
	numericType dEpsilon2_i, dEpsilon2_j, dEpsilon2_g;
	{
		numericType denInv = 1.0 / (1.0 - sqr(chi * cos_theta_i_menos_j));
		epsilon2 = sqrt(denInv);

		dEpsilon2_g = 0.0;
		dEpsilon2_i = -epsilon2 * sqr(chi) * cos_theta_i_menos_j
				* sin_theta_i_menos_j * denInv;
		dEpsilon2_j = -dEpsilon2_i;
	}
	numericType epsilon;
	numericType dEpsilon_i, dEpsilon_j, dEpsilon_g;
	{
		numericType epsilon2_nu_menos_1 = pow(epsilon2, nu - 1);
		numericType epsilonLinha_nu_menos_1 = pow(epsilonLinha, mu - 1);
		numericType prodComum = epsilon0 * epsilon2_nu_menos_1
				* epsilonLinha_nu_menos_1;

		epsilon = prodComum * epsilon2 * epsilonLinha;

		dEpsilon_i = prodComum * (nu * epsilonLinha * dEpsilon2_i + mu
				* epsilon2 * dEpsilonLinha_i);
		dEpsilon_j = prodComum * (nu * epsilonLinha * dEpsilon2_j + mu
				* epsilon2 * dEpsilonLinha_j);
		dEpsilon_g = prodComum * (nu * epsilonLinha * dEpsilon2_g + mu
				* epsilon2 * dEpsilonLinha_g);
	}
	numericType fra, fra6;
	fra = sigma0 / (dij - sigma + sigma0);
	fra6 = sqr(cube(fra));

	numericType dUi, dUj, dUr, dUg, dU_den;

	*U = 4 * epsilon * fra6 * (fra6 - 1.0);

	dU_den = 48.0 * epsilon / sigma0 * fra * fra6 * (fra6 - 0.5);

	dUr = -dU_den;
	dUi = 4 * fra6 * (fra6 - 1.0) * dEpsilon_i + dU_den * dSigma_i;
	dUj = 4 * fra6 * (fra6 - 1.0) * dEpsilon_j + dU_den * dSigma_j;
	dUg = 4 * fra6 * (fra6 - 1.0) * dEpsilon_g + dU_den * dSigma_g;

	numericType f_x = -dUg * eij.y / dij + dUr * eij.x;
	numericType f_y = +dUg * eij.x / dij + dUr * eij.y;

	*forca_i = newVetor(f_x, f_y);
	*torque_i = -dUi;
	*profPoco = -epsilon;
}

//
//inline void interacaoNova(const vetor rij, const numericType theta_i,
//		const numericType theta_j, const numericType chi,
//		const numericType chi_, const numericType sigma0,
//		const numericType epsilon0, const double nu, const double mu,
//		numericType *U, vetor *forca_i, vetor *forca_j, numericType *torque_i,
//		numericType *torque_j) {
//
//	numericType chiHalf = 0.5 * chi;
//	numericType chi_Half = 0.5 * chi_;
//	numericType dij = length(rij);
//	vetor eij = divEscalar(dij, rij);
//	numericType gamma = atan2(rij.y, rij.x);
//
//	//numericType f_ = 0.06;
//	numericType cos_theta_i_menos_j = cos(theta_i - theta_j);
//	numericType cos_gamma_menos_theta_i = cos(gamma - theta_i);
//	numericType cos_gamma_menos_theta_j = cos(gamma - theta_j);
//
//	numericType sin_theta_i_menos_j = sin(theta_i - theta_j);
//	numericType sin_gamma_menos_theta_i = sin(gamma - theta_i);
//	numericType sin_gamma_menos_theta_j = sin(gamma - theta_j);
//
//	numericType somaCos2 = sqr(
//			cos_gamma_menos_theta_i + cos_gamma_menos_theta_j);
//	numericType difCos2 =
//			sqr(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j);
//
//	numericType sigma;
//	numericType dSigma_i, dSigma_j, dSigma_g;
//	{
//		numericType denSomaChi = 1 + chi * cos_theta_i_menos_j;
//		numericType denDifChi = 1 - chi * cos_theta_i_menos_j;
//		numericType denSomaChi_inv = 1.0 / denSomaChi;
//		numericType denDifChi_inv = 1.0 / denDifChi;
//		numericType denSomaChi_inv2 = sqr(denSomaChi_inv);
//		numericType denDifChi_inv2 = sqr(denDifChi_inv);
//		numericType R_chi = 1 - chiHalf * (somaCos2 * denSomaChi_inv + difCos2
//				* denDifChi_inv);
//		sigma = sigma0 / sqrt(R_chi);
//		//		if ((dij - sigma) > (4.5 * sigma0)) {
//		//			continue;
//		//		}
//		numericType dRi, dRj, dRg;
//		dRi = -chiHalf * ( //
//				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//						* sin_gamma_menos_theta_i * denSomaChi + chi * somaCos2
//						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
//						(2
//								* (cos_gamma_menos_theta_i
//										- cos_gamma_menos_theta_j)
//								* sin_gamma_menos_theta_i * denDifChi - chi
//								* difCos2 * sin_theta_i_menos_j)
//								* denDifChi_inv2);
//		dRj = -chiHalf * ( //
//				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//						* sin_gamma_menos_theta_j * denSomaChi - chi * somaCos2
//						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
//						(-2 * (cos_gamma_menos_theta_i
//								- cos_gamma_menos_theta_j)
//								* sin_gamma_menos_theta_j * denDifChi + chi
//								* difCos2 * sin_theta_i_menos_j)
//								* denDifChi_inv2);
//		dRg = chi * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
//				* denSomaChi_inv + //
//				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
//						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
//						* denDifChi_inv);
//
//		dSigma_i = -0.5 * (sigma / R_chi) * dRi;
//		dSigma_j = -0.5 * (sigma / R_chi) * dRj;
//		dSigma_g = -0.5 * (sigma / R_chi) * dRg;
//	}
//
//	numericType dist = dij - sigma;
//	if (dist < 0.0) {
//
//		numericType kappa = 1.0E-2;
//		numericType d_new = (dij + 0.01 /*+ sigma - sigma0*/);
//
//		numericType dij_menos1 = 1.0 / d_new;
//		numericType e_kappa_r = 100.0 * exp(-kappa * d_new);
//
//		*U = e_kappa_r * dij_menos1;
//		numericType f = (1.0 + kappa * d_new) * e_kappa_r * dij_menos1
//				* dij_menos1;
//		double theta = 2E2;
//		f += theta * uniformDev(&seedLong);
//		//		printf("dij = %Lf, sigma = %Lf, sigma0 = %Lf, f=%LF\n", dij, sigma,
//		//				sigma0, f);
//
//		*forca_i = multEscalar(-f, eij);
//		*forca_j = multEscalar(f, eij);
//
//		//		*U = 0;
//		//
//		//		*forca_i = newVetor(0.0, 0.0);
//		//		*forca_j = newVetor(0.0, 0.0);
//		//
//		*torque_i = 0.0;
//		*torque_j = 0.0;
//		return;
//	}
//	if (dist == 0.0) {
//		*U = 0;
//
//		*forca_i = newVetor(0.0, 0.0);
//		*forca_j = newVetor(0.0, 0.0);
//
//		*torque_i = 0.0;
//		*torque_j = 0.0;
//		return;
//	}
//
//	numericType epsilonLinha;
//	numericType dEpsilonLinha_i, dEpsilonLinha_j, dEpsilonLinha_g;
//	{
//
//		numericType denSomaChi = 1 + chi_ * cos_theta_i_menos_j;
//		numericType denDifChi = 1 - chi_ * cos_theta_i_menos_j;
//		numericType denSomaChi_inv = 1.0 / denSomaChi;
//		numericType denDifChi_inv = 1.0 / denDifChi;
//		numericType denSomaChi_inv2 = sqr(denSomaChi_inv);
//		numericType denDifChi_inv2 = sqr(denDifChi_inv);
//
//		numericType R_chi_ = 1 - chi_Half * (somaCos2 * denSomaChi_inv
//				+ difCos2 * denDifChi_inv);
//		epsilonLinha = R_chi_;
//
//		numericType dRi, dRj, dRg;
//		dRi = -chi_Half * ( //
//				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//						* sin_gamma_menos_theta_i * denSomaChi + chi_
//						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
//						(2
//								* (cos_gamma_menos_theta_i
//										- cos_gamma_menos_theta_j)
//								* sin_gamma_menos_theta_i * denDifChi - chi_
//								* difCos2 * sin_theta_i_menos_j)
//								* denDifChi_inv2);
//		dRj = -chi_Half * ( //
//				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//						* sin_gamma_menos_theta_j * denSomaChi - chi_
//						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
//						(-2 * (cos_gamma_menos_theta_i
//								- cos_gamma_menos_theta_j)
//								* sin_gamma_menos_theta_j * denDifChi + chi_
//								* difCos2 * sin_theta_i_menos_j)
//								* denDifChi_inv2);
//		dRg = chi_ * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
//				* denSomaChi_inv + //
//				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
//						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
//						* denDifChi_inv);
//		dEpsilonLinha_i = dRi;
//		dEpsilonLinha_j = dRj;
//		dEpsilonLinha_g = dRg;
//	}
//	numericType epsilon2;
//	numericType dEpsilon2_i, dEpsilon2_j, dEpsilon2_g;
//	{
//		numericType denInv = 1.0 / (1.0 - sqr(chi * cos_theta_i_menos_j));
//		epsilon2 = sqrt(denInv);
//
//		dEpsilon2_g = 0.0;
//		dEpsilon2_i = -epsilon2 * sqr(chi) * cos_theta_i_menos_j
//				* sin_theta_i_menos_j * denInv;
//		dEpsilon2_j = -dEpsilon2_i;
//	}
//	numericType epsilon;
//	numericType dEpsilon_i, dEpsilon_j, dEpsilon_g;
//	{
//		numericType epsilon2_nu_menos_1 = pow(epsilon2, nu - 1);
//		numericType epsilonLinha_nu_menos_1 = pow(epsilonLinha, mu - 1);
//		numericType prodComum = epsilon0 * epsilon2_nu_menos_1
//				* epsilonLinha_nu_menos_1;
//
//		epsilon = prodComum * epsilon2 * epsilonLinha;
//
//		dEpsilon_i = prodComum * (nu * epsilonLinha * dEpsilon2_i + mu
//				* epsilon2 * dEpsilonLinha_i);
//		dEpsilon_j = prodComum * (nu * epsilonLinha * dEpsilon2_j + mu
//				* epsilon2 * dEpsilonLinha_j);
//		dEpsilon_g = prodComum * (nu * epsilonLinha * dEpsilon2_g + mu
//				* epsilon2 * dEpsilonLinha_g);
//	}
//
//	numericType B = 2.0;
//
//	numericType dUi, dUj, dUr, dUg, dU_den;
//	numericType f = pow(dist, -B);
//
//	*U = epsilon * f;
//
//	dU_den = epsilon * B * pow(dist, -B - 1.0);
//
//	dUr = -dU_den;
//	dUi = dU_den * dSigma_i + f * dEpsilon_i;
//	dUj = dU_den * dSigma_j + f * dEpsilon_j;
//	dUg = dU_den * dSigma_g + f * dEpsilon_g;
//
//	numericType f_x = -dUg * eij.y / dij + dUr * eij.x;
//	numericType f_y = +dUg * eij.x / dij + dUr * eij.y;
//
//	*forca_i = newVetor(f_x, f_y);
//	*forca_j = newVetor(-f_x, -f_y);
//
//	*torque_i = -dUi;
//	*torque_j = -dUj;
//}
//
//inline void interacaoNova_2(const vetor rij, const numericType theta_i,
//		const numericType theta_j, const numericType chi,
//		const numericType chi_, const numericType sigma0,
//		const numericType epsilon0, const double nu, const double mu,
//		numericType *U, vetor *forca_i, vetor *forca_j, numericType *torque_i,
//		numericType *torque_j) {
//
//	numericType chiHalf = 0.5 * chi;
//	numericType chi_Half = 0.5 * chi_;
//	numericType dij = length(rij);
//	vetor eij = divEscalar(dij, rij);
//	numericType gamma = atan2(rij.y, rij.x);
//
//	//numericType f_ = 0.06;
//	numericType cos_theta_i_menos_j = cos(theta_i - theta_j);
//	numericType cos_gamma_menos_theta_i = cos(gamma - theta_i);
//	numericType cos_gamma_menos_theta_j = cos(gamma - theta_j);
//
//	numericType sin_theta_i_menos_j = sin(theta_i - theta_j);
//	numericType sin_gamma_menos_theta_i = sin(gamma - theta_i);
//	numericType sin_gamma_menos_theta_j = sin(gamma - theta_j);
//
//	numericType somaCos2 = sqr(
//			cos_gamma_menos_theta_i + cos_gamma_menos_theta_j);
//	numericType difCos2 =
//			sqr(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j);
//
//	numericType sigma;
//	numericType dSigma_i, dSigma_j, dSigma_g;
//	{
//		numericType denSomaChi = 1 + chi * cos_theta_i_menos_j;
//		numericType denDifChi = 1 - chi * cos_theta_i_menos_j;
//		numericType denSomaChi_inv = 1.0 / denSomaChi;
//		numericType denDifChi_inv = 1.0 / denDifChi;
//		numericType denSomaChi_inv2 = sqr(denSomaChi_inv);
//		numericType denDifChi_inv2 = sqr(denDifChi_inv);
//		numericType R_chi = 1 - chiHalf * (somaCos2 * denSomaChi_inv + difCos2
//				* denDifChi_inv);
//		sigma = sigma0 / sqrt(R_chi);
//		//		if ((dij - sigma) > (4.5 * sigma0)) {
//		//			continue;
//		//		}
//		numericType dRi, dRj, dRg;
//		dRi = -chiHalf * ( //
//				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//						* sin_gamma_menos_theta_i * denSomaChi + chi * somaCos2
//						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
//						(2
//								* (cos_gamma_menos_theta_i
//										- cos_gamma_menos_theta_j)
//								* sin_gamma_menos_theta_i * denDifChi - chi
//								* difCos2 * sin_theta_i_menos_j)
//								* denDifChi_inv2);
//		dRj = -chiHalf * ( //
//				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//						* sin_gamma_menos_theta_j * denSomaChi - chi * somaCos2
//						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
//						(-2 * (cos_gamma_menos_theta_i
//								- cos_gamma_menos_theta_j)
//								* sin_gamma_menos_theta_j * denDifChi + chi
//								* difCos2 * sin_theta_i_menos_j)
//								* denDifChi_inv2);
//		dRg = chi * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
//				* denSomaChi_inv + //
//				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
//						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
//						* denDifChi_inv);
//
//		dSigma_i = -0.5 * (sigma / R_chi) * dRi;
//		dSigma_j = -0.5 * (sigma / R_chi) * dRj;
//		dSigma_g = -0.5 * (sigma / R_chi) * dRg;
//	}
//	numericType epsilonLinha;
//	numericType dEpsilonLinha_i, dEpsilonLinha_j, dEpsilonLinha_g;
//	{
//
//		numericType denSomaChi = 1 + chi_ * cos_theta_i_menos_j;
//		numericType denDifChi = 1 - chi_ * cos_theta_i_menos_j;
//		numericType denSomaChi_inv = 1.0 / denSomaChi;
//		numericType denDifChi_inv = 1.0 / denDifChi;
//		numericType denSomaChi_inv2 = sqr(denSomaChi_inv);
//		numericType denDifChi_inv2 = sqr(denDifChi_inv);
//
//		numericType R_chi_ = 1 - chi_Half * (somaCos2 * denSomaChi_inv
//				+ difCos2 * denDifChi_inv);
//		epsilonLinha = R_chi_;
//
//		numericType dRi, dRj, dRg;
//		dRi = -chi_Half * ( //
//				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//						* sin_gamma_menos_theta_i * denSomaChi + chi_
//						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
//						(2
//								* (cos_gamma_menos_theta_i
//										- cos_gamma_menos_theta_j)
//								* sin_gamma_menos_theta_i * denDifChi - chi_
//								* difCos2 * sin_theta_i_menos_j)
//								* denDifChi_inv2);
//		dRj = -chi_Half * ( //
//				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//						* sin_gamma_menos_theta_j * denSomaChi - chi_
//						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
//						(-2 * (cos_gamma_menos_theta_i
//								- cos_gamma_menos_theta_j)
//								* sin_gamma_menos_theta_j * denDifChi + chi_
//								* difCos2 * sin_theta_i_menos_j)
//								* denDifChi_inv2);
//		dRg = chi_ * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
//				* denSomaChi_inv + //
//				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
//						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
//						* denDifChi_inv);
//		dEpsilonLinha_i = dRi;
//		dEpsilonLinha_j = dRj;
//		dEpsilonLinha_g = dRg;
//	}
//	numericType epsilon2;
//	numericType dEpsilon2_i, dEpsilon2_j, dEpsilon2_g;
//	{
//		numericType denInv = 1.0 / (1.0 - sqr(chi * cos_theta_i_menos_j));
//		epsilon2 = sqrt(denInv);
//
//		dEpsilon2_g = 0.0;
//		dEpsilon2_i = -epsilon2 * sqr(chi) * cos_theta_i_menos_j
//				* sin_theta_i_menos_j * denInv;
//		dEpsilon2_j = -dEpsilon2_i;
//	}
//	numericType epsilon;
//	numericType dEpsilon_i, dEpsilon_j, dEpsilon_g;
//	{
//		numericType epsilon2_nu_menos_1 = pow(epsilon2, nu - 1);
//		numericType epsilonLinha_nu_menos_1 = pow(epsilonLinha, mu - 1);
//		numericType prodComum = epsilon0 * epsilon2_nu_menos_1
//				* epsilonLinha_nu_menos_1;
//
//		epsilon = prodComum * epsilon2 * epsilonLinha;
//
//		dEpsilon_i = prodComum * (nu * epsilonLinha * dEpsilon2_i + mu
//				* epsilon2 * dEpsilonLinha_i);
//		dEpsilon_j = prodComum * (nu * epsilonLinha * dEpsilon2_j + mu
//				* epsilon2 * dEpsilonLinha_j);
//		dEpsilon_g = prodComum * (nu * epsilonLinha * dEpsilon2_g + mu
//				* epsilon2 * dEpsilonLinha_g);
//	}
//
//	numericType B = 2.0;
//	numericType dist = dij - sigma /*+ sigma0*/;
//	if (!(dist > 0)) {
//		*U = 0;
//
//		*forca_i = newVetor(0.0, 0.0);
//		*forca_j = newVetor(0.0, 0.0);
//
//		*torque_i = 0.0;
//		*torque_j = 0.0;
//		return;
//	}
//	numericType dUi, dUj, dUr, dUg, dU_den;
//	numericType f = pow(dist, -B);
//
//	*U = epsilon * f;
//
//	dU_den = epsilon * B * pow(dist, -B - 1.0);
//
//	dUr = -dU_den;
//	dUi = dU_den * dSigma_i + f * dEpsilon_i;
//	dUj = dU_den * dSigma_j + f * dEpsilon_j;
//	dUg = dU_den * dSigma_g + f * dEpsilon_g;
//
//	numericType f_x = -dUg * eij.y / dij + dUr * eij.x;
//	numericType f_y = +dUg * eij.x / dij + dUr * eij.y;
//
//	*forca_i = newVetor(f_x, f_y);
//	*forca_j = newVetor(-f_x, -f_y);
//
//	*torque_i = -dUi;
//	*torque_j = -dUj;
//}
//
//inline void interacaoGayBerneV_parede(const vetor rij, const numericType theta_i,
//		const numericType theta_j, const numericType chi,
//		const numericType chi_, const numericType sigma0,
//		const numericType epsilon0, const double nu, const double mu,
//		numericType *U, vetor *forca_i, vetor *forca_j, numericType *torque_i,
//		numericType *torque_j) {
//
//	numericType chiHalf = 0.5 * chi;
//	numericType chi_Half = 0.5 * chi_;
//	numericType dij = length(rij);
//	vetor eij = divEscalar(dij, rij);
//	numericType gamma = atan2(rij.y, rij.x);
//
//	//numericType f_ = 0.06;
//	numericType cos_theta_i_menos_j = cos(theta_i - theta_j);
//	numericType cos_gamma_menos_theta_i = cos(gamma - theta_i);
//	numericType cos_gamma_menos_theta_j = cos(gamma - theta_j);
//
//	numericType sin_theta_i_menos_j = sin(theta_i - theta_j);
//	numericType sin_gamma_menos_theta_i = sin(gamma - theta_i);
//	numericType sin_gamma_menos_theta_j = sin(gamma - theta_j);
//
//	numericType somaCos2 = sqr(
//			cos_gamma_menos_theta_i + cos_gamma_menos_theta_j);
//	numericType difCos2 =
//			sqr(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j);
//
//	numericType sigma;
//	numericType dSigma_i, dSigma_j, dSigma_g;
//	{
//		numericType denSomaChi = 1 + chi * cos_theta_i_menos_j;
//		numericType denDifChi = 1 - chi * cos_theta_i_menos_j;
//		numericType denSomaChi_inv = 1.0 / denSomaChi;
//		numericType denDifChi_inv = 1.0 / denDifChi;
//		numericType denSomaChi_inv2 = sqr(denSomaChi_inv);
//		numericType denDifChi_inv2 = sqr(denDifChi_inv);
//		numericType R_chi = 1 - chiHalf * (somaCos2 * denSomaChi_inv + difCos2
//				* denDifChi_inv);
//		sigma = sigma0 / sqrt(R_chi);
//		//		if ((dij - sigma) > (4.5 * sigma0)) {
//		//			continue;
//		//		}
//		numericType dRi, dRj, dRg;
//		dRi = -chiHalf * ( //
//				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//						* sin_gamma_menos_theta_i * denSomaChi + chi * somaCos2
//						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
//						(2
//								* (cos_gamma_menos_theta_i
//										- cos_gamma_menos_theta_j)
//								* sin_gamma_menos_theta_i * denDifChi - chi
//								* difCos2 * sin_theta_i_menos_j)
//								* denDifChi_inv2);
//		dRj = -chiHalf * ( //
//				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//						* sin_gamma_menos_theta_j * denSomaChi - chi * somaCos2
//						* sin_theta_i_menos_j) * denSomaChi_inv2 + //
//						(-2 * (cos_gamma_menos_theta_i
//								- cos_gamma_menos_theta_j)
//								* sin_gamma_menos_theta_j * denDifChi + chi
//								* difCos2 * sin_theta_i_menos_j)
//								* denDifChi_inv2);
//		dRg = chi * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
//				* denSomaChi_inv + //
//				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
//						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
//						* denDifChi_inv);
//
//		dSigma_i = -0.5 * (sigma / R_chi) * dRi;
//		dSigma_j = -0.5 * (sigma / R_chi) * dRj;
//		dSigma_g = -0.5 * (sigma / R_chi) * dRg;
//	}
//	numericType epsilonLinha;
//	numericType dEpsilonLinha_i, dEpsilonLinha_j, dEpsilonLinha_g;
//	{
//
//		numericType denSomaChi = 1 + chi_ * cos_theta_i_menos_j;
//		numericType denDifChi = 1 - chi_ * cos_theta_i_menos_j;
//		numericType denSomaChi_inv = 1.0 / denSomaChi;
//		numericType denDifChi_inv = 1.0 / denDifChi;
//		numericType denSomaChi_inv2 = sqr(denSomaChi_inv);
//		numericType denDifChi_inv2 = sqr(denDifChi_inv);
//
//		numericType R_chi_ = 1 - chi_Half * (somaCos2 * denSomaChi_inv
//				+ difCos2 * denDifChi_inv);
//		epsilonLinha = R_chi_;
//
//		numericType dRi, dRj, dRg;
//		dRi = -chi_Half * ( //
//				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//						* sin_gamma_menos_theta_i * denSomaChi + chi_
//						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
//						(2
//								* (cos_gamma_menos_theta_i
//										- cos_gamma_menos_theta_j)
//								* sin_gamma_menos_theta_i * denDifChi - chi_
//								* difCos2 * sin_theta_i_menos_j)
//								* denDifChi_inv2);
//		dRj = -chi_Half * ( //
//				(2 * (cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//						* sin_gamma_menos_theta_j * denSomaChi - chi_
//						* somaCos2 * sin_theta_i_menos_j) * denSomaChi_inv2 + //
//						(-2 * (cos_gamma_menos_theta_i
//								- cos_gamma_menos_theta_j)
//								* sin_gamma_menos_theta_j * denDifChi + chi_
//								* difCos2 * sin_theta_i_menos_j)
//								* denDifChi_inv2);
//		dRg = chi_ * ((cos_gamma_menos_theta_i + cos_gamma_menos_theta_j)
//				* (sin_gamma_menos_theta_i + sin_gamma_menos_theta_j)
//				* denSomaChi_inv + //
//				(cos_gamma_menos_theta_i - cos_gamma_menos_theta_j)
//						* (sin_gamma_menos_theta_i - sin_gamma_menos_theta_j)
//						* denDifChi_inv);
//		dEpsilonLinha_i = dRi;
//		dEpsilonLinha_j = dRj;
//		dEpsilonLinha_g = dRg;
//	}
//	numericType epsilon2;
//	numericType dEpsilon2_i, dEpsilon2_j, dEpsilon2_g;
//	{
//		numericType denInv = 1.0 / (1.0 - sqr(chi * cos_theta_i_menos_j));
//		epsilon2 = sqrt(denInv);
//
//		dEpsilon2_g = 0.0;
//		dEpsilon2_i = -epsilon2 * sqr(chi) * cos_theta_i_menos_j
//				* sin_theta_i_menos_j * denInv;
//		dEpsilon2_j = -dEpsilon2_i;
//	}
//	numericType epsilon;
//	numericType dEpsilon_i, dEpsilon_j, dEpsilon_g;
//	{
//		numericType epsilon2_nu_menos_1 = pow(epsilon2, nu - 1);
//		numericType epsilonLinha_nu_menos_1 = pow(epsilonLinha, mu - 1);
//		numericType prodComum = epsilon0 * epsilon2_nu_menos_1
//				* epsilonLinha_nu_menos_1;
//
//		epsilon = prodComum * epsilon2 * epsilonLinha;
//
//		dEpsilon_i = prodComum * (nu * epsilonLinha * dEpsilon2_i + mu
//				* epsilon2 * dEpsilonLinha_i);
//		dEpsilon_j = prodComum * (nu * epsilonLinha * dEpsilon2_j + mu
//				* epsilon2 * dEpsilonLinha_j);
//		dEpsilon_g = prodComum * (nu * epsilonLinha * dEpsilon2_g + mu
//				* epsilon2 * dEpsilonLinha_g);
//	}
//	numericType fra, fra6;
//	fra = sigma0 / (dij - sigma + sigma0);
//	fra6 = sqr(cube(fra));
//
//	numericType dUi, dUj, dUr, dUg, dU_den;
//
//	*U = 4 * epsilon * fra6 * (fra6 - 1.0);
//
//	dU_den = 48.0 * epsilon / sigma0 * fra * fra6 * (fra6 - 0.5);
//
//	dUr = -dU_den;
//	dUi = 4 * fra6 * (fra6 - 1.0) * dEpsilon_i + dU_den * dSigma_i;
//	dUj = 4 * fra6 * (fra6 - 1.0) * dEpsilon_j + dU_den * dSigma_j;
//	dUg = 4 * fra6 * (fra6 - 1.0) * dEpsilon_g + dU_den * dSigma_g;
//
//	numericType f_x = -dUg * eij.y / dij + dUr * eij.x;
//	numericType f_y = +dUg * eij.x / dij + dUr * eij.y;
//
//	*forca_i = newVetor(f_x, f_y);
//	*forca_j = newVetor(-f_x, -f_y);
//
//	*torque_i = -dUi;
//	*torque_j = -dUj;
//}
