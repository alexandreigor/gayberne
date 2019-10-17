#include <math.h>
#include "paredes.h"
#include "vetor.h"
#include "math_comp.h"
#include "gayberne.h"
#include <stdio.h>

extern int N, N_2;
extern double Lx, Ly, wTop, wBot;
extern double Lx_h, Ly_h;
extern double Req, proporcao_extr_lado;
extern const double mu, nu;
static const boolean REFLETIR = false;
extern numericType *rug_px, *rug_pyc, *rug_pyb;
extern int rug_n;
extern numericType *rug_ang_seg_c, *rug_ang_seg_b;
extern numericType *rug_px_m, *rug_pyc_m, *rug_pyb_m;

inline void paredesCorredor(particle part, numericType *u, vetor *forca,
		numericType *ta, double sigma0, double epsilon0) {
	numericType dist = distanciaEquilibrioGBV(newVetor(0, 1), part.theta, 0,
			chi, sigma0);
	numericType r_y_1 = part.r.y - (wBot - 0.5 * sigma0);
	numericType r_y_2 = part.r.y - (wTop + 0.5 * sigma0);
	numericType r_y = (r_y_1 > r_y_2) ? r_y_2 : r_y_1;
	r_y = part.r.y - (wBot - sigma0);
	if (r_y < dist) {
		if (REFLETIR) {
			if (part.v.y < 0) {
				part.v.y = -part.v.y;
			}
		} else {
			numericType _tau_i = 0.0, _tau_j = 0.0;
			vetor _f_i, _f_j;
			interacaoGayBerneV(newVetor(0, -r_y), part.theta, 0, chi, chi_,
					sigma0, epsilon0, nu, mu, u, &_f_i, &_f_j, &_tau_i, &_tau_j);
			//	interacaoNova(newVetor(0, -r_y), part.theta, 0, chi, chi_, sigma0, 0.2, nu,
			//			mu, u, &_f_i, &_f_j, &_tau_i, &_tau_j);
			*forca = _f_i;
			*ta = _tau_i;
		}

		return;
	}
	r_y = part.r.y - (wTop + sigma0);
	if (r_y > -dist) {
		if (REFLETIR) {
			if (part.v.y > 0) {
				part.v.y = -part.v.y;
			}
		} else {
			numericType _tau_i = 0.0, _tau_j = 0.0;
			vetor _f_i, _f_j;
			interacaoGayBerneV(newVetor(0, -r_y), part.theta, 0, chi, chi_,
					sigma0, epsilon0, nu, mu, u, &_f_i, &_f_j, &_tau_i, &_tau_j);
			//			interacaoNova_2(newVetor(0, -r_y), part.theta, 0, chi, chi_,
			//					sigma0, 0.2, nu, mu, u, &_f_i, &_f_j, &_tau_i, &_tau_j);
			*forca = _f_i;
			*ta = _tau_i;
		}
		return;
	}
	*forca = newVetor(0.0, 0.0); //inc(&forca, _f_i);
	*ta = 0.0;
	*u = 0.0;
}

inline void paredesTodas(particle part, numericType *u, vetor *forca,
		numericType *ta, double sigma0, double epsilon0) {

	numericType profPoco = 0.0;
	vetor forcaV, forcaH;
	numericType uV, uH;
	numericType taV, taH;
	interacaoGayBerneV_WallV(part.r.y, part.theta, chi, chi_, sigma0, epsilon0,
			nu, mu, &uV, &forcaV, &taV, &profPoco, Ly, 0.0);
	interacaoGayBerneV_WallH(part.r.x, part.theta, chi, chi_, sigma0, epsilon0,
			nu, mu, &uH, &forcaH, &taH, &profPoco, Lx, 0.0);
	*forca = add(forcaV, forcaH);
	*ta = taV + taH;
	*u = uV + uH;
}

inline void paredeNenhuma(particle part, numericType *u, vetor *forca,
		numericType *ta, double sigma0, double epsilon0) {
	*forca = newVetor(0.0, 0.0); //inc(&forca, _f_i);
	*ta = 0.0;
	*u = 0.0;
}

inline void paredesRugosas(particle part, numericType *u, vetor *forca,
		numericType *ta, double sigma0, double epsilon0) {
	numericType _u = 0.0, _tau_i = 0.0, _tau_j = 0.0;
	*forca = newVetor(0.0, 0.0);
	*ta = 0.0;
	*u = 0.0;

	register int i, j;

	numericType menor_r_len = fabs(300 * Lx);
	numericType menor_ang = 0.0;
	vetor menor_r = newVetor(0.0, 0.0);

	for (i = 1; i < rug_n; i++) {
		int i_menos_1 = i - 1;
		numericType pxi = rug_px[i], pxo = rug_px[i_menos_1];
		numericType pxm = rug_px_m[i_menos_1];//= 0.5 * (pxi + pxo);
		numericType dx = pxm - part.r.x;
		if (dx > Lx_h) {
			pxi -= Lx_h;
			pxo -= Lx_h;
		} else if (dx < -Lx_h) {
			pxi += Lx_h;
			pxo += Lx_h;
		}

		numericType py_i[2] = { rug_pyc[i], rug_pyb[i] };
		numericType py_o[2] = { rug_pyc[i_menos_1], rug_pyb[i_menos_1] };
		numericType ang[2] = { rug_ang_seg_c[i_menos_1],
				rug_ang_seg_b[i_menos_1] };
		for (j = 0; j < 2; j++) {
			numericType pyi = py_i[j], pyo = py_o[j];
			numericType angSegmento = ang[j];//atan2(pyi - pyo, pxi - pxo);
			numericType pxNear, pyNear;
			distanciaDePontoAteReta(part.r.x, part.r.y, pxo, pyo, pxi, pyi,
					true, true, &pxNear, &pyNear);
			vetor r = newVetor(pxNear - part.r.x, pyNear - part.r.y);
			numericType r_wall = length(r);

			if (r_wall < menor_r_len) {
				menor_r_len = r_wall;
				menor_ang = angSegmento;
				menor_r = r;
			}

		}
	}

	{
		vetor _f_i, _f_j;
		numericType profPoco = 0.0;
		interacaoGayBerneV_Rep(menor_r, part.theta, menor_ang, chi, chi_,
				sigma0, epsilon0, nu, mu, &_u, &_f_i, &_f_j, &_tau_i, &_tau_j,
				&profPoco, 1);
		inc(forca, _f_i);
		// TODO: Ignorando força "rotacional" por causa de giro
		// incontrolado quando o sistema trava e uma parede continua
		// dando torque a uma partícula
		//			*ta += _tau_i;
		*u += _u;
	}
}

inline void paredesCorredor2(particle part, numericType *u, vetor *forca,
		numericType *ta, double sigma0, double epsilon0) {
	numericType profPoco = 0.0;
	interacaoGayBerneV_WallV(part.r.y, part.theta, chi, chi_, sigma0, epsilon0,
			nu, mu, u, forca, ta, &profPoco, wTop, wBot);
	*u -= profPoco;
}
