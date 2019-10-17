#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>
#include <GL/freeglut.h>
#include <sys/stat.h>
#include "gb_main.h"
#include "vetor.h"					// Definição e rotinas de manipulação do tipo vetor
#include "math_comp.h"
#include "gayberne.h"
#include "gasdev.h"
#include "paredes.h"
#include "glutDraw.h"
#include "SO_specific_func.h"

#define SEED_LONG_PAD -982217308L

static funcaoInteracaoParede interacaoComParedes;
const boolean MED_CONSERV = false;
const double MED_CONSERV_LIM = 5E3;

boolean PERFECT_RUN = true;
boolean USE_GLUT = true;
boolean TESTE_INVERSAO_VEL = false;

// ======= TODO: Parâmetros principais aqui
long seedLong0 = SEED_LONG_PAD - 96778698798;
double dt_scale = 1E-3;
int N = 32;
double densid = 0.31415;
double proporcao_extr_lado = 3.0;
double noiseStr = 0.0;
double tc = 0.2;
double t_relax_ang = 0.2;
// ======= FIM: Parâmetros principais

wallType paredes = WALLS_CORREDOR;
sppType usedSPP = SPP_OBJETIVO;
potencial_range LIMITE_POTENCIAL = POTENCIAL_EQUIL_RANGE;
boolean SPP_NOS_ANGULOS = false;
boolean TRUNCAR_RUIDO = true;

const boolean PLOT_ON_SCREEN = false;
const String base = "runs/gb_";

boolean measure = true;
double measureTimeStart = 4E3;
double measureTimeEnd = 5E3;
boolean FAZER_ARQUIVO = true;

static boolean loadFile = false;
static String fileToLoadName = NULL;

long seedLong;
double L = 32;
double L_h = 10;
double Lx_h, Ly_h;
double v0 = 1.0;
double tau = 1.0;
double tau_inv;

double dt, dt2, dt_2, dt2_2;

double noiseStrSqrt, noiseMaxAmplitude, noiseMaxAmplitude_2;

double t = 0.0;

double fe_fs = 0.3;

int N_2;
particle* p;

double Lx = 20;
double Ly = 10;
double wTop = 7.5;
double wBot = 2.50;

boolean running = false;
boolean RUN_CONT = true;
double Req = 1.0;

numericType I, I_inv;
numericType vMax;

double *e;
FILE *outputFile;

const double epsilon0 = 3.0;
const double mu = 2.0;
const double nu = 1.0;

static unsigned long refreshTime = 500, sysOutTime = 500, step = 0;
static double e_accum = 0, e_accumG = 0, e_measure = 0;

static String nome_arq_saida;
static String nome_arq_saida_old;
static unsigned int plot_num_cols;
static unsigned int* plot_cols;
static String* plot_titles;
// ------Randon gen.------------
static int iset;
static float gset;
static long iy;
static long *iv;
//------------------------------
//-------Rugosidades------------
int rug_k = 5;
int rug_n;
double rug_b_f = 0.1;
double rug_h = 2.0;
numericType *rug_px, *rug_pyc, *rug_pyb;
numericType *rug_ang_seg_c, *rug_ang_seg_b;
numericType *rug_px_m, *rug_pyc_m, *rug_pyb_m;
boolean rug_inverteCima = false, rug_inverteBaixo = true;
//------------------------------

inline String trimString(String string);
static inline void singleStep();
static inline void finish();
static inline void interacoes(const int rug_n, const particle partic[rug_n],
		vetor forca[rug_n], numericType ta[rug_n], numericType *U);
static inline void preinit();
static inline void init();
static inline void reduzirDistancias(vetor *r);
static inline void initFileParams(char mode);
static inline void plotFile(const String inputFileName,
		const unsigned int col_x, const unsigned int rug_n,
		const unsigned int col_y[rug_n], const String title[rug_n]);
void setDt(double _dt);
void setNoiseStr(const double _noiseStr);
void setN(int _N);
void setDimensoes(const double _Lx, const double _Ly);
void setReq(const double _Req);
void setProporcao(const double _proporcao_extr_lado);
void setRugosidade(const unsigned int rug_k, const double rug_b_f,
		const double rug_h, const boolean inverteCima,
		const boolean inverteBaixo);
void setPerfect(const boolean _perfect);
void processarEntrada(int argc, char **argv);
inline void print_usage(String exeName, unsigned int num_options,
		const String options[num_options],
		const unsigned char num_option_params[num_options],
		const char *types[num_options], const String *paramsDefinition);
static inline void deallocAll();
static inline String getSavFileName(const long step);
static void pontosRugosidade();
inline void carregarDados(String nomeArq);
inline void salvarDados(const String nomeArq, boolean comBackupDoArquivo);
inline boolean fileExists(String nomeArq);
void avaliaRandom(long nQtde, float media, float desvio, long *semente,
		int numClasses, const char distrib);
inline vetor generateNoise(long *seed);

inline String newString(unsigned int size) {
	return (String) calloc(size, sizeof(char));
}

inline String getTimeStamp() {
	String timestamp = newString(strlen("yyyy-MM-dd_hhmmss"));
	time_t t = time(0);
	struct tm * now = localtime(&t);
	sprintf(timestamp, "%04i-%02i-%02i_%02i%02i%02i", 1900 + (now->tm_year),
			now->tm_mon + 1, now->tm_mday, now->tm_hour, now->tm_min,
			now->tm_sec);
	return timestamp;
}

int main(int argc, String *argv) {
	atexit(deallocAll);

	if (false) {
		float med = 1.7777, std = 1.0;
		long n = 1000000;
		long seed = -111689900L;
		char distrib = 'U';
		int nClasses = 31;
		avaliaRandom(n, med, std, &seed, nClasses, distrib);
		if (distrib == 'U' || distrib == 'u') {
			double vl = std * sqrt(3.0);
			printf("\n\nUniform(min=%lf, max=%lf)\n", med - vl, med + vl);
		}
		n = 1000000;
		numericType x_ac = 0, y_ac = 0, xy_ac = 0, xx_ac = 0.0, yy_ac = 0.0;
		numericType xM, yM, xxM, xyM, yyM;
		std = 1;
		setNoiseStr(sqr(std));
		long contM = 0;
		FILE* f = fopen("/home/igor/Desktop/te.txt", "w");
		register long i;
		for (i = 0; i < n; i++) {
			vetor v;
			numericType len;
			do {
				v = generateNoise(&seed);
				len = sqrLength(v);
			} while (TRUNCAR_RUIDO && len > noiseMaxAmplitude_2);
			numericType le = length(v);
			if (le > (3 * std)) {
				contM++;
				//				v = newVetorByDirection(3 * std, v);
			}
			fprintf(f, "%Lf\t%Lf\n", v.x, v.y);
			x_ac += v.x;
			y_ac += v.y;
			xx_ac += v.x * v.x;
			xy_ac += v.x * v.y;
			yy_ac += v.y * v.y;

			if (!(i % 1000)) {
				int i_ = i + 1;
				xM = x_ac / i_;
				yM = y_ac / i_;
				xxM = xx_ac / i_;
				xyM = xy_ac / i_;
				yyM = yy_ac / i_;
				printf("At %.2f%%:\n", (100.0 * i_) / n);
				printf("cov_xx = %Lf\n", xxM - xM * xM);
				printf("cov_xy = %Lf\n", xyM - xM * yM);
				printf("cov_yy = %Lf\n", yyM - yM * yM);
			}
		}
		fclose(f);

		xM = x_ac / n;
		yM = y_ac / n;
		xxM = xx_ac / n;
		xyM = xy_ac / n;
		yyM = yy_ac / n;
		printf("At %.2f%%:\n", 100.0);
		printf("cov_xx = %Lf\n", xxM - xM * xM);
		printf("cov_xy = %Lf\n", xyM - xM * yM);
		printf("cov_yy = %Lf\n", yyM - yM * yM);

		printf("contM = %.3f %%\n", 100.0 * contM / n);
		exit(0);
	}

	preinit();
	processarEntrada(argc, argv);
	init();
	if (false) {
		double max = M_2PI, d;
		setProporcao(3.0);
		setReq(1.0);
		for (d = 0.1; d <= max; d += 0.01) {
			//setProporcao(d);
			numericType x1 = distanciaEquilibrioGBV(newVetor(0, 1), 0, d, chi,
					sigma0);
			numericType x2 = distanciaContatoGBV(newVetor(0, 1), 0, d, chi,
					sigma0);
			printf("%f\t%Lf\t%Lf\n", d, x1, x2);

			//			numericType u, ti, tj;
			//			vetor fi, fj;
			//			interacaoGayBerneV(newVetor(d, 0.0), 0, M_PI_2, chi, chi_, sigma0,
			//					epsilon0, nu, mu, &u, &fi, &fj, &ti, &tj);
			//			printf("%f\t%Lf\n", d, u);

		}
		exit(0);

	}
	if (loadFile) {
		printf("Carregando dados de arquivo \'%s\'.\n", fileToLoadName);
		carregarDados(fileToLoadName);
	}
	printf("dt=%f, noiseStr=%f, densidade=%f, sigma0 = %lf\n", dt, noiseStr,
			densid, sigma0);
	initFileParams(PERFECT_RUN ? 'P' : 'R');
	printf("%s\n", nome_arq_saida);
	if (FAZER_ARQUIVO) {
		struct stat _st;
		if (stat("runs", &_st) != 0) {
			system("mkdir runs 2>/dev/null");
			printf("Diretório \"runs\" criado...\n");
		}

		printf("loadFile=%i\tnome_ant=%s\n", loadFile, nome_arq_saida_old);
		if (loadFile && strlen(nome_arq_saida_old)) {
			int resp = copyFile(nome_arq_saida_old, nome_arq_saida);
			printf("Cópia: %i\n", resp);
		}
		if (fileExists(nome_arq_saida)) {
			outputFile = fopen(nome_arq_saida, "a");
		} else {
			outputFile = fopen(nome_arq_saida, "w");
			fprintf(outputFile, "t\teff\tE\tU\tK\t<eff>\t<effG>\n");
		}
	}
	if (loadFile) {
		reinitGasDev(iset, gset);
		reinitRan1(iy, iv);
		reinitGD();
		free(iv);
	}
	printf("chi=%Le\tchi'=%Le\n", chi, chi_);
	if (USE_GLUT) {
		/* Não usar atexit() e glutCloseFunc simultaneamente. */
		glutInit(&argc, argv);
		glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
				GLUT_ACTION_CONTINUE_EXECUTION);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		{
			int res_x = glutGet(GLUT_SCREEN_WIDTH);
			int res_y = glutGet(GLUT_SCREEN_HEIGHT);
			int pad = 700;
			int x_size = pad, y_size = pad;
			if (Lx > Ly) {
				y_size = pad * (Ly / Lx);
			} else {
				x_size = pad * (Lx / Ly);
			}
			glutInitWindowSize(x_size, y_size);
			glutInitWindowPosition((res_x - x_size) / 2, (res_y - y_size) / 2);
		}
		glutCreateWindow("Potencial elíptico de Gay-Berne");

		glClearColor(1.0, 1.0, 1.0, 0.0);
		glMatrixMode(GL_PROJECTION);
		gluOrtho2D(0, Lx, 0, Ly);

		glutDisplayFunc(desenha);
		glutKeyboardFunc(keyboard);
		glutSpecialFunc(specialKeyboard);
		glutCloseFunc(finish);
		if (running) {
			glutIdleFunc(singleStep);
		} else {
			glutIdleFunc(NULL);
		}

		glutMainLoop();
		printf("Glut encerrado.\nAplicação fechando...");
	} else {
		atexit(finish);
		while (true) {
			singleStep();
		}
	}
	return 0;
}

void setDt(const double _dt) {
	dt = _dt;
	dt2 = sqr(dt);
	dt_2 = 0.5 * dt;
	dt2_2 = 0.5 * dt2;
}

void setNoiseStr(const double _noiseStr) {
	noiseStr = _noiseStr;
	noiseStrSqrt = sqrt(_noiseStr);
	noiseMaxAmplitude = 3.0 * noiseStrSqrt;
	noiseMaxAmplitude_2 = 9.0 * noiseStr;
	printf("noiseStr=%lf\t(noiseStrSqrt=%lf)\n", noiseStr, noiseStrSqrt);
}

void setN(const int _N) {
	register int i;
	N = _N;
	N_2 = N / 2;

	free(p);
	free(e);
	p = (particle*) calloc(N, sizeof(particle));
	e = (double*) calloc(N, sizeof(double));

	for (i = 0; i < N_2; i++) {
		e[i] = 1.0;
	}
	for (i = N_2; i < N; i++) {
		e[i] = -1.0;
	}

}

void setDimensoes(const double _Lx, const double _Ly) {
	Lx = _Lx;
	Ly = _Ly;
	Lx_h = 0.5 * Lx;
	Ly_h = 0.5 * Ly;
}

void setProporcao(const double _proporcao_extr_lado) {
	proporcao_extr_lado = _proporcao_extr_lado;
	chi = (sqr(proporcao_extr_lado) - 1) / (sqr(proporcao_extr_lado) + 1);
	fe_fs = (0.4 + 0.6 / proporcao_extr_lado) / proporcao_extr_lado;
	chi_ = (1 - pow(fe_fs, 1.0 / mu)) / (1 + pow(fe_fs, 1.0 / mu));
	I = (M_PI * 1 * proporcao_extr_lado * (1 + sqr(proporcao_extr_lado)) * sqr(
			sqr(Req))) / 4;
	I_inv = 1.0 / I;
}

void setReq(const double _Req) {
	Req = _Req;
	sigma0 = Req;
	setProporcao(proporcao_extr_lado);
}

inline void preinit() {

	seedLong = seedLong0;

	setPerfect(false);

	setDt(dt_scale);
	setN(N);
	setNoiseStr(noiseStr);
	setDimensoes(Lx, Ly);
	setRugosidade(rug_k, rug_b_f, rug_h, rug_inverteCima, rug_inverteBaixo);

	tau_inv = 1.0 / tau;
}

inline void init() {

	printf("N=%i\tN/2=%i\n", N, N_2);

	if (PERFECT_RUN) {
		setNoiseStr(0.0);
	}

	switch (paredes) {
	case WALLS_CORREDOR:
		interacaoComParedes = &paredesCorredor2;
		break;
	case WALLS_CORREDOR_RUGOSO:
		if (rug_h) {
			interacaoComParedes = &paredesRugosas;
		} else {
			//			interacaoComParedes = &paredesCorredor;
			interacaoComParedes = &paredesRugosas;
		}
		break;
	case WALLS_ALL:
		interacaoComParedes = &paredesTodas;
		break;
	case WALLS_NONE:
		interacaoComParedes = &paredeNenhuma;
		break;
	default:
		printf("Tipo de parede desconhecido!");
		exit(1);
		break;
	}

	dx_PBC[0] = dx_PBC[3] = dx_PBC[5] = -Lx;
	dx_PBC[1] = dx_PBC[6] = 0.0;
	dx_PBC[2] = dx_PBC[4] = dx_PBC[7] = Lx;

	dy_PBC[0] = dy_PBC[1] = dy_PBC[2] = -Ly;
	dy_PBC[3] = dy_PBC[4] = 0.0;
	dy_PBC[5] = dy_PBC[6] = dy_PBC[7] = Ly;

	double dy;
	double y0 = 0;
	double wBot_eff = wBot, wTop_eff = wTop;
	{
		double areaDisponivel;//
		if (paredes == WALLS_ALL || paredes == WALLS_NONE) {
			dy = Ly;
			y0 = 0.0;
		} else {
			if (paredes == WALLS_CORREDOR_RUGOSO) {
				wBot_eff += (rug_inverteBaixo) ? rug_h : 0.0;
				wTop_eff -= (rug_inverteCima) ? rug_h : 0.0;
				printf("wBot = %f, wTop = %f\n", wBot, wTop);
				printf("wBot_eff = %f, wTop_eff = %f\n", wBot_eff, wTop_eff);
			}
			dy = wTop_eff - wBot_eff - 2 * EQUILIBRIO_UTIL * sigma0;
			y0 = wBot_eff + EQUILIBRIO_UTIL * sigma0;

			printf("wBot = %f, wTop = %f\n", wBot, wTop);
			printf("wBot_eff = %f, wTop_eff = %f\n", wBot_eff, wTop_eff);
			printf("dy = %f, y0 = %f\n", dy, y0);
		}
		areaDisponivel = Lx * dy;
		double areaParticula = densid * areaDisponivel / ((double) N);
		Req = sqrt(4 * areaParticula / proporcao_extr_lado / M_PI);
		densid = (N * areaParticula) / areaDisponivel;
	}

	setReq(Req);
	setProporcao(proporcao_extr_lado);

	register int i, j;
	long numSort = 0;
	long limSort = 30 * N * N;
	int k = 0;
	while (k < 1) {
		numSort = 0;
		for (i = 0; i < N; i++) {
			boolean ok = false;
			while (!ok) {
				p[i].theta = M_2PI * uniformDev(&seedLong);
				p[i].omega = 0.0;
				p[i].tau = 0.0;
				numSort++;

				if (numSort > limSort) {
					break;
				}

				numericType contato_x = distanciaEquilibrioGBV(newVetor(1, 0),
						p[i].theta, -p[i].theta, chi, sigma0);
				numericType contato_y = distanciaEquilibrioGBV(newVetor(0, 1),
						p[i].theta, -p[i].theta, chi, sigma0);

				switch (paredes) {
				case WALLS_CORREDOR:
				case WALLS_CORREDOR_RUGOSO: {
					p[i].r.x = Lx * uniformDev(&seedLong);
					p[i].r.y = (dy - contato_y) * uniformDev(&seedLong) + y0
							+ contato_y / 2.0;
					break;
				}
				case WALLS_ALL: {
					p[i].r.x = (Lx - contato_x) * uniformDev(&seedLong)
							+ contato_x / 2;
					p[i].r.y = (Ly - contato_y) * uniformDev(&seedLong)
							+ contato_y / 2.0;
					break;
				}
				default: {
					p[i].r.x = Lx * uniformDev(&seedLong);
					p[i].r.y = Ly * uniformDev(&seedLong);
					break;
				}
				}

				ok = true;
				for (j = 0; j < i; j++) {
					vetor r = sub(p[i].r, p[j].r);
					reduzirDistancias(&r);
					if (overlapGBV(r, p[i].theta, p[j].theta, chi, sigma0)) {
						ok = false;
						break;
					}
				}
			}
			//double ang = M_2PI * uniformDev(&seedLong);
			p[i].v = newVetor(gasdev(&seedLong), gasdev(&seedLong));//multEscalar(v0, newVetor(cos(ang), sin(ang)));
			p[i].a = newVetor(0.0, 0.0);
		}
		k++;
	}
	if (numSort >= limSort) {
		double w = Lx;
		double d = dy;
		double l = (1 + EQUILIBRIO_UTIL) * sigma0;
		double k = (1 + EQUILIBRIO_UTIL) * sigma0 * proporcao_extr_lado;
		int nx, ny;
		printf("l=%lf k=%lf\n", l, k);
		printf("d=%lf w=%lf\n", d, w);

		nx = w / k;
		ny = d / l;

		printf("nx=%i, ny=%i\n", nx, ny);
		int b = 0;
		if (nx * ny >= N) {
			nx = ceil((double) N / (double) ny);
			for (i = 0; i < nx; i++) {
				double px = i * (w / nx) + (1 + EQUILIBRIO_UTIL) * sigma0 * 0.5;
				double py0 = y0;// + 2.0 * au * (1 + EQUILIBRIO_UTIL) * sigma0;
				for (j = 0; j < ny; j++) {
					double py = j * (d / ny) + py0;
					p[b].r = newVetor(px, py);
					p[b].v = newVetor(0.0, 0.0);
					p[b].a = newVetor(0.0, 0.0);
					p[b].theta = 0.0;
					p[b].omega = 0.0;
					p[b].tau = 0.0;
					if (b < N_2) {
						e[b] = 1.0;
					} else {
						e[b] = -1.0;
					}
					b++;
					if (b >= N) {
						break;
					}
				}
				if (b >= N) {
					break;
				}
			}
			for (i = 0; i < N; i++) {
				int ind = (int) (N * uniformDev(&seedLong));
				particle aux = p[i];
				p[i] = p[ind];
				p[ind] = aux;
			}
		} else {
			printf("Sorteio abortado e quadro falhou.");
			exit(1);
		}
	}
	vMax = 0.0;
	numericType vMax2 = sqr(vMax);
	for (i = 0; i < N; i++) {
		double v2 = sqrLength(p[i].v);
		if (v2 > vMax2) {
			vMax = sqrt(v2);
		}
	}
}

inline vetor generateNoise(long *seed) {
	double vx = gasDev(0.0, noiseStrSqrt, seed);
	double vy = gasDev(0.0, noiseStrSqrt, seed);
	return newVetor(vx, vy);
}

inline void interacoes(const int rug_n, const particle partic[rug_n],
		vetor forca[rug_n], numericType ta[rug_n], numericType *U) {
	register int i, j;

	*U = 0;
	double Uwall = 0;
	for (i = 0; i < rug_n; i++) {
		forca[i] = newVetor(0.0, 0.0);
		ta[i] = 0.0;
	}

	for (i = 0; i < rug_n; i++) {

		for (j = 0; j < i; j++) {
			vetor rij = sub(partic[j].r, partic[i].r);
			reduzirDistancias(&rij);
			numericType _u = 0.0, _tau_i = 0.0, _tau_j = 0.0;
			numericType profPoco = 0.0;
			vetor _f_i, _f_j;

			char teste = false;
			if (!teste) {
				//TODO: interação
				interacaoGayBerneV_Rep(rij, partic[i].theta, partic[j].theta,
						chi, chi_, sigma0, epsilon0, nu, mu, &_u, &_f_i, &_f_j,
						&_tau_i, &_tau_j, &profPoco, LIMITE_POTENCIAL);

			} else {
				numericType sigma = distanciaEquilibrioGBV(rij,
						partic[i].theta, partic[j].theta, chi, sigma0);
				long double dij = length(rij);
				vetor eij = divEscalar(dij, rij);
				long double dij_menos1 = 1.0 / dij;
				double theta = 2E2;
				if (dij < sigma) {
					double kappa = 1.0E-0;
					double e_kappa_r = exp(-kappa * dij);
					_u = e_kappa_r * dij_menos1;
					double f = (1.0 + kappa * dij) * e_kappa_r * dij_menos1
							* dij_menos1;

					f += theta * uniformDev(&seedLong);

					_f_i = multEscalar(-f, eij);
					_f_j = multEscalar(f, eij);

					//					printf("f=%f (%f)\n", f, (float) length(_f_i));
				} else {
					_f_j = _f_i = newVetor(0, 0);
					_u = 0;
				}
			}

			profPoco = 0.0;
			*U += (_u - profPoco);
			inc(&forca[i], _f_i);
			inc(&forca[j], _f_j);
			ta[i] += _tau_i;
			ta[j] += _tau_j;

			if (noiseStr) {
				vetor v;
				numericType len;
				do {
					v = generateNoise(&seedLong);
					len = sqrLength(v);
				} while (TRUNCAR_RUIDO && len > noiseMaxAmplitude_2);
				inc(&forca[i], v);
			}

		}
		numericType _u = 0.0;
		vetor _forca = newVetor(0.0, 0.0);
		numericType _ta = 0.0;
		// TODO: Testando as paredes diferentes
		if (!TESTE_INVERSAO_VEL) {
			interacaoComParedes(partic[i], &_u, &_forca, &_ta, sigma0, epsilon0);
		}

		inc(&forca[i], _forca);
		*ta += _ta;
		Uwall += _u;
	}
	for (i = 0; i < rug_n; i++) {
		ta[i] *= I_inv;
	}
	//	*U += Uwall;

}

inline double eficienciaAtual() {
	register int i;
	double summation = 0.0;
	for (i = 0; i < N; i++) {
		summation += p[i].v.x * e[i];
	}
	return (summation / N / v0);
}

inline numericType forcaReduzirAngulo(const particle p) {

	// Quadrante do ângulo
	int Q = (int) (p.theta / M_PI_2) + 1;

	double d_tau;
	switch (Q) {
	case 1:
		d_tau = -p.theta;
		break;
	case 4:
		d_tau = M_2PI - p.theta;
		break;
	default:
		d_tau = M_PI - p.theta;
		break;
	}
	return (+d_tau - cube(p.omega)) / t_relax_ang;
}

inline void gravarMedidas(double forma, double x, double y1, double y2,
		String comment) {
	struct stat st;
	if (stat("measure", &st) != 0) {
		system("mkdir measure 2>/dev/null");
		printf("Diretório \"measure\" criado...\n");
	}

	String nomeArquivo = newString(200);
	sprintf(nomeArquivo, "measure/gb_f%g.dat", forma);
	trimString(nomeArquivo);
	FILE *measureFile = fopen(nomeArquivo, "a+");
	fprintf(measureFile, "%f\t%f\t%f\t# %s\n", x, y1, y2, comment);
	fclose(measureFile);
}

inline void singleStep() {

	vetor f[N];
	numericType ta[N];

	register int i;

	step++;
	t += dt;

	for (i = 0; i < N; i++) {
		inc(&(p[i].r), add(multEscalar(dt, p[i].v), multEscalar(dt2_2, p[i].a)));
		p[i].theta += dt * p[i].omega + dt2_2 * p[i].tau;

		p[i].theta = fmod(p[i].theta, M_2PI);

		if (paredes == WALLS_NONE) {
			if (p[i].r.x < 0) {
				p[i].r.x += Lx;
			} else if (p[i].r.x > Lx) {
				p[i].r.x -= Lx;
			}
			if (p[i].r.y < 0) {
				p[i].r.y += Ly;
			} else if (p[i].r.y > Ly) {
				p[i].r.y -= Ly;
			}
		} else if (paredes == WALLS_ALL) {
			//			if (p[i].r.x < 0) {
			//				p[i].r.x = dL;
			//				p[i].v.x = -p[i].v.x;
			//			} else if (p[i].r.x > Lx) {
			//				p[i].r.x = Lx - dL;
			//				p[i].v.x = -p[i].v.x;
			//			}
			//			if (p[i].r.y < 0) {
			//				p[i].r.y = dL;
			//				p[i].v.y = -p[i].v.y;
			//			} else if (p[i].r.y > Ly) {
			//				p[i].r.y = Ly - dL;
			//				p[i].v.y = -p[i].v.y;
			//			}
		} else {
			if (p[i].r.x < 0) {
				p[i].r.x += Lx;
			} else if (p[i].r.x > Lx) {
				p[i].r.x -= Lx;
			}
		}
	}
	numericType U = 0.0, K = 0.0;
	interacoes(N, p, f, ta, &U);
	vMax = 0.0;
	for (i = 0; i < N; i++) {

		numericType vx_est = p[i].v.x + (f[i].x + p[i].a.x) * dt_2;
		numericType vy_est = p[i].v.y + (f[i].y + p[i].a.y) * dt_2;

		switch (usedSPP) {
		case SPP_OBJETIVO: {
			f[i].x += (v0 * e[i] - vx_est) / tc;
			f[i].y -= vy_est / tc;
			if (SPP_NOS_ANGULOS) {
				ta[i] += forcaReduzirAngulo(p[i]);
			}
			break;
		}
		case SPP_MOVIMENTO: {
			double V2 = vx_est * vx_est + vy_est * vy_est;
			f[i].x += ((1.0 - V2) * vx_est);
			f[i].y += ((1.0 - V2) * vy_est);
			break;
		}
		default:
			break;
		}

		inc(&(p[i].v), multEscalar(dt_2, add(p[i].a, f[i])));
		//TODO: teste novas paredes
		if (TESTE_INVERSAO_VEL) {
			switch (paredes) {
			case WALLS_CORREDOR: {
				numericType dist_lim = 0.5 * distanciaEquilibrioGBV(
						newVetor(0, 1), p[i].theta, -p[i].theta, chi, sigma0);//*/
				numericType limite_expandido = 1.1 * dist_lim;
				if (p[i].r.y > (wTop - limite_expandido)) {
					p[i].v.y = -fabs(p[i].v.y);
				} else if (p[i].r.y < (wBot + limite_expandido)) {
					p[i].v.y = fabs(p[i].v.y);
				}
				break;
			}
			case WALLS_ALL: {
				numericType dist_lim = 0.5 * distanciaEquilibrioGBV(
						newVetor(0, 1), p[i].theta, -p[i].theta, chi, sigma0);//*/
				numericType limite_expandido = 1.1 * dist_lim;
				if (p[i].r.y > (Ly - limite_expandido)) {
					p[i].v.y = -fabs(p[i].v.y);
				} else if (p[i].r.y < (0.0 + limite_expandido)) {
					p[i].v.y = fabs(p[i].v.y);
				}
				dist_lim = 0.5 * distanciaEquilibrioGBV(newVetor(1, 0),
						p[i].theta, -p[i].theta, chi, sigma0);//*/
				limite_expandido = 1.1 * dist_lim;
				if (p[i].r.x > (Lx - limite_expandido)) {
					p[i].v.x = -fabs(p[i].v.x);
				} else if (p[i].r.x < (0.0 + limite_expandido)) {
					p[i].v.x = fabs(p[i].v.x);
				}
				break;
			}
			default:
				break;
			}
		}
		p[i].omega += dt_2 * (p[i].tau + ta[i]);

		p[i].a = f[i];
		p[i].tau = ta[i];

		double v_ = length(p[i].v);
		K += 0.5 * sqr(v_) + 0.5 * I * sqr(p[i].omega);
		if (v_ > vMax) {
			vMax = v_;
		}

	}

	double x = eficienciaAtual();
	e_accum += x;
	if (step % refreshTime == 0 || !RUN_CONT) {
		desenha();
	}
	if (step % sysOutTime == 0 || !RUN_CONT) {
		numericType E = U + K;
		desenha();
		e_accumG += e_accum;
		e_accum /= sysOutTime;
		printf("%f\t%lf\t%Lf\t%Lf\t%Lf\t%lf\t%lf\n", t, x, E, U, K, e_accum,
				e_accumG / step);
		if (FAZER_ARQUIVO) {
			fprintf(outputFile, "%f\t%lf\t%Lf\t%Lf\t%Lf\t%lf\t%lf\n", t, x, E,
					U, K, e_accum, e_accumG / step);
		}
		e_accum = 0;
	}

	if (measure) {
		if (t > measureTimeStart) {
			e_measure += x;
		}
		if (t > measureTimeEnd) {
			double time_diff = (measureTimeEnd - measureTimeStart) / dt;
			e_measure /= time_diff;
			double e_accumG_meas = e_accumG / step;
			gravarMedidas(proporcao_extr_lado, noiseStr, e_measure,
					e_accumG_meas, nome_arq_saida);

			String name2 = getSavFileName(step);
			salvarDados(name2, true);
			free(name2);

			if (USE_GLUT) {
				glutDisplayFunc(NULL);
				glutLeaveMainLoop();
			} else {
				exit(0);
			}
			//finish();
		}

	}

	if (!RUN_CONT) {
		keyboard(' ', 0, 0);
	}

	if (MED_CONSERV && (t > MED_CONSERV_LIM)) {
		String name2 = getSavFileName(step);
		salvarDados(name2, true);
		free(name2);
		if (USE_GLUT) {
			glutDisplayFunc(NULL);
			glutLeaveMainLoop();
		} else {
			exit(0);
		}
	}

	//	printf("step=%li t=%f p[0].r=%s p[0].v=%s p[0].theta=%f\n", step, t,
	//			toString(p[0].r), toString(multEscalar(dt, p[0].v)), p[0].theta);

}

inline void reduzirDistancias(vetor *r) {
	switch (paredes) {
	case WALLS_NONE: {
		if (r->x > Lx_h) {
			r->x -= Lx;
		} else if (r->x < -Lx_h) {
			r->x += Lx;
		}
		if (r->y > Ly_h) {
			r->y -= Ly;
		} else if (r->y < -Ly_h) {
			r->y += Ly;
		}
		break;
	}
	case WALLS_CORREDOR:
	case WALLS_CORREDOR_RUGOSO: {
		if (r->x > Lx_h) {
			r->x -= Lx;
		} else if (r->x < -Lx_h) {
			r->x += Lx;
		}
		break;
	}
	default:
		break;
	}

}

inline String trimString(String string) {
	int n = strlen(string);
	//printf("trimString(\'%s\' /* len: %i */);\n", string, n);
	string = (String) realloc(string, (n + 1) * sizeof(char));
	return string;
}
// mode:
//	'P': perfect - medir energias
//	'R': run - medições ordinárias
inline void initFileParams(char mode) {
	String st = newString(300);
	String ts = getTimeStamp();
	switch (mode) {
	case 'P': {
		plot_num_cols = 3;
		plot_cols = (unsigned int*) calloc(plot_num_cols, sizeof(unsigned int));
		plot_titles = (String*) calloc(plot_num_cols, sizeof(String));

		plot_cols[0] = 3;
		plot_titles[0] = "E";

		plot_cols[1] = 4;
		plot_titles[1] = "U";

		plot_cols[2] = 5;
		plot_titles[2] = "K";

		sprintf(st, "%sP_dt%f_ts%s.dat", base, dt, ts);
		break;
	}
	case 'R': {
		plot_num_cols = 1;
		plot_cols = (unsigned int*) calloc(plot_num_cols, sizeof(unsigned int));
		plot_titles = (String*) calloc(plot_num_cols, sizeof(String));

		plot_cols[0] = 2;
		plot_titles[0] = "Eff";

		//		sprintf(st, "%sR_dt%g_N%i_r%g_d%g_f%g_s%li-%hu-%hu-%hu_ts%s.dat", base,
		//				dt, N, noiseStr, densid, proporcao_extr_lado, seedLong0,
		//				seed0[0], seed0[1], seed0[2], ts);
		sprintf(st, "%s%c_dt%g_N%i_r%g_d%g_f%g_s%li_ts%s.dat", base, paredes,
				dt, N, noiseStr, densid, proporcao_extr_lado, seedLong0, ts);
		break;
	}
	}
	st = trimString(st);
	nome_arq_saida = st;
}

static inline void deallocAll() {
	printf("Desalocando memória utilizada...\n");
	free(p);
	free(plot_cols);
	free(plot_titles);
	free(nome_arq_saida);
	free(nome_arq_saida_old);
	free(e);
	free(rug_px);
	free(rug_pyc);
	free(rug_pyb);
	free(rug_ang_seg_c);
	free(rug_ang_seg_b);
	free(rug_px_m);
	free(rug_pyc_m);
	free(rug_pyb_m);
}

inline void finish() {
	if (FAZER_ARQUIVO) {
		fclose(outputFile);
		plotFile(nome_arq_saida, 1, plot_num_cols, plot_cols, plot_titles);
	}

	printf("finish();\n");
	exit(0);
}

void specialKeyboard(int key, int x, int y) {
	static double fatorAumentoVel = 2;
	switch (key) {
	case GLUT_KEY_PAGE_DOWN: {
		if (noiseStr < 10) {
			setNoiseStr(0);
		} else {
			setNoiseStr(noiseStr - 10);
		}
		break;
	}
	case GLUT_KEY_PAGE_UP: {
		setNoiseStr(noiseStr + 10);
		break;
	}
	case GLUT_KEY_LEFT: {
		setDt(dt * fatorAumentoVel);
		printf("dt=%lf\n", dt);
		break;
	}
	case GLUT_KEY_RIGHT: {
		setDt(dt / fatorAumentoVel);
		printf("dt=%lf\n", dt);
		break;
	}
	default:
		break;
	}

}
void keyboard(unsigned char key, int x, int y) {
	switch (key) {
	case 27: {
		if (USE_GLUT) {
			glutIdleFunc(NULL);
		}
		finish();
		return;
	}
	case 'd': {
		RUN_CONT = !RUN_CONT;
		return;
	}
	case ' ': {
		if (running) {
			glutIdleFunc(NULL);
			running = false;
		} else {
			glutIdleFunc(singleStep);
			running = true;
		}
		return;
	}
	case '+': {
		setNoiseStr(noiseStr + 1);
		return;
	}
	case '-': {
		if (noiseStr < 1) {
			setNoiseStr(0);
		} else {
			setNoiseStr(noiseStr - 1);
		}
		return;
	}
	case 's': {
		String name2 = getSavFileName(step);
		salvarDados(name2, true);
		free(name2);
		return;
	}
	case 't': {
		TRUNCAR_RUIDO = !TRUNCAR_RUIDO;
		printf("TRUNCAR_RUIDO=%s\n", (TRUNCAR_RUIDO) ? "true" : "false");
		return;
	}
	}
	glutPostRedisplay();
}

static inline void plotFile(const String inputFileName,
		const unsigned int col_x, const unsigned int rug_n,
		const unsigned int col_y[rug_n], const String title[rug_n]) {
	FILE *gnuplotScript;
	String comando = newString(500);
	String nomeArquivoScript = "gnuplot.gnu";
	String comando2 = newString(9 + strlen(nomeArquivoScript)); //"rm -f '...'"
	register int i;
	int codSaida;

	gnuplotScript = fopen(nomeArquivoScript, "w");

	fprintf(gnuplotScript, "set encoding iso_8859_1\n");
	fprintf(gnuplotScript, "plot ");
	if (rug_n) {
		fprintf(gnuplotScript, " \\\n '%s' using %i:%i title '%s' with lines",
				inputFileName, col_x, col_y[0], title[0]);
	}
	for (i = 1; i < rug_n; i++) {
		fprintf(gnuplotScript, ", \\\n '%s' using %i:%i title '%s' with lines",
				inputFileName, col_x, col_y[i], title[i]);
	}
	fprintf(gnuplotScript, "\n");
	fprintf(gnuplotScript, "set terminal png size 800,600 \n");
	fprintf(gnuplotScript, "set output '%s.png'\n", inputFileName);
	fprintf(gnuplotScript, "replot\n");
	fclose(gnuplotScript);

	if (PLOT_ON_SCREEN) {
		sprintf(comando, "gnome-terminal -e \"gnuplot '%s' - \" & ",
				nomeArquivoScript);
	} else {
		sprintf(comando, "gnuplot '%s' 2>/dev/null", nomeArquivoScript);
	}
	comando = trimString(comando);
	sprintf(comando2, "rm -f '%s'", nomeArquivoScript);
	codSaida = system(comando);
	if (codSaida) {
		printf("Comando retornou mensagem de erro: %i\n%s\n", codSaida, comando);
	}
	//	system(comando2);

}

void processarEntrada(int argc, char **argv) {
	const unsigned int num_flags = 14;
	const String flags[] = { "-seed", "-load", "-nograph", "-noise", "-shape",
			"-dens", "-measure", "-N", "-repOnly", "-rug", "-perfect", "-help",
			"-dt", "-sppAngle" };
	const unsigned char nexts[] = { 1, 1, 0, 1, 1, 1, 2, 1, 1, 4, 0, 0, 1, 1 };
	const char *types[] = { "l", "s", "", "d", "d", "d", "dd", "l", "d",
			"ldds", "", "", "d", "l" };
	const String
			paramsDefinition[] =
					{
							"Semente (long) para o gerador gasdev",
							"Arquivo a ser carregado",
							"Intensidade do ruído",
							"Forma das partículas",
							"Densidade da distribuição das partículas",
							"Início do intervalo de medição",
							"Final do intervalo de medição",
							"Número de partículas",
							"Usar o potencial completo (0) ou "
								"apenas a parte repulsiva do mesmo com corte em sigma (1) ou com corte no equilíbrio (2)",
							"Número (inteiro) de rugosidades durante o comprimento do corredor, que determina \"l\" da rugosidade",
							"Tamanho (fracionado) b_f = (b/l) que determina \"b\" da rugosidade",
							"Profundidade da rugosidade \"h\"",
							"Algum dos valores (00, 01, 10, 11) representando se a rugosidade é para \"fora\" (0) "
								"ou para \"dentro\" (1) do corredor, em cima e em baixo, respectivamente",
							"Intervalo de tempo de evolução do sistema",
							"Indica se haverá (1) ou não (0) propulsão angular" };
	int i, j, k_;
	char perfect = 0;

	for (i = 1; i < argc; i++) {
		char *arg = argv[i];
		char **nextArg;

		for (j = 0; j < num_flags; j++) {
			if (!strcasecmp(arg, flags[j])) {
				break;
			}
		}

		if (j >= num_flags) {
			printf("ERRO: <%s> não é opção válida!\n", arg);
			print_usage(argv[0], num_flags, flags, nexts, types,
					paramsDefinition);
			exit(1);
		}

		int numArgs = nexts[j];
		nextArg = (String*) calloc(numArgs, sizeof(String));

		for (k_ = 0; k_ < numArgs; k_++) {
			int index = i + k_ + 1;
			if (index >= argc) {
				free(nextArg);
				printf("Argumentos insuficientes para opção <%s>!", flags[j]);
				print_usage(argv[0], num_flags, flags, nexts, types,
						paramsDefinition);
				exit(1);
			}
			nextArg[k_] = argv[index];
		}

		switch (j) {
		case 0: { // Semente do gasdev
			char *end = 0;
			long lida = strtol(nextArg[0], &end, 10);
			if (strlen(end)) {
				printf("Formato numérico inesperado em <%s %s>.\n", flags[j],
						nextArg[0]);
				exit(1);
			}
			printf("Seed lida: %li\n", lida);
			if (lida > 0) {
				lida *= -1;
			}
			seedLong = lida;
			seedLong0 = lida;
			break;
		}
		case 1: { // Arquivo a ser carregado
			loadFile = true;
			fileToLoadName = nextArg[0];
			break;
		}
		case 2: { // No-Graphics
			USE_GLUT = false;
			break;
		}
		case 3: { // Ruído
			char *end = 0;
			noiseStr = strtod(nextArg[0], &end);
			if (strlen(end)) {
				printf("Formato numérico inesperado em <%s %s>.\n", flags[j],
						nextArg[0]);
				exit(1);
			}
			setNoiseStr(noiseStr);
			break;
		}
		case 4: { // forma das partículas
			char *end = 0;
			proporcao_extr_lado = strtod(nextArg[0], &end);
			if (strlen(end)) {
				printf("Formato numérico inesperado em <%s %s>.\n", flags[j],
						nextArg[0]);
				exit(1);
			}
			break;
		}
		case 5: { // Densidade
			char *end = 0;
			densid = strtod(nextArg[0], &end);
			if (strlen(end)) {
				printf("Formato numérico inesperado em <%s %s>.\n", flags[j],
						nextArg[0]);
				exit(1);
			}
			break;
		}

		case 6: { // Intervalo de tempo de medição
			double ini, fim;
			char *end1, *end2;
			ini = strtod(nextArg[0], &end1);
			fim = strtod(nextArg[1], &end2);
			if (strlen(end1) + strlen(end2)) {
				printf("Formato numérico inesperado em <%s %s %s>.\n",
						flags[j], nextArg[0], nextArg[1]);
				exit(1);
			}
			if (fim <= ini) {
				printf(
						"Erro em <%s %s %s>: Fim (%s) é menor ou igual ao começo (%s)!\n",
						flags[j], nextArg[0], nextArg[1], nextArg[1],
						nextArg[0]);
				exit(1);
			}
			measure = true;
			measureTimeStart = ini;
			measureTimeEnd = fim;
			FAZER_ARQUIVO = true;
			break;
		}
		case 7: { // Número de partículas
			char *end = 0;
			N = strtol(nextArg[0], &end, 10);
			if (strlen(end)) {
				printf("Formato numérico inesperado em <%s %s>.\n", flags[j],
						nextArg[0]);
				exit(1);
			}
			setN(N);
			break;
		}
		case 8: { // Usar apenas a parte repulsiva do potencial
			char *end = 0;
			long val = strtol(nextArg[0], &end, 10);
			if (strlen(end) || (val != 0 && val != 1 && val != 2)) {
				printf(
						"Formato numérico inesperado em <%s %s>. Apenas '0', '1' e '2' são aceitos.\n",
						flags[j], nextArg[0]);
				exit(1);
			}
			LIMITE_POTENCIAL = val;
			break;
		}
		case 9: {
			char *end1 = 0, *end2 = 0, *end3 = 0;
			char cim = 'e', bai = 'e';
			boolean erroModos = false;
			long n_rug = strtol(nextArg[0], &end1, 10);
			double bf_rug = strtod(nextArg[1], &end2);
			double h_rug = strtod(nextArg[2], &end3);
			if (strlen(nextArg[3]) != 2) {
				erroModos = true;
			} else {
				cim = nextArg[3][0];
				bai = nextArg[3][1];
				erroModos = (cim != '0' && cim != '1') || (bai != '0' && bai
						!= '1');
			}
			if (strlen(end1) + strlen(end2) + strlen(end3)) {
				printf("Formato numérico inesperado em <%s %s %s %s %s>.\n",
						flags[j], nextArg[0], nextArg[1], nextArg[2],
						nextArg[3]);
				exit(1);
			}
			if (erroModos) {
				printf(
						"Valor \"%s\" inválido em <%s %s %s %s %s>. Apenas \"00\", \"01\",\"10\",\"11\" são aceitos.\n",
						nextArg[3], flags[j], nextArg[0], nextArg[1],
						nextArg[2], nextArg[3]);
				exit(1);
			}
			paredes = WALLS_CORREDOR_RUGOSO;
			setRugosidade(n_rug, bf_rug, h_rug, cim == '1', bai == '1');
			break;
		}
		case 10: {
			perfect = 1;
			break;
		}
		case 11: {
			print_usage(argv[0], num_flags, flags, nexts, types,
					paramsDefinition);
			exit(0);
		}
		case 12: { // dt
			char *end = 0;
			double dt = strtod(nextArg[0], &end);
			if (strlen(end)) {
				printf("Formato numérico inesperado em <%s %s>.\n", flags[j],
						nextArg[0]);
				exit(1);
			}
			setDt(dt);
			break;
		}
		case 13: {
			char *end = 0;
			int b = strtol(nextArg[0], &end, 10);
			if (strlen(end)) {
				printf("Formato numérico inesperado em <%s %s>.\n", flags[j],
						nextArg[0]);
				exit(1);
			}
			SPP_NOS_ANGULOS = !(b == 0);
			break;
		}
		default: {
			printf(
					"Não devia chegar aqui!\n\tprocessarEntrada(...) { j = %i; }\n",
					j);
			print_usage(argv[0], num_flags, flags, nexts, types,
					paramsDefinition);
			exit(1);
			break;
		}
		}
		i += numArgs;
		free(nextArg);
	}
	// TODO: trocar aqui
	setPerfect(false);
}

inline void print_usage(String exeName, unsigned int num_options,
		const String options[num_options],
		const unsigned char num_option_params[num_options],
		const char *types[num_options], const String *paramsDefinition) {
	int i;
	int argCount = 0;
	int totalArgs = 0;
	for (i = 0; i < num_options; i++) {
		totalArgs += num_option_params[i];
	}
	char *paramsDefs[totalArgs];
	printf("Usage:\n\t%s", exeName);
	for (i = 0; i < num_options; i++) {
		int j;
		unsigned int numArgs = num_option_params[i];
		printf("\n\t\t%s", options[i]);
		for (j = 0; j < numArgs; j++) {
			char type = types[i][j];
			String argStr = newString(20);
			switch (type) {
			case 'd': {
				sprintf(argStr, "<double_%i>", argCount);
				break;
			}
			case 'l': {
				sprintf(argStr, "<long_%i>", argCount);
				break;
			}
			case 's': {
				sprintf(argStr, "<string_%i>", argCount);
				break;
			}
			default: {
				sprintf(argStr, "<unk_%i>", argCount);
				break;
			}
			}
			argStr = trimString(argStr);
			printf(" %s", argStr);
			String argDef = newString(
					(strlen(argStr) + strlen(paramsDefinition[argCount]) + 3));
			sprintf(argDef, "%s:\t%s", argStr, paramsDefinition[argCount]);
			paramsDefs[argCount] = argDef;
			argCount++;
		}
	}
	printf("\n");
	if (totalArgs) {
		printf("Descrição dos parâmetros:\n");
	}
	for (i = 0; i < totalArgs; i++) {
		printf("\t%s\n", paramsDefs[i]);
	}

}

void setPerfect(const boolean _perfect) {
	PERFECT_RUN = _perfect;
	//	LIMITE_POTENCIAL = PERFECT_RUN ? POTENCIAL_SIGMA_RANGE
	//			: POTENCIAL_SIGMA_RANGE;
	if (PERFECT_RUN) {
		paredes = WALLS_NONE;
		usedSPP = SPP_NONE;
		SPP_NOS_ANGULOS = false;

		setNoiseStr(0.0);
	}
}

static inline String getSavFileName(const long step) {
	int len = strlen(nome_arq_saida);
	int size2 = len + 15;
	String name2 = newString(size2);
	strncpy(name2, nome_arq_saida, len - 4);
	sprintf(name2, "%s_%li.sav", name2, step);
	trimString(name2);
	return name2;
}

void setRugosidade(const unsigned int k, const double b_f, const double h,
		const boolean inverteCima, const boolean inverteBaixo) {
	rug_k = k;
	rug_n = 2 * k + 1;
	rug_b_f = b_f;
	rug_h = h;

	free(rug_px);
	free(rug_pyc);
	free(rug_pyb);
	free(rug_ang_seg_c);
	free(rug_ang_seg_b);
	free(rug_px_m);
	free(rug_pyc_m);
	free(rug_pyb_m);

	rug_px = (numericType*) calloc(rug_n, sizeof(numericType));
	rug_pyc = (numericType*) calloc(rug_n, sizeof(numericType));
	rug_pyb = (numericType*) calloc(rug_n, sizeof(numericType));
	rug_ang_seg_c = (numericType*) calloc(rug_n - 1, sizeof(numericType));
	rug_ang_seg_b = (numericType*) calloc(rug_n - 1, sizeof(numericType));
	rug_px_m = (numericType*) calloc(rug_n - 1, sizeof(numericType));
	rug_pyc_m = (numericType*) calloc(rug_n - 1, sizeof(numericType));
	rug_pyb_m = (numericType*) calloc(rug_n - 1, sizeof(numericType));

	rug_inverteCima = inverteCima;
	rug_inverteBaixo = inverteBaixo;

	pontosRugosidade();

	//	if (!h && paredes == WALLS_CORREDOR_RUGOSO) {
	//		printf("Pronfundidade das rugosidades = 0.0; "
	//			"Assumindo paredes lisas como otimização...\n");
	//		interacaoComParedes = &paredesCorredor;
	//	}
}

static void pontosRugosidade() {
	numericType l = Lx / rug_k;
	numericType b = l * rug_b_f;
	numericType r = l - b;

	numericType dyc = (rug_inverteCima) ? -rug_h : rug_h;
	numericType dyb = (rug_inverteBaixo) ? rug_h : -rug_h;

	int i;
	rug_px[0] = 0;
	rug_pyc[0] = wTop;
	rug_pyb[0] = wBot;
	for (i = 1; i < rug_n; i++) {
		int i_menos_1 = i - 1;
		if (i % 2) {
			rug_px[i] = rug_px[i_menos_1] + b;
			rug_pyc[i] = wTop + dyc;
			rug_pyb[i] = wBot + dyb;
		} else {
			rug_px[i] = rug_px[i_menos_1] + r;
			rug_pyc[i] = wTop;
			rug_pyb[i] = wBot;
		}
		rug_ang_seg_c[i_menos_1] = atan2(rug_pyc[i] - rug_pyc[i_menos_1],
				rug_px[i] - rug_px[i_menos_1]);
		rug_ang_seg_b[i_menos_1] = atan2(rug_pyb[i] - rug_pyb[i_menos_1],
				rug_px[i] - rug_px[i_menos_1]);
		rug_px_m[i_menos_1] = 0.5 * (rug_px[i] + rug_px[i_menos_1]);
		rug_pyc_m[i_menos_1] = 0.5 * (rug_pyc[i] + rug_pyc[i_menos_1]);
		rug_pyb_m[i_menos_1] = 0.5 * (rug_pyb[i] + rug_pyb[i_menos_1]);
	}
	return;
}

inline boolean fileExists(String nomeArq) {
	FILE *f;
	if ((f = fopen(nomeArq, "r"))) {
		fclose(f);
		return true;
	}
	return false;
}

void carregarDados(String nomeArq) {
	FILE *fileToLoad = fopen(nomeArq, "r");
	if (fileToLoad == NULL) {
		printf("Arquivo \'%s\' não existe!\n", nomeArq);
		exit(1);
	}
	double dt, noiseStr, Lx, Ly, Req, proporcao_extr_lado;
	int N;
	char par;
	int repu;
	register int i;
	double rug_h_local;
	fscanf(fileToLoad, "%i\n%le\n%le\n%le\n"
		"%li\n"
		"%li\t%i\t%e\n" //
			"%li", //
			&N, &noiseStr, &proporcao_extr_lado, &dt, //
			&seedLong0, //
			&seedLong, &iset, &gset,//
			&iy);
	int nTab = getNTab();
	iv = (long*) calloc(nTab, sizeof(long));
	for (i = 0; i < nTab; i++) {
		fscanf(fileToLoad, "\t%li", &(iv[i]));
	}
	fscanf(fileToLoad, "\n%le\n%le\n%le\n%le\n%i\n%c\n", &Lx, &Ly, &densid,
			&Req, &repu, &par);
	LIMITE_POTENCIAL = (char) repu;
	paredes = par;
	switch (par) {
	case WALLS_CORREDOR:
		interacaoComParedes = &paredesCorredor2;
		fscanf(fileToLoad, "%le\t%le\n", &wTop, &wBot);
		break;
	case WALLS_CORREDOR_RUGOSO:
		interacaoComParedes = &paredesRugosas;
		int invCima, invBaixo;
		fscanf(fileToLoad, "%le\t%le\t%i\t%le\t%le\t%i\t%i\n", &wTop, &wBot,
				&rug_k, &rug_b_f, &rug_h_local, &invCima, &invBaixo);
		rug_inverteCima = (char) invCima;
		rug_inverteBaixo = (char) invBaixo;
		break;
	case WALLS_ALL:
		interacaoComParedes = &paredesTodas;
		break;
	case WALLS_NONE:
		interacaoComParedes = &paredeNenhuma;
		break;
	default:
		printf("Não devia ter chegado aqui!");
		exit(1);
		break;
	}
	char spp;
	fscanf(fileToLoad, "%c\n", &spp);
	usedSPP = spp;
	switch (spp) {
	case SPP_OBJETIVO: {
		int sppAngulo;
		fscanf(fileToLoad, "%i\n", &sppAngulo);
		SPP_NOS_ANGULOS = sppAngulo;
		break;
	}
	case SPP_MOVIMENTO:
		break;
	case SPP_NONE:
		break;
	default:
		printf("Não devia ter chegado aqui!");
		exit(1);
		break;
	}
	nome_arq_saida_old = newString(500);
	int k = fscanf(fileToLoad, "$name$ %s\n", nome_arq_saida_old);
	nome_arq_saida_old = trimString(nome_arq_saida_old);
	printf("SPP_NOS_ANGULOS=%i\tk=%i\n", SPP_NOS_ANGULOS, k);
	fscanf(fileToLoad, "%le\t%lu\t%le\t%le\t%le\t%lu\t%lu\t%le\n", &t, &step,
			&e_accum, &e_accumG, &e_measure, &refreshTime, &sysOutTime,
			&t_relax_ang);
	printf("t_relax=%lf\tsysOutTime=%lu\n", t_relax_ang, sysOutTime);
	setDt(dt);
	setN(N);
	setNoiseStr(noiseStr);
	setDimensoes(Lx, Ly);
	setProporcao(proporcao_extr_lado);
	setReq(Req);
	if (par == WALLS_CORREDOR_RUGOSO) {
		setRugosidade(rug_k, rug_b_f, rug_h_local, rug_inverteCima,
				rug_inverteBaixo);
	}
	for (i = 0; i < N; i++) {
		fscanf(fileToLoad, "%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\t%Le\n",
				&p[i].r.x, &p[i].r.y, &p[i].v.x, &p[i].v.y, &p[i].a.x,
				&p[i].a.y, &p[i].theta, &p[i].omega, &p[i].tau);
	}
	fclose(fileToLoad);
	if (FAZER_ARQUIVO) {
		//fclose();
	}
}

void salvarDados(const String nomeArq, boolean comBackupDoArquivo) {
	String n2 = newString(strlen(nomeArq) + 3);
	sprintf(n2, "%s.1", nomeArq);
	FILE *fileToWrite = fopen(n2, "w");
	register int i;
	int nTab = getNTab();
	iv = (long*) calloc(nTab, sizeof(long));
	getRan1Status(&iy, iv);
	getGasDevStatus(&iset, &gset);
	fprintf(fileToWrite, "%i\n%.15le\n%.15le\n%.15le\n"
		"%li\n"
		"%li\t%i\t%.15e\n"//
			"%li", //
			N, noiseStr, proporcao_extr_lado, dt, //
			seedLong0, //
			seedLong, iset, gset,//
			iy);//
	for (i = 0; i < nTab; i++) {
		fprintf(fileToWrite, "\t%li", iv[i]);
	}
	free(iv);
	fprintf(fileToWrite, "\n%.15le\n%.15le\n%.15le\n%.15le\n%i\n%c\n", Lx, Ly,
			densid, Req, LIMITE_POTENCIAL, paredes);
	switch (paredes) {
	case WALLS_CORREDOR:
		fprintf(fileToWrite, "%.15le\t%.15le\n", wTop, wBot);
		break;
	case WALLS_CORREDOR_RUGOSO:
		fprintf(fileToWrite, "%.15le\t%.15le\t%i\t%.15le\t%.15le\t%i\t%i\n",
				wTop, wBot, rug_k, rug_b_f, rug_h, rug_inverteCima,
				rug_inverteBaixo);
		break;
	default:
		fprintf(fileToWrite, "\n");
		break;
	}

	fprintf(fileToWrite, "%c\n", usedSPP);
	switch (usedSPP) {
	case SPP_OBJETIVO:
		fprintf(fileToWrite, "%i\n", SPP_NOS_ANGULOS);
		break;
	default:
		fprintf(fileToWrite, "\n");
		break;
	}

	String nome;
	if (comBackupDoArquivo) {
		fflush(outputFile);
		String backupName = newString(strlen(nomeArq) + 5);
		sprintf(backupName, "%s.dat", nomeArq);
		nome = backupName;
		copyFile(nome_arq_saida, backupName);
	} else {
		nome = nome_arq_saida;
	}
	fprintf(fileToWrite, "$name$ %s\n", nome);
	fprintf(fileToWrite,
			"%.15le\t%lu\t%.15le\t%.15le\t%.15le\t%lu\t%lu\t%.15le\n", t, step,
			e_accum, e_accumG, e_measure, refreshTime, sysOutTime, t_relax_ang);

	for (i = 0; i < N; i++) {
		fprintf(
				fileToWrite,
				"%.19Le\t%.19Le\t%.19Le\t%.19Le\t%.19Le\t%.19Le\t%.19Le\t%.19Le\t%.19Le\n",
				p[i].r.x, p[i].r.y, p[i].v.x, p[i].v.y, p[i].a.x, p[i].a.y,
				p[i].theta, p[i].omega, p[i].tau);
	}
	fclose(fileToWrite);
	moveFile(n2, nomeArq);
	free(n2);
	if (comBackupDoArquivo) {
		free(nome);
	}
}

void avaliaRandom(long nQtde, float media, float desvio, long *semente,
		int numClasses, const char distrib) {
	register long i;
	float x[nQtde];
	float sum = 0.0;
	float min, max;
	min = media + 9 * desvio;
	max = media - 9 * desvio;
	for (i = 0; i < nQtde; i++) {
		switch (distrib) {
		case 'N':
		case 'n':
			x[i] = gasDev(media, desvio, semente);
			break;
		case 'U':
		case 'u':
			x[i] = uniformDev2(media, desvio, semente);
			break;
		default:
			printf("Distribuição desconhecida.");
			return;
		}

		sum += x[i];
		if (x[i] < min) {
			min = x[i];
		}
		if (x[i] > max) {
			max = x[i];
		}
	}

	float diff = max - min;
	float tam = diff / numClasses;
	long dist[numClasses];
	float mean = sum / nQtde;
	float var = 0.0;
	for (i = 0; i < numClasses; i++) {
		dist[i] = 0;
	}
	for (i = 0; i < nQtde; i++) {
		float sub = x[i] - mean;
		var += sub * sub;
		int index = (int) ((x[i] - min) / tam);
		if (index > numClasses) {
			printf("Erro\n");
			exit(0);
		}
		dist[index]++;
	}

	var /= nQtde;
	float dev = sqrt(var);
	printf("Média: %f\nDesvio: %f\nVariância: %f\n", mean, dev, var);
	for (i = 0; i < numClasses; i++) {
		printf("[%f,%f]\t%li\n", min + i * tam, min + (i + 1) * tam, dist[i]);
	}
}
