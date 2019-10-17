/*
 * glutDraw.c
 *
 *  Created on: 01/04/2012
 *      Author: igor
 */
#include <stdlib.h>
#include <math.h>
#include <GL/freeglut.h>
#include "glutDraw.h"
#include "gb_main.h"
#include "math_comp.h"

extern int N, N_2;
extern double Lx, Ly, wTop, wBot;
extern double Req, proporcao_extr_lado;
extern const char USE_GLUT;
extern const wallType paredes;
extern particle* p;
extern const double epsilon0, mu, nu;
static boolean init = false;
extern numericType *rug_px, *rug_pyc, *rug_pyb;
extern int rug_n;

inline void reinitGD() {
	init = false;
}

inline void desenhaCirculo(double raio) {
	static int numPassosDesenhoCirculo = 32;
	static double *cs = NULL, *ss = NULL;

	register int i;
	if (cs == NULL) {
		cs = (double*) malloc(numPassosDesenhoCirculo * sizeof(double));
		ss = (double*) malloc(numPassosDesenhoCirculo * sizeof(double));
		for (i = 0; i < numPassosDesenhoCirculo; i++) {
			double ang = i * M_2PI / numPassosDesenhoCirculo;
			cs[i] = cos(ang);
			ss[i] = sin(ang);
		}
	}
	glPushMatrix();
	glBegin(GL_POLYGON);
	{
		for (i = 0; i < numPassosDesenhoCirculo; ++i) {
			glVertex2d(raio * cs[i], raio * ss[i]);
		}
	}
	glEnd();
	glPopMatrix();
}

inline void desenhaParticula(const particle part) {
	static char j0, jMax;
	if (!init) {
		j0
				= (paredes == WALLS_CORREDOR || paredes
						== WALLS_CORREDOR_RUGOSO) ? 3 : 0;
		jMax
				= (paredes == WALLS_CORREDOR || paredes
						== WALLS_CORREDOR_RUGOSO) ? 5
						: ((paredes == WALLS_NONE) ? 8 : 0);
		init = true;
	}

	double angGrausOrientacao = rad2grau(part.theta);
	register int j;
	glPushMatrix();
	{
		glTranslated(part.r.x, part.r.y, 0.0);
		glPushMatrix();
		glRotated(angGrausOrientacao, 0, 0, 1);
		glScaled(proporcao_extr_lado, 1, 1);
		desenhaCirculo(Req / 2.0);
		glPopMatrix();
	}
	glPopMatrix();

	for (j = j0; j < jMax; j++) {
		glPushMatrix();
		{
			glTranslated(part.r.x + dx_PBC[j], part.r.y + dy_PBC[j], 0.0);
			glRotated(angGrausOrientacao, 0, 0, 1);
			glScaled(proporcao_extr_lado, 1, 1);
			desenhaCirculo(Req / 2.0);
		}
		glPopMatrix();
	}
}

inline void desenha() {
	static float cor1[] = { 0.5, 0.5, 0.5 };
	static float cor2[] = { 1.0, 1.0, 1.0 };
	static float corLinha[] = { 0.0, 0.0, 0.0 };
	static float corCorredor[] = { 0.8, 0.8, 0.8 };

	if (!USE_GLUT) {
		return;
	}
	register int i;

	glClear(GL_COLOR_BUFFER_BIT);

	glPointSize(3.0);
	glPushMatrix();

	glLineWidth(3.0);
	glColor3fv(corLinha);
	if (paredes == WALLS_CORREDOR) {
		glColor3fv(corCorredor);
		glBegin(GL_POLYGON);
		{
			glVertex2d(0.0, wBot);
			glVertex2d(Lx, wBot);
			glVertex2d(Lx, wTop);
			glVertex2d(0.0, wTop);
		}
		glEnd();
		glColor3fv(corLinha);
		glBegin(GL_LINE);
		{
			glVertex2d(0.0, wBot);
			glVertex2d(Lx, wBot);
			glVertex2d(0.0, wTop);
			glVertex2d(Lx, wTop);
		}
		glEnd();
	} else if (paredes == WALLS_ALL) {
		glLineWidth(15);
		glBegin(GL_LINE_LOOP);
		{
			glVertex2d(0.0, 0.0);
			glVertex2d(0.0, Ly);
			glVertex2d(Lx, Ly);
			glVertex2d(Lx, 0.0);
		}
		glEnd();
	} else if (paredes == WALLS_CORREDOR_RUGOSO) {
		register int i;
		glBegin(GL_LINE_STRIP);
		{
			for (i = 0; i < rug_n; i++) {
				glVertex2d(rug_px[i], rug_pyc[i]);
			}
		}
		glEnd();
		glBegin(GL_LINE_STRIP);
		{
			for (i = 0; i < rug_n; i++) {
				glVertex2d(rug_px[i], rug_pyb[i]);
			}
		}
		glEnd();
	}
	glLineWidth(1);

	glColor3fv(cor1);
	for (i = 0; i < N_2; ++i) {
		desenhaParticula(p[i]);
	}
	glColor3fv(cor2);
	for (i = N_2; i < N; ++i) {
		desenhaParticula(p[i]);
	}
	glPopMatrix();
	glFlush();
}
