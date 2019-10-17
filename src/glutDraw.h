/*
 * glutDraw.h
 *
 *  Created on: 01/04/2012
 *      Author: igor
 */

#ifndef GLUTDRAW_H_
#define GLUTDRAW_H_
#include "gb_main.h"

double dx_PBC[8], dy_PBC[8];

void reinitGD();
void desenhaCirculo(double raio);
void desenhaParticula(const particle part);
void desenha();
void specialKeyboard(int key, int x, int y);
void keyboard(unsigned char key, int x, int y);

#endif /* GLUTDRAW_H_ */
