/*
 * gb_main.h
 *
 *  Created on: 31/03/2012
 *      Author: igor
 */

#ifndef GB_MAIN_H_
#define GB_MAIN_H_
#include "vetor.h"

typedef enum boolean {
	false = 0, true = 1
} boolean;
typedef char* String;
typedef enum wallType {
	WALLS_NONE = 'N',
	WALLS_ALL = 'A',
	WALLS_CORREDOR = 'C',
	WALLS_CORREDOR_RUGOSO = 'R'
} wallType;
typedef enum sppType {
	SPP_OBJETIVO = 'O', SPP_MOVIMENTO = 'M', SPP_NONE = 'N'
} sppType;

typedef struct _partic particle; // Tipo de dado que representa uma partícula
// Especificação da estrutura _partic
struct _partic {
	vetor r; // Posição
	vetor v; // Velocidade
	vetor a;
	numericType theta; // Orientação
	numericType omega; //
	numericType tau; //
};

numericType chi;
numericType chi_;
double sigma0;

#endif /* GB_MAIN_H_ */
