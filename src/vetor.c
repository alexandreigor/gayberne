#include "vetor.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

inline vetor newVetor(const numericType x, const numericType y) {
	vetor _vetor;
	_vetor.x = x;
	_vetor.y = y;

	return _vetor;
}

inline vetor newVetorByDirection(const numericType modulo, const vetor direcao) {
	numericType l = length(direcao);
	if (l == 0.0 || modulo == 1.0) {
		return direcao;
	} else {
		double fat = modulo / l;
		return multEscalar(fat, direcao);
	}
}

inline vetor add(const vetor vetor1, const vetor vetor2) {
	vetor vetorResp;
	vetorResp.x = vetor1.x + vetor2.x;
	vetorResp.y = vetor1.y + vetor2.y;

	return vetorResp;
}

inline void inc(vetor *vector, const vetor incremento) {
	vector->x += incremento.x;
	vector->y += incremento.y;
}

inline vetor sub(const vetor vetor1, const vetor vetor2) {
	vetor vetorResp;
	vetorResp.x = vetor1.x - vetor2.x;
	vetorResp.y = vetor1.y - vetor2.y;

	return vetorResp;
}

inline void dec(vetor *vector, const vetor decremento) {
	vector->x -= decremento.x;
	vector->y -= decremento.y;
}

inline numericType sqrLength(const vetor _vetor) {
	return prodEscalar(_vetor, _vetor);
}

inline numericType length(const vetor _vetor) {
	return sqrt(prodEscalar(_vetor, _vetor));
}

inline numericType prodVetorial(const vetor vetor1, const vetor vetor2) {
	return (vetor1.x * vetor2.y) - (vetor2.x * vetor1.y);
}

inline vetor multEscalar(const numericType k, const vetor _vetor) {
	vetor vetorResp;
	vetorResp.x = k * _vetor.x;
	vetorResp.y = k * _vetor.y;
	return vetorResp;
}
inline vetor divEscalar(const numericType k, const vetor _vetor) {
	vetor vetorResp;
	vetorResp.x = _vetor.x / k;
	vetorResp.y = _vetor.y / k;
	return vetorResp;
}

inline vetor rescale(const numericType k, vetor *vector) {
	vector->x *= k;
	vector->y *= k;
	return *vector;
}

inline numericType prodEscalar(const vetor vetor1, const vetor vetor2) {
	return vetor1.x * vetor2.x + vetor1.y * vetor2.y;
}

inline void setZero(vetor *vetor1) {
	vetor1->x = 0.0;
	vetor1->y = 0.0;
}

inline void vprint(char *nome, const vetor _vetor) {
	printf("%s=( %Lf , %Lf )\n", nome, _vetor.x, _vetor.y);
}

inline char* toString(const vetor _vetor) {
	char *str;
	str = (char*) malloc(sizeof(char[50]));
	sprintf(str, "(%03.5Lf, %03.5Lf)", _vetor.x, _vetor.y);
	return str;
}

