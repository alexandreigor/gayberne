#ifndef VETOR_H_
#define VETOR_H_

#define numericType long double

typedef struct _vetor vetor;

struct _vetor {
	numericType x;
	numericType y;
};

vetor newVetor(const numericType x, const numericType y);
vetor newVetorByDirection(const numericType modulo, const vetor direcao);
vetor add(const vetor vetor1, const vetor vetor2);
void inc(vetor *vector, const vetor incremento);
vetor sub(const vetor vetor1, const vetor vetor2);
void dec(vetor *vector, const vetor decremento);
vetor multEscalar(const numericType k, const vetor _vetor);
vetor divEscalar(const numericType k, const vetor _vetor);
vetor rescale(const numericType k, vetor *vector);
numericType prodEscalar(const vetor vetor1, const vetor vetor2);
numericType sqrLength(const vetor vetor);
numericType length(const vetor vetor);
numericType prodVetorial(const vetor vetor1, const vetor vetor2);
void setZero(vetor *_vetor);
void vprint(char *nome, const vetor _vetor);
char* toString(const vetor _vetor);

#endif /* VETOR_H_ */
