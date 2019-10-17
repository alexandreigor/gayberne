#include <math.h>
#include "math_comp.h"
#include "ran1.h"
#include "gasdev.h"

inline double uniformDev(long *seed) {
	return ran1(seed);
}

inline double uniformDev2(double ave, double stdDev, long *semente) {
	static long double M_SQRT3 = 1.7320508075688772935274463415058723L;
	return (ran1(semente) - 0.5) * (2 * stdDev * M_SQRT3) + ave;
}

// Números aleatórios com distribuição normal
// com média ave e desvio padrão stdDev.
inline double gasDev(double ave, double stdDev, long *semente) {
	return stdDev * gasdev(semente) + ave;
}

inline long double sqr(long double x) {
	return x * x;
}

inline long double cube(long double x) {
	return x * x * x;
}

inline long double powInt(const long double x, const unsigned int y) {
	long double ret = x;
	if (y == 0) {
		return 1;
	}
	if (y == 1) {
		return x;
	}

	ret = sqr(powInt(x, y / 2));
	if (y % 2) {
		ret *= x;
	}

	return ret;
}

inline double rad2grau(double rad) {
	return FATOR_RADIANO_PARA_GRAU * rad;
}

inline double grau2rad(double grau) {
	return FATOR_GRAU_PARA_RADIANO * grau;
}

inline long double lineLength(const long double x, const long double y,
		const long double x0, const long double y0) {
	const long double dx = (x - x0), dy = (y - y0);
	return sqrt(dx * dx + dy * dy);
}

inline void pontoMaisProximo(const long double x, const long double y,
		const long double x0, const long double y0, const long double x1,
		const long double y1, const char segmento, long double *px,
		long double *py) {
	long double ox, oy;
	if (!(x1 - x0)) {
		ox = x0;
		oy = y;
	} else if (!(y1 - y0)) {
		ox = x;
		oy = y0;
	} else {
		double left, tg = -1 / ((y1 - y0) / (x1 - x0));
		left = (x1 * (x * tg - y + y0) + x0 * (x * -tg + y - y1)) / (tg * (x1
				- x0) + y0 - y1);
		ox = left;
		oy = tg * left - tg * x + y;
	}
	if (segmento && !(ox >= fmin(x0, x1) && ox <= fmax(x0, x1) && oy >= fmin(
			y0, y1) && oy <= fmax(y0, y1))) {
		long double l1 = lineLength(x, y, x0, y0), l2 =
				lineLength(x, y, x1, y1);
		if (l1 > l2) {
			*px = x1;
			*py = y1;
		} else {
			*px = x0;
			*py = y0;
		}
	} else {
		*px = ox;
		*py = oy;
	}
}

inline long double distanciaDePontoAteReta(const long double x,
		const long double y, const long double x0, const long double y0,
		const long double x1, const long double y1, const char segmento,
		const char retornarPontoMaisProximo, long double *px, long double *py) {
	if (retornarPontoMaisProximo || segmento) {
		long double ox, oy;
		if (!(x1 - x0)) {
			ox = x0;
			oy = y;
		} else if (!(y1 - y0)) {
			ox = x;
			oy = y0;
		} else {
			double left, tg = -1 / ((y1 - y0) / (x1 - x0));
			left = (x1 * (x * tg - y + y0) + x0 * (x * -tg + y - y1)) / (tg
					* (x1 - x0) + y0 - y1);
			ox = left;
			oy = tg * left - tg * x + y;
		}
		*px = ox;
		*py = oy;
		if (segmento && !(ox >= fmin(x0, x1) && ox <= fmax(x0, x1) && oy
				>= fmin(y0, y1) && oy <= fmax(y0, y1))) {
			long double l1 = lineLength(x, y, x0, y0), l2 = lineLength(x, y,
					x1, y1);
			if (l1 > l2) {
				*px = x1;
				*py = y1;
				return l2;
			} else {
				*px = x0;
				*py = y0;
				return l1;
			}
		} else {
			return lineLength(x, y, ox, oy);
		}
	}
	double a = y0 - y1, b = x1 - x0, c = x0 * y1 - y0 * x1;
	return fabs(a * x + b * y + c) / sqrt(a * a + b * b);
}

