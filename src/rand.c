#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**
 * Get a random float in [0, 1].
 */
float rand_f(void) {
	return ((float) rand() / RAND_MAX);
}

/**
 * Get a random number following normal distribution
 */
float normal(const int mean, const int s) {
	float w = 2, x1, x2, y1, y2;
	while (w >= 1) {
		x1 = (2 * ((double) rand() / ((double) (RAND_MAX)))) - 1;
		x2 = (2 * ((double) rand() / ((double) (RAND_MAX)))) - 1;
		w = x1 * x1 + x2 * x2;
	}

	w = sqrt((-2 * log(w)) / w);
	y1 = (x1 * w * s) + mean;
	y2 = (x2 * w * s) + mean;
	return y1;
}

//int main() {
//	int i = 0;
//	for (i = 0; i < 100; i++) {
//		printf("%d = %.2f\n", i, normal(0, 1));
//	}
//	return 0;
//}
