#include <math.h>
#include "md.h"

/*
 * Returns distance between particles i and j in each component
 */
void dist(struct particle i, struct particle j, double box, double s[3]) {
	for (int q = 0; q < 3; q++) {
		// distance between i and j
		s[q] = i.r[q][0] - j.r[q][0];

		// apply periodic boundary condition
		s[q] = s[q] - box * round(s[q] / box);
	}
}
