#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include "md.h"

void integrate(int key, const gsl_rng *rng, unsigned long int N, double temp,
		double nu, double dt, struct particle p[N], double *pe, double
		*ke, double *etot, double *inst_temp) {

	double sumv2;
	double m = 1; // particles mass, change this later
	double sigma = sqrt(temp);
	double randunif;

	if (key == 1) {
		for (unsigned long int i = 0; i < N; i++) {
			for (int q = 0; q < 3; q++) {
				p[i].r[q][0] = p[i].r[q][0] + dt * p[i].r[q][1] + dt * dt * p[i].r[q][2] / 2;
				p[i].r[q][1] = p[i].r[q][1] + dt * p[i].r[q][2] / 2;
			}
		}
	} else if (key == 2) {
		sumv2 = 0;
		for (unsigned long int i = 0; i < N; i++) {
			for (int q = 0; q < 3; q++) {
				p[i].r[q][1] = p[i].r[q][1] + dt * p[i].r[q][2] / 2;

				sumv2 += p[i].r[q][1] * p[i].r[q][1];
			}
		}

		*inst_temp = m * sumv2 / (3 * N);

		for (unsigned long int i = 0; i < N; i++) {
			randunif = gsl_ran_flat(rng, 0, 1);
			if (randunif < (nu * dt)) {
				for (int q = 0; q < 3; q++) {
					p[i].r[1][1] = gsl_ran_gaussian(rng, sigma);
				}
			}
		}

		*ke = 0.5 * sumv2;
		*etot = *pe + *ke;
	}
}

