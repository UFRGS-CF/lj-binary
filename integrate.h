#ifndef INTEGRATE_H
#define INTEGRATE_H

void integrate(int key, const gsl_rng *rng, unsigned long int N, double
		box_length, double temp, double nu, double dt, struct particle
		p[N], double *pe, double *ke, double *etot);

#endif
