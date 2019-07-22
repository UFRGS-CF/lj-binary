/**
 * Molecular Dynamics in canonical (NVT) ensemble with Andersen thermostat and
 * Velocity Verlet integrator
 */

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "md.h"
#include "integrate.h"
#include "io.h"

/*
 * Returns smallest integer a such that a^3 > N
 */
int cubic(unsigned long int N) {
	int a = pow(N, 1.0/3);
	if (a*a*a == N)
		return a;
	else
		return a+1;
}

void lattice(unsigned long int N, double box, struct particle p[N]) {
	int max = cubic(N); // maximum box side index
	unsigned long int n = 0; // how many particles were placed in lattice

	for (int i = 0; i < max; i++) {
		for (int j = 0; j < max; j++) {
			for (int k = 0; k < max; k++) {
				if (n < N) {
					p[n].r[0][0] = box * ((float)k) / max;
					p[n].r[1][0] = box * ((float)j) / max;
					p[n].r[2][0] = box * ((float)i) / max;
					n++;
				}
			}
		}
	}

}

void init(unsigned long int N, const gsl_rng *rng, double temp, double box,
		struct particle p[N], double ca) {

	double sumv[3] = {0, 0, 0}, vcm[3] = {0, 0, 0};
	double sumv2 = 0, v2[3];
	double mv2, fs;

	// setup particles in cubic lattice
	lattice(N, box, p);

	// setup velocities
	for (int q = 0; q < 3; q++) {
		for (unsigned long int i = 0; i < N; i++) {
			p[i].r[q][1] = gsl_ran_flat(rng, -1, 1);

			sumv[q] += p[i].r[q][1];

			v2[q] = p[i].r[q][1] * p[i].r[q][1];
			sumv2 += v2[q];
		}

		vcm[q] = sumv[q] / N; // velocity center of mass
	}

	mv2 = sumv2 / N; // mean square velocity

	// match system KE with temp, remove net CM drift
	fs = sqrt(3 * temp / mv2); // calculate scale factor
	for (int q = 0; q < 3; q++) {
		for (unsigned long int i = 0; i < N; i++) {
			p[i].r[q][1] = (p[i].r[q][1] - vcm[q]) * fs;
		}
	}

	for (unsigned long int i = 0; i < N; i++) {
		if(gsl_ran_flat(rng, 0, 1) < ca)
			p[i].type = 0;
		else
			p[i].type = 1;
	}
}

void forces(unsigned long int N, double box, struct particle p[N], double rcut,
		double *pe, double *virial, struct interaction inter[2][2]) {

	double eps, sig;
	double r2, r6i, rcut2, ecut, ff;
	double s[3]; // separation between particles

	*pe = 0;
	*virial = 0;
	rcut2 = rcut * rcut;

	for (unsigned long int i = 0; i < N; i++)
		p[i].r[0][2] = p[i].r[1][2] = p[i].r[2][2] = 0;

	for (unsigned long int i = 0; i < (N - 1); i++) {
		for (unsigned long int j = i + 1; j < N; j++) {

			dist(p[i], p[j], box, s);
			r2 = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];

			if (r2 <= rcut2) {
				eps = inter[(int) p[i].type][(int) p[j].type].epsilon;
				sig = inter[(int) p[i].type][(int) p[j].type].sigma;

				r6i = 1.0 / (r2 * r2 * r2);
				ff = 48 * eps * r6i * (pow(sig, 12) * r6i - pow(sig, 6) * 0.5);

				for (int q = 0; q < 3; q++) {
					p[i].r[q][2] += ff * s[q] / r2;
					p[j].r[q][2] -= ff * s[q] / r2;
				}

				ecut = 4 * eps * (pow(sig/rcut, 12) - pow(sig/rcut, 6));
				*pe += 4 * eps * r6i * (pow(sig, 12) * r6i - pow(sig, 6)) - ecut;
				*virial += ff;
			}
		}
	}
}

int main(int argc, char *argv[]) {

	const gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

	/** Values that can be defined in command line arguments **/
	unsigned long int sample_frequency = 100; // snapshot writing frequency
	unsigned long int N = 2000;               // max number of particles
	unsigned long int n_steps = 1;            // integration steps
	double T = 1;                             // system temperature
	double nu = 2;                            // heat bath coupling frequency
	double rho = 0.85;                        // density
	double dt = 1E-3;                         // time step
	double rc = 2.5;                          // cutoff radius
	double ca = 0.5;                          // species A concentration
	double epsilon = 1;                       // epsilon AA
	double sigma = 1;                         // sigma AA
	double alpha = 1;                         // epsilon AB / epsilon AA
	double beta = 1;                          // sigma BB / sigma AA
	double delta = 1;                         // sigma BB / sigma AA
	double gamma = 1;                         // sigma AB / sigma AA

	/** Simulation variables */
	struct particle p[N]; // particles
	double pe;            // potential energy
	double ke;            // potential energy
	double etot;          // total energy
	double inst_temp;     // instantaneous system temperature
	double temp0, drift = 0;
	double virial, pressure;
	double box_volume, box_length;
	double cb; // particle species B concentration

	unsigned long int step_count = 0;
	FILE *out;
	char filename[50];

	// for now let's consider two particle species
	struct interaction inter[2][2];

	/* Parse command line arguments */
	read_flags(argc, argv, &N, &n_steps, &T, &nu, &dt, &rho, &rc,
			&sample_frequency, &alpha, &beta, &delta, &gamma, &ca);

	cb = 1 - ca; // B species concentration

	/* Set parameters */
	inter[0][0].sigma = sigma;
	inter[0][0].epsilon = epsilon;

	// Symmetric matrix
	inter[0][1].sigma = inter[1][0].sigma = gamma * sigma;
	inter[0][1].epsilon = inter[1][0].epsilon = alpha * epsilon;

	inter[1][1].sigma = beta * sigma;
	inter[1][1].epsilon = delta * epsilon;

	/** Initialization calls */
	box_volume = ((double)N) / rho;
	box_length = pow(box_volume, 1.0 / 3);
	init(N, rng, T, box_length, p, ca);
	forces(N, box_length, p, rc, &pe, &virial, inter);

	/** MD loop */
	printf("# step\ttemp\ttemp drift\tpressure\n");
	do {
		integrate(1, rng, N, T, nu, dt, p, &pe, &ke, &etot);
		forces(N, box_length, p, rc, &pe, &virial, inter);
		integrate(2, rng, N, T, nu, dt, p, &pe, &ke, &etot);

		inst_temp = ke * 2.0 / (3 * N);
		pressure = rho * inst_temp + virial / (3 * box_volume);

		if (step_count == 0)
			temp0 = inst_temp;
		else
			drift = (inst_temp - temp0) / temp0;

		printf("%lu\t%lf\t%E\t%lf\n", step_count, inst_temp, drift, pressure);

		step_count++;

		// Save simulation snapshot every sample_frequency steps
		if(!(step_count % sample_frequency)) {
			sprintf(filename, "snapshot_%lu.dat", step_count);
			out = fopen(filename, "w");
			snapshot(out, N, p);
			fclose(out);
		}
	} while (step_count <= n_steps);

	return 0;
}
