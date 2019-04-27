/**
 * NVT molecular dynamics with Andersen thermostat and vel. verlet integrator
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "md.h"

/*
 * Write simultion coordinates and velocities to output stream fp
 */
void snapshot(FILE *fp, unsigned long int N, struct particle p[N]) {

	fprintf(fp, "rx\try\trz\tvx\tvy\tvz\n");
	for (unsigned long int i = 0; i < N; ++i) {
		fprintf(fp, "%lu\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n",
				p[i].type,
				p[i].r.x, p[i].r.y, p[i].r.z, p[i].v.x,
				p[i].v.y, p[i].v.z);
	}
}

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
					p[n].r.x = box * ((float)k) / max;
					p[n].r.y = box * ((float)j) / max;
					p[n].r.z = box * ((float)i) / max;
					n++;
				}
			}
		}
	}

}

void init(unsigned long int N, const gsl_rng *rng, double temp, double box,
		struct particle p[N], double ca) {

	struct cartesian sumv, vcm;
	double sumv2 = 0;
	double mv2, fs;

	sumv.x = sumv.y = sumv.z = 0;

	// setup particles in cubic lattice
	lattice(N, box, p);

	// generate random velocities in range [-1, 1]
	for (unsigned long int i = 0; i < N; i++) {
		p[i].v.x = gsl_ran_flat(rng, -1, 1);
		p[i].v.y = gsl_ran_flat(rng, -1, 1);
		p[i].v.z = gsl_ran_flat(rng, -1, 1);

		sumv.x += p[i].v.x;
		sumv.y += p[i].v.x;
		sumv.z += p[i].v.z;

		sumv2 += pow(p[i].v.x, 2) + pow(p[i].v.y, 2) + pow(p[i].v.z, 2);
	}

	vcm.x = sumv.x / N; // velocity center of mass
	vcm.y = sumv.y / N; // velocity center of mass
	vcm.z = sumv.z / N; // velocity center of mass

	mv2 = sumv2 / N; // mean square velocity

	// match system KE with temp, remove net CM drift
	fs = sqrt(3 * temp / mv2); // calculate scale factor
	for (unsigned long int i = 0; i < N; i++) {
		p[i].v.x = (p[i].v.x - vcm.x) * fs;
		p[i].v.y = (p[i].v.y - vcm.y) * fs;
		p[i].v.z = (p[i].v.z - vcm.z) * fs;
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
	double r2, r2i, r6i, rcut2, ecut, ff;
	struct cartesian s;

	*pe = 0;
	*virial = 0;
	rcut2 = rcut * rcut;

	for (unsigned long int i = 0; i < N; i++)
		p[i].f.x = p[i].f.y = p[i].f.z = 0;

	for (unsigned long int i = 0; i < (N - 1); i++) {
		for (unsigned long int j = i + 1; j < N; j++) {

			s = dist(p[i], p[j], box);
			r2 = pow(s.x, 2) + pow(s.y, 2) + pow(s.z, 2);

			if (r2 <= rcut2) {
				eps = inter[(int) p[i].type][(int) p[j].type].epsilon;
				sig = inter[(int) p[i].type][(int) p[j].type].sigma;

				r2i = 1.0 / r2;
				r6i = pow(r2i, 3);
				ff = 48 * eps * r2i * r6i * (pow(sig, 12) * r6i - pow(sig, 6) * 0.5);

				p[i].f.x += ff * s.x;
				p[j].f.x -= ff * s.x;

				p[i].f.y += ff * s.y;
				p[j].f.y -= ff * s.y;

				p[i].f.z += ff * s.z;
				p[j].f.z -= ff * s.z;

				ecut = 4 * eps * (pow(sig/rcut, 12) - pow(sig/rcut, 6));
				*pe += 4 * eps * r6i * (pow(sig, 12) * r6i - pow(sig, 6)) - ecut;
				*virial += ff;
			}
		}
	}
}

void integrate(int key, const gsl_rng *rng, unsigned long int N, double temp,
		double nu, double dt, struct particle p[N], double *pe, double
		*ke, double *etot, double *inst_temp) {

	double sumv2;
	double m = 1; // particles mass, change this later
	double sigma = sqrt(temp);
	double randunif;

	if (key == 1) {
		for (unsigned long int i = 0; i < N; i++) {
			p[i].r.x = p[i].r.x + dt * p[i].v.x + dt * dt * p[i].f.x / 2;
			p[i].r.y = p[i].r.y + dt * p[i].v.y + dt * dt * p[i].f.y / 2;
			p[i].r.z = p[i].r.z + dt * p[i].v.z + dt * dt * p[i].f.z / 2;

			p[i].v.x = p[i].v.x + dt * p[i].f.x / 2;
			p[i].v.y = p[i].v.y + dt * p[i].f.y / 2;
			p[i].v.z = p[i].v.z + dt * p[i].f.z / 2;
		}
	} else if (key == 2) {
		sumv2 = 0;

		for (unsigned long int i = 0; i < N; i++) {
			p[i].v.x = p[i].v.x + dt * p[i].f.x / 2;
			p[i].v.y = p[i].v.y + dt * p[i].f.y / 2;
			p[i].v.z = p[i].v.z + dt * p[i].f.z / 2;

			sumv2 += pow(p[i].v.x, 2) + pow(p[i].v.y, 2) + pow(p[i].v.z, 2);
		}

		*inst_temp = m * sumv2 / (3 * N);

		for (unsigned long int i = 0; i < N; i++) {
			randunif = gsl_ran_flat(rng, 0, 1);
			if (randunif < (nu * dt)) {
				p[i].v.x = gsl_ran_gaussian(rng, sigma);
				p[i].v.y = gsl_ran_gaussian(rng, sigma);
				p[i].v.z = gsl_ran_gaussian(rng, sigma);
			}
		}

		*ke = 0.5 * sumv2;
		*etot = *pe + *ke;
	}
}

void help() {
	printf("NVT Lennard Jones simulation software.\n");
	printf("Configuration options:\n");
	printf("\t -N \t\tNumber of particles\n");
	printf("\t -T \t\tSystem temperature\n");
	printf("\t -nu \t\tHeath bath collision frequency\n");
	printf("\t -rho \t\tDensity\n");
	printf("\t -ns \t\tNumber of integration steps\n");
	printf("\t -dt \t\tTime step\n");
	printf("\t -rc \t\tCutoff radius\n");
	printf("\t -fs \t\tSnapshot sample frequency\n");
	printf("\t -epsilon \tEpsilon for type A particle\n");
	printf("\t -sigma \tSigma for type A particle\n");
	printf("\t -alpha \tepsilon AB / epsilon AA\n");
	printf("\t -beta \t\tepsilon BB / epsilon AA\n");
	printf("\t -delta \tsigma BB / sigma AA\n");
	printf("\t -gamma \t\tsigma AB / sigma AA\n");
	printf("\t -ca \tParticle A concentration\n");
	printf("\t -h  \t\tPrint this message\n");
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
	// Start at 1 because 0 is program name
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-N"))
			N = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-ns"))
			n_steps = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-T"))
			T = atof(argv[++i]);
		else if (!strcmp(argv[i], "-nu"))
			nu = atof(argv[++i]);
		else if (!strcmp(argv[i], "-rho"))
			rho = atof(argv[++i]);
		else if (!strcmp(argv[i], "-dt"))
			dt = atof(argv[++i]);
		else if (!strcmp(argv[i], "-rc"))
			rc = atof(argv[++i]);
		else if (!strcmp(argv[i], "-sf"))
			sample_frequency = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-alpha"))
			alpha = atof(argv[++i]);
		else if (!strcmp(argv[i], "-beta"))
			beta = atof(argv[++i]);
		else if (!strcmp(argv[i], "-delta"))
			delta = atof(argv[++i]);
		else if (!strcmp(argv[i], "-gamma"))
			gamma = atof(argv[++i]);
		else if (!strcmp(argv[i], "-ca"))
			ca = atof(argv[++i]);
		else if (!strcmp(argv[i], "-h")) {
			help();
			exit(0);
		} else {
			fprintf(stderr, "Error: '%s' not recognized.\n", argv[i]);
			exit(1);
		}
	}

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
		integrate(1, rng, N, T, nu, dt, p, &pe, &ke, &etot, &inst_temp);
		forces(N, box_length, p, rc, &pe, &virial, inter);
		integrate(2, rng, N, T, nu, dt, p, &pe, &ke, &etot, &inst_temp);

		if (step_count == 0)
			temp0 = inst_temp;
		else
			drift = (inst_temp - temp0) / temp0;

		pressure = rho * inst_temp + virial / box_volume;
		printf("%lu\t%lf\t%E\t%lf\n", step_count, inst_temp, drift,
				pressure);

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
