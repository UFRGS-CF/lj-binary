#include "md.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*
 * Write simultion coordinates and velocities to output stream fp
 */
void snapshot(FILE *fp, unsigned long int N, struct particle p[N]) {

	fprintf(fp, "rx\try\trz\tvx\tvy\tvz\n");
	for (unsigned long int i = 0; i < N; ++i) {
		fprintf(fp, "%lu\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n",
				p[i].type,
				p[i].r[0][0], p[i].r[1][0], p[i].r[2][0],
				p[i].r[0][1], p[i].r[1][1], p[i].r[2][1]);
	}
}

void help() {
	printf("NVT Lennard Jones simulation software.\n");
	printf("Configuration options:\n");
	printf("\t -N \t\tNumber of particles\n");
	printf("\t -T \t\tSystem temperature\n");
	printf("\t -nu \t\tHeat bath collision frequency\n");
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

void read_flags(int argc, char *argv[], unsigned long int *N, unsigned long int
		*n_steps, double *T, double *nu, double *dt, double *rho,
		double *rc, unsigned long int *sample_frequency, double *alpha,
		double *beta, double *delta, double *gamma, double *ca) {

	// Start at 1 because 0 is program name
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-N"))
			*N = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-ns"))
			*n_steps = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-T"))
			*T = atof(argv[++i]);
		else if (!strcmp(argv[i], "-nu"))
			*nu = atof(argv[++i]);
		else if (!strcmp(argv[i], "-rho"))
			*rho = atof(argv[++i]);
		else if (!strcmp(argv[i], "-dt"))
			*dt = atof(argv[++i]);
		else if (!strcmp(argv[i], "-rc"))
			*rc = atof(argv[++i]);
		else if (!strcmp(argv[i], "-sf"))
			*sample_frequency = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-alpha"))
			*alpha = atof(argv[++i]);
		else if (!strcmp(argv[i], "-beta"))
			*beta = atof(argv[++i]);
		else if (!strcmp(argv[i], "-delta"))
			*delta = atof(argv[++i]);
		else if (!strcmp(argv[i], "-gamma"))
			*gamma = atof(argv[++i]);
		else if (!strcmp(argv[i], "-ca"))
			*ca = atof(argv[++i]);
		else if (!strcmp(argv[i], "-h")) {
			help();
			exit(0);
		} else {
			fprintf(stderr, "Error: '%s' not recognized.\n", argv[i]);
			exit(1);
		}
	}
}
