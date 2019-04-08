/*
 *
 * Reads any number of simulation snapshots and computes g(r)
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "md.h"

void snapshot_in(FILE *snapshot, unsigned long int N, struct particle p[N]) {

	double vx, vy, vz;

	fscanf(snapshot, "%*[^\n]\n", NULL); // ignore file header

	for (unsigned long int i = 0; i < N; i++) {
		fscanf(snapshot, "%lu\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",
				&(p[i].type),
				&(p[i].r.x), &(p[i].r.y), &(p[i].r.z),
				&vx, &vy, &vz);
	}
}

void update_histogram(const unsigned long int N, const double box_lenght,
		struct particle p[N], unsigned long int *H, const double rcut,
		const double bin_size) {

	struct cartesian s;
	double r2;

	int histogram_index;
	const double rcut2 = rcut * rcut;

	for (unsigned long int i = 0; i < (N - 1); i++) {
		for (unsigned long int j = i + 1; j < N; j++) {

			if(p[i].type == 0 && p[j].type == 1) {

				s = dist(p[i], p[j], box_lenght);
				r2 = pow(s.x, 2) + pow(s.y, 2) + pow(s.z, 2);

				if (r2 <= rcut2) {
					histogram_index = (int)(sqrt(r2) / bin_size);
					H[histogram_index] += 2;
				}
			}

		}
	}
}

void help() {
	fprintf(stdout, "Options:\n");
	fprintf(stdout, "\t -rho [real]\t\tDensity\n");
	fprintf(stdout, "\t -dr [real]\t\tBin size\n");
	fprintf(stdout, "\t -rc [real]\t\tCutoff radius\n");
	fprintf(stdout, "\t -start [integer]\tFirst snapshot\n");
	fprintf(stdout, "\t -step [integer]\tSnapshot step\n");
	fprintf(stdout, "\t -stop [integer]\tFinal snapshot\n");
	fprintf(stdout, "\t -h           \t\tPrint help\n");
}

int main(int argc, char *argv[]) {

	/* Command line arguments */
	unsigned long int start = 10000, step = 10, stop = 100000;
	double rho = 0.85;
	double rcut = 2.5;
	double bin_size = 0.02;

	/* Sampling variables */
	struct particle *p; //particles
	unsigned long int N = 0; // number of particles
	unsigned long int *H;    // histogram
	unsigned long int n_bins;
	unsigned long int rdf_sample; // number of snapshots samples in *H
	double *g;                    // rdf
	FILE *snapshot;
	char filename[100];            // store snapshot file name
	double box_volume, box_lenght; // simulation box properties
	double bin_volume, N_ideal_gas, dist;

	/* Parse command line arguments */
	// Start at 1 because 0 is program name
	for (int i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-start"))
			start = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-stop"))
			stop = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-step"))
			step = atof(argv[++i]);
		else if (!strcmp(argv[i], "-rho"))
			rho = atof(argv[++i]);
		else if (!strcmp(argv[i], "-dr"))
			bin_size = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-rc"))
			rcut = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-h")) {
			help();
			exit(0);
		} else {
			fprintf(stderr, "Error: '%s' not recognized.\n", argv[i]);
			exit(1);
		}
	}

	/* Init calls */
	n_bins = (unsigned long int)(rcut / bin_size) + 1;

	rdf_sample = 0;
	sprintf(filename, "snapshot_%lu.dat", start);
	snapshot = fopen(filename, "r");
	if (!snapshot) {
		fprintf(stdout, "Error: could not read %s\n", filename);
		return 1;
	} else {
		// Count number of lines in the file
		while (!feof(snapshot)) {
			char ch = fgetc(snapshot);
			if (ch == '\n')
				N++;
		}
	}
	fclose(snapshot);
	// There is a trailing line, so number of particles is
	// number of lines minus 1
	N = N - 1;

	H = calloc(n_bins, sizeof(unsigned long int));
	g = calloc(n_bins, sizeof(double));
	p = calloc(N, sizeof(struct particle));

	box_volume = ((double)N) / rho;
	box_lenght = pow(box_volume, 1.0 / 3);

	/* Calculate rdf */
	for (unsigned long int i = start; i <= stop; i += step) {
		sprintf(filename, "snapshot_%lu.dat", i);
		snapshot = fopen(filename, "r");
		if (!snapshot) {
			fprintf(stdout, "Error: could not read %s\n", filename);
			return 1;
		} else {
			snapshot_in(snapshot, N, p);
			update_histogram(N, box_lenght, p, H, rcut, bin_size);
			rdf_sample++;
		}
		fclose(snapshot);
	}

	/* Normalize and output rdf */
	for (unsigned long int i = 0; i < n_bins; i++) {
		// assume particle is halfway in the bin
		dist = bin_size * (i + 0.5);
		// volume between bins i+1 and i
		bin_volume = (pow(i + 1, 3) - pow(i, 3)) * pow(bin_size, 3.0);
		// number of ideal gas particles in bin_volume
		// 0.5 because concentration is 1/2
		N_ideal_gas = (4.0 / 3) * M_PI * bin_volume * (rho / 2);
		// normalize
		g[i] = ((double)H[i]) / (N_ideal_gas * (N/2) * rdf_sample);

		fprintf(stdout, "%.4lf\t%.8lf\n", dist, g[i]);
	}

	return 0;
}
