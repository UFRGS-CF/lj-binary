/**
 * NVT molecular dynamics with Andersen thermostat and vel. verlet integrator
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

struct cartesian {
	double x, y, z;
};

/**
 * Return random number in [0, 1]
 */
double ranf() { return (double)rand() / RAND_MAX; }

/**
 * Return value from gaussian distribution with l0 mean and stdev sigma
 */
double gauss(double sigma, double l0) {

	double l;
	double v1, v2;
	double r = 2;

	do {
		v1 = 2.0 * ranf() - 1;
		v2 = 2.0 * ranf() - 1;
		r = pow(v1, 2) + pow(v2, 2);
	} while (r >= 1);

	l = l0 + sigma * v1 * sqrt(-2.0 * log(r) / r);

	return l;
}

/*
 * Write simultion coordinates and velocities to output stream fp
 */
void snapshot(FILE *fp, unsigned long int N, struct cartesian r[N], struct
		cartesian v[N]) {

	fprintf(fp, "rx\try\trz\tvx\tvy\tvz\n");
	for (unsigned long int i = 0; i < N; ++i) {
		fprintf(fp, "%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n",
				r[i].x, r[i].y, r[i].z, v[i].x, v[i].y,
				v[i].z);
	}
}

void init(unsigned long int N, double temp, double box, struct cartesian r[N],
		struct cartesian v[N]) {

	struct cartesian sumv, vcm;
	double sumv2 = 0;
	double mv2, fs;
	int a, max;
	unsigned long int n = 0;

	sumv.x = sumv.y = sumv.z = 0;

	a = pow(N, 1.0 / 3);

	if (pow(a, 3) != N)
		max = a + 1;
	else
		max = a;

	for (int i = 0; i < max; i++) {
		for (int j = 0; j < max; j++) {
			for (int k = 0; k < max; k++) {
				if (n < N) {
					r[n].x = box * ((float)k) / max;
					r[n].y = box * ((float)j) / max;
					r[n].z = box * ((float)i) / max;
					n++;
				}
			}
		}
	}

	// generate random velocities in range [-1, 1]
	for (unsigned long int i = 0; i < N; i++) {
		v[i].x = 2 * ((double)rand() / RAND_MAX) - 1;
		v[i].y = 2 * ((double)rand() / RAND_MAX) - 1;
		v[i].z = 2 * ((double)rand() / RAND_MAX) - 1;

		sumv.x += v[i].x;
		sumv.y += v[i].x;
		sumv.z += v[i].z;

		sumv2 += pow(v[i].x, 2) + pow(v[i].y, 2) + pow(v[i].z, 2);
	}

	vcm.x = sumv.x / N; // velocity center of mass
	vcm.y = sumv.y / N; // velocity center of mass
	vcm.z = sumv.z / N; // velocity center of mass

	mv2 = sumv2 / N; // mean square velocity

	// match system KE with temp, remove net CM drift
	fs = sqrt(3 * temp / mv2); // calculate scale factor
	for (unsigned long int i = 0; i < N; i++) {
		v[i].x = (v[i].x - vcm.x) * fs;
		v[i].y = (v[i].y - vcm.y) * fs;
		v[i].z = (v[i].z - vcm.z) * fs;
	}
}

void forces(unsigned long int N, double box, struct cartesian r[N], double rcut,
		struct cartesian f[N], double *pe, double *virial) {

	double r2, r2i, r6i, rcut2, ecut, ff;
	struct cartesian s;

	*pe = 0;
	*virial = 0;

	rcut2 = rcut * rcut;
	ecut = 4 * (1 / pow(rcut, 12) - 1 / pow(rcut, 6));

	for (unsigned long int i = 0; i < N; i++)
		f[i].x = f[i].y = f[i].z = 0;

	for (unsigned long int i = 0; i < (N - 1); i++) {
		for (unsigned long int j = i + 1; j < N; j++) {

			// distance between i and j
			s.x = r[i].x - r[j].x;
			s.y = r[i].y - r[j].y;
			s.z = r[i].z - r[j].z;

			// periodic boundary condition
			s.x = s.x - box * round(s.x / box);
			s.y = s.y - box * round(s.y / box);
			s.z = s.z - box * round(s.z / box);

			r2 = pow(s.x, 2) + pow(s.y, 2) + pow(s.z, 2);

			if (r2 <= rcut2) {
				// calculate forces
				r2i = 1 / r2;
				r6i = pow(r2i, 3);
				ff = 48 * r2i * r6i * (r6i - 0.5);

				f[i].x += ff * s.x;
				f[j].x -= ff * s.x;

				f[i].y += ff * s.y;
				f[j].y -= ff * s.y;

				f[i].z += ff * s.z;
				f[j].z -= ff * s.z;

				*pe += 4 * r6i * (r6i - 1) - ecut;
				*virial += ff;
			}
		}
	}
}

void integrate(int key, unsigned long int N, double temp, double nu, double dt,
		struct cartesian r[N], struct cartesian v[N],
		struct cartesian f[N], double *pe, double *ke, double *etot,
		double *inst_temp) {

	double sumv2;
	double m = 1; // particles mass, change this later
	double sigma = sqrt(temp);
	double randunif;

	if (key == 1) {
		for (unsigned long int i = 0; i < N; i++) {
			r[i].x = r[i].x + dt * v[i].x + dt * dt * f[i].x / 2;
			r[i].y = r[i].y + dt * v[i].y + dt * dt * f[i].y / 2;
			r[i].z = r[i].z + dt * v[i].z + dt * dt * f[i].z / 2;

			v[i].x = v[i].x + dt * f[i].x / 2;
			v[i].y = v[i].y + dt * f[i].y / 2;
			v[i].z = v[i].z + dt * f[i].z / 2;
		}
	} else if (key == 2) {
		sumv2 = 0;

		for (unsigned long int i = 0; i < N; i++) {
			v[i].x = v[i].x + dt * f[i].x / 2;
			v[i].y = v[i].y + dt * f[i].y / 2;
			v[i].z = v[i].z + dt * f[i].z / 2;

			sumv2 += pow(v[i].x, 2) + pow(v[i].y, 2) + pow(v[i].z, 2);
		}

		*inst_temp = m * sumv2 / (3 * N);

		for (unsigned long int i = 0; i < N; i++) {
			randunif = ((double)rand()) / RAND_MAX;
			if (randunif < (nu * dt)) {
				v[i].x = gauss(sigma, 0);
				v[i].y = gauss(sigma, 0);
				v[i].z = gauss(sigma, 0);
			}
		}

		*ke = 0.5 * sumv2;
		*etot = *pe + *ke;
	}
}

void help() {
	printf("NVT Lennard Jones simulation software.\n");
	printf("Configuration options:\n");
	printf("\t -N \tNumber of particles\n");
	printf("\t -T \tSystem temperature\n");
	printf("\t -nu \tHeath bath collision frequency\n");
	printf("\t -rho \tDensity\n");
	printf("\t -ns \tNumber of integration steps\n");
	printf("\t -dt \tTime step\n");
	printf("\t -rc \tCutoff radius\n");
	printf("\t -fs \tSnapshot sample frequency\n");
	printf("\t -h  \tPrint this message\n");
}

int main(int argc, char *argv[]) {

	/** Values that can be defined in command line arguments **/
	unsigned long int sample_frequency = 100; // snapshot writing frequency
	unsigned long int N = 1000;    // max number of particles
	unsigned long int n_steps = 1; // integration steps
	double T = 1;                  // system temperature
	double nu = 2;                 // heat bath coupling frequency
	double rho = 0.85;             // density
	double dt = 1E-3;              // time step
	double rc = 2.5;               // cutoff radius

	/** Simulation variables */
	struct cartesian r[N]; // particle positions
	struct cartesian v[N]; // particle velocities
	struct cartesian f[N]; // forces on particles
	double pe;             // potential energy
	double ke;             // potential energy
	double etot;           // total energy
	double inst_temp;      // instantaneous system temperature
	double temp0, drift = 0;
	double virial, pressure;
	double box_volume, box_length;

	unsigned long int step_count = 0;
	FILE *out;
	char filename[50];

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
		else if (!strcmp(argv[i], "-h")) {
			help();
			exit(0);
		} else {
			fprintf(stderr, "Error: '%s' not recognized.\n", argv[i]);
			exit(1);
		}
	}

	/** Initialization calls */
	srand(time(NULL));
	box_volume = ((double)N) / rho;
	box_length = pow(box_volume, 1.0 / 3);
	init(N, T, box_length, r, v);
	forces(N, box_length, r, rc, f, &pe, &virial);

	/** MD loop */
	printf("# step\ttemp\ttemp drift\tpressure\n");
	do {
		integrate(1, N, T, nu, dt, r, v, f, &pe, &ke, &etot, &inst_temp);
		forces(N, box_length, r, rc, f, &pe, &virial);
		integrate(2, N, T, nu, dt, r, v, f, &pe, &ke, &etot, &inst_temp);

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
			snapshot(out, N, r, v);
			fclose(out);
		}
	} while (step_count <= n_steps);

	return 0;
}
