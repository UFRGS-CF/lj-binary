/*
 * md.h
 *
 * Functions and structures useful to simulations and data post-processing
 *
 */

struct cartesian {
	double x, y, z;
};

struct particle {
	struct cartesian r;
	struct cartesian v;
	struct cartesian f;
	unsigned long int type;
};

struct interaction {
	double epsilon;
	double sigma;
};

struct cartesian dist(struct particle i, struct particle j, double box);
