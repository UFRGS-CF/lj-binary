/*
 * md.h
 *
 * Functions and structures useful to simulations and data post-processing
 *
 */

struct particle {
	// i: 0, 1, 2 (spatial coordinates)
	// j: 0 (position), 1 (velocity), 2 (acceleration)
	double r[3][3];
	unsigned long int type;
};

struct interaction {
	double epsilon;
	double sigma;
};

void dist(struct particle i, struct particle j, double box, double s[3]);
