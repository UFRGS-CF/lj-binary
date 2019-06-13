void snapshot(FILE *fp, unsigned long int N, struct particle p[N]);
void help();
void read_flags(int argc, char *argv[], unsigned long int *N, unsigned long int
		*n_steps, double *T, double *nu, double *dt, double *rho,
		double *rc, unsigned long int *sample_frequency, double *alpha,
		double *beta, double *delta, double *gamma, double *ca);
