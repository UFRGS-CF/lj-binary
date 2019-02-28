struct cartesian {
	double x, y, z;
};

struct particle {
	struct cartesian r;
	struct cartesian v;
	struct cartesian f;
};
