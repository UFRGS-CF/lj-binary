#include <math.h>
#include "md.h"

/*
 * Returns distance between particles i and j in each component
 */
struct cartesian dist(struct particle i, struct particle j, double box){

	struct cartesian s;

	// distance between i and j
	s.x = i.r.x - j.r.x;
	s.y = i.r.y - j.r.y;
	s.z = i.r.z - j.r.z;

	// apply periodic boundary condition
	s.x = s.x - box * round(s.x / box);
	s.y = s.y - box * round(s.y / box);
	s.z = s.z - box * round(s.z / box);

	return s;
}
