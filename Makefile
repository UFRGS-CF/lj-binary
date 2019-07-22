## Insert GSL path here if static compiling is desired
#LIBGSL=gsl/.libs/libgsl.a gsl/cblas/.libs/libgslcblas.a
LIBGSL=-lgsl -lgslcblas

WARN=-Wall -Wextra

all: md rdf

md: main.c
	gcc main.c integrate.c io.c md.c -lm ${LIBGSL} -O2 ${WARN} -o md.out

rdf: rdf.c
	gcc rdf.c md.c -lm -O2 ${WARN} -o rdf.out
