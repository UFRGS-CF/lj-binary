nvt: nvt.c
	gcc -lm -lgsl -lgslcblas -O2 -Wall -Wextra nvt.c -o md.out

rdf: rdf.c
	gcc -lm -O2 -Wall -Wextra rdf.c -o rdf.out
