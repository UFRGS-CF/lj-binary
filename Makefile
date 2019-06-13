all: md rdf

md: main.c
	gcc -lm -lgsl -lgslcblas -O2 -Wall -Wextra main.c integrate.c io.c md.c -o md.out

rdf: rdf.c
	gcc -lm -O2 -Wall -Wextra rdf.c md.c -o rdf.out
