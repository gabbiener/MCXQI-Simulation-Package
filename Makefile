all: make_source SPDC propagate_to aperture print_rays design_object absorption_object

print_rays : vector.o ray.o
	gcc -lm -Wall print_rays.c vector.o ray.o -o print_rays

aperture : vector.o ray.o
	gcc -lm -Wall aperture.c vector.o ray.o -o aperture

propagate_to: vector.o ray.o
	gcc -lm -Wall propagate_to.c vector.o ray.o -o propagate_to

make_source: vector.o ray.o rng-double.o
	gcc -lm -Wall make_source.c rng-double.o vector.o ray.o -o make_source

SPDC: vector.o ray.o rng-double.o
	gcc -lm -Wall SPDC.c rng-double.o vector.o ray.o -o SPDC

design_object: vector.o ray.o rng-double.o
	gcc -lm -Wall design_object.c rng-double.o vector.o ray.o -o design_object

absorption_object: vector.o ray.o rng-double.o
	gcc -lm -Wall absorption_object.c rng-double.o vector.o ray.o -o absorption_object

random:
	gcc -lm -Wall -c rng-double.c rng-double.o

ray: vector.o
	gcc -lm -Wall -c ray.c vector.o

vector: 
	gcc -lm -Wall -c vector.c
clean:
	rm -rf *.o 
	rm -rf make_source aperture print_rays propagate_to SPDC design_object absorption_object