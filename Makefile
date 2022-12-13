all: openMPI MPI Imperative

openMPI:
	mpic++ -fopenmp -O5 parallelOMPI.cpp -o parallelO.o
MPI:
	mpicc -O5 parallelMPI.cpp -o parallel.o
Imperative:
	mpicc -o imperative.o imperative.c -O3

clean:
	rm parallelO.o parallel.o imperative.o