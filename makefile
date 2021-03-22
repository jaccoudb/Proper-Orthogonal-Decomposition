FC = gfortran
FLIB = -llapack -lblas
FOPT = -O3

all: pod_svd_program.o runcode
pod_svd_program.o: pod_svd_program.f90
	$(FC) $(FOPT) pod_svd_program.f90 -c
runcode: pod_svd_program.o
	$(FC) $(FOPT) pod_svd_program.o -o runcode $(FLIB)
clean:
	rm *.o *.mod runcode
