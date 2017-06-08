# Makefile for NEW inversion code
#-------------------------------------------
# JM Borrero March 21, 2013
#FC = mpif90
#CC = mpicxx
FC = gfortran
CC = gcc
F90F = -c -O3
F90F1 = -c -O3
F77F1 = -c -O3
F77F2 = -c -O3 -r8
CF = -c
# For Laptop
LIBS = -L /opt/cfitsio -lcfitsio \
	-L /opt/lapack-3.5.0 -llapack -lblas -ltmglib \
	-L /usr/lib64 -lpthread

LDFLAGS =

OBJ = cons_param.o \
	derivvar.o \
	checknan.o \
	precision.o \
	voigt.o \
	splines.o \
	time_param.o \
	strings.o \
	log.o \
	lines_database.o \
	atom_database.o \
	damping.o \
	magsplit.o \
	absorption_matrix.o \
	chemical_equilibrium.o \
	simul_box.o \
	background.o \
	rtesolver.o \
	deriv_test.o \
	test.o

cons_param.o: cons_param.f90
	$(FC) $(F90F) cons_param.f90

derivvar.o: derivvar.f90  cons_param.o
	$(FC) $(F90F) derivvar.f90

checknan.o: checknan.f90
	$(FC) $(F90F) checknan.f90

precision.o: precision.f90
	$(FC) $(F90F) precision.f90

voigt.o: voigt.f
	$(FC) $(F77F1) voigt.f

splines.o: splines.f90 cons_param.o
	$(FC) $(F90F) splines.f90

time_param.o: time_param.f90 cons_param.o
	$(FC) $(F90F) time_param.f90

strings.o: strings.f90 precision.o
	$(FC) $(F90F) strings.f90

log.o: log.f90 cons_param.o
	$(FC) $(F90F) log.f90

damping.o: damping.f90 cons_param.o atom_database.o \
	log.o splines.o derivvar.o
	$(FC) $(F90F) damping.f90

lines_database.o: lines_database.f90 cons_param.o log.o damping.o
	$(FC) $(F90F) lines_database.f90

atom_database.o: atom_database.f90 cons_param.o
	$(FC) $(F90F) atom_database.f90

magsplit.o: magsplit.f90 cons_param.o lines_database.o simul_box.o
	$(FC) $(F90F) magsplit.f90

absorption_matrix.o: absorption_matrix.f90 cons_param.o lines_database.o simul_box.o \
	magsplit.o atom_database.o voigt.o checknan.o
	$(FC) $(F90F) absorption_matrix.f90

chemical_equilibrium.o: chemical_equilibrium.f90 cons_param.o lines_database.o \
	atom_database.o derivvar.o
	$(FC) $(F90F) chemical_equilibrium.f90

simul_box.o: simul_box.f90 cons_param.o strings.o log.o \
	lines_database.o chemical_equilibrium.o
	$(FC) $(F90F) simul_box.f90

background.o: background.f90
	$(FC) $(F90F) background.f90

rtesolver.o: rtesolver.f90 cons_param.o simul_box.o absorption_matrix.o checknan.o
	$(FC) $(F90F) rtesolver.f90

deriv_test.o: deriv_test.f90 cons_param.o lines_database.o atom_database.o simul_box.o \
	chemical_equilibrium.o
	$(FC) $(F90F) deriv_test.f90

test.o: test.f90 log.o time_param.o cons_param.o lines_database.o simul_box.o \
	atom_database.o chemical_equilibrium.o damping.o absorption_matrix.o rtesolver.o background.o \
	checknan.o deriv_test.o derivvar.o
	$(FC) $(F90F) test.f90

test: $(OBJ)
	$(FC) $(OBJ) $(LDFLAGS) $(LIBS) -o test.x

clean:
	rm -f *.o *.a *.x *~ *.mod
