PROG =	dmat2

SRCS =	main.f90 matrix.f90 matrix_mapping.f90 numbers.f90 proc_mapping.f90 \
	scalapack_interfaces.f90

OBJS =	main.o matrix.o matrix_mapping.o numbers.o proc_mapping.o \
	scalapack_interfaces.o

LIBS =-lblacs-openmpi -lblacsF77init-openmpi -lscalapack-openmpi -llapack -lblas	

CC = cc
CFLAGS = -O
FC = f77
FFLAGS = -O
F90 = mpif90
F90FLAGS = -O -g -std=f2008 -Wall -Wextra -Wimplicit-interface -Wuse-without-only -fcheck=all -finit-real=snan
LDFLAGS = 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

main.o: matrix.o numbers.o
matrix.o: matrix_mapping.o numbers.o scalapack_interfaces.o
matrix_mapping.o: numbers.o proc_mapping.o
proc_mapping.o: 
scalapack_interfaces.o: numbers.o
