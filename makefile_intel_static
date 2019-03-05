#
EXT=.exe
NAME1=AlloModes
DIR=./
NAMEFUL1=$(DIR)/$(NAME1)$(EXT)

FC = ifort
#FFLAGS = -c -u -O
FFLAGS = -c -u -O3
CC = icc
CPP = icpc
CFLAGS = -c -O -Dintel
CPPFLAGS = -c -O
LDFLAGS = -O
#LIBRARIES = -lpthread -lifcore -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -lgmp -liomp5
LIBRARIES = -lpthread -lgmp -qopenmp -qopenmp-link=static -mkl -static-intel
LIBS = -L/opt/intel/lib/ -L/opt/intel/mkl/lib/ -L/usr/local/lib/ -L/opt/intel/mkl/lib/intel64

.c.o :
	$(CC) $(CFLAGS) $<

.cpp.o :
	$(CPP) $(CPPFLAGS) $<

.f.o :
	$(FC) $(FFLAGS) $<

OBJECTS1 = \
$(NAME1).o \
protein_topology.o geometry_tools.o \
readpdb.o select_variables.o \
sorting_tools.o basic_arpack.o \
bestfitm.o \
hessian.o \
listcontact_cutoff.o \
listcontact_del.o delaunay.o ran2.o sos_minor_gmp.o \
normal_modes.o rescale_eigvect.o \
compute_bfact.o overlap.o \
select_atoms.o eigen_deriv.o score_allostery.o lbfgsb.o linpack.o \
write_results.o write_pdb.o

$(NAMEFUL1) : $(OBJECTS1)
	$(FC) -o $(NAMEFUL1) $(LDFLAGS) $(OBJECTS1) $(LIBS) $(LIBRARIES)

all: $(OBJECTS1)
	$(FC) -o $(NAMEFUL1) $(LDFLAGS) $(OBJECTS1) $(LIBS) $(LIBRARIES)

clean:
	touch junk.o; rm -f *.o $(NAMEFUL1)

$(OBJECTS1) :
