#Makefile
CMBONLYLIB = libactpol.a
OBJS = ACTPol_s2_like.o

RM = rm -f

LIBS ?= $(LAPACKL)

F90C    ?= mpif90 -f90=ifort
FFLAGS  = -mkl -openmp -O3 -no-prec-div -fpp -DMPI -g -traceback
LIBS =  -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
#------------------


#                       Rules.

all: $(CMBONLYLIB) test

$(CMBONLYLIB): $(OBJS)
	ar r $@ $^
 
test: $(OBJS) test.o
	$(F90C) $(FFLAGS) -o $@ test.o $(OBJS) $(LIBS)

%.o: %.f90
	$(F90C) $(FFLAGS) $(INCS) -c -o $@ $<

%.o: %.F90
	$(F90C) $(FFLAGS) $(INCS) -c -o $@ $<

clean:
	$(RM) *.o *.mod *.log *~ *.a test

distclean: clean
	$(RM) test

