ARCH =intel
# compiler, can be intel, pgf90, gfortran ....
# OMP = on
OMP =
# for OpenMP option
#NC_LIB=/usr/local/netcdf/intel/4.1.1/
NC_LIB=/share/apps/netcdf/intel/4.1.1/
#NC_LIB2=-I/usr/local/netcdf/4.1.1-intel/include

ifeq ($(ARCH),intel)
F90= ifort
LD= ifort
LDFLAGS = -L$(NC_LIB)lib -I$(NC_LIB)include -lnetcdff   #-pg
FFLAGS = -L$(NC_LIB)lib -I$(NC_LIB)include -lnetcdff -nowarn
ifdef OMP
	LDFLAGS += -openmp
  FFLAGS  += -openmp
endif
#-check all
endif

ifeq ($(ARCH),pgf90)

F90= pgf90
LD= pgf90
FFLAGS = -O -Mbackslash
ifdef OMP
	LDFLAGS += -mp
        FFLAGS  += -mp
endif
endif



ifeq ($(ARCH),gfortran)
F90= gfortran
LD= gfortran
FFLAGS =  -O -ffree-line-length-none
ifdef OMP
	LDFLAGS += -fopenmp
        FFLAGS  += -fopenmp
endif
endif


CMD = fclimdex.exe

SRCS = modules.f90 main.f90 utility.f90 CDD.f90 FD.f90 GSL.f90 QC.f90 R95p.f90 Rnnmm.f90 RX5day.f90 TX10p.f90 TXX.f90 inp_nc_f90.f90 out_nc_f90.f90
 
OBJECTS = $(SRCS:%.f90=%.o)

.SUFFIXES:
.SUFFIXES: .o .f90

.PHONY: clean 
#del_data

.f90.o:
	$(F90) -c $(FFLAGS) $<

fclimdex.exe: $(OBJECTS)
	$(LD) $(LDFLAGS) -o $(CMD) $(OBJECTS)

clean:
	rm -f *.o *.mod core

#del_data:
#	@echo   Warning: all data will be deleted !
#	rm -f ./data/*
