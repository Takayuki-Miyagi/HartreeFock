#--------------------------------------------------
# Make file for TTFcalc code
#--------------------------------------------------
TARGET=HartreeFock
INSTLDIR= $(HOME)/bin
Host= $(shell if hostname|grep -q apt1; then echo apt; \
  elif hostname|grep -q oak; then echo oak; \
  elif hostname|grep -q cedar; then echo cedar; \
  else echo other; fi)
HOST=$(strip $(Host))
$(info Host:$(HOST))

#--------------------------------------------------
# Default Parameters
#--------------------------------------------------
ifeq ($(HOST),other)
  FDEP=makedepf90
  CC=gcc
  FC=gfortran
  LDFLAGS= -I/usr/local/include -L/usr/local/lib -llapack -lblas -lgsl
  OMP = -fopenmp
  FFLAGS=-O3 -Dsingle_precision
  CFLAGS=-O3
  FF2C= -ff2c
  FDFLAGS=
  #FDFLAGS+=-DModelSpaceDebug
  #FDFLAGS+=-DNOperatorsDebug
  FDFLAGS+=-fbounds-check -Wall -fbacktrace -O -Wuninitialized
endif

ifeq ($(HOST),apt)
  #FDEP=
  CC=gcc
  FC=ifort
  LDFLAGS=-mkl -lgsl
  OMP = -openmp
  FFLAGS=-O3 -no-ipo -static -Dsingle_precision
  CFLAGS=-O3
  FF2C=
  #FDFLAGS+=-check-all
endif

#-----------------------------
# oak (oak.arc.ubc.ca)
# -heap-arrays option is needed for built in transpose function with
#  large dimension matrix
#-----------------------------
ifeq ($(strip $(HOST)),oak)
  FDEP=
  CC=icc
  FC=ifort
  MKL=-L$(MKLROOT)/lib/ -L$(MKLROOT)/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -Wl,-rpath,$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/../compiler/lib/
  LDFLAGS=$(MKL) -lgsl
  OMP = -qopenmp
  FFLAGS=-O3 -heap-arrays
  CFLAGS=-O3
  LINT = -i8
  FF2C=
  #FDFLAGS =-check all
  FFLAGS+= -Dsingle_precision
endif

#--------------------------------------------------
# Source Files
#--------------------------------------------------
SRCDIR = src
DEPDIR = .
SRCC = $(wildcard $(SRCDIR)/*.c)
SRCF77 = $(wildcard $(SRCDIR)/*.f)
SRCF90 = $(wildcard $(SRCDIR)/*.f90)
SRCF95 = $(wildcard $(SRCDIR)/*.F90)

OBJDIR = obj
OBJC = $(SRCC:$(SRCDIR)/%.c=$(OBJDIR)/%.o)
OBJF77 = $(SRCF77:$(SRCDIR)/%.f=$(OBJDIR)/%.o)
OBJF90 = $(SRCF90:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
OBJF95 = $(SRCF95:$(SRCDIR)/%.F90=$(OBJDIR)/%.o)

MODDIR = src
MODF90 = $(SRCF90:$(SRCDIR)/%.f90=$(MODDIR)/%.mod)
MODF95 = $(SRCF95:$(SRCDIR)/%.F90=$(MODDIR)/%.mod)

#--------------------------------------------------
# Source Files (LinAlgf90)
#--------------------------------------------------
LINSRCDIR = LinAlgf90/src
LINSRC = $(wildcard $(LINSRCDIR)/*.f90)
LINOBJ = $(LINSRC:$(LINSRCDIR)/%.f90=$(OBJDIR)/%.o)
LINMOD = $(LINSRC:$(LINSRCDIR)/%.f90=$(MODDIR)/%.mod)

#--------------------------------------------------
SRCS = $(SRCC) $(SRCF77) $(SRCF90) $(SRCF95) $(LINSRC)
OBJS = $(OBJC) $(OBJF77) $(OBJF90) $(OBJF95) $(LINOBJ)
MODSUP = $(MODF90) $(MODF95) $(LINMOD)
MODS = $(shell echo $(MODSUP) | tr A-Z a-z)
#$(info "soruce files:" $(SRCS))
#$(info "object files:" $(OBJS))
#$(info "module files:" $(MODS))

MODOUT=
ifeq ($(HOST),other)
	MODOUT=-J$(MODDIR)
endif

ifeq ($(HOST),apt)
	MODOUT=-module $(MODDIR)
endif

ifeq ($(HOST),oak)
	MODOUT=-module $(MODDIR)
endif

#--------------------------------------------------
# Rules
#--------------------------------------------------
all: dirs $(TARGET)
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -o $(TARGET).exe $^ $(LDFLAGS)
	if test -d $(INSTLDIR); then \
		: ; \
	else \
		mkdir $(INSTLDIR); \
	fi
	ln -sf $(PWD)/$(TARGET).exe $(INSTLDIR)

$(OBJDIR)/%.o:$(SRCDIR)/%.c
	$(CC) $(CFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(LINSRCDIR)/%.f90
	$(FC) $(FF2C) $(FFLAGS) $(OMP) $(FDFLAGS) $(MODOUT) -o $@ -c $<
#
dirs:
	if test -d $(OBJDIR); then \
		: ; \
	else \
		mkdir $(OBJDIR); \
	fi
	if test -d $(MODDIR); then \
		: ; \
	else \
		mkdir $(MODDIR); \
	fi


dep:
	$(FDEP) $(SRCS) -b $(OBJDIR)/ > $(DEPDIR)/makefile.d

clean:
	rm -f $(MODDIR)/*.mod
	rm -f $(OBJS)

#--------------------------------------------------
#-include $(wildcard $(DEPDIR)/*.d)
obj/ClassSys.o : src/ClassSys.f90
obj/CommonLibrary.o : src/CommonLibrary.f90 obj/ClassSys.o
obj/DefineOperators.o : src/DefineOperators.F90 obj/CommonLibrary.o
obj/HFInput.o : src/HFInput.F90 obj/ClassSys.o
obj/HFMain.o : src/HFMain.F90 obj/WriteOperator.o obj/MBPT.o obj/HartreeFock.o obj/Operators.o obj/ModelSpace.o obj/HFInput.o obj/CommonLibrary.o obj/Profiler.o
obj/HartreeFock.o : src/HartreeFock.F90 obj/CommonLibrary.o obj/Profiler.o obj/Operators.o obj/LinAlgLib.o
obj/MBPT.o : src/MBPT.F90 obj/CommonLibrary.o obj/Profiler.o obj/Operators.o obj/ModelSpace.o
obj/MPIFunction.o : src/MPIFunction.F90
obj/ModelSpace.o : src/ModelSpace.F90 obj/LinAlgLib.o obj/CommonLibrary.o obj/ClassSys.o obj/Profiler.o obj/SingleParticleState.o
obj/NOperators.o : src/NOperators.F90 obj/Profiler.o obj/ClassSys.o obj/DefineOperators.o obj/CommonLibrary.o obj/ModelSpace.o obj/LinAlgLib.o
obj/Operators.o : src/Operators.F90 obj/Profiler.o obj/DefineOperators.o obj/NOperators.o obj/ModelSpace.o
obj/Profiler.o : src/Profiler.F90 obj/MPIFunction.o obj/ClassSys.o
obj/SingleParticleState.o : src/SingleParticleState.F90 obj/ClassSys.o
obj/WriteOperator.o : src/WriteOperator.F90 obj/HFInput.o obj/Profiler.o obj/ClassSys.o obj/Operators.o
obj/LinAlgLib.o : LinAlgf90/src/LinAlgLib.f90 obj/MatVecComplex.o obj/MatVecDouble.o obj/MatVecSingle.o obj/MatrixComplex.o obj/MatrixDouble.o obj/MatrixSingle.o obj/VectorComplex.o obj/VectorDouble.o obj/VectorSingle.o obj/SingleDouble.o obj/LinAlgParameters.o
obj/LinAlgParameters.o : LinAlgf90/src/LinAlgParameters.f90
obj/MatVecComplex.o : LinAlgf90/src/MatVecComplex.f90 obj/MatrixComplex.o obj/VectorComplex.o obj/LinAlgParameters.o
obj/MatVecDouble.o : LinAlgf90/src/MatVecDouble.f90 obj/MatrixDouble.o obj/VectorDouble.o obj/LinAlgParameters.o
obj/MatVecSingle.o : LinAlgf90/src/MatVecSingle.f90 obj/MatrixSingle.o obj/VectorSingle.o obj/LinAlgParameters.o
obj/MatrixComplex.o : LinAlgf90/src/MatrixComplex.f90 obj/VectorComplex.o obj/LinAlgParameters.o
obj/MatrixDouble.o : LinAlgf90/src/MatrixDouble.f90 obj/VectorDouble.o obj/LinAlgParameters.o
obj/MatrixSingle.o : LinAlgf90/src/MatrixSingle.f90 obj/VectorSingle.o obj/LinAlgParameters.o
obj/SingleDouble.o : LinAlgf90/src/SingleDouble.f90 obj/MatrixDouble.o obj/VectorDouble.o obj/MatrixSingle.o obj/VectorSingle.o
obj/VectorComplex.o : LinAlgf90/src/VectorComplex.f90 obj/LinAlgParameters.o
obj/VectorDouble.o : LinAlgf90/src/VectorDouble.f90 obj/LinAlgParameters.o
obj/VectorSingle.o : LinAlgf90/src/VectorSingle.f90 obj/LinAlgParameters.o
