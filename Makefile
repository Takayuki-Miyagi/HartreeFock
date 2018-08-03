#--------------------------------------------------
# Make file for TTFcalc code
#--------------------------------------------------
TARGET=HartreeFock
INSTLDIR=~/Desktop/HFtest/
HOST= $(shell if hostname|grep -q apt1; then echo apt; else echo other; fi)
$(info Host:$(HOST))

#--------------------------------------------------
# Default Parameters
#--------------------------------------------------
ifeq ($(HOST),other)
#	FDEP=makedepf90
	CC=gcc
	FC=gfortran
	LDFLAGS=-llapack -lblas
	OMP = -fopenmp
	FFLAGS=-O3 -Dsingle_precision
	CFLAGS=-O3
	FF2C= -ff2c
	#FDFLAGS=-Ddebug
	#FDFLAGS+=-fbounds-check -Wall -fbacktrace -O -Wuninitialized
endif

ifeq ($(HOST),apt)
#	FDEP=
	CC=gcc
	FC=ifort
	LDFLAGS=-mkl
	OMP = -openmp
	FFLAGS=-O3 -no-ipo -static -Dsingle_precision
	CFLAGS=-O3
	FF2C=
	#FDFLAGS=-Ddebug
	#FDFLAGS+=-check-all
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

MODDIR = mod
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
	MODOUT=-module$(MODDIR)
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
	cp $(TARGET).exe $(INSTLDIR)
	cp *.py $(INSTLDIR)

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
$(OBJDIR)/wallclock.o : $(SRCDIR)/wallclock.c
$(OBJDIR)/common_library.o : $(SRCDIR)/common_library.f90 $(OBJDIR)/class_sys.o
$(OBJDIR)/HFSolver.o : $(SRCDIR)/HFSolver.f90 $(OBJDIR)/Optimizer.o $(OBJDIR)/LinAlgLib.o $(OBJDIR)/VectorDouble.o $(OBJDIR)/MatrixDouble.o $(OBJDIR)/NormalOrdering.o $(OBJDIR)/ScalarOperator.o $(OBJDIR)/Read3BME.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/HFInput.o $(OBJDIR)/MPIFunction.o
$(OBJDIR)/class_sys.o : $(SRCDIR)/class_sys.f90
$(OBJDIR)/Optimizer.o : $(SRCDIR)/Optimizer.f90 $(OBJDIR)/LinAlgLib.o $(OBJDIR)/VectorDouble.o $(OBJDIR)/MatrixDouble.o
$(OBJDIR)/MBPT3.o : $(SRCDIR)/MBPT3.F90 $(OBJDIR)/ScalarOperator.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/HFInput.o $(OBJDIR)/common_library.o
$(OBJDIR)/NormalOrdering.o : $(SRCDIR)/NormalOrdering.F90 $(OBJDIR)/MPIFunction.o $(OBJDIR)/Read3BME.o $(OBJDIR)/common_library.o $(OBJDIR)/LinAlgLib.o $(OBJDIR)/VectorDouble.o $(OBJDIR)/MatrixDouble.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/ScalarOperator.o $(OBJDIR)/HFInput.o
$(OBJDIR)/WriteOperator.o : $(SRCDIR)/WriteOperator.F90 $(OBJDIR)/ScalarOperator.o $(OBJDIR)/Read3BME.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/HFInput.o $(OBJDIR)/MPIFunction.o $(OBJDIR)/class_stopwatch.o $(OBJDIR)/class_sys.o
$(OBJDIR)/MPIFunction.o : $(SRCDIR)/MPIFunction.F90
$(OBJDIR)/ModelSpace.o : $(SRCDIR)/ModelSpace.F90 $(OBJDIR)/LinAlgLib.o $(OBJDIR)/VectorDouble.o $(OBJDIR)/MatrixDouble.o $(OBJDIR)/MPIFunction.o $(OBJDIR)/HFInput.o $(OBJDIR)/common_library.o
$(OBJDIR)/ScalarOperator.o : $(SRCDIR)/ScalarOperator.F90 $(OBJDIR)/class_stopwatch.o $(OBJDIR)/common_library.o $(OBJDIR)/Read3BME.o $(OBJDIR)/LinAlgLib.o $(OBJDIR)/VectorDouble.o $(OBJDIR)/MatrixDouble.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/HFInput.o $(OBJDIR)/class_sys.o
$(OBJDIR)/HFInput.o : $(SRCDIR)/HFInput.F90 $(OBJDIR)/class_sys.o $(OBJDIR)/MPIFunction.o
$(OBJDIR)/class_stopwatch.o : $(SRCDIR)/class_stopwatch.F90 $(OBJDIR)/MPIFunction.o
$(OBJDIR)/Read3BME.o : $(SRCDIR)/Read3BME.F90 $(OBJDIR)/class_sys.o $(OBJDIR)/common_library.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/MPIFunction.o $(OBJDIR)/HFInput.o
$(OBJDIR)/HFMain.o : $(SRCDIR)/HFMain.F90 $(OBJDIR)/class_stopwatch.o $(OBJDIR)/MPIFunction.o $(OBJDIR)/MBPT3.o $(OBJDIR)/WriteOperator.o $(OBJDIR)/HFSolver.o $(OBJDIR)/NormalOrdering.o $(OBJDIR)/ScalarOperator.o $(OBJDIR)/Read3BME.o $(OBJDIR)/ModelSpace.o $(OBJDIR)/HFInput.o $(OBJDIR)/common_library.o
$(OBJDIR)/LinAlgLib.o : $(LINSRCDIR)/LinAlgLib.f90 $(OBJDIR)/Parameters.o $(OBJDIR)/MatVecComplex.o $(OBJDIR)/MatVecDouble.o $(OBJDIR)/MatrixComplex.o $(OBJDIR)/MatrixDouble.o $(OBJDIR)/VectorComplex.o $(OBJDIR)/VectorDouble.o
$(OBJDIR)/MatVecDouble.o : $(LINSRCDIR)/MatVecDouble.f90 $(OBJDIR)/MatrixDouble.o $(OBJDIR)/VectorDouble.o
$(OBJDIR)/VectorComplex.o : $(LINSRCDIR)/VectorComplex.f90 $(OBJDIR)/Parameters.o
$(OBJDIR)/MatrixDouble.o : $(LINSRCDIR)/MatrixDouble.f90 $(OBJDIR)/VectorDouble.o
$(OBJDIR)/VectorDouble.o : $(LINSRCDIR)/VectorDouble.f90 $(OBJDIR)/Parameters.o
$(OBJDIR)/MatrixComplex.o : $(LINSRCDIR)/MatrixComplex.f90 $(OBJDIR)/VectorComplex.o
$(OBJDIR)/Parameters.o : $(LINSRCDIR)/Parameters.f90
$(OBJDIR)/MatVecComplex.o : $(LINSRCDIR)/MatVecComplex.f90 $(OBJDIR)/MatrixComplex.o $(OBJDIR)/VectorComplex.o
