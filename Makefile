#--------------------------------------------------
# Make file for TTFcalc code
#--------------------------------------------------
TARGET=HartreeFock
INSTLDIR= $(HOME)/bin
EXEDIR=$(PWD)/exe
DEBUG_MODE=off
MODDIR = mod
Host= $(shell if hostname|grep -q apt1; then echo apt; \
  elif hostname|grep -q oak; then echo oak; \
  elif hostname|grep -q cedar; then echo cedar; \
  else echo other; fi)
HOST=$(strip $(Host))
$(info Host:$(HOST))

FDEP=
FC=
LFLAGS=
FFLAGS=
DFLAGS=
#--------------------------------------------------
# Default Parameters
#--------------------------------------------------
ifeq ($(HOST),other)
  FDEP=makedepf90
  FC=gfortran
  LFLAGS= -I/usr/local/include -L/usr/local/lib -llapack -lblas -lgsl -lz
  FFLAGS = -O3 -Dsingle_precision_three_body_file -fopenmp
  #FFLAGS += -ff2c
  DFLAGS=
  ifeq ($(DEBUG_MODE),on)
    DFLAGS+=-DModelSpaceDebug
    #FDFLAGS+=-DTwoBodyOperatorDebug
    DFLAGS+=-fbounds-check -Wall -fbacktrace -O -Wuninitialized
  endif
endif

ifeq ($(HOST),apt)
  #FDEP=
  FC=ifort
  LFLAGS=-mkl -lgsl -lz
  FFLAGS=-O3 -no-ipo -static -Dsingle_precision_three_body_file -openmp
  ifeq ($(DEBUG_MODE),on)
    #DFLAGS+=-DModelSpaceDebug
    #FDFLAGS+=-DTwoBodyOperatorDebug
    DFLAGS+=-check all
  endif
endif

#-----------------------------
# oak (oak.arc.ubc.ca)
# -heap-arrays option is needed for built in transpose function with
#  large dimension matrix
#-----------------------------
ifeq ($(strip $(HOST)),oak)
  FDEP=
  FC=ifort
  MKL=-L$(MKLROOT)/lib/ -L$(MKLROOT)/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -Wl,-rpath,$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/../compiler/lib/
  LFLAGS=$(MKL) -lgsl -lz
  FFLAGS=-O3 -heap-arrays -qopenmp
  FFLAGS += -Dsingle_precision_three_body_file
  ifeq ($(DEBUG_MODE),on)
    DFLAGS+=-check all
  endif
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
	$(FC) $(FFLAGS) $(DFLAGS) -o $(TARGET).exe $^ $(LFLAGS)
	if test -d $(EXEDIR); then \
		: ; \
	else \
		mkdir -p $(EXEDIR); \
	fi
	mv $(TARGET).exe $(EXEDIR)
	@echo "#####################################################################################"
	@echo "To complete the installation, do 'make install'."
	@echo "Edit '$(PWD)/exe/run_hf_mbpt.py' and excecute."
	@echo "#####################################################################################"

$(OBJDIR)/%.o:$(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(LINSRCDIR)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<
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

install:
	ln -sf $(EXEDIR)/$(TARGET).exe $(INSTLDIR)

clean:
	rm -f $(MODDIR)/*.mod
	rm -f $(OBJS)

#--------------------------------------------------
-include $(wildcard $(DEPDIR)/*.d)
