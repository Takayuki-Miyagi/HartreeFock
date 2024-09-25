#--------------------------------------------------
# Make file for HartreeFock code
#
# When we use ifort,
# -heap-arrays option is needed for built in transpose function with
#  large dimension matrix
#--------------------------------------------------
TARGET=HartreeFock
INSTLDIR=$(HOME)/bin
EXEDIR=$(PWD)/exe
MODDIR = mod
MPI=off
Host= $(shell if hostname|grep -q apt1; then echo apt; \
  elif hostname|grep -q oak; then echo oak; \
  elif hostname|grep -q cedar; then echo cedar; \
  elif hostname|grep -q strongint; then echo strongint; \
  else echo other; fi)
HOST=$(strip $(Host))
DEBUG_MODE=off

OS = Linux
ifneq (,$(findstring arwin,$(shell uname)))
  OS = OSX
endif
$(info HOST is $(HOST).)
$(info OS is $(OS).)
$(info Debug mode is $(DEBUG_MODE).)

FDEP=
FC=
LFLAGS=  # library
FFLAGS=  # option
DFLAGS=  # debug
LINT=    # 8-byte integer

#--------------------------------------------------
# Default Parameters
#--------------------------------------------------
ifeq ($(strip $(HOST)),other)
  FDEP=makedepf90
  FC=gfortran -fallow-argument-mismatch # if GCC version > 10.xx
  ifeq ($(OS), OSX)
    LFLAGS+= -I/usr/local/include -L/usr/local/lib
    LFLAGS+= -I/opt/homebrew/include -L/opt/homebrew/lib
    LFLAGS+= -L/opt/homebrew/opt/openblas/lib
  endif
  LFLAGS+= -lblas -llapack -lgsl -lz
  FFLAGS=-O3
  CFLAGS=-O3
  FFLAGS+= -fopenmp
  FFLAGS+= -Dsingle_precision_three_body_file
  #FFLAGS+= -Dhalf_precision_three_body_file
  #FFLAGS+= -ff2c # for dot product (LinAlgf90)
  ifeq ($(DEBUG_MODE),on)
    DFLAGS+=-Wall -pedantic -fbounds-check -O -Wuninitialized -fbacktrace
    #FDFLAGS+=-ffpe-trap=invalid,zero,overflow # Note: gsl larguerre signal
    ifneq ($(OS), OSX)
      DFLAGS+= -pg -g
    endif
  endif
endif

#-----------------------------
# apt1
#-----------------------------
ifeq ($(strip $(HOST)),apt)
  FC=ifort
  LFLAGS+= -mkl -lgsl -lz
  FFLAGS=-O3 -heap-arrays -static
  FFLAGS+= -openmp
  FFLAGS+= -Dsingle_precision_three_body_file
  ifeq ($(DEBUG_MODE),on)
    DFLAGS+=-check all
  endif
endif

#-----------------------------
# oak (oak.arc.ubc.ca)
#-----------------------------
ifeq ($(strip $(HOST)),oak)
  FC=ifort
  MKL=-L$(MKLROOT)/lib/ -L$(MKLROOT)/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -Wl,-rpath,$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/../compiler/lib/
  LFLAGS+= $(MKL) -lgsl -lz
  FFLAGS=-O3 -heap-arrays
  FFLAGS+= -qopenmp
  FFLAGS+= -Dsingle_precision_three_body_file
  ifeq ($(DEBUG_MODE),on)
    DFLAGS+=-check all
  endif
endif

#-----------------------------
# cedar
#-----------------------------
ifeq ($(strip $(HOST)),cedar)
  FC=ifort
  LFLAGS+= -lmkl -lgsl -lz
  FFLAGS=-O3 -heap-arrays
  FFLAGS+= -qopenmp
  FFLAGS+= -Dsingle_precision_three_body_file
  ifeq ($(DEBUG_MODE),on)
    DFLAGS+=-check all
  endif
endif

ifeq ($(DEBUG_MODE),on)
  #DFLAGS+=-DModelSpaceDebug
  #DFLAGS+=-DTwoBodyOperatorDebug
endif

#--------------------------------------------------
# Strongint cluster
#--------------------------------------------------
ifeq ($(strip $(HOST)),strongint)
  FC=gfortran
  #LFLAGS+= -I/usr/local/include -L/usr/local/lib -I/usr/lib64/gfortran/modules
  LFLAGS+= -lblas -llapack -lz -lgsl -lgslcblas
  FFLAGS= -O3
  FFLAGS+= -fopenmp -fdec-math
  ifeq ($(strip $(TARGET)),HartreeFock_half)
    FFLAGS+= -Dhalf_precision_three_body_file
  endif
  DFLAGS+= -Wall
  ifeq ($(DEBUG_MODE),on)
    DFLAGS+= -pedantic -fbounds-check -O -Wuninitialized -fbacktrace
    #FDFLAGS+=-ffpe-trap=invalid,zero,overflow # Note: gsl larguerre signal
    ifneq ($(OS), OSX)
      DFLAGS+= -pg -g
    endif
  endif
  LINT= -fdefault-integer-8
endif


#--------------------------------------------------
# Source Files
#--------------------------------------------------

SRCDIR = src
SRCDIR_myfort = submodule/myfort/src
SRCDIR_HF = main
DEPDIR = .
OBJDIR = obj

SRCS=
OBJS=
MODS=

SRCC:=$(wildcard $(SRCDIR)/*.c)
SRCF77:=$(wildcard $(SRCDIR)/*.f)
SRCF90:=$(wildcard $(SRCDIR)/*.f90)
SRCF95:=$(wildcard $(SRCDIR)/*.F90)
OBJC:=$(addprefix $(OBJDIR)/, $(patsubst %c, %o, $(notdir $(SRCC))))
OBJF77:=$(addprefix $(OBJDIR)/, $(patsubst %f, %o, $(notdir $(SRCF77))))
OBJF90:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90))))
OBJF95:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95))))
SRCS= $(SRCC) $(SRCF77) $(SRCF90) $(SRCF95)
OBJS= $(OBJC) $(OBJF77) $(OBJF90) $(OBJF95)

SRCC_myfort:=$(wildcard $(SRCDIR_myfort)/*.c)
SRCF77_myfort:=$(wildcard $(SRCDIR_myfort)/*.f)
SRCF90_myfort:=$(wildcard $(SRCDIR_myfort)/*.f90)
SRCF95_myfort:=$(wildcard $(SRCDIR_myfort)/*.F90)
OBJC_myfort:=$(addprefix $(OBJDIR)/, $(patsubst %c, %o, $(notdir $(SRCC_myfort))))
OBJF77_myfort:=$(addprefix $(OBJDIR)/, $(patsubst %f, %o, $(notdir $(SRCF77_myfort))))
OBJF90_myfort:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90_myfort))))
OBJF95_myfort:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95_myfort))))
SRCS_myfort= $(SRCC_myfort) $(SRCF77_myfort) $(SRCF90_myfort) $(SRCF95_myfort)
OBJS_myfort= $(OBJC_myfort) $(OBJF77_myfort) $(OBJF90_myfort) $(OBJF95_myfort)

SRCC_HF:=$(wildcard $(SRCDIR_HF)/*.c)
SRCF77_HF:=$(wildcard $(SRCDIR_HF)/*.f)
SRCF90_HF:=$(wildcard $(SRCDIR_HF)/*.f90)
SRCF95_HF:=$(wildcard $(SRCDIR_HF)/*.F90)
OBJC_HF:=$(addprefix $(OBJDIR)/, $(patsubst %c, %o, $(notdir $(SRCC_HF))))
OBJF77_HF:=$(addprefix $(OBJDIR)/, $(patsubst %f, %o, $(notdir $(SRCF77_HF))))
OBJF90_HF:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90_HF))))
OBJF95_HF:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95_HF))))
SRCS_HF= $(SRCC_HF) $(SRCF77_HF) $(SRCF90_HF) $(SRCF95_HF)
OBJS_HF= $(OBJC_HF) $(OBJF77_HF) $(OBJF90_HF) $(OBJF95_HF)


SRCS_ALL = $(SRCS) $(SRCS_myfort) $(SRCS_HF)
OBJS_ALL = $(OBJS) $(OBJS_myfort) $(OBJS_HF)

MODOUT=
ifeq ($(strip $(HOST)),other)
  MODOUT=-J$(MODDIR)
endif

ifeq ($(strip $(HOST)),strongint)
  MODOUT=-J$(MODDIR)
endif

ifeq ($(strip $(HOST)),apt)
  MODOUT=-module $(MODDIR)
endif

ifeq ($(strip $(HOST)),oak)
  MODOUT=-module $(MODDIR)
endif

ifeq ($(strip $(HOST)),cedar)
  MODOUT=-module $(MODDIR)
endif

#--------------------------------------------------
# Rules
#--------------------------------------------------
all: dirs $(TARGET)
$(TARGET): $(OBJS_ALL)
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


$(OBJDIR)/%.o:$(SRCDIR)/%.c
	$(CC) $(CFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f
	$(FC) $(FFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<

$(OBJDIR)/%.o:$(SRCDIR_myfort)/%.c
	$(CC) $(CFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_myfort)/%.f
	$(FC) $(FFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_myfort)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_myfort)/%.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<

$(OBJDIR)/%.o:$(SRCDIR_HF)/%.c
	$(CC) $(CFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_HF)/%.f
	$(FC) $(FFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_HF)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_HF)/%.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<

dep:
	$(FDEP) $(SRCS_ALL) -b $(OBJDIR)/ > $(DEPDIR)/makefile.d

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

install:
	ln -sf $(EXEDIR)/$(TARGET).exe $(INSTLDIR)
	@printf "####################################################################\n"
	@printf " make sure that '$(INSTLDIR)' is included in PATH \n"
	@printf "####################################################################\n"

clean:
	rm -f $(TARGET).exe
	rm -rf $(OBJDIR)
	rm -rf $(MODDIR)

#--------------------------------------------------
-include $(wildcard $(DEPDIR)/*.d)
