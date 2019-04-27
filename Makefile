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
MODDIR = src
MPI=off
Host= $(shell if hostname|grep -q apt1; then echo apt; \
  elif hostname|grep -q oak; then echo oak; \
  elif hostname|grep -q cedar; then echo cedar; \
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
  FC=gfortran
  ifeq ($(OS), OSX)
    LFLAGS+= -I/usr/local/include -L/usr/local/lib
  endif
  LFLAGS+= -lblas -llapack -lgsl -lz
  FFLAGS=-O3
  CFLAGS=-O3
  FFLAGS+= -fopenmp
  FFLAGS+= -Dsingle_precision_three_body_file
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
  DFLAGS+=-DModelSpaceDebug
  #DFLAGS+=-DTwoBodyOperatorDebug
endif



#--------------------------------------------------
# Source Files
#--------------------------------------------------

SRCDIR = src
SRCDIR_LinAlg = submodule/LinAlgf90/src
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

SRCC_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.c)
SRCF77_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.f)
SRCF90_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.f90)
SRCF95_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.F90)
OBJC_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %c, %o, $(notdir $(SRCC_LinAlg))))
OBJF77_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %f, %o, $(notdir $(SRCF77_LinAlg))))
OBJF90_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90_LinAlg))))
OBJF95_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95_LinAlg))))
SRCS_LinAlg= $(SRCC_LinAlg) $(SRCF77_LinAlg) $(SRCF90_LinAlg) $(SRCF95_LinAlg)
OBJS_LinAlg= $(OBJC_LinAlg) $(OBJF77_LinAlg) $(OBJF90_LinAlg) $(OBJF95_LinAlg)

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


SRCS_ALL = $(SRCS) $(SRCS_LinAlg) $(SRCS_HF)
OBJS_ALL = $(OBJS) $(OBJS_LinAlg) $(OBJS_HF)

MODOUT=
ifeq ($(strip $(HOST)),other)
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

$(OBJDIR)/%.o:$(SRCDIR_LinAlg)/%.c
	$(CC) $(CFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_LinAlg)/%.f
	$(FC) $(FFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_LinAlg)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_LinAlg)/%.F90
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
	rm -f $(MODDIR)/*.mod
	rm -f $(OBJS_ALL)

#--------------------------------------------------
-include $(wildcard $(DEPDIR)/*.d)
