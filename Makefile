#--------------------------------------------------
# Make file for TTFcalc code
#--------------------------------------------------
TARGET=HartreeFock
INSTLDIR=~/Desktop/HFtest/

#--------------------------------------------------
# Default Parameters
#--------------------------------------------------
FDEP=makedepf90
CC=gcc
FC=gfortran
LDFLAGS=-llapack -lblas
OMP = -fopenmp
FFLAGS=-O3
CFLAGS=-O3
FDFLAGS=-Ddebug
#FDFLAGS+=-fbounds-check -Wall -fbacktrace -O -Wuninitialized

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

$(OBJDIR)/%.o:$(SRCDIR)/%.c
	$(CC) $(CFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -J$(SRCDIR) -o $@ -c $<  # for debug
#	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -J$(MODDIR) -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -J$(SRCDIR) -o $@ -c $< # for debug
#	$(FC) $(FFLAGS) $(OMP) $(FDFLAGS) -J$(MODDIR) -o $@ -c $<
$(OBJDIR)/%.o:$(LINSRCDIR)/%.f90
	$(FC) -ff2c $(FFLAGS) $(OMP) $(FDFLAGS) -J$(SRCDIR) -o $@ -c $< # for debug
#	$(FC) -ff2c $(FFLAGS) $(OMP) $(FDFLAGS) -J$(MODDIR) -o $@ -c $<

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
	rm -f $(TARGET).exe
	rm -f $(MODS)
	rm -f $(OBJS)

#--------------------------------------------------
-include $(wildcard $(DEPDIR)/*.d)
