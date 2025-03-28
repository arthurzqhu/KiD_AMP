#
# Makefile for KiD model
#
#
#

MODEL=KiD

# compiler options
include compiler_options.inc
include defines.inc

F90D=$(F90) $(FPP) $(DEFS) -DCOMPILER=$(COMPILER) # Pass through defines including compiler

# The dependencies are autogenerated in mkdepend.py
include .includes

ifdef DEF_UM
VPATHO=UM_source/*.o
OBJSF90_UM=$(patsubst %.f90,$(OBJDIR)/%.o,$(wildcard UM_source/*.f90)) 
endif

ifdef DEF_4A
OBJSF90_4A=$(patsubst %.f90,$(OBJDIR)/%.o,$(wildcard 4A_code/*.f90)) 
endif

OBJSF=$(patsubst %.f,$(OBJDIR)/%.o,$(wildcard *.f))
OBJSF90=$(patsubst %.f90,$(OBJDIR)/%.o,$(wildcard *.f90))
OBJS=$(OBJSF) $(OBJSF90)


all: defs $(EXECUTABLE)

defs:
	touch case.used
	[ -d $(OBJDIR)/UM_source ] || [ "$(OBJSF90_UM)" == " " ] || mkdir -p "$(OBJDIR)/UM_source";
	[ -d $(OBJDIR)/4A_code ] || [ "$(OBJSF90_4A)" == " " ] || mkdir -p "$(OBJDIR)/4A_code";
	if [[ "$(CASE)" != `cat case.used` ]]; then \
	touch parameters.f90; \
	echo $(CASE) > case.used; \
	fi

$(OBJSF90_UM): $(OBJDIR)/%.o: %.f90 
	[ -d $(OBJDIR)/UM_source ] || mkdir -p "$(OBJDIR)/UM_source"; \
	$(F90D) -c $< -o $@ $(INCDF)

$(OBJSF90_4A): $(OBJDIR)/%.o: %.f90 
	[ -d $(OBJDIR)/4A_code ] || mkdir -p "$(OBJDIR)/4A_code"; \
	$(F90D) -c $< -o $@ $(INCDF)

$(OBJSF90): $(OBJDIR)/%.o: %.f90 
	$(F90D) -cpp -c $< -o $@ $(INCDF) $(FKBINC)

$(OBJSF): $(OBJDIR)/%.o: %.f 
	$(F90D)  -c $< -o $@ $(INCDF)

$(EXECUTABLE): $(OBJS) $(OBJSF90_UM) $(OBJSF90_4A)
	$(LD) -o $(EXECDIR)/$(EXECUTABLE) $(OBJDIR)/*.o $(OBJSF90_4A) $(OBJSF90_UM) $(LDFLAGS)
