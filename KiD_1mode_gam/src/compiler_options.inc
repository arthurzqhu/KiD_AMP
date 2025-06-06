
# COMPILER=/usr/local/bin/gfortran #Gnu 
COMPILER=gfortran# Gnu 
#COMPILER=pgf90#    Portland
#COMPILER=ifort#     Intel
#COMPILER=nagfor#      NAG 

VPATH = ./:./UM_source

#NCPATH=/usr/local/netcdf/$(COMPILER)/
#NCPATH=/home/h01/lemdev/netcdf/$(COMPILER)/
# NCPATH=/share/apps/netcdf-c-4.7.4-intel/
# FKBPATH=/Users/arthurhu/github/FKB/build
FKBPATH=/home/arthurhu/FKB/build
NCPATH=/home/arthurhu/.local/Cellar/netcdf-fortran/4.6.1_2
NCPATHC=/home/arthurhu/.local/Cellar/netcdf/4.9.3
# NCPATH=/usr/local/Cellar/netcdf-fortran/4.6.1
# NCPATHC=/usr/local/Cellar/netcdf/4.9.2_1

TESTING=

ifeq ($(COMPILER), gfortran)
FPP= -x f95-cpp-input
ifdef TESTING
OPTS=-g -O0 -W -Wall -fbounds-check -fbacktrace \
-ffpe-trap=invalid,zero,overflow -frange-check #-Werror -fsanitize=address
else
OPTS=-O3 -pg -g -fallow-argument-mismatch -fno-backtrace
endif
endif

ifeq ($(COMPILER), nagfor)
FPP= -fpp
ifdef TESTING
OPTS=-C=all  -info -g -gline
else
OPTS=-O3 -pg -g
endif
endif

ifeq ($(COMPILER), pgf90)
FPP= -Mpreprocess
ifdef TESTING
OPTS=-g -Mbounds -Mchkstk -Mchkptr -Mchkfpstk -Minform=inform -C 
else
OPTS=-O3 -pg -g
endif
endif

ifeq ($(COMPILER), ifort)
FPP= -fpp
ifdef TESTING
OPTS=-g -O0 -check all -warn all,nodec,interfaces \
 -traceback -fpe0
else
OPTS=-O3 -traceback -pg -g
endif
endif

F90=$(COMPILER) $(OPTS)
LD=$(F90)   # Linker

INCDF=-I$(NCPATH)/include
LNCDF=-L$(NCPATH)/lib -lnetcdff -L$(NCPATHC)/lib -lnetcdf 
FKBINC=-I$(FKBPATH)/include
FKBLIB=-L$(FKBPATH)/lib -lneural
LDFLAGS=$(FKBINC) $(FKBLIB) $(INCDF) $(LNCDF) 
#LDFLAGS=-lnetcdf -lnetcdff

EXECUTABLE=$(MODEL)_$(CASE).exe

F2PY=/usr/bin/f2py 
F2PY_OPTS=--fcompiler-exec=$(COMPILER) $(LDFLAGS) 
