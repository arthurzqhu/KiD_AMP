
#COMPILER=gfortran# Gnu 
#COMPILER=pgf90#    Portland
COMPILER=ifort#     Intel
#COMPILER=nagfor#      NAG 

VPATH = ./:./UM_source

#NCPATH=/usr/local/netcdf/$(COMPILER)/
#NCPATH=/home/h01/lemdev/netcdf/$(COMPILER)/
NCPATH=/share/apps/netcdf-c-4.7.4-intel/

TESTING=

ifeq ($(COMPILER), gfortran)
FPP= -x f95-cpp-input
ifdef TESTING
OPTS=-g -O0 -W -Wall -pedantic-errors -fbounds-check -fbacktrace \
-ffpe-trap=invalid,zero,overflow -frange-check #-Werror
else
OPTS=-O3 -pg -g
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
LNCDF=-L$(NCPATH)/lib -lnetcdf -lnetcdff
LDFLAGS=$(INCDF) $(LNCDF)
#LDFLAGS=-lnetcdf -lnetcdff

EXECUTABLE=$(MODEL)_$(CASE).exe

F2PY=/usr/bin/f2py 
F2PY_OPTS=--fcompiler-exec=$(COMPILER) $(LDFLAGS) 
