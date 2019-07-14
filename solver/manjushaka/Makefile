##
##
## Introduction
## ============
##
## Makefile for iqist/src/ctqmc/manjushaka directory.
##
## Usage
## =====
##
## 'make', 'make lib', 'make clean', 'make clean-dat' or 'make clean-all'
##
## Author
## ======
##
## This building system is designed, created, implemented, and maintained by
##
## Li Huang // email: lihuang.dmft@gmail.com
##
## History
## =======
##
## 02/28/2015 by li huang (created)
## 08/17/2015 by li huang (modified)
## 07/13/2019 by Jianzhou Zhao (last modified)
##
##

.SUFFIXES: .f90

include ../../make.inc

modc = ctqmc_control.o ctqmc_context.o
dmft = ctqmc_dmft.o
core = ctqmc_solver.o
lev1 = ctqmc_flavor.o ctqmc_update.o
lev2 = ctqmc_record.o ctqmc_status.o ctqmc_stream.o ctqmc_util.o
lev3 = ctqmc_dump.o ctqmc_print.o
main = ctqmc_main.o
open = ctqmc_open.o

CAPI = ctqmc_api.o
CLIB = libMM.a

objects = $(modc) $(dmft) $(core) $(lev1) $(lev2) $(lev3)
exe_obj = $(objects) $(main) $(CLIB)
lib_obj = $(objects) $(open) $(CAPI)

default: all

all: exe

api:
	@echo [CP] $@
	@cp ../capi/capi.mod .
	@cp ../capi/ctqmc_api.o .

cssl:
	@echo [CP] $@
	@cp ../base/libMM.a .

csml:
	@echo [CP] $@
	@cp ../base/constants.mod .
	@cp ../base/mmpi.mod .
	@cp ../base/spring.mod .
	@cp ../base/stack.mod .
	@cp ../base/parser.mod .
	@cp ../base/linkedlist.mod .

exe: api cssl csml ctqmc-exe

ctqmc-exe: $(exe_obj)
	@echo [LD] $@
	@$(LD) $(exe_obj) -o ctqmc $(LDFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS)

lib: api cssl csml ctqmc-lib
ctqmc-lib: $(lib_obj)
	@echo [AR] $@
	@$(AR) libctqmc.a $(lib_obj)

pylib: lib
	@echo [F2PY] $@
	@cpp -E -DPYAPI ctqmc_open.f90 > ctqmc_open_py.f90
	@$(F2PY) -c --f90exec=$(F2PYC) -m pyiqist -L. -lctqmc -lMM $(F2PYL) ctqmc_open_py.f90
	@rm -fr *.dSYM ctqmc_open_py.f90

.f90.o:
	@echo [MPIF90] $@
	@$(CPP) $(CPPFLAGS) $*.f90 > $*_cpp.f90
	@$(MPIF90) $(F90FLAGS) -o $*.o -c $*_cpp.f90
	@rm -f $*_cpp.f90

clean:
	@echo "Cleaning solver gardenia"
	@rm -f *.mod
	@rm -f *.o
	@rm -f *.so
	@rm -f ctqmc
	@rm -f libctqmc.a
	@rm -f libMM.a
