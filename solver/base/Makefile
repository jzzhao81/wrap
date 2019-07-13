##
##
## Introduction
## ============
##
## Makefile for iqist/src/base directory.
##
## Usage
## =====
##
## 'make' or 'make clean'
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
## 09/18/2014 by li huang (created)
## 08/17/2015 by li huang (modified)
## 07/13/2019 by Jianzhou Zhao (last modified)
##
##

.SUFFIXES: .f90

include ../../make.inc

mod = m_constants.o m_leja.o m_linkedlist.o m_mpi.o m_parser.o m_skynet.o m_sparse.o m_spring.o m_stack.o
sub = s_error.o s_fft.o s_integrator.o s_matrix.o s_spline.o s_util.o s_vector.o
objects = $(mod) $(sub)

default: all

all: build-mod build-sub build-lib

build-mod: $(mod)
build-sub: $(sub)
build-lib: $(objects)
	@echo [AR] $@
	@$(AR) $(ARFLAGS) libMM.a $(objects)

.f90.o:
	@echo [MPIF90] $@
	@$(CPP) $(CPPFLAGS) $*.f90 > $*_cpp.f90
	@$(MPIF90) $(F90FLAGS) -o $*.o -c $*_cpp.f90
	@rm -f $*_cpp.f90

clean:
	@echo "Cleaning solver base"
	@rm -f *.mod
	@rm -f *.o
	@rm -f libMM.a

