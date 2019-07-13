##
##
## Introduction
## ============
##
## Makefile for iqist/src/capi directory.
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
## 05/10/2015 by li huang (created)
## 08/17/2015 by li huang (modified)
## 07/13/2019 by Jianzhou Zhao (last modified)
##
##

.SUFFIXES: .f90

include ../../make.inc

mod = ctqmc_api.o atomic_api.o

default: all

all: build_mod

build_mod: $(mod)

.f90.o:
	@echo [MPIF90] $@
	@$(MPIF90) $(F90FLAGS) -c $*.f90

clean:
	@echo "Cleaning solver capi"
	@rm -f *.mod
	@rm -f *.o

