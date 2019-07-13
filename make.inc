# make.inc.

# compilation rules

.SUFFIXES :
.SUFFIXES : .o .c .f .f90

MAKE           = make

# IFLAGS = how to locate directories with *.h or *.f90 file to be included
#          typically -I$(MKLROOT)/include -I/some/other/directory/
IFLAGS         = -I$(MKLROOT)/include -I/opt/local/include

MPIF90         = mpif90
F90            = gfortran

CPP            = cpp
CPPFLAGS       = -P -traditional $(DFLAGS) $(IFLAGS)

F90FLAGS       = -g -Wall $(IFLAGS)

LD             = mpif90
LDFLAGS        = -g
BLAS_LIBS      =
LAPACK_LIBS    = /opt/local/lib/libopenblas.a

AR             = ar
ARFLAGS        = ruv

F2PY           = f2py
F2PYC          = $(F90)
F2PYL          = $(LAPACK_LIBS) $(BLAS_LIBS)
