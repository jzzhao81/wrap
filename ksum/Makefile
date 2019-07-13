# @Author: Jianzhou Zhao <jzzhao>
# @Date:   2018-11-02T16:20:00+01:00
# @Email:  jianzhou.zhao@gmail.com
# @Last modified by:   jzzhao
# @Last modified time: 2019-03-15T16:29:43+01:00

include ../make.inc

OBJ=    mod_param.o mod_dmft.o mod_eigen.o mod_mpi.o mod_mtrx.o dmft_read.o \
        dmft_finalize.o dmft_mu.o dmft_fermi.o dmft_green.o dmft_proj.o \
        dmft_init.o dmft_dm.o math_spline.o dmft_sechk.o dmft_output.o \
	      lda_occ.o dmft_mix.o dmft_main.o dmft_occup.o

ksum: $(OBJ)
	@echo [LD] $@
	@rm -f ksum
	@$(LD) $(LDFLAGS) -o ksum $(OBJ) $(LAPACK_LIBS) $(BLAS_LIBS)

clean:
	@echo "Cleaning ksum"
	@rm -f *.o *.mod

$(OBJ):%.o:%.f90
	@echo [MPIF90] $@
	@$(CPP) $(CPPFLAGS) $*.f90 > $*_cpp.f90
	@$(MPIF90) $(F90FLAGS) -o $*.o -c $*_cpp.f90
	@rm -f $*_cpp.f90
