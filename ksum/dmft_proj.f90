!=========>=========>=========>=========>=========>=========>=========!
! project self-energy to bloch state
!=========>=========>=========>=========>=========>=========>=========!

SUBROUTINE dmft_proj_ks(nbnd,norb,mtrx_in,mtrx_out,proj)

  USE mod_param, ONLY : dp, zzero
  USE mod_mtrx

  IMPLICIT NONE

  INTEGER,INTENT(in) :: nbnd, norb
  COMPLEX(kind=dp),INTENT(in) :: mtrx_in(norb,norb)
  COMPLEX(kind=dp),INTENT(in) :: proj(nbnd,norb)
  COMPLEX(kind=dp),INTENT(out) :: mtrx_out(nbnd,nbnd)
  COMPLEX(kind=dp), ALLOCATABLE :: tmp(:,:)

  ALLOCATE(tmp(nbnd,norb))

  mtrx_out=zzero
  tmp=zzero

  CALL mtrx_mult(proj,mtrx_in,tmp,'N','N')

  CALL mtrx_mult(tmp,proj,mtrx_out,'N','C')

  DEALLOCATE(tmp)

  RETURN

END SUBROUTINE dmft_proj_ks


!=========>=========>=========>=========>=========>=========>=========!
! project Green's function to local basis
!=========>=========>=========>=========>=========>=========>=========!
SUBROUTINE dmft_proj_loc(norb,nbnd,mtrx_in,mtrx_out,proj)

  USE mod_param, ONLY : dp, zzero
  USE mod_mtrx

  IMPLICIT NONE

  INTEGER,INTENT(in) :: norb, nbnd
  COMPLEX(kind=dp),INTENT(in) :: mtrx_in(nbnd,nbnd)
  COMPLEX(kind=dp),INTENT(in) :: proj(nbnd,norb)
  COMPLEX(kind=dp),INTENT(out) :: mtrx_out(norb,norb)
  COMPLEX(kind=dp),ALLOCATABLE :: tmp(:,:)

  ALLOCATE(tmp(norb,nbnd))

  mtrx_out=zzero
  tmp=zzero

  CALL mtrx_mult(proj,mtrx_in,tmp,'C','N')

  CALL mtrx_mult(tmp,proj,mtrx_out,'N','N')

  DEALLOCATE(tmp)

  RETURN

END SUBROUTINE dmft_proj_loc
