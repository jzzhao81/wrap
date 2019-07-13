!=========>=========>=========>=========>=========>=========>=========!
! compute fermi distribution
!=========>=========>=========>=========>=========>=========>=========!

REAL(kind=8) FUNCTION fermi(freq,beta)

  USE mod_param,ONLY : dp,dzero

  IMPLICIT NONE

  REAL(kind=dp), INTENT(in) :: freq,beta

  fermi=dzero

  IF(freq .LE. dzero) THEN
     fermi=1.d0/(dexp(freq*beta)+1.d0)
  ELSE
     fermi=dexp(-freq*beta)/(1.d0+dexp(-freq*beta))
  END IF

END FUNCTION fermi

! REAL(kind=8) FUNCTION gauss(ene, sgm)
!
!   USE mod_param, ONLY : dp, dzero, pi
!
!   IMPLICIT NONE
!   REAL(kind=dp), INTENT(in) :: ene, sgm
!   INTEGER, PARAMETER :: ntot = 5001
!   REAL(kind=dp), PARAMETER :: emin = -100.d0
!   REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: grid, val
!   INTEGER :: ii
!   REAL(kind=dp) :: dis
!
!   ALLOCATE(grid(ntot), val(ntot))
!
!   gauss = dzero
!   grid  = dzero
!
!   dis = (ene - emin)/REAL(ntot, dp)
!   DO ii = 1, ntot
!      grid(ii) = emin + REAL(ii-1,dp)*dis
!      val(ii)  = EXP(-grid(ii)**2.d0/2.d0)
!   ENDDO
!
!   DO ii = 1, ntot-1
!      gauss = gauss + (val(ii)+val(ii+1))*5.d-1*dis
!   ENDDO
!   gauss = 1.d0 - gauss / SQRT(2.d0*pi)
!
! END FUNCTION gauss


COMPLEX(kind=8) FUNCTION cfermi(freq,beta)

  USE mod_param,ONLY : dp,zzero,dzero

  IMPLICIT NONE

  COMPLEX(kind=dp), INTENT(in) :: freq
  REAL(kind=dp) , INTENT(in) :: beta

  cfermi=zzero

  IF(REAL(freq) .LE. dzero) THEN
     cfermi=1.d0/(EXP(freq*beta)+1.d0)
  ELSE
     cfermi=EXP(-freq*beta)/(1.d0+EXP(-freq*beta))
  END IF

END FUNCTION cfermi
