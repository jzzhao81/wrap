!=========>=========>=========>=========>=========>=========>=========!
! calculate Total occupation
!=========>=========>=========>=========>=========>=========>=========!

SUBROUTINE dmft_occup(mu,eigs,occ)

  USE mod_param,ONLY : dp,dzero,zzero
  USE mod_dmft, ONLY : beta,wk,freq_log,nlog,num_bnd
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  !variable definition
  !===================================================================!
  ! input chemical potential
  REAL(kind=dp),INTENT(in) :: mu
  ! input eigen value of ham
  COMPLEX(kind=dp), INTENT(in) :: eigs(num_bnd,nlog,nkpp(iproc))
  ! output total electron
  COMPLEX(kind=dp), INTENT(out) :: occ

  ! loop parameter
  INTEGER :: ikp, ibnd, ifreq
  ! iwn = iw
  COMPLEX(kind=dp) :: iwn
  ! fermi distribution
  REAL(kind=dp) :: dfermi, ene
  COMPLEX(kind=dp) :: domega

  ! external function of fermi distribution
  REAL(kind=dp), EXTERNAL :: fermi
  !===================================================================!

  ! initialize variable
  occ = zzero
  ene = dzero

  ! calculate density matrix
  !===================================================================!
  DO ikp=1,nkpp(iproc)

     dfermi=dzero
     domega=zzero
     ! Fermi distribution is frequence independent
     DO ibnd=1,num_bnd
        ene = REAL(eigs(ibnd,nlog,ikp)) - mu
        dfermi = dfermi + fermi(ene,beta)
     END DO

     DO ifreq = 1,nlog
        iwn=CMPLX(mu,freq_log(ifreq),dp)
        DO ibnd = 1,num_bnd
           domega = domega+(1.d0/(iwn-eigs(ibnd,ifreq,ikp)))
           domega = domega-(1.d0/(iwn-eigs(ibnd,nlog,ikp)))
        END DO
     END DO ! end freq loop

     occ=occ+wk(disp(iproc)+ikp)*(dfermi+REAL(domega,dp)/beta)

  END DO ! loop over num_kp
  !===================================================================!

  RETURN

END SUBROUTINE dmft_occup
