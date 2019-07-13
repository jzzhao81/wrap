! SUBROUTINE INDEX
! dmft_init : main dmft_init subroutine
! ekk_init : get ekk from ek
! snorm   : normalize S-matrix
SUBROUTINE dmft_init()

  IMPLICIT NONE

  !transform arrays to log mesh point
  CALL dmft_logmesh()

  ! mix old and new self-energy
  CALL dmft_mix()

  RETURN

END SUBROUTINE dmft_init
!=====================================================================!

! calclate eimp from ek & smtrx
!=====================================================================!
SUBROUTINE cal_eimp()

  USE mod_param,ONLY : dp,zzero
  USE mod_dmft,ONLY : num_kp,num_orb,num_bnd,ek,wk,eimp,smtrx
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: ikp, ibnd
  COMPLEX(kind=dp),ALLOCATABLE :: empi(:,:)
  COMPLEX(kind=dp),ALLOCATABLE :: etmp(:,:)
  COMPLEX(kind=dp),ALLOCATABLE :: smat(:,:)
  COMPLEX(kind=dp),ALLOCATABLE :: emat(:,:)
  COMPLEX(kind=dp), EXTERNAL :: mtrace

  ALLOCATE(empi(num_orb,num_orb))
  ALLOCATE(etmp(num_orb,num_orb))
  ALLOCATE(smat(num_bnd,num_orb))
  ALLOCATE(emat(num_bnd,num_bnd))

  smat = zzero; eimp = zzero
  empi = zzero; etmp = zzero

  ! we compute eimp from ek & s-matrix
  !===================================================================!
  DO ikp=1+iproc,num_kp, nproc

     smat = smtrx(:,:,ikp)
     emat = zzero
     DO ibnd = 1, num_bnd
        emat(ibnd,ibnd)=ek(ibnd,ikp)
     ENDDO

     CALL dmft_proj_loc(num_orb,num_bnd,emat,etmp,smat)

     empi=empi+etmp*wk(ikp)

  END DO

  CALL MPI_ALLREDUCE(empi,eimp,SIZE(eimp),MPI_DOUBLE_COMPLEX,MPI_&
       &SUM,MPI_COMM_WORLD,ierr)

  DEALLOCATE(smat,emat)
  DEALLOCATE(empi,etmp)

  RETURN

END SUBROUTINE cal_eimp
!=====================================================================!


! for the first DMFT loop, initialize frequence, self-energy
!=====================================================================!
SUBROUTINE dmft_mkmesh()

  USE mod_param,ONLY : dp, pi
  USE mod_dmft,ONLY : num_freq,freq,beta
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: itmp

  ! build matsubara frequency mesh: rmesh
  DO itmp=1,num_freq
     freq(itmp)=(2.d0*REAL(itmp-1,dp)+1.d0)*(pi/beta)
  ENDDO ! over j={1,mfreq} loop

  RETURN

END SUBROUTINE dmft_mkmesh
!=====================================================================!

!=====================================================================!
! log mesh point
!=====================================================================!
SUBROUTINE dmft_logmesh()

  USE mod_param,ONLY : dp, zzero, dzero
  USE mod_dmft, ONLY : num_freq,nlog,nhead,ntail,sgm,freq,sgm_log,freq_log

  IMPLICIT NONE

  INTEGER :: index_log(nlog)

  INTEGER :: ipt, tmp, ifreq
  REAL(kind=dp) :: eta

  ! this is a dense parameter, determined by nhead & ntail
  eta=dlog(DBLE(num_freq)/DBLE(nhead))/DBLE(ntail)

  ! we keep first nhead point as they were
  DO ipt=1,nhead
     index_log(ipt)=ipt
  END DO

  ! we keep the last ntail point as log mesh
  DO ipt=1,ntail
     tmp=INT(nhead*dexp(eta*ipt)+0.5d0)
     index_log(ipt+nhead)=tmp
  END DO

  freq_log = dzero
  sgm_log  = zzero
  DO ifreq=1,nlog
     freq_log(ifreq)=freq(index_log(ifreq))
     sgm_log(:,ifreq)=sgm(:,index_log(ifreq))
  END DO

  RETURN

END SUBROUTINE dmft_logmesh
