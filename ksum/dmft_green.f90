!=========>=========>=========>=========>=========>=========>=========!
! Get Green's function on local state
! Get hybridization function
!=========>=========>=========>=========>=========>=========>=========!

SUBROUTINE dmft_green()

  USE mod_param,ONLY : dp,zzero,dzero
  USE mod_dmft,ONLY : nlog,num_kp,num_bnd,num_orb,gloc,hyb_log,sgm_log,ek,&
       &smtrx,freq_log,wk,eimp,mu_dmft,wk
  USE mod_mtrx
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: ifreq, ikp, iorb, ibnd
  ! green's function in bloch state
  COMPLEX(kind=dp),ALLOCATABLE :: glat(:,:)
  !self-energy in bloch state
  COMPLEX(kind=dp),ALLOCATABLE :: sgmk(:,:)
  ! tmp local green's function
  COMPLEX(kind=dp),ALLOCATABLE :: gtmp(:,:)
  COMPLEX(kind=dp),ALLOCATABLE :: gmpi(:,:,:)
  ! tmp hybrid function
  COMPLEX(kind=dp),ALLOCATABLE :: hmpi(:,:)
  COMPLEX(kind=dp),ALLOCATABLE :: htmp(:,:)
  COMPLEX(kind=dp),ALLOCATABLE :: liwn(:,:),biwn(:,:)
  COMPLEX(kind=dp),ALLOCATABLE :: invg(:,:)
  COMPLEX(kind=dp),ALLOCATABLE :: stmp(:,:)
  ! temp smtrx
  COMPLEX(kind=dp),ALLOCATABLE :: smat(:,:)
  ! temp ek
  COMPLEX(kind=dp),ALLOCATABLE :: emat(:,:)

  ! allocate iw for local basis & Bloch basis
  ! and zero them
  !===================================================================!
  ALLOCATE(glat(num_bnd,num_bnd))
  ALLOCATE(sgmk(num_bnd,num_bnd))
  ALLOCATE(biwn(num_bnd,num_bnd))
  ALLOCATE(liwn(num_orb,num_orb))
  ALLOCATE(gtmp(num_orb,num_orb))
  ALLOCATE(gmpi(num_orb,num_orb,nlog))
  ALLOCATE(stmp(num_orb,num_orb))
  ALLOCATE(smat(num_bnd,num_orb))
  ALLOCATE(emat(num_bnd,num_bnd))

  glat = zzero; sgmk = zzero
  biwn = zzero; liwn = zzero
  gtmp = zzero; gmpi = zzero
  stmp = zzero; smat = zzero

  !get green's function on local state
  !===================================================================!
  DO ifreq=1+iproc,nlog,nproc

     gtmp=zzero

     !mu+iw, sigma
     biwn=zzero; stmp=zzero
     DO ibnd=1,num_bnd
        biwn(ibnd,ibnd)=CMPLX(mu_dmft,freq_log(ifreq),dp)
     END DO
     DO iorb=1,num_orb
        stmp(iorb,iorb)=sgm_log(iorb,ifreq)
     END DO

     ! calculate Green's function on local basis
     !================================================================!
     DO ikp=1,num_kp

        !project self-energy on KS state
        smat = smtrx(:,:,ikp)
        emat = zzero
        DO ibnd = 1, num_bnd
           emat(ibnd,ibnd) = ek(ibnd,ikp)
        ENDDO

        CALL dmft_proj_ks(num_bnd,num_orb,stmp,sgmk,smat)

        ! lattice Green's function
        glat=biwn-emat-sgmk

        CALL mtrx_invs(num_bnd,glat)

        ! project lattice Green's function on local state
        CALL dmft_proj_loc(num_orb,num_bnd,glat,gtmp,smat)

        gmpi(:,:,ifreq) = gmpi(:,:,ifreq)+gtmp*wk(ikp)

     END DO ! loop over k-points
     !================================================================!

  END DO ! loop over freq

  CALL MPI_ALLREDUCE(gmpi,gloc,SIZE(gloc),MPI_DOUBLE_COMPLEX,MPI_SUM,&
       &MPI_COMM_WORLD,ierr)

  ! deallocate tmp gloc array
  DEALLOCATE(glat,sgmk)
  DEALLOCATE(gtmp,gmpi)
  DEALLOCATE(smat)

  !get hybrid function on local state
  !===================================================================!
  ALLOCATE(hmpi(num_orb,nlog))
  ALLOCATE(htmp(num_orb,num_orb))
  ALLOCATE(invg(num_orb,num_orb))

  invg = zzero; hmpi = zzero; htmp = zzero

  DO ifreq=1+iproc,nlog,nproc

     ! mu+iw, sigma
     liwn=zzero; stmp=zzero
     DO iorb=1,num_orb
        liwn(iorb,iorb)=CMPLX(mu_dmft,freq_log(ifreq),dp)
     END DO
     DO iorb=1,num_orb
        stmp(iorb,iorb)=sgm_log(iorb,ifreq)
     END DO

     ! get inverse of local green's function
     invg=gloc(:,:,ifreq)

     CALL mtrx_invs(num_orb,invg)

     htmp=liwn-eimp-invg-stmp

     DO iorb = 1, num_orb
        hmpi(iorb,ifreq) = htmp(iorb,iorb)
     END DO

  END DO ! loop over freq
  !===================================================================!

  CALL MPI_ALLREDUCE(hmpi,hyb_log,SIZE(hyb_log),MPI_DOUBLE_COMPLEX,&
       &MPI_SUM,MPI_COMM_WORLD,ierr)

  ! deallocate invg, no use after now
  DEALLOCATE(invg)
  DEALLOCATE(htmp)
  DEALLOCATE(hmpi)
  DEALLOCATE(stmp)
  DEALLOCATE(biwn)
  DEALLOCATE(liwn)
  !===================================================================!

  RETURN

END SUBROUTINE dmft_green
