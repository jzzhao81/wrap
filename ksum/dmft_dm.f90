!=========>=========>=========>=========>=========>=========>=========!
! calculate density matrix
! We have to diagonal the non-Hermitian Hamiltonian again
! Because we cannot save large lvec and rvec
!=========>=========>=========>=========>=========>=========>=========!

SUBROUTINE dmft_dm()

  USE mod_param,ONLY : dp,dzero,zzero
  USE mod_dmft,ONLY : num_bnd,num_kp,num_orb,num_aorb,sgm_log,dmft1,mu_dmft,&
       &nlog,smtrx,ek,beta,wk,freq_log,rho,occ_dmft
  USE mod_mtrx
  USE mod_eig
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  !variable definition
  !===================================================================!
  ! input chemical potential
  REAL(kind=dp) :: mu
  COMPLEX(kind=dp) :: occ
  COMPLEX(kind=dp) :: iwn
  ! loop parameter
  INTEGER :: ikp, ibnd, ifreq, iorb
  ! a tmp energy
  REAL(kind=dp) :: ene
  ! fermi distribution
  COMPLEX(kind=dp),ALLOCATABLE :: vtmp(:,:), stmp(:,:)
  COMPLEX(kind=dp),ALLOCATABLE :: dfermi(:,:), domega(:,:)
  ! ew_min : highest freq eigenvalue
  COMPLEX(kind=dp),ALLOCATABLE :: evk_inf(:),evk(:)
  COMPLEX(kind=dp),ALLOCATABLE :: vl_inf(:,:),vr_inf(:,:)
  COMPLEX(kind=dp),ALLOCATABLE :: vl(:,:),vr(:,:)
  ! self-energy in bloch state
  COMPLEX(kind=dp),ALLOCATABLE :: sgmk(:,:)
  ! hamiltonian, LDA eigenvalue plus self-energy
  COMPLEX(kind=dp),ALLOCATABLE :: hmat(:,:)
  ! external function of fermi distribution
  REAL(kind=dp), EXTERNAL :: fermi
  ! temp rho matrx
  COMPLEX(kind=dp),ALLOCATABLE :: rho_mpi(:,:,:)
  ! occ_dmft for mpi
  COMPLEX(kind=dp),ALLOCATABLE :: occ_mpi(:)
  ! rho on local basis
  COMPLEX(kind=dp),ALLOCATABLE :: lrho(:,:)
  ! temp smtrx
  COMPLEX(kind=dp),ALLOCATABLE :: smat(:,:)
  ! temp Ek
  COMPLEX(kind=dp),ALLOCATABLE :: emat(:,:)


  !===================================================================!
  mu=mu_dmft

  ! mpi begin
  CALL mpi_begin(num_kp, num_bnd*num_bnd)

  ! allocate fermi array
  ALLOCATE(dfermi(num_bnd,num_bnd),domega(num_bnd,num_bnd))
  ALLOCATE(evk(num_bnd),evk_inf(num_bnd))
  ALLOCATE(vl(num_bnd,num_bnd),vl_inf(num_bnd,num_bnd))
  ALLOCATE(vr(num_bnd,num_bnd),vr_inf(num_bnd,num_bnd))
  ALLOCATE(hmat(num_bnd,num_bnd), vtmp(num_bnd,num_bnd))
  ALLOCATE(stmp(num_orb,num_orb), smat(num_bnd,num_orb))
  ALLOCATE(emat(num_bnd,num_bnd))
  ALLOCATE(rho_mpi(num_bnd,num_bnd,nkpp(iproc)))
  ALLOCATE(sgmk(num_bnd,num_bnd))

  ! initialize variable
  vl      = zzero; vr      = zzero
  vl_inf  = zzero; vr_inf  = zzero
  evk     = zzero; evk_inf = zzero
  dfermi  = zzero; domega  = zzero
  rho_mpi = zzero; sgmk    = zzero
  hmat    = zzero; vtmp    = zzero

  ! calculate density matrix
  !===================================================================!
  DO ikp=1, nkpp(iproc)

     stmp=zzero; emat = zzero
     DO iorb=1,num_orb
        stmp(iorb,iorb)=sgm_log(iorb,nlog)
     END DO
     DO ibnd=1,num_bnd
        emat(ibnd,ibnd)=ek(ibnd,disp(iproc)+ikp)
     ENDDO

     smat = smtrx(:,:,disp(iproc)+ikp)

     ! calculate eigenvalue of highest freq
     CALL dmft_proj_ks(num_bnd,num_orb,stmp,sgmk,smat)
     hmat=emat+sgmk
     CALL eigsys(num_bnd,hmat,evk_inf,vl_inf,vr_inf)

     ! get fermi contribution to density matrix
     vtmp = zzero
     DO ibnd=1,num_bnd
        ene = REAL(evk_inf(ibnd)) - mu
        vtmp(ibnd,:) = vl_inf(ibnd,:)
        ! Computes the product of a vector by a scalar
        CALL zdscal(num_bnd,fermi(ene,beta),vtmp(ibnd,:),1)
     END DO
     CALL mtrx_mult(vr_inf,vtmp,dfermi,'N','N')

     ! Omega contribution to density matrix
     domega =zzero
     IF(.NOT. dmft1) THEN ! for dmft loop after first loop

        DO ifreq=1,nlog

           stmp=zzero
           DO iorb=1,num_orb
              stmp(iorb,iorb)=sgm_log(iorb,ifreq)
           END DO

           ! calculate eigenvalue of omega dependent part
           CALL dmft_proj_ks(num_bnd,num_orb,stmp,sgmk,smat)
           hmat=emat+sgmk
           CALL eigsys(num_bnd,hmat,evk,vl,vr)

           DO ibnd=1,num_bnd
              iwn=1.d0/(CMPLX(mu,freq_log(ifreq),dp)-evk(ibnd))
              vtmp(ibnd,:) = vl(ibnd,:)
              CALL zscal(num_bnd,iwn,vtmp(ibnd,:),1)
           END DO
           CALL mtrx_mult(vr,vtmp,hmat,'N','N')
           domega = domega + hmat

           DO ibnd=1,num_bnd
              iwn=1.d0/(CMPLX(mu,freq_log(ifreq),dp)-evk_inf(ibnd))
              vtmp(ibnd,:) = vl_inf(ibnd,:)
              CALL zscal(num_bnd,iwn,vtmp(ibnd,:),1)
           END DO
           CALL mtrx_mult(vr_inf,vtmp,hmat,'N','N')
           domega = domega - hmat

        END DO ! end freq loop

     END IF ! end if dmft1

     rho_mpi(:,:,ikp) = dfermi + domega/beta

  END DO ! loop over num_kp

  CALL mpi_allgatherv(rho_mpi,SIZE(rho_mpi),mpi_double_complex,rho,ndpp,ndsp,&
       mpi_double_complex,mpi_comm_world,ierr)
  !===================================================================!

  ! deallocate temp array
  DEALLOCATE(dfermi,domega)
  DEALLOCATE(evk_inf,evk)
  DEALLOCATE(vl,vl_inf)
  DEALLOCATE(vr,vr_inf)
  DEALLOCATE(stmp,vtmp)
  DEALLOCATE(rho_mpi)
  DEALLOCATE(sgmk,emat)

  ! get local occ
  !==================================================================!
  ! allocate lrho
  ALLOCATE(occ_mpi(num_orb))
  ALLOCATE(lrho(num_orb,num_orb))

  occ_mpi=zzero
  occ_dmft=zzero
  ! 1. project rho onto local basis
  ! 2. sum for each orbital
  DO ikp=1+iproc,num_kp,nproc

     smat = smtrx(:,:,ikp)
     hmat = rho(:,:,ikp)
     lrho=zzero

     CALL dmft_proj_loc(num_orb,num_bnd,hmat,lrho,smat)

     DO iorb=1,num_orb
        occ_mpi(iorb)=occ_mpi(iorb)+wk(ikp)*lrho(iorb,iorb)
     END DO

  END DO ! loop over num_kp

  ! reduce variable
  CALL MPI_ALLREDUCE(occ_mpi,occ_dmft,SIZE(occ_dmft),MPI_DOUBLE_COMPLEX,&
       &MPI_SUM,MPI_COMM_WORLD,ierr)
  !===================================================================!
  DEALLOCATE(lrho)

  ! evaluate valence electron occupation
  !===================================================================!
  occ_mpi=zzero
  occ=zzero
  DO ikp=1+iproc,num_kp,nproc
     hmat = rho(:,:,ikp)
     DO ibnd=1,num_bnd
        occ_mpi(1)=occ_mpi(1)+wk(ikp)*hmat(ibnd,ibnd)
     END DO
  END DO

  CALL mpi_allreduce(occ_mpi(1), occ, 1, mpi_double_complex,mpi_sum,&
       &mpi_comm_world,ierr)
  !===================================================================!

  DEALLOCATE(occ_mpi)
  DEALLOCATE(smat)
  DEALLOCATE(hmat)

  IF(iproc==master) THEN
     WRITE(*,*)
     WRITE(*,"(' Occupation from DMFT :',2X,F7.3,3X,F7.3)")REAL(occ),&
          &REAL(SUM(occ_dmft(:num_aorb)),dp)
  END IF

  RETURN

END SUBROUTINE dmft_dm


COMPLEX(kind=8) FUNCTION mtrace(n,mtrx)

  USE mod_param, ONLY : zzero

  IMPLICIT NONE

  INTEGER, INTENT(in) :: n
  COMPLEX(kind=8), INTENT(in) :: mtrx(n,n)

  INTEGER :: ii

  mtrace=zzero

  DO ii=1,n
     mtrace=mtrace+mtrx(ii,ii)
  END DO

  RETURN

END FUNCTION mtrace
