!=========>=========>=========>=========>=========>=========>=========!
! calculate density matrix
!=========>=========>=========>=========>=========>=========>=========!
SUBROUTINE lda_occ()

  USE mod_param,ONLY : dp,dzero,zzero,dhalf
  USE mod_dmft,ONLY : num_bnd,num_kp,num_orb,num_aorb,smtrx,ek,occ_lda,wk,beta
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  !variable definition
  !===================================================================!
  ! input chemical potential
  REAL(kind=dp) :: occ,occ_mpi
  ! loop parameter
  INTEGER :: ikp, ibnd
  ! temp rho matrx
  REAL(kind=dp),ALLOCATABLE :: rho_mpi(:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: fmat(:,:)
  ! lda rho matrx
  REAL(kind=dp),ALLOCATABLE :: rho_lda(:,:)
  ! rho on local basis
  COMPLEX(kind=dp),ALLOCATABLE :: lrho(:,:)
  ! nflda for mpi
  COMPLEX(kind=dp),ALLOCATABLE :: nflda(:,:)
  COMPLEX(kind=dp),ALLOCATABLE :: nfmpi(:,:)
  ! temp smtrx
  COMPLEX(kind=dp),ALLOCATABLE :: smat(:,:)
  ! external function of fermi distribution
  COMPLEX(kind=dp), EXTERNAL :: mtrace
  REAL(kind=dp), EXTERNAL :: fermi
  !===================================================================!

  ! allocate fermi array
  ALLOCATE(rho_mpi(num_bnd,num_kp))
  ALLOCATE(rho_lda(num_bnd,num_kp))
  ALLOCATE(smat(num_bnd,num_orb))

  rho_mpi = dzero; rho_lda = dzero

  ! calculate density matrix
  !===================================================================!
  DO ikp=1+iproc,num_kp,nproc
     ! get fermi part
     DO ibnd=1,num_bnd
        rho_mpi(ibnd,ikp)=fermi(REAL(ek(ibnd,ikp)),beta)
     END DO

  END DO ! loop over num_kp

  CALL MPI_ALLREDUCE(rho_mpi,rho_lda,SIZE(rho_lda),MPI_DOUBLE_precision,&
       &MPI_SUM,MPI_COMM_WORLD,ierr)
  !===================================================================!
  ! deallocate temp array
  DEALLOCATE(rho_mpi)

  ! get local occ
  !==================================================================!
  ! allocate lrho
  ALLOCATE(lrho(num_orb,num_orb))
  ALLOCATE(fmat(num_bnd,num_bnd))
  ALLOCATE(nfmpi(num_orb,num_orb))
  ALLOCATE(nflda(num_orb,num_orb))

  nfmpi=zzero; nflda=zzero

  ! 1. project rho onto local basis
  ! 2. sum for each orbital
  DO ikp=1+iproc,num_kp,nproc

     smat = smtrx(:,:,ikp)

     fmat=zzero
     DO ibnd=1,num_bnd
        fmat(ibnd,ibnd)=CMPLX(rho_lda(ibnd,ikp),dzero,dp)
     ENDDO

     CALL dmft_proj_loc(num_orb,num_bnd,fmat,lrho,smat)

     nfmpi=nfmpi+lrho*wk(ikp)

  END DO ! loop over num_kp

  ! reduce variable
  CALL MPI_ALLREDUCE(nfmpi,nflda,SIZE(nflda),MPI_DOUBLE_COMPLEX,MPI_SUM,&
       &MPI_COMM_WORLD,ierr)
  !===================================================================!
  DEALLOCATE(lrho)
  DEALLOCATE(nfmpi)
  DEALLOCATE(smat)
  DEALLOCATE(fmat)

  occ_lda = REAL(mtrace(num_aorb,nflda(1:num_aorb,1:num_aorb)))

  occ_mpi=zzero
  ! evaluate valence electron occupation
  !===================================================================!
  DO ikp=1+iproc,num_kp,nproc
     DO ibnd=1,num_bnd
        occ_mpi=occ_mpi+rho_lda(ibnd,ikp)*wk(ikp)
     END DO
  END DO

  CALL mpi_allreduce(occ_mpi,occ,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
  !===================================================================!
  DEALLOCATE(rho_lda)
  DEALLOCATE(nflda)

  IF(iproc==master) THEN
     WRITE(*,"(' Occupation from LDA :',2X,F7.3,3X,F7.3)") occ,occ_lda
  END IF

  RETURN

END SUBROUTINE lda_occ
