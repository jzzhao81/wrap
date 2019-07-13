!=========>=========>=========>=========>=========>=========>=========!
! calculate occupation, chemical potential
!
! In this subroutine, we use Van Wijngaarden–Dekker–Brent Method to
! find chemical potential, generally this is the recomment method
! in numerical recipes, reference NUMERICAL RECIPES IN FORTRAN 77: P352
! or
! http://www.mapleprimes.com/posts/38757-Brents-Method-For-Root-Finding
!=========>=========>=========>=========>=========>=========>=========!
SUBROUTINE dmft_mu()

  USE mod_param,ONLY : dp,dzero,max_mu_step,occ_eps,ha2ev,zzero
  USE mod_dmft,ONLY : num_ele,mu_dmft,ek,nlog,num_kp,smtrx,num_bnd,&
       sgm_log,num_orb
  USE mod_eig
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: istep,icount,ifreq,ikp,iorb,ibnd
  ! zer0
  REAL(kind=dp), PARAMETER :: mu_broad=1.d1
  ! machine floating-point precision ()
  REAL(kind=dp), PARAMETER :: eps=1.d-10

  REAL(kind=dp) :: lmu,umu,mmu
  REAL(kind=dp) :: locc,uocc,mocc
  REAL(kind=dp) :: p,q,t1,t2,tol1,judge_ib
  REAL(kind=dp) :: pstep,nstep

  REAL(kind=dp) :: time1, time2
  COMPLEX(kind=dp) :: vocc, vmpi

  ! eigenvalue
  COMPLEX(kind=dp), ALLOCATABLE :: eval(:,:,:),eigs(:)
  ! sgm_tmp
  COMPLEX(kind=dp),ALLOCATABLE :: stmp(:,:)
  ! self-energy in bloch state
  COMPLEX(kind=dp),ALLOCATABLE :: sgmk(:,:)
  ! hamiltonian, LDA eigenvalue plus self-energy
  COMPLEX(kind=dp),ALLOCATABLE :: hams(:,:)
  ! temp smtrx
  COMPLEX(kind=dp),ALLOCATABLE :: smat(:,:)
  ! temp eigenmatrix
  COMPLEX(kind=dp),ALLOCATABLE :: emat(:,:)

  LOGICAL :: cflag

  ! goto 1001

  IF(iproc==master) THEN
     WRITE(6,*)
     WRITE(6,"(' Seaching chemical potential !')")
  ENDIF
  icount=0

  CALL mpi_begin(num_kp, num_bnd*nlog)

  ALLOCATE(eigs(num_bnd))
  ALLOCATE(eval(num_bnd,nlog,nkpp(iproc)))
  ALLOCATE(hams(num_bnd,num_bnd))
  ALLOCATE(stmp(num_orb,num_orb))
  ALLOCATE(sgmk(num_bnd,num_bnd))
  ALLOCATE(smat(num_bnd,num_orb))
  ALLOCATE(emat(num_bnd,num_bnd))

  eval=zzero

  ! Diagonal Non-Hermitian Hamiltonian
  ! Save eigenvalues, left and right eigenvectors
  ! perform diagonalization for each k-point
  DO ikp = 1, nkpp(iproc)
     DO ifreq = 1, nlog

        ! prepare self-energy
        stmp=zzero; emat=zzero
        DO iorb=1,num_orb
           stmp(iorb,iorb)=sgm_log(iorb,ifreq)
        END DO
        DO ibnd=1,num_bnd
           emat(ibnd,ibnd)=ek(ibnd,disp(iproc)+ikp)
        ENDDO

        ! initial
        smat = zzero; hams = zzero; eigs = zzero
        ! project self-energy to bloch basis
        smat = smtrx(:,:,disp(iproc)+ikp)
        CALL dmft_proj_ks(num_bnd,num_orb,stmp,sgmk,smat)

        ! build hams and diagonal
        hams = emat+sgmk
        CALL eigsys(num_bnd,hams,eigs)
        ! empi(:,isymop,ikp,ifreq) = eigs
        eval(:,ifreq,ikp) = eigs

     END DO ! loop over freq
  END DO ! loop over num_kp

  DEALLOCATE(eigs)
  DEALLOCATE(hams)
  DEALLOCATE(stmp)
  DEALLOCATE(sgmk)
  DEALLOCATE(smat)

  !initialize variable
  !===================================================================!
  mmu=mu_dmft
  cflag=.FALSE.
  ! get up and down boundary
  DO WHILE(.NOT. cflag)
     lmu=mmu-mu_broad
     umu=mmu+mu_broad
     CALL dmft_occup(lmu,eval,vmpi)
     CALL mpi_allreduce(vmpi,vocc,1,mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
     locc=REAL(vocc)
     CALL dmft_occup(umu,eval,vmpi)
     CALL mpi_allreduce(vmpi,vocc,1,mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
     uocc=REAL(vocc)
     IF((uocc-num_ele)*(locc-num_ele) > dzero) THEN
        IF(ABS(uocc-num_ele) < ABS(locc-num_ele))THEN
           mmu=umu
        ELSE
           mmu=lmu
        END IF
        IF(mmu<MINVAL(REAL(ek)) .OR. mmu>MAXVAL(REAL(ek))) THEN
           CALL mpi_fnl(1, 'Failed in Search Mu !')
        END IF
        cflag=.FALSE.
     ELSE
        cflag=.TRUE.
     END IF
  END DO
  !===================================================================!

  mmu=lmu
  mocc=locc

  ! information title
  IF(iproc==master) THEN
     WRITE(6,"(' Step',4X,' CP ', 6X, ' ntot ', 5X, ' nlda')")
     WRITE(6,"(' !======================================!')")
  END IF

  !search chemical potential using bisection method
  !===================================================================!
  DO istep=1,max_mu_step

     time1=MPI_WTIME()

     pstep=umu-lmu

     IF(ABS(mocc-num_ele) < ABS(uocc-num_ele)) THEN
        lmu=umu
        umu=mmu
        mmu=lmu
        locc=uocc
        uocc=mocc
        mocc=locc
     END IF

     ! check tol1 if converged return
     tol1=2.d0*eps*ABS(umu)+5.d-1*occ_eps
     nstep=5.d-1*(mmu-umu)

     IF(ABS(nstep) < tol1 .OR.(ABS(uocc-num_ele)<occ_eps) ) THEN
        IF(iproc==master)THEN
           WRITE(6,"(' !======================================!')")
           WRITE(6,"(' Total interpolation count : ',I3)")icount
        END IF
        mu_dmft=umu
        RETURN
     END IF

     ! if posible, use inverse quadratic interpolation
     IF(ABS(nstep)>tol1 .AND. ABS(locc-num_ele)>ABS(uocc-num_ele)) THEN

        t1=(uocc-num_ele)/(locc-num_ele)

        IF(lmu == mmu) THEN
           p=2.d0*nstep*t1
           q=1.d0-t1
        ELSE
           q  = (locc-num_ele)/(mocc-num_ele)
           t2 = (uocc-num_ele)/(mocc-num_ele)
           p  = t1*(2.d0*nstep*q*(q-t2)-(umu-lmu)*(t2-1.d0))
           q  = (q-1.d0)*(t2-1.d0)*(t1-1.d0)
        END IF

        IF(p>0.d0) q=-q

        p=ABS(p)

        judge_ib = MIN(3.d0*nstep*q-5.d-1*ABS(tol1*q),ABS(pstep*q))
        judge_ib = judge_ib/2.d0

        IF(p < judge_ib) THEN
           nstep  = p/q
           icount = icount+1
        END IF

     END IF ! endif use inverse quadratic interpolation

     lmu  = umu
     locc = uocc

     IF(ABS(nstep) < tol1) THEN
        IF(nstep>0) THEN
           nstep=tol1
        ELSE
           nstep=-tol1
        END IF
     END IF

     umu=umu+nstep
     CALL dmft_occup(umu,eval,vmpi)
     CALL mpi_allreduce(vmpi,vocc,1,mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
     uocc=REAL(vocc)

     IF((uocc-num_ele)*(mocc-num_ele) > dzero) THEN
        mocc = locc
        mmu  = lmu
     END IF

     time2=MPI_WTIME()

     IF(iproc==master) THEN
        WRITE(6,"(1X,I3,1X,F10.6,2X,F10.6,2X,F6.2)")istep,umu,uocc,num_ele
     ENDIF

  END DO ! loop over istep
  !===================================================================!

  IF(iproc==master) THEN
     WRITE(6,"(' !======================================!')")
     WRITE(6,"(' Total interpolation count : ',I3)")icount
     WRITE(6,"(' Caution : exceeding maximum iterations !')")
  END IF
  ! give results to global variable

  mu_dmft = umu

  DEALLOCATE(eval)

  RETURN

END SUBROUTINE dmft_mu
