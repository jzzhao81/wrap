!=========>=========>=========>=========>=========>=========>=========!
! This subroutine is to read dmft input file
! 1. read lda input : dmft_read_lda
! 2. read ctqmc input : dmft_read_inp
! 3. read self-energy : dmft_read_sgm
!=========>=========>=========>=========>=========>=========>=========!
SUBROUTINE dmft_read

  USE mod_param,ONLY : dp,dzero,zzero

  IMPLICIT NONE

  ! read dmft_inf.dat
  CALL dmft_read_lda()
  !===================================================================!

  ! read solver.ctqmc.in
  CALL dmft_read_inp()
  !===================================================================!

  ! evaluate local occ from LDA results
  CALL lda_occ()
  !===================================================================!

  ! if edc.dat exist, then read edc from file
  !===================================================================!
  CALL dmft_read_edc()

  ! if solver.sgm.dat exist then read
  !===================================================================!
  CALL dmft_read_sgm()

  RETURN

END SUBROUTINE dmft_read

! real parameter from lda
! Unit : eV
!=====================================================================!
SUBROUTINE dmft_read_lda

  USE mod_param,ONLY : zzero,dzero,dp
  USE mod_dmft,ONLY : num_kp,num_bnd,num_ele,num_atm,lmoment,ek,wk,eimp,&
       & smtrx,rho,occ_dmft,num_aorb,num_orb,nspin,nemin,nemax,mu_lda
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: file_id=1001
  INTEGER :: ikp,ibnd,iorb,ii
  REAL(kind=dp) :: rpart, ipart
  LOGICAL :: alive

  ! CARUTION !!!!!!!!
  ! In this code, we modified for TWO EQUAL correlated atoms
  ! num_atm=3

  IF(iproc==master) THEN
     INQUIRE(file='eig.dat', exist=alive)
     IF(.NOT. alive) THEN
        CALL mpi_fnl(1, 'File : eig.dat does not exist !')
     END IF
     OPEN(file_id,file='eig.dat',form='formatted',status='old')
     READ(file_id,*) num_kp,num_atm,num_aorb,lmoment,nspin
     READ(file_id,*) nemin, nemax, num_bnd,num_ele,mu_lda
     READ(file_id,*)
     num_orb = num_aorb * num_atm
  END IF

  CALL mpi_barrier(mpi_comm_world,ierr)
  CALL mpi_bcast(num_kp,1,mpi_integer,master,mpi_comm_world,ierr)
  CALL mpi_bcast(num_aorb,1,mpi_integer,master,mpi_comm_world,ierr)
  CALL mpi_bcast(num_orb,1,mpi_integer,master,mpi_comm_world,ierr)
  CALL mpi_bcast(num_atm,1,mpi_integer,master,mpi_comm_world,ierr)
  CALL mpi_bcast(num_bnd,1,mpi_integer,master,mpi_comm_world,ierr)
  CALL mpi_bcast(lmoment,1,mpi_integer,master,mpi_comm_world,ierr)
  CALL mpi_bcast(nspin,1,mpi_integer,master,mpi_comm_world,ierr)
  CALL mpi_bcast(nemin,1,mpi_integer,master,mpi_comm_world,ierr)
  CALL mpi_bcast(nemax,1,mpi_integer,master,mpi_comm_world,ierr)
  CALL mpi_bcast(num_ele,1,mpi_double_precision,master,mpi_comm_world,ierr)
  CALL mpi_bcast(mu_lda,1,mpi_double_precision,master,mpi_comm_world,ierr)

  ! allocate LDA variables
  ALLOCATE(wk(num_kp))
  ALLOCATE(ek(num_bnd,num_kp))
  ALLOCATE(eimp(num_orb,num_orb))
  ALLOCATE(smtrx(num_bnd,num_orb,num_kp))
  ALLOCATE(rho(num_bnd,num_bnd,num_kp))
  ALLOCATE(occ_dmft(num_orb))

  ek    = zzero
  smtrx = zzero
  eimp  = zzero
  rho   = zzero
  occ_dmft = zzero

  IF(iproc==master) THEN
     DO ikp=1,num_kp
        READ(file_id,*) ii, wk(ikp)
        DO ibnd = 1,num_bnd
           READ(file_id,*) ii, rpart
           ek(ibnd,ikp) = CMPLX(rpart-mu_lda,dzero,dp)
        END DO
        READ(file_id,*)
        READ(file_id,*)
     END DO
     CLOSE(file_id)
  END IF

  CALL mpi_barrier(mpi_comm_world,ierr)
  CALL mpi_bcast(ek,SIZE(ek),mpi_double_complex,master,mpi_comm_world,&
       &ierr)
  CALL mpi_bcast(wk,SIZE(wk),mpi_double_precision,master,mpi_comm_world,&
       &ierr)

  IF(iproc==master) THEN

     INQUIRE(file='udmft.dat', exist=alive)
     IF(.NOT. alive) THEN
        CALL mpi_fnl(1, 'File : udmft.dat does not exist !')
     END IF
     OPEN(file_id,file='udmft.dat',form='formatted',status='old')

     READ(file_id,*)
     READ(file_id,*)
     READ(file_id,*)

     DO ikp=1,num_kp
        READ(file_id,*)
        DO iorb=1,num_orb
           DO ibnd = 1, num_bnd
              READ(file_id,*) ii, ii, rpart, ipart
              smtrx(ibnd,iorb,ikp)=CMPLX(rpart,ipart,dp)
           END DO
        END DO
        READ(file_id,*)
        READ(file_id,*)
     ENDDO

     !close dmft_inf.dat
     CLOSE(file_id)

  END IF

  ! broadcast varaible
  CALL mpi_barrier(mpi_comm_world,ierr)
  CALL mpi_bcast(smtrx,SIZE(smtrx),mpi_double_complex,master,mpi_comm_world,&
       &ierr)

  ! get eimp from ek & smtrx
  CALL cal_eimp()
  !===================================================================!

  RETURN

END SUBROUTINE dmft_read_lda

!read CTQMC(azalea) solver input : solver.ctqmc.in
!=====================================================================!
SUBROUTINE dmft_read_inp()

  USE mod_param,ONLY : dp, dzero, zzero
  USE mod_dmft,ONLY : file_inp,uj,num_freq,mu_dmft,beta,nhead,ntail,&
       &nlog,alpha,sgm,freq,gloc,hyb,nlog,num_orb,sgm_log,freq_log,&
       &hyb_log,sgm_old
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: max_tail = 512
  INTEGER,PARAMETER :: file_id=1001
  INTEGER :: iline
  LOGICAL :: alive

  IF(iproc==master) THEN

     INQUIRE(file=file_inp, exist=alive)
     IF(.NOT. alive) THEN
        CALL mpi_fnl(1,'File : dmft.main.in does not exist !')
     END IF

     ! open dmft.main.in
     OPEN(file_id, file=file_inp, form='formatted', status='old')

     DO iline=1, 3
        READ(file_id,*)
     END DO

     !read U and J
     READ(file_id,*) uj(1)
     READ(file_id,*) uj(2)
     READ(file_id,*)

     ! read mu_dmft from ctqmc solver
     READ(file_id,*) mu_dmft
     ! read beta
     READ(file_id,*) beta
     ! read mixing parameter
     READ(file_id,*) alpha
     READ(file_id,*)

     !read freq number
     READ(file_id,*) num_freq
     ! read nhead, number of data sampling
     READ(file_id,*) nhead

     !close file
     CLOSE(file_id)

  END IF

  ! broadcast varaible
  CALL mpi_barrier(mpi_comm_world,ierr)
  CALL mpi_bcast(uj,SIZE(uj),mpi_double_precision,master,mpi_comm_world,ierr)
  CALL mpi_bcast(num_freq,1,mpi_integer,master,mpi_comm_world,ierr)
  CALL mpi_bcast(mu_dmft,1,mpi_double_precision,master,mpi_comm_world,ierr)
  CALL mpi_bcast(beta,1,mpi_double_precision,master,mpi_comm_world,ierr)
  CALL mpi_bcast(nhead,1,mpi_integer,master,mpi_comm_world,ierr)
  CALL mpi_bcast(alpha,1,mpi_double_precision,master,mpi_comm_world,ierr)

  IF (nhead < max_tail) THEN
     ntail = nhead
  ELSE
     ntail = max_tail
  ENDIF
  nlog=ntail+nhead

  ! allocate self-energe, and freqence
  ALLOCATE(sgm(num_orb,num_freq))
  ALLOCATE(freq(num_freq))
  ALLOCATE(hyb(num_orb,num_freq))
  ! allocate those in log mesh
  ALLOCATE(sgm_log(num_orb,nlog))
  ALLOCATE(sgm_old(num_orb,nlog))
  ALLOCATE(freq_log(nlog))
  ALLOCATE(gloc(num_orb,num_orb,nlog))
  ALLOCATE(hyb_log(num_orb,nlog))

  sgm  = zzero; freq = dzero
  hyb  = zzero; gloc = zzero
  sgm_log  = zzero; freq_log = dzero; hyb_log  = zzero

  ! generate frequency mesh
  CALL dmft_mkmesh()

  RETURN

END SUBROUTINE dmft_read_inp

! read self-energy
!=====================================================================!
SUBROUTINE dmft_read_sgm()

  USE mod_param,ONLY : dp,dzero, zzero
  USE mod_dmft,ONLY : file_sgm,num_aorb,num_freq,sgm,edc,num_atm,dmft1
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: file_id=1001
  INTEGER :: iorb, ifreq, otmp, iatm
  REAL(kind=dp) :: r1, i1, rtmp
  LOGICAL :: alive

  IF(iproc==master) THEN

     alive = .FALSE.
     INQUIRE(file=file_sgm, exist=alive)
     IF (alive) THEN
        !open solver.sgm.dat
        OPEN(file_id, file=file_sgm, form='formatted', status='old')
        ! read self-energe for each orbital
        DO iorb=1, num_aorb, 2
           ! read for each freqence
           DO ifreq = 1, num_freq
              READ(file_id,*) otmp, rtmp, r1, i1
              DO iatm = 1, num_atm
                 sgm(iorb+(iatm-1)*num_aorb,ifreq)=CMPLX(r1-edc,i1,dp)
              ENDDO
           END DO
           ! jump two empty lines
           READ(file_id,*)
           READ(file_id,*)
        END DO
        DO iorb=2, num_aorb, 2
           ! read for each freqence
           DO ifreq = 1, num_freq
              READ(file_id,*) otmp, rtmp, r1, i1
              DO iatm = 1, num_atm
                 sgm(iorb+(iatm-1)*num_aorb,ifreq)=CMPLX(r1-edc,i1,dp)
              ENDDO
           END DO
           ! jump two empty lines
           READ(file_id,*)
           READ(file_id,*)
        END DO
        ! close solver.sgm.dat
        CLOSE(file_id)
        dmft1 = .FALSE.
     ELSE
        sgm = zzero
        dmft1 = .TRUE.
     ENDIF

  END IF

  ! broadcast varaible
  CALL mpi_barrier(mpi_comm_world,ierr)
  CALL mpi_bcast(sgm,SIZE(sgm),mpi_double_complex,master,mpi_comm_world,&
       &ierr)

  RETURN

END SUBROUTINE dmft_read_sgm


SUBROUTINE dmft_read_edc()

  USE mod_param, ONLY : dp
  USE mod_dmft, ONLY : edc,uj,occ_lda
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(kind=dp) :: occ_read
  LOGICAL :: alive

  IF(iproc==master) THEN
     INQUIRE(file='edc.dat',exist=alive)
     IF(alive) THEN
        OPEN(5001,file='edc.dat',status='old',action='read')
        READ(5001,*) occ_read, edc
        CLOSE(5001)
        WRITE(*,"(' Occ from LDA input :', F8.4,', Expected Double-Couting :', F8.4)")&
             & occ_lda, uj(1)*(occ_lda-5.d-1)-5.d-1*uj(2)*(occ_lda-1.d0)
        WRITE(*,"(' Occ from EDC input :', F8.4,', Expected Double-Couting :', F8.4)")&
             & occ_read, uj(1)*(occ_read-5.d-1)-5.d-1*uj(2)*(occ_read-1.d0)
        WRITE(*,"(' Double-Counting from input :',F8.4)") edc
     ELSE
        edc=uj(1)*(occ_lda-5.d-1)-5.d-1*uj(2)*(occ_lda-1.d0)
     END IF
  ENDIF
  CALL mpi_bcast(edc,1,mpi_double_precision,master,mpi_comm_world,ierr)
  !===================================================================!

  RETURN

END SUBROUTINE dmft_read_edc
