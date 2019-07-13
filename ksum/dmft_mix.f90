!=====================================================================!
! mix with the old and new self-energy
!=====================================================================!

SUBROUTINE dmft_mix

  USE mod_param,ONLY : dp,dzero,zzero
  USE mod_dmft, ONLY : nlog,num_orb,alpha,sgm_old,sgm_log
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: file_id=1001
  INTEGER :: ifreq, iorb

  CALL dmft_read_sgm_old()

  ! check self-energy from CTQMC, calculate differece between last two steps before mix
  IF(iproc==master) CALL dmft_sechk()

  ! only mix diagonal part of self-energy
  DO ifreq = 1, nlog
     DO iorb = 1, num_orb
        sgm_log(iorb,ifreq)=sgm_log(iorb,ifreq)*alpha+sgm_old(iorb,ifreq)*(1.d0-alpha)
     END DO
  END DO

  RETURN

END SUBROUTINE dmft_mix

! read self-energy
!=====================================================================!
SUBROUTINE dmft_read_sgm_old()

  USE mod_param,ONLY : dp,dzero
  USE mod_dmft,ONLY  : file_osgm,num_orb,nlog,sgm_old,sgm_log,dmft1
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: file_id=1001
  INTEGER :: iorb, ifreq, otmp
  REAL(kind=dp) :: r1, i1, ow
  LOGICAL :: alive

  IF(iproc==master) THEN

     ! check if file_osgm alive or not
     INQUIRE(file=file_osgm, exist=alive)

     IF( dmft1 .OR. (.NOT. alive) ) THEN
        ! give SE this loop as oold
        sgm_old=sgm_log
        ! if exist, read hyb last loop
     ELSE
        !open solver.sgm.dat
        OPEN(file_id, file=file_osgm, form='formatted', status='old')
        ! read self-energe for each orbital
        DO iorb=1, num_orb
           ! read for each freqence
           DO ifreq = 1, nlog
              READ(file_id,*) otmp, ow, r1, i1
              sgm_old(iorb,ifreq)=CMPLX(r1,i1,dp)
           END DO
           ! jump two empty lines
           READ(file_id,*)
           READ(file_id,*)
        END DO
        ! close solver.sgm.dat
        CLOSE(file_id)
     ENDIF

  END IF

  ! broadcast varaible
  CALL mpi_barrier(mpi_comm_world,ierr)
  CALL mpi_bcast(sgm_old,SIZE(sgm_old),mpi_double_complex,master,mpi_comm_world,&
       &ierr)

  RETURN

END SUBROUTINE dmft_read_sgm_old
