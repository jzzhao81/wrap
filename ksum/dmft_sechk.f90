!=====================================================================!
! check self-energy & write difference in dmft.step.dat
!=====================================================================!

SUBROUTINE dmft_sechk()

  USE mod_param, ONLY : dp,zzero
  USE mod_dmft, ONLY : sgm_old, sgm_log, dmft1, num_orb, nloop, nlog, num_atm

  IMPLICIT NONE

  INTEGER, PARAMETER :: file = 1001
  CHARACTER(len=13),PARAMETER :: file_step='dmft.step.dat'
  LOGICAL :: alive
  INTEGER :: iorb, ifreq
  REAL(kind=dp) :: rpart, ipart
  ! sum of self-energy for new & old
  COMPLEX(kind=dp) :: snew, sdiff

  ! sum new self-energy
  snew = zzero
  DO ifreq = 1, nlog
     DO iorb = 1, num_orb/num_atm
        rpart = ABS(REAL(sgm_log(iorb,ifreq)))
        ipart = ABS(AIMAG(sgm_log(iorb,ifreq)))
        snew = snew + CMPLX( rpart, ipart, dp )
     ENDDO
  ENDDO
  snew = snew/REAL(num_orb,dp)/REAL(nlog,dp)

  ! diff new & old
  ! Sum self-energy difference
  sdiff = zzero
  DO ifreq = 1, nlog
     DO iorb = 1, num_orb/num_atm
        rpart = ABS(REAL(sgm_log(iorb,ifreq))-REAL(sgm_old(iorb,ifreq)))
        ipart = ABS(AIMAG(sgm_log(iorb,ifreq))-AIMAG(sgm_old(iorb,ifreq)))
        sdiff = sdiff + CMPLX( rpart, ipart, dp )
     ENDDO
  ENDDO
  sdiff = sdiff/REAL(num_orb/num_atm,dp)/REAL(nlog,dp)

  ! check if file exist
  INQUIRE(file=file_step,exist=alive)

  ! if file exist, count number of loop before
  !===================================================================!
  IF(dmft1) THEN
     IF(alive) THEN
        OPEN(file,file=file_step,status='old')
        CLOSE(file,status='DELETE')
     ENDIF
     OPEN(file,file=file_step,form='formatted',status='new')
     WRITE(file,"(' loop       Self-energy ',11X,'  diff')")
     CLOSE(file)
  ENDIF
  !===================================================================!

  ! Determine the loop number
  !===================================================================!
  nloop = -1
  OPEN (file, file = file_step)
  DO
     READ (file,*, END=2001)
     nloop = nloop + 1
  END DO
2001 CLOSE (file)
  !===================================================================!

  ! Write self-energy difference
  !===================================================================!
  OPEN(file,file=file_step,form='formatted',status='old',position='APPEND')
  WRITE(file,"(1X,I3,4F13.3)")nloop,REAL(snew),AIMAG(snew),REAL(sdiff),&
       &AIMAG(sdiff)
  CLOSE(file)
  !===================================================================!

  RETURN

END SUBROUTINE dmft_sechk
