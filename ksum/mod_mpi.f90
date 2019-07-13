!define some mpi global viables

MODULE mod_mpi

  ! ierr  : Error message
  ! iproc : Process index
  ! nproc : number of processes
  INTEGER, PUBLIC :: ierr, iproc, nproc, errorcode

  ! master process
  INTEGER, PARAMETER, PUBLIC :: master=0

  ! Number of K-Points per-nodes
  INTEGER, DIMENSION(:), ALLOCATABLE :: nkpp

  ! K-Points displs per-nodes
  INTEGER, DIMENSION(:), ALLOCATABLE :: disp

  ! Number data per node
  INTEGER, DIMENSION(:), ALLOCATABLE :: ndpp

  ! Total displs
  INTEGER, DIMENSION(:), ALLOCATABLE :: ndsp

  ! Here are some subroutines for mpi

CONTAINS

  ! mpi initialize, called at the beginning of the program

  SUBROUTINE mpi_ini()

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    CHARACTER(len=3) :: MON(12)
    DATA mon/'Jan','Feb','Mar','Apr','May','Jun', &
         &         'Jul','Aug','Sep','Oct','Nov','Dec'/
    CHARACTER(len=8)      :: date
    CHARACTER(len=10)     :: time
    CHARACTER(len=5)      :: zone
    INTEGER, DIMENSION(8) :: value

    !initialize MPI
    !===================================================================!
    CALL mpi_init(ierr)
    CALL mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)
    CALL mpi_comm_rank(MPI_COMM_WORLD,iproc,ierr)
    CALL DATE_AND_TIME(date,time,zone,value)
    IF(iproc==master) THEN
       WRITE(6,*)
       ! Print Start time
       WRITE(6,201) mon(value(2)),value(3),value(1),value(5:7)
       WRITE(6,"(' Total Number of CPUs : ', I4)")nproc
       WRITE(6,*) 'The more CPUs, the faster !'
       WRITE(6,*) 'Good luck !'
       WRITE(6,*)
    END IF
    !===================================================================!

    ! allocate Number of K-points per-node
    ALLOCATE(nkpp(0:nproc-1), disp(0:nproc-1))
    ALLOCATE(ndpp(0:nproc-1), ndsp(0:nproc-1))

    nkpp = 0
    disp = 0
    ndpp = 0
    ndsp = 0

    RETURN

201 FORMAT(' KSUM start at :  ', A3,1X,I2,1X,I4,',',2X,2(I2.2,':'),I2.2)

  END SUBROUTINE mpi_ini

  ! mpi finalization, called at the end of the program

  SUBROUTINE mpi_fnl(ecode, msg)

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    INTEGER, INTENT(in) :: ecode
    CHARACTER(*), INTENT(in), OPTIONAL :: msg

    DEALLOCATE(nkpp, disp)
    DEALLOCATE(ndpp, ndsp)

    IF ( iproc == master ) THEN
       IF ( ecode /= 0 ) THEN
          WRITE(*, "(' KSUM STOP with Errors ! ')")
          IF (PRESENT(msg)) THEN
             WRITE(*,"(' Error message : ')")
             WRITE(*,*) TRIM(ADJUSTL(msg))
          END IF
          WRITE(*,*)

          CALL mpi_abort(mpi_comm_world, ierr)

       ELSE
          WRITE(*, "(' KSUM STOP Normally ! ')")
          WRITE(*,*)

       ENDIF
    ENDIF

    CALL mpi_finalize(ierr)
    STOP

    RETURN

  END SUBROUTINE mpi_fnl

  ! Give value for nkpp & disp

  SUBROUTINE mpi_begin(ntot,ndat)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ntot, ndat
    INTEGER :: icpu

    nkpp = 0
    disp = 0
    ndpp = 0
    ndsp = 0

    DO icpu = 0, nproc-1
       IF(icpu < MOD(ntot,nproc)) THEN
          nkpp(icpu) = ntot / nproc + 1
       ELSE
          nkpp(icpu) = ntot / nproc
       END IF
       ndpp(icpu) = nkpp(icpu) * ndat
       IF(icpu == 0 ) THEN
          disp(icpu) = 0
       ELSE
          disp(icpu) = disp(icpu-1) + nkpp(icpu-1)
       END IF
       ndsp(icpu) = disp(icpu) * ndat
    END DO

    RETURN

  END SUBROUTINE mpi_begin




END MODULE mod_mpi
