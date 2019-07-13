!=========>=========>=========>=========>=========>=========>=========!
! This is ksum program for DMFT
!=========>=========>=========>=========>=========>=========>=========!
PROGRAM dmft_main

  USE mod_param, ONLY : dp
  USE mod_mpi

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  REAL(kind=dp) :: tstart,tfinish, time1, time2

  ! MPI initialize
  !===================================================================!
  CALL mpi_ini()

  ! program start time
  tstart=MPI_WTIME()
  !===================================================================!

  ! read variable
  ! initialize variable
  !===================================================================!
  time1=MPI_WTIME()

  ! read input file and initialize
  CALL dmft_read()

  time2=MPI_WTIME()

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  IF(iproc==0) THEN
     CALL dmft_time(time2-time1,'dmft_read')
  END IF

  ! initialize viarable
  CALL dmft_init()

  ! get dmft chemical potential
  !===================================================================!
  time1=MPI_WTIME()

  CALL dmft_mu()

  CALL dmft_dm()

  time2=MPI_WTIME()

  IF(iproc==master) THEN
     CALL dmft_time(time2-time1,'dmft_mu')
  END IF

  !===================================================================!

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! get Green's function on local state
  !===================================================================!

  time1=MPI_WTIME()

  CALL dmft_green()

  time2=MPI_WTIME()

  IF(iproc==master) THEN
     WRITE(6,*)
     CALL dmft_time(time2-time1,'dmft_green')
  END IF
  !===================================================================!

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! write output file, and deallocate arrays
  !===================================================================!
  IF(iproc==master) THEN
     CALL dmft_output()
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  CALL dmft_finalize()

  ! program finish time
  tfinish=MPI_WTIME()
  !===================================================================!

  ! write time
  IF(iproc==master) THEN
     WRITE(6,*)
     CALL dmft_time(tfinish-tstart,'total')
     WRITE(6,*)
  END IF

  CALL mpi_fnl(0)

END PROGRAM dmft_main
