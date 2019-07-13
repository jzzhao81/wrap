
!=========>=========>=========>=========>=========>=========>=========!
! dmft main output subroutine
!=========>=========>=========>=========>=========>=========>=========!
SUBROUTINE dmft_output()

  IMPLICIT NONE

  ! write out occupation for local orbital
  !===================================================================!
  CALL dmft_output_occ()

  ! open main output file
  !===================================================================!
  ! CALL dmft_output_main()

  !write hybrid function back to ctqmc solver
  !===================================================================!
  CALL dmft_output_hyb()

  ! write dmft.gloc.dat, here we both give real and imag part
  !===================================================================!
  CALL dmft_output_gloc()

  ! write eimp to CTQMC solver
  !===================================================================!
  CALL dmft_output_eimp()

  !write density matrix back to bstate
  !===================================================================!
  ! CALL dmft_output_state

  ! CALL dmft_output_qmc()

  CALL dmft_output_edc()

  CALL dmft_output_sgm()

  RETURN

END SUBROUTINE dmft_output

! write time to standard output
!=====================================================================!
SUBROUTINE dmft_time(rtime,name)

  USE mod_param,ONLY : dp
  IMPLICIT NONE

  REAL(kind=dp),INTENT(in) :: rtime
  CHARACTER(len=*),INTENT(in) :: name
  INTEGER :: hr, min
  REAL(kind=dp) :: sec
  CHARACTER(len=20) :: fmt

  hr=INT(rtime/36.d2)
  min=INT((rtime-DBLE(hr)*36.d2)/6.d1)
  sec=MOD(rtime,6.d1)

  WRITE(fmt,*) LEN(name)

  IF(hr/=0) THEN
     WRITE(6,"(' Time for ',A"//ADJUSTL(fmt)//",':',I3,' hr ',I3,' min&
          &',F6.2,' sec !')") name,hr,min,sec
  ELSEIF(min/=0) THEN
     WRITE(6,"(' Time for ',A"//ADJUSTL(fmt)//",':',I3,' min ',F5.2,' &
          &sec !')") name,min,sec
  ELSE
     WRITE(6,"(' Time for ',A"//ADJUSTL(fmt)//",':',F6.2,' sec !')") n&
          &ame,sec
  END IF

  RETURN

END SUBROUTINE dmft_time


! write input file for cqtmc solver
! only change the chemical potential
! $$$$   CAUTION   $$$$
!--------------------------------------------
! In Li Huang's code
! version         Mu
!--------------------------------------------
! azalea          contains in Eimp
! lavender        contains in solver.ctqmc.in
!=====================================================================!
SUBROUTINE dmft_output_qmc()

  USE mod_dmft, ONLY : file_inp

  IMPLICIT NONE

  INTEGER, PARAMETER :: file_id=1001
  INTEGER, PARAMETER :: muline=22
  INTEGER, PARAMETER :: maxline=39

  INTEGER :: iline
  CHARACTER(len=10) :: mu
  CHARACTER(len=75), ALLOCATABLE :: inline(:)

  ! orginal input
  ALLOCATE(inline(maxline))

  ! write(mu,"(F10.6)") mu_dmft
  WRITE(mu,"(F10.6)") 0.d0

  OPEN(file_id,file=file_inp)

  ! read original input
  DO iline=1,maxline
     READ(file_id,"(A75)") inline(iline)
  END DO

  ! delete orginal input file
  CLOSE(file_id,status='delete')

  ! open new input file
  OPEN(file_id,file=file_inp,status='new')

  ! write new input
  DO iline=1,maxline

     ! write chemical potentital keep other unchange
     IF(iline==muline) THEN
        WRITE(file_id,"(A10,3X,'! Chemical Potential')") ADJUSTL(mu)
     ELSE
        WRITE(file_id,"(A75)") inline(iline)
     END IF

  END DO

  ! close input file
  CLOSE(file_id)

  ! remove input
  DEALLOCATE(inline)

  RETURN

END SUBROUTINE dmft_output_qmc

! write eimp to solver.eimp.in
!=====================================================================!
SUBROUTINE dmft_output_eimp()

  USE mod_param,ONLY : gen_std,dp
  USE mod_dmft,ONLY : file_eimp,num_aorb,num_orb,eimp,edc,mu_dmft

  IMPLICIT NONE

  INTEGER,ALLOCATABLE :: symm(:)
  INTEGER,PARAMETER :: file_id=1001
  INTEGER :: inum, iorb, jorb, gen_flag=0

  ALLOCATE(symm(num_orb))

  ! determine generated energy level
  !===================================================================!
  symm(1)=1
  DO iorb=2,num_aorb
     gen_flag=0
     DO jorb=1,iorb-1
        IF(ABS(REAL(eimp(iorb,iorb))-REAL(eimp(jorb,jorb))).LT.gen_std)&
             &THEN
           symm(iorb)=symm(jorb)
           gen_flag=1
           EXIT
        END IF
     END DO
     IF(gen_flag==0) THEN
        symm(iorb)=MAXVAL(symm(1:iorb-1))+1
     END IF
  END DO
  !===================================================================!

  !write energy level back to ctqmc solver
  !===================================================================!
  OPEN(file_id,file=file_eimp,form='formatted',status='unknown')
  inum=1
  DO iorb=1,num_aorb,2
     WRITE(file_id,'(i5,f16.8,i5)') inum, REAL(eimp(iorb,iorb))-edc,&
          &symm(iorb)
     inum=inum+1
  END DO
  DO iorb=2,num_aorb,2
     WRITE(file_id,'(i5,f16.8,i5)') inum, REAL(eimp(iorb,iorb))-edc,&
          &symm(iorb)
     inum=inum+1
  END DO
  CLOSE(file_id)
  DEALLOCATE(symm)

  !write energy level back to ctqmc solver
  !===================================================================!
  OPEN(file_id,file='dmft.eimp.dat',form='formatted',status='unknown')
  DO iorb=1,num_orb,2
     WRITE(file_id,'(i5,6f16.8)') (iorb+1)/2, edc, mu_dmft, &
          &REAL(eimp(iorb,iorb)),REAL(eimp(iorb+1,iorb+1)), &
          &REAL(eimp(iorb,iorb))-edc-mu_dmft,REAL(eimp(iorb+1,iorb+1))-edc-mu_dmft
  END DO
  CLOSE(file_id)

  RETURN

END SUBROUTINE dmft_output_eimp

! write dmft main output
!=====================================================================!
! SUBROUTINE dmft_output_main
!
!   USE mod_param,ONLY : ev2k
!   USE mod_dmft,ONLY : file_out,lmoment,num_kp,num_bnd,num_ele,num_orb,&
!        &occ_lda,uj,edc,mu_dmft,num_freq,alpha,beta,eip,etot
!
!   IMPLICIT NONE
!
!   INTEGER,PARAMETER :: file_id=1001
!   CHARACTER(len=1) :: name_orb
!
!   OPEN(file_id,file=file_out,status='unknown')
!
!   ! determin orbital name, just for fun, useless
!   SELECT CASE(lmoment)
!   CASE(0)
!      name_orb='s'
!   CASE(1)
!      name_orb='p'
!   CASE(2)
!      name_orb='d'
!   CASE(3)
!      name_orb='f'
!   CASE default
!      WRITE(6,*) 'Wrong angle moment ! '
!      WRITE(6,"(' L = ',I1)")lmoment
!      STOP
!   END SELECT
!
!   WRITE(file_id,*)
!   WRITE(file_id,"(' Target orbital : ',a1)")name_orb
!
!   WRITE(file_id,"(' Number of k-points : ',i6)") num_kp
!   WRITE(file_id,"(' Number of bands from LDA : ',i3)") num_bnd
!   WRITE(file_id,"(' Number of total valence electron : ',f6.1)") num_ele
!   WRITE(file_id,"(' Number of correlated orbitals : ',i2)") num_orb
!
!   WRITE(file_id,"(' Occupation from LDA : ',f8.4)") occ_lda
!
!   WRITE(file_id,"(' Number of matsubara frequency :',I5)") num_freq
!   WRITE(file_id,"(' Inverse of temperature (beta) :',F6.2)") beta
!   WRITE(file_id,"(' Temperature : ',F6.1,' K')") ev2k/beta
!   WRITE(file_id,"(' Mixing parameter :',F6.2)") alpha
!   WRITE(file_id,"(' Chemical potential (mu) = ',F12.8)") mu_dmft
!
!   WRITE(file_id,*)
!   WRITE(file_id,"(' Energy information:')")
!   WRITE(file_id,"(' U & J : ',f8.4, f8.4)") uj(1), uj(2)
!   WRITE(file_id,"(' Double-couning :'F8.4,' eV')") edc
!   WRITE(file_id,"(' First moment of Green''s function : ',F10.4,' eV')") eip
!   WRITE(file_id,"(' Total energy from DMFT : ',F10.4,' eV')") etot
!
!   !close output file
!   CLOSE(file_id)
!
!   RETURN
!
! END SUBROUTINE dmft_output_main

! output occupation info
!=====================================================================!
SUBROUTINE dmft_output_occ

  USE mod_param,ONLY : dp,dzero
  USE mod_dmft,ONLY : file_occ,occ_dmft,num_atm,num_aorb

  IMPLICIT NONE

  INTEGER, PARAMETER :: file_id=1001
  INTEGER :: iorb, iatm
  REAL(kind=dp) :: oup,odn, otot

  otot=dzero

  ! open occ output file
  OPEN(file_id,file=file_occ,status='unknown')

  DO iatm = 1, num_atm
     oup=dzero; odn=dzero
     ! write occ for each orbital
     WRITE(file_id,"(1X,'Orbital ',5X,' Up ',10X,' Down')")
     DO iorb=1,num_aorb,2
        WRITE(file_id,"(3X,I2,4X,F12.8,3X,F12.8)")(iorb+1)/2,&
             &REAL(occ_dmft((iatm-1)*num_aorb+iorb)),&
             &REAL(occ_dmft((iatm-1)*num_aorb+iorb+1))
        oup=oup+REAL(occ_dmft((iatm-1)*num_aorb+iorb)); odn=odn+REAL(occ_dmft((iatm-1)*num_aorb+iorb+1))
     END DO
     ! write total up & down occ
     WRITE(file_id,"(2X,'total',2X,F12.8,3X,F12.8)")oup,odn
  ENDDO

  ! occ total
  otot=REAL(SUM(occ_dmft))

  ! write occ total
  WRITE(file_id,*)
  WRITE(file_id,"(' Occupation in DMFT (occ_dmft) : ' F12.8)") otot

  ! close file dmft.occ.dat
  CLOSE(file_id)

  RETURN

END SUBROUTINE dmft_output_occ

! write hybrid function
!=====================================================================!
SUBROUTINE dmft_output_hyb()

  USE mod_param,ONLY : dzero
  USE mod_dmft,ONLY : file_ihyb,file_ohyb,num_aorb,num_orb,num_freq,hyb,freq,&
       &nlog,hyb_log,freq_log

  IMPLICIT NONE

  INTEGER,PARAMETER :: file_id=1001
  INTEGER :: iorb,ifreq

  ! write solver.hyb.in
  !===================================================================!
  CALL dmft_interpolate_hyb()

  ! open data file: solver.hyb.in, we set real part to zero
  OPEN(file_id, file=file_ihyb, form='formatted', status='unknown')

  ! write it
  DO iorb=1,num_aorb,2
     DO ifreq=1,num_freq
        WRITE(file_id,'(i5,5e18.10)') iorb, freq(ifreq), &
             REAL(hyb(iorb,ifreq)),AIMAG(hyb(iorb,ifreq)),dzero,dzero
     ENDDO ! over j={1,num_freq} loop
     WRITE(file_id,*) ! write empty lines
     WRITE(file_id,*)
  ENDDO ! over i={1,num_orb} loop
  ! write it
  DO iorb=2,num_aorb,2
     DO ifreq=1,num_freq
        WRITE(file_id,'(i5,5e18.10)') iorb, freq(ifreq), &
             REAL(hyb(iorb,ifreq)),AIMAG(hyb(iorb,ifreq)),dzero,dzero
     ENDDO ! over j={1,num_freq} loop
     WRITE(file_id,*) ! write empty lines
     WRITE(file_id,*)
  ENDDO ! over i={1,num_orb} loop

  ! close data file
  CLOSE(file_id)
  !===================================================================!

  ! write dmft.hyb.dat
  !===================================================================!
  ! write dmft.hyb.dat, here we both give real and imag part
  OPEN(file_id, file=file_ohyb, form='formatted', status='unknown')

  ! write it
  DO iorb=1,num_orb,2
     DO ifreq=1,nlog
        WRITE(file_id,'(i5,5e18.10)') (iorb+1)/2, freq_log(ifreq), &
             REAL(hyb_log(iorb,ifreq)),AIMAG(hyb_log(iorb,ifreq)), &
             REAL(hyb_log(iorb+1,ifreq)), AIMAG(hyb_log(iorb+1,ifreq))
     ENDDO ! over j={1,num_freq} loop
     WRITE(file_id,*) ! write empty lines
     WRITE(file_id,*)
  ENDDO ! over i={1,num_orb} loop

  ! close data file
  CLOSE(file_id)
  !===================================================================!

  RETURN

END SUBROUTINE dmft_output_hyb

! write self-energy to file, for mix in next step.
! sgm = sgm - edc
!=====================================================================!
SUBROUTINE dmft_output_sgm()

  USE mod_dmft,ONLY : file_osgm,num_orb,nlog,sgm_log,freq_log

  IMPLICIT NONE

  INTEGER,PARAMETER :: file_id=1001
  INTEGER :: iorb,ifreq

  ! write dmft.sgm.dat
  !===================================================================!
  ! write dmft.sgm.dat, here we both give real and imag part
  OPEN(file_id, file=file_osgm, form='formatted', status='unknown')
  ! write it
  DO iorb=1,num_orb
     DO ifreq=1,nlog
        WRITE(file_id,'(i5,5e18.10)') iorb, freq_log(ifreq), &
             REAL(sgm_log(iorb,ifreq)), AIMAG(sgm_log(iorb,ifreq))
     ENDDO ! over j={1,num_freq} loop
     WRITE(file_id,*) ! write empty lines
     WRITE(file_id,*)
  ENDDO ! over i={1,num_orb} loop
  ! close data file
  CLOSE(file_id)
  !===================================================================!

  RETURN

END SUBROUTINE dmft_output_sgm

! write gloc to dmft.gloc.dat
!=====================================================================!
SUBROUTINE dmft_output_gloc()

  USE mod_dmft,ONLY : file_gloc,num_orb,nlog,freq_log,gloc,num_atm

  IMPLICIT NONE

  INTEGER,PARAMETER :: file_id=1001
  INTEGER :: iorb,ifreq

  OPEN(file_id, file=file_gloc, form='formatted', status='unknown')

  ! write it
  DO iorb=1,num_orb/num_atm,2
     DO ifreq=1,nlog
        WRITE(file_id,'(i5,5e18.10)') (iorb+1)/2, freq_log(ifreq), &
             REAL(gloc(iorb,iorb,ifreq)), &
             AIMAG(gloc(iorb,iorb,ifreq)), &
             REAL(gloc(iorb+1,iorb+1,ifreq)), &
             AIMAG(gloc(iorb+1,iorb+1,ifreq))
     ENDDO ! over j={1,num_freq} loop
     WRITE(file_id,*) ! write empty lines
     WRITE(file_id,*)
  ENDDO ! over i={1,num_orb} loop

  ! close data file
  CLOSE(file_id)

  RETURN

END SUBROUTINE dmft_output_gloc

! write bstate input
!=====================================================================!
SUBROUTINE dmft_output_state

  USE mod_param,ONLY : dp
  USE mod_dmft,ONLY : file_dm,num_kp,rho,mu_dmft,nemin,nemax,mu_lda

  IMPLICIT NONE

  INTEGER,PARAMETER :: file_id=1101

  ! integer  :: ikp,isymop,ibnd,jbnd
  ! real(kind=8) :: eps

  ! INTEGER    ::  lwork, lrwork, liwork
  ! complex*16 :: vec(num_bnd,num_bnd)
  ! real*8     :: val(num_bnd)
  ! complex*16, allocatable :: work(:)
  ! real*8,     allocatable :: rwork(:)
  ! integer,    allocatable :: iwork(:)
  ! integer :: info

  ! write to file dmft_dm.dat
  OPEN(file_id,file=file_dm,form='unformatted')
  WRITE(file_id) num_kp,nemin,nemax
  WRITE(file_id) mu_lda+mu_dmft
  WRITE(file_id) rho
  CLOSE(file_id)

  RETURN

END SUBROUTINE dmft_output_state

SUBROUTINE dmft_output_edc

  USE mod_param,ONLY : dp
  USE mod_dmft,ONLY : edc, occ_lda

  IMPLICIT NONE

  LOGICAL :: alive

  ! write double counting to file
  INQUIRE(file='edc.dat',exist=alive)
  IF(.NOT. alive) THEN
     OPEN(5001,file='edc.dat',status='new')
     WRITE(5001,"(2F10.5)") occ_lda, edc
     CLOSE(5001)
  ENDIF

  RETURN

END SUBROUTINE dmft_output_edc


SUBROUTINE dmft_interpolate_hyb()

  USE mod_param, ONLY : dzero, dp
  USE mod_dmft , ONLY : num_freq,num_orb,nlog,ntail,hyb,hyb_log,freq,freq_log

  IMPLICIT NONE

  INTEGER :: iorb
  INTEGER :: istart,nint    ! number for interpolation
  REAL(kind=dp),ALLOCATABLE :: real_array(:), imag_array(:)

  nint=num_freq-nlog/2
  istart=1+nlog/2

  ! allocate a temp array
  ALLOCATE(real_array(istart:num_freq))
  ALLOCATE(imag_array(istart:num_freq))
  real_array = dzero
  imag_array = dzero

  hyb(:,1:nlog/2)=hyb_log(:,1:nlog/2)

  DO iorb=1,num_orb
     CALL math_spline(ntail,freq_log(istart:nlog),REAL(hyb_log(iorb,istart:nlog)),&
          &nint,freq(istart:num_freq),real_array(istart:num_freq))
     CALL math_spline(ntail,freq_log(istart:nlog),AIMAG(hyb_log(iorb,istart:nlog)),&
          &nint,freq(istart:num_freq),imag_array(istart:num_freq))
     hyb(iorb,istart:num_freq)=CMPLX(real_array(istart:num_freq),&
          &imag_array(istart:num_freq),dp)
  END DO

  ! deallocate the temp-array
  DEALLOCATE(real_array)
  DEALLOCATE(imag_array)

END SUBROUTINE dmft_interpolate_hyb
