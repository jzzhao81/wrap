!=====================================================================!
! definition
!=====================================================================!
MODULE mod_dmft

  USE mod_param,ONLY : dp

  ! the following variables used for log mesh
  ! nhead is the number should keep original point
  ! ntail is the number should keep in log mesh
  INTEGER  :: nlog, nhead, ntail

  ! numbers
  INTEGER :: num_freq, num_kp, num_bnd, num_orb, num_atm, num_aorb, lmoment
  INTEGER :: nemin, nemax ! band range from LDA, this is determined by energy cut off

  ! number of spin
  INTEGER :: nspin

  ! logical parameter first dmft loop
  LOGICAL :: dmft1
  INTEGER :: nloop

  ! num_ele = number of valence eletron
  ! uj(2) = U & J
  ! occ_lda = real(sum(nflda))
  ! occ_dmft = local occupation in dmft
  ! edc = double counting
  ! mu_dmft = chemical potential
  REAL(kind=dp) :: num_ele, uj(2), occ_lda, mu_dmft, mu_lda, edc
  COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: occ_dmft

  ! beta=1/Kb*T
  REAL(kind=dp) :: beta

  ! self-energy in local orbital
  COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE :: sgm
  COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE :: sgm_log
  COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE :: sgm_old

  ! array for freq from CTQMC
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: freq
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: freq_log

  ! eigenvalue from LDA, here we assume it's a complex number
  COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE :: ek

  ! k-point weight
  REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: wk

  ! energy level
  COMPLEX(kind=dp), DIMENSION(:,:), ALLOCATABLE :: eimp

  ! overlap matrix <i|k>
  COMPLEX(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: smtrx

  ! Local Green's function, We only store the log mesh
  COMPLEX(kind=dp),DIMENSION(:,:,:),ALLOCATABLE :: gloc

  ! hybrid's function
  COMPLEX(kind=dp),DIMENSION(:,:),ALLOCATABLE :: hyb
  COMPLEX(kind=dp),DIMENSION(:,:),ALLOCATABLE :: hyb_log

  ! density matrix
  COMPLEX(kind=dp),DIMENSION(:,:,:),ALLOCATABLE :: rho

  ! mix parameter
  REAL(kind=dp) :: alpha

  ! file name parameter
  !===================================================================!
  ! lda input
  CHARACTER(len=12), PARAMETER :: file_lda='dmft_inf.dat'

  ! dmft input
  CHARACTER(len=15), PARAMETER :: file_inp='dmft.main.in'

  ! solver input
  CHARACTER(len=15), PARAMETER :: solver_inp='solver.ctqmc.in'

  ! self-energy input
  CHARACTER(len=14), PARAMETER :: file_sgm='solver.sgm.dat'

  ! self-energy output, for mix purpose
  CHARACTER(len=12), PARAMETER :: file_osgm='dmft.sgm.dat'

  ! eimp input
  CHARACTER(len=14), PARAMETER :: file_eimp='solver.eimp.in'

  ! main output
  CHARACTER(len=13), PARAMETER :: file_out='dmft.main.out'

  ! occ output
  CHARACTER(len=12), PARAMETER :: file_occ='dmft.occ.dat'

  ! hyb output
  CHARACTER(len=13), PARAMETER :: file_ihyb='solver.hyb.in'
  CHARACTER(len=12), PARAMETER :: file_ohyb='dmft.hyb.dat'

  ! gloc output
  CHARACTER(len=13), PARAMETER :: file_gloc='dmft.gloc.dat'

  ! bstate input
  CHARACTER(len=11), PARAMETER :: file_dm='dmft_dm.dat'
  !===================================================================!

  ! energy variable
  !===================================================================!
  ! impurity potential energy
  REAL(kind=dp) :: eip
  REAL(kind=dp) :: tdc
  REAL(kind=dp) :: etot
  !===================================================================!

  ! occ from solver
  REAL(kind=dp) :: socc

  ! if perform dmft band
  LOGICAL :: ifband = .FALSE.

END MODULE mod_dmft
