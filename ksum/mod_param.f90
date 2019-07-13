!=====================================================================!
! Parameter module
!=====================================================================!
module mod_param

  ! Double precision
  integer, public, parameter :: dp = kind(1.0d0)

  ! stanard output
  integer, public, parameter :: stdout = 6

  ! maximum search chemical potential steps
  integer, public, parameter :: max_mu_step=50

  ! criteria for search mu
  real(kind=dp),public,parameter :: occ_eps=1.d-8

  ! Hartree to eV
  real(kind=dp),public,parameter :: ha2ev=27.21138386_dp

  ! Pi
  real(kind=dp),public,parameter :: pi=3.1415926535898_dp
  real(kind=dp),public,parameter :: srpi=1.772453850905516_dp

  ! eV 2 K
  real(kind=dp), public, parameter :: ev2k=11604.505008098d0

  ! zero in complex number
  real(kind=dp), public,parameter :: dzero=0.d0
  complex(kind=dp),public,parameter :: zzero=cmplx(dzero,dzero,dp)
  real(kind=dp), public,parameter :: dhalf=5.d-1

  ! generate standard
  real(kind=dp),parameter :: gen_std=1.d-4

  ! Euler's number 
  real(kind=dp) :: lne=2.718281828d0

end module mod_param
