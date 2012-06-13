module fastcube_mod
  ! fastcube.f90 itself does not contain modules, since f2py does not work very well with them
  ! So, here is one that provides explicit interfaces to its routines
  interface
     subroutine calculate_cubes(R, Mu, V, Ems, validmask, NI, NJ, NL, cubes, NU, NY, NX,&
          & NK, inc_degrees, xmin, ymin, umin, dx, dy, du)
       implicit none
!!! Input Arguments
       integer, intent(in) :: NI     ! Size of input radius array
       integer, intent(in) :: NJ     ! Size of input angle array
       integer, intent(in) :: NL     ! Number of emission lines
       real, intent(in), dimension(NI) :: R ! Input radius array
       real, intent(in), dimension(NJ) :: Mu ! Input angle (cosine) array
       real, intent(in), dimension(NJ, NI) :: V ! Input velocity field
       real, intent(in), dimension(NL, NJ, NI) :: Ems ! Input line emissivities
       ! which combinations of R and mu have valid data
       logical, intent(in), dimension(NJ, NI) :: validmask 
       integer, intent(in) :: NX, NY, NU ! Dimensions of output cubes
!!! Output Argument
       real, intent(out), dimension(NL, NU, NY, NX) :: cubes
!!! Auxiliary input arguments
       integer :: NK                 ! Number of phi values to use
       real :: inc_degrees           ! Inclination of proplyd axis to line of sight
       ! Limits of PPV cube
       real :: xmin, ymin, umin
       ! Note that [xyu]max are not used
       real :: dx, dy, du
     end subroutine calculate_cubes
  end interface
end module fastcube_mod



program test_fastcube
  ! Pure Fortran test of routines in fastcube.f90
  use fastcube_mod, only: calculate_cubes
  implicit none
  integer, parameter :: N = 50       ! All arrays the same size
  integer, parameter :: NI = N   ! Size of input radius array
  integer, parameter :: NJ = N   ! Size of input angle array
  integer, parameter :: NL = 5   ! Number of emission lines
  real, dimension(NI) :: R ! Input radius array
  real, dimension(NJ) :: Mu ! Input angle (cosine) array
  real, dimension(NJ, NI) :: V ! Input velocity field
  real, dimension(NL, NJ, NI) :: Ems ! Input line emissivities
  ! which combinations of R and mu have valid data
  logical, dimension(NJ, NI) :: validmask 
  integer, parameter :: NX = N, NY = N, NU = N ! Dimensions of output cubes
!!! Output Argument
  real, dimension(NL, NU, NY, NX) :: cubes
!!! Auxiliary input arguments
  integer :: NK = N                ! Number of phi values to use
  real :: inc_degrees = 30.0           ! Inclination of proplyd axis to line of sight
  real, parameter :: Rmax = 9.0
  ! Limits of PPV cube
  real :: xmin = -1.0, ymin = -1.0, umin = -1.0
  ! Note that [xyu]max are not used
  real :: dx = 0.01, dy = 0.01, du = 0.01
  integer :: i, j, k
  

  R = (/(1.0 + (Rmax - 1.0)*real(i-1)/real(NI), i = 1, NI)/)
  Mu = (/(real(j-1)/real(NJ), j = 1, NJ)/)
  V = 1.0
  Ems = 1.0
  validmask = .true. 
  call calculate_cubes(R, Mu, V, Ems, validmask, NI, NJ, NL, cubes, NU, NY, NX,&
       & NK, inc_degrees, xmin, ymin, umin, dx, dy, du)

  print *, "Sum of cubes: ", sum(cubes)

end program test_fastcube
