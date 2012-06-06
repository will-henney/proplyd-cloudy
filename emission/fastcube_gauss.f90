! Wrapper routine to be called from python
! Compiled with f2py
! Even though it is a subroutine, it appears a function in python
! 
! Note that R must be sent in dimensionless form (units of i-front radius)
subroutine calculate_cubes(R, Mu, V, T, Ems, awt, validmask, NI, NJ, NL, cubes, NU, NY, NX,&
     & NK, inc_degrees, xmin, ymin, umin, dx, dy, du)
  implicit none
!!! Input Arguments
  integer, intent(in) :: NI     ! Size of input radius array
  integer, intent(in) :: NJ     ! Size of input angle array
  integer, intent(in) :: NL     ! Number of emission lines
  real, intent(in), dimension(NI) :: R ! Input radius array
  real, intent(in), dimension(NJ) :: Mu ! Input angle (cosine) array
  real, intent(in), dimension(NJ, NI) :: V, T ! Input velocity field and temperature
  real, intent(in), dimension(NL, NJ, NI) :: Ems ! Input line emissivities
  real, intent(in), dimension(NL) :: awt         ! Atomic weight of emitting species for each line
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

  print '(a)',  "NK, inc_degrees, xmin, ymin, umin, dx, dy, du"
  print *,  NK, inc_degrees, xmin, ymin, umin, dx, dy, du
  print *, "sum(Ems) = ", sum(Ems)

  call calculate_cubes_internal

contains

  subroutine calculate_cubes_internal
    ! Actually do the work
    real, parameter :: PI = 3.1415926535897932385
    ! Local variables
    real :: u0, x, y
    real, allocatable, dimension(:) :: Phi, Ugrid, Doppler
    integer :: i, j, k
    real :: dphi, dmu, dr, dvol
    real :: cphi, sphi, ctheta, stheta, cosi, sini
    integer :: jneg, jpos, ineg, ipos
    integer :: ix, iy
    logical :: is_outside_X, is_outside_Y
    real, parameter :: xi0 = 0.128486646972 ! Doppler broadening parameter for H at 1 K
    integer :: iline
    real :: xi

    allocate( Phi(NK), Ugrid(NU), Doppler(NU) )
    cubes = 0.0
    dphi = 2.0*PI/real(NK)
    Phi = (/ (real(k)*dphi, k = 1, NK ) /)
    cosi = cos(inc_degrees*PI/180.0)
    sini = sin(inc_degrees*PI/180.0)

    Ugrid = (/(umin + real(i-1)*du, i = 1, NU)/)

    philoop: do k = 1, NK
       cphi = cos(Phi(k))
       sphi = sin(Phi(k))
       thetaloop: do j = 1, NJ-1  
          ctheta = Mu(j)
          stheta = sqrt(1.0 - ctheta**2)
          ! Trapezium rule requires half-sized dmu at end points
          jneg = max(1, j - 1)
          jpos = min(NJ, j + 1)
          dmu = -0.5*(mu(jpos) - mu(jneg))
          radiusloop: do i = 1, NI
             if (.not. validmask(j,i)) cycle radiusloop ! invalid radius for this mu
             ineg = max(1, i - 1)
             ipos = min(NI, i + 1)
             dr = 0.5*(R(ipos) - R(ineg))
             dvol =  abs(dphi * dmu * (R(i)**2) * dr)
             u0 = -V(j,i)*(sini*stheta*cphi + cosi*ctheta)
             x = R(i)*(-cosi*stheta*cphi + sini*ctheta)
             y = R(i)*stheta*sphi

             ix = 1 + int(floor((x-xmin)/dx))
             iy = 1 + int(floor((y-ymin)/dy))
             is_outside_X = ix < 1 .or. ix > NX
             is_outside_Y = iy < 1 .or. iy > NY
             ! Skip this volume element if we are outside the box
             if (is_outside_X .or. is_outside_Y) cycle radiusloop

             do iline = 1, NL
                ! Calculate Doppler broadening parameter
                xi = xi0 * sqrt(T(j, i)/awt(iline))
                ! calculate the Doppler profile
                doppler = exp( -((Ugrid-u0)/xi)**2 )
                ! Note that this should really be integrated over the width of the velocity bin
                ! Normalize the profile
                Doppler = Doppler / sum(Doppler)
                ! Add the emissions into all velocity bins of this cube
                cubes(iline, :, iy, ix) = cubes(iline, :, iy, ix) + dvol*Ems(iline, j, i)*Doppler
             end do
          end do radiusloop
       end do thetaloop
    end do philoop

    print *, "sum(cubes) = ", sum(cubes)
  end subroutine calculate_cubes_internal

end subroutine calculate_cubes

