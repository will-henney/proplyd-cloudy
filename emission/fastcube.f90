! Wrapper routine to be called from python
! Compiled with f2py
! Even though it is a subroutine, it appears a function in python
! 
! Note that R must be sent in dimensionless form (units of i-front radius)
subroutine calculate_cubes(R, Mu, V, Ems, NI, NJ, NL, cubes, NX, NY, NU)
  implicit none
!!! Input Arguments
  integer, intent(in) :: NI     ! Size of input radius array
  integer, intent(in) :: NJ     ! Size of input angle array
  integer, intent(in) :: NL     ! Number of emission lines
  real, intent(in), dimension(NI) :: R ! Input radius array
  real, intent(in), dimension(NJ) :: Mu ! Input angle (cosine) array
  real, intent(in), dimension(NI, NJ) :: V ! Input velocity field
  real, intent(in), dimension(NL, NI, NJ) :: Ems ! Input line emissivities
  integer, intent(in) :: NX, NY, NU ! Dimensions of output cubes
!!! Output Argument
  real, intent(out), dimension(NL, NX, NY, NU) :: cubes

  call calculate_cubes_internal

contains

  subroutine calculate_cubes_internal
    ! Actually do the work
    real, parameter :: PI = 3.1415926535897932385
    ! FIXME: Hardwire these parameters for now
    ! FIXME: eventually they should be passed as arguments
    integer, parameter :: NK = 400
    real, parameter :: xmin = -1.0, xmax = 3.0, ymin = -2.0, ymax = 2.0
    real, parameter :: umin = -50.0, umax = 50.0
    real, parameter :: inclination = 60.0
    ! Local variables
    real :: dx, dy, du, u, x, y
    real, dimension(NK) :: Phi
    integer :: i, j, k
    real :: dphi, dmu, dr, dvol
    real :: cphi, sphi, ctheta, stheta, cosi, sini
    integer :: jneg, jpos, ineg, ipos
    integer :: iu, ix, iy
    logical :: is_outside_V, is_outside_X, is_outside_Y

    cubes = 0.0
    dphi = 2.0*PI/real(NK)
    Phi = (/ (real(k)*dphi, k = 1, NK ) /)
    cosi = cos(inclination*PI/180.0)
    sini = sin(inclination*PI/180.0)
    dx = (xmax - xmin)/real(NX-1)
    dy = (ymax - ymin)/real(NY-1)
    du = (umax - umin)/real(NU-1)

    philoop: do k = 1, NK
       cphi = cos(Phi(k))
       sphi = sin(Phi(k))
       thetaloop: do j = 1, NJ
          ctheta = Mu(j)
          stheta = sqrt(1.0 - ctheta**2)
          ! Trapezium rule requires half-sized dmu at end points
          jneg = max(1, j - 1)
          jpos = min(NJ, j + 1)
          dmu = -0.5*(mu(jpos) - mu(jneg))
          radiusloop: do i = 1, NI
             if (V(i,j) == -1.0) cycle radiusloop
             ineg = max(1, i - 1)
             ipos = min(NI, i + 1)
             dr = 0.5*(R(ipos) - R(ineg))
             dvol =  abs(dphi * dmu * (R(i)**2) * dr)
             u = -V(i,j)*(sini*stheta*cphi + cosi*ctheta)
             iu = int(floor((u-umin)/du))
             x = R(i)*(cosi*stheta*cphi + sini*ctheta)
             y = R(i)*stheta*sphi
             ix = int(floor((x-xmin)/dx))
             iy = int(floor((y-ymin)/dy))

             is_outside_V = iu < 1 .or. iu > NU
             is_outside_X = ix < 1 .or. ix > NX
             is_outside_Y = iy < 1 .or. iy > NY
             ! Skip this volume element if we are outside the box
             if (is_outside_V .or. is_outside_X .or. is_outside_Y) cycle radiusloop

             ! Otherwise, add the emissions into the cubes
             cubes(:, ix, iy, iu) = cubes(:, ix, iy, iu) + dvol*Ems(:, i, j)
          end do radiusloop
       end do thetaloop
    end do philoop
    
  end subroutine calculate_cubes_internal

end subroutine calculate_cubes

