! Wrapper routine to be called from python
! Compiled with f2py
! Even though it is a subroutine, it appears a function in python
! 
! Note that R must be sent in dimensionless form (units of i-front radius)
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
  ! local vars
  integer :: idebug = 1


  if (idebug >= 1) then
     print '(a)',  "NK, inc_degrees, xmin, ymin, umin, dx, dy, du"
     print *,  NK, inc_degrees, xmin, ymin, umin, dx, dy, du
     print *, "sum(Ems) = ", sum(Ems)
  end if
  
  call calculate_cubes_internal

contains

  subroutine calculate_cubes_internal
    ! Actually do the work
    real, parameter :: PI = 3.1415926535897932385
    ! Local variables
    real :: u, x, y
    real, allocatable, dimension(:) :: Phi
    integer :: i, j, k
    real :: dphi, dmu, dr, dvol
    real, save :: cphi, sphi, ctheta, stheta, cosi, sini
    integer, save :: jneg, jpos, ineg, ipos
    integer, save :: iu, ix, iy
    logical, save :: is_outside_V, is_outside_X, is_outside_Y
    character(len=*), parameter :: interpolation = "linear"
    ! Variables for linear interpolation
    real, save :: uu, xx, yy, au, ax, ay, bu, bx, by
    real, save :: c000, c001, c010, c011, c100, c101, c110, c111
    !$OMP THREADPRIVATE(cphi, sphi, ctheta, stheta, jneg, jpos, ineg, ipos)
    !$OMP THREADPRIVATE(iu, ix, iy, is_outside_V, is_outside_X, is_outside_Y)
    !$OMP THREADPRIVATE(uu, xx, yy, au, ax, ay, bu, bx, by)
    !$OMP THREADPRIVATE(c000, c001, c010, c011, c100, c101, c110, c111)



    allocate( Phi(NK) )
    cubes = 0.0
    dphi = 2.0*PI/real(NK)
    Phi = (/ (real(k)*dphi, k = 1, NK ) /)
    cosi = cos(inc_degrees*PI/180.0)
    sini = sin(inc_degrees*PI/180.0)

    !$OMP PARALLEL DO
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
             u = -V(j,i)*(sini*stheta*cphi + cosi*ctheta)
             x = R(i)*(-cosi*stheta*cphi + sini*ctheta)
             y = R(i)*stheta*sphi

             if (interpolation == "nearest") then
                iu = 1 + int(floor((u-umin)/du))
                ix = 1 + int(floor((x-xmin)/dx))
                iy = 1 + int(floor((y-ymin)/dy))
                is_outside_V = iu < 1 .or. iu > NU
                is_outside_X = ix < 1 .or. ix > NX
                is_outside_Y = iy < 1 .or. iy > NY
                ! Skip this volume element if we are outside the box
                if (is_outside_V .or. is_outside_X .or. is_outside_Y) cycle radiusloop
                ! Otherwise, add the emissions into the cubes
                cubes(:, iu, iy, ix) = cubes(:, iu, iy, ix) + dvol*Ems(:, j, i)
             else if (interpolation == "linear") then
                ! Linear interpolation in 3D
                !
                ! map to real scale from 0.0 -> N
                uu = (u-umin)/du
                xx = (x-xmin)/dx
                yy = (y-ymin)/dy
                ! nint gives the "left" cell 
                iu = nint(uu)
                ix = nint(xx)
                iy = nint(yy)
                ! have to be slightly more picky in this case
                is_outside_V = iu < 1 .or. iu + 1 > NU
                is_outside_X = ix < 1 .or. ix + 1 > NX
                is_outside_Y = iy < 1 .or. iy + 1 > NY
                ! Skip this volume element if we are outside the box
                if (is_outside_V .or. is_outside_X .or. is_outside_Y) cycle radiusloop
                ! find distance in fractional cell widths from center of left cell
                au = uu - real(iu) + 0.5
                ax = xx - real(ix) + 0.5
                ay = yy - real(iy) + 0.5
                ! same for distance from center of right cell
                bu = 1.0 - au
                bx = 1.0 - ax
                by = 1.0 - ay
                ! the eight interpolation coefficients
                c000 = bu*bx*by
                c001 = bu*bx*ay
                c010 = bu*ax*by
                c011 = bu*ax*ay
                c100 = au*bx*by
                c101 = au*bx*ay
                c110 = au*ax*by
                c111 = au*ax*ay
                ! add into the eight cells
                cubes(:, iu  , iy  , ix  ) = cubes(:, iu  , iy  , ix  ) + c000*dvol*Ems(:, j, i)
                cubes(:, iu  , iy  , ix+1) = cubes(:, iu  , iy  , ix+1) + c001*dvol*Ems(:, j, i)
                cubes(:, iu  , iy+1, ix  ) = cubes(:, iu  , iy+1, ix  ) + c010*dvol*Ems(:, j, i)
                cubes(:, iu  , iy+1, ix+1) = cubes(:, iu  , iy+1, ix+1) + c011*dvol*Ems(:, j, i)
                cubes(:, iu+1, iy  , ix  ) = cubes(:, iu+1, iy  , ix  ) + c100*dvol*Ems(:, j, i)
                cubes(:, iu+1, iy  , ix+1) = cubes(:, iu+1, iy  , ix+1) + c101*dvol*Ems(:, j, i)
                cubes(:, iu+1, iy+1, ix  ) = cubes(:, iu+1, iy+1, ix  ) + c110*dvol*Ems(:, j, i)
                cubes(:, iu+1, iy+1, ix+1) = cubes(:, iu+1, iy+1, ix+1) + c111*dvol*Ems(:, j, i)
             else
                print '(2a)', "Unknown interpolation method: ", interpolation
                stop
             end if

             if (dvol*Ems(1, j, i) >= 1.e-12 .and. idebug >= 2) then
                print *, "Large value warning:"
                print *, "iu, iy, ix, cubes(1, iu, iy, ix )"
                print *, iu, iy, ix, cubes(1, iu, iy, ix )
                print *, "j, i, dvol, Ems(1, j, i)"
                print *, j, i, dvol, Ems(1, j, i)
                print *, "interpolation coefficients:"
                print *, au, bu, ay, by, ax, bx
             end if

             

          end do radiusloop
       end do thetaloop
    end do philoop
    !$OMP END PARALLEL DO

    print *, "sum(cubes) = ", sum(cubes)
  end subroutine calculate_cubes_internal

end subroutine calculate_cubes


subroutine weighted_temperature(R, Mu, T, Ems, validmask, NI, NJ, NL, Tmaps, NY, NX,&
     & NK, inc_degrees, xmin, ymin, dx, dy)
  implicit none
!!! Input Arguments
  integer, intent(in) :: NI     ! Size of input radius array
  integer, intent(in) :: NJ     ! Size of input angle array
  integer, intent(in) :: NL     ! Number of emission lines
  real, intent(in), dimension(NI) :: R ! Input radius array
  real, intent(in), dimension(NJ) :: Mu ! Input angle (cosine) array
  real, intent(in), dimension(NJ, NI) :: T ! Input temperature field
  real, intent(in), dimension(NL, NJ, NI) :: Ems ! Input line emissivities
  ! which combinations of R and mu have valid data
  logical, intent(in), dimension(NJ, NI) :: validmask 
  integer, intent(in) :: NX, NY ! Dimensions of output cubes
!!! Output Argument
  real, intent(out), dimension(NL, NY, NX) :: Tmaps
!!! Auxiliary input arguments
  integer :: NK                 ! Number of phi values to use
  real :: inc_degrees           ! Inclination of proplyd axis to line of sight
  ! Limits of PPV cube
  real :: xmin, ymin
  ! Note that [xyu]max are not used
  real :: dx, dy

  call weighted_temperature_internal

contains

  subroutine weighted_temperature_internal
    ! Actually do the work
    real, parameter :: PI = 3.1415926535897932385
    ! Local variables
    real :: x, y
    real, allocatable, dimension(:) :: Phi
    integer :: i, j, k
    real :: dphi, dmu, dr, dvol
    real :: cphi, sphi, ctheta, stheta, cosi, sini
    integer :: jneg, jpos, ineg, ipos
    integer :: ix, iy
    logical :: is_outside_X, is_outside_Y
    character(len=*), parameter :: interpolation = "linear"
    ! Variables for linear interpolation
    real :: xx, yy, ax, ay, bx, by
    real :: c00, c01, c10, c11
    ! Surface brightness maps 
    real, allocatable, dimension(:,:,:) :: Smaps
    real, dimension(NL) :: product
    
    ! Smaps is the same size as Tmaps and is used to hold the surface brightness so that 
    ! we can normalize by it at the end
    allocate( Phi(NK), Smaps(NL, NY, NX) )
    Tmaps = 0.0
    Smaps = 0.0
    dphi = 2.0*PI/real(NK)
    Phi = (/ (real(k)*dphi, k = 1, NK ) /)
    cosi = cos(inc_degrees*PI/180.0)
    sini = sin(inc_degrees*PI/180.0)

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
             x = R(i)*(-cosi*stheta*cphi + sini*ctheta)
             y = R(i)*stheta*sphi

             if (interpolation == "nearest") then
                ix = 1 + int(floor((x-xmin)/dx))
                iy = 1 + int(floor((y-ymin)/dy))
                is_outside_X = ix < 1 .or. ix > NX
                is_outside_Y = iy < 1 .or. iy > NY
                ! Skip this volume element if we are outside the box
                if (is_outside_X .or. is_outside_Y) cycle radiusloop
                ! Otherwise, add the emissions into the cubes
                Tmaps(:, iy, ix) = Tmaps(:, iy, ix) + dvol*Ems(:, j, i)*T(j, i)
                Smaps(:, iy, ix) = Smaps(:, iy, ix) + dvol*Ems(:, j, i)
             else if (interpolation == "linear") then
                ! Linear interpolation in 3D
                !
                ! map to real scale from 0.0 -> N
                xx = (x-xmin)/dx
                yy = (y-ymin)/dy
                ! nint gives the "left" cell 
                ix = nint(xx)
                iy = nint(yy)
                ! have to be slightly more picky in this case
                is_outside_X = ix < 1 .or. ix + 1 > NX
                is_outside_Y = iy < 1 .or. iy + 1 > NY
                ! Skip this volume element if we are outside the box
                if (is_outside_X .or. is_outside_Y) cycle radiusloop
                ! find distance in fractional cell widths from center of left cell
                ax = xx - real(ix) + 0.5
                ay = yy - real(iy) + 0.5
                ! same for distance from center of right cell
                bx = 1.0 - ax
                by = 1.0 - ay
                ! the eight interpolation coefficients
                c00 = bx*by
                c01 = bx*ay
                c10 = ax*by
                c11 = ax*ay
                ! add into the eight cells
                product = dvol*Ems(:, j, i)*T(j, i)
                Tmaps(:, iy  , ix  ) = Tmaps(:, iy  , ix  ) + c00*product
                Tmaps(:, iy  , ix+1) = Tmaps(:, iy  , ix+1) + c01*product
                Tmaps(:, iy+1, ix  ) = Tmaps(:, iy+1, ix  ) + c10*product
                Tmaps(:, iy+1, ix+1) = Tmaps(:, iy+1, ix+1) + c11*product
                product = dvol*Ems(:, j, i)
                Smaps(:, iy  , ix  ) = Smaps(:, iy  , ix  ) + c00*product
                Smaps(:, iy  , ix+1) = Smaps(:, iy  , ix+1) + c01*product
                Smaps(:, iy+1, ix  ) = Smaps(:, iy+1, ix  ) + c10*product
                Smaps(:, iy+1, ix+1) = Smaps(:, iy+1, ix+1) + c11*product
             else
                print '(2a)', "Unknown interpolation method: ", interpolation
                stop
             end if

          end do radiusloop
       end do thetaloop
    end do philoop

    Tmaps = Tmaps / Smaps

  end subroutine weighted_temperature_internal

end subroutine weighted_temperature

