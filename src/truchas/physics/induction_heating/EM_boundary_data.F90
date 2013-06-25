!!
!!  The EM_BOUNDARY_DATA Module
!!
!!  Neil N. Carlson <nnc@newmexico.com>
!!  Last revised 23 Aug 2003
!!

#include "f90_assert.fpp"

module EM_boundary_data

  use kind_module, only: rk => real_kind
  use MaxwellBoundaryData
  implicit none
  private
  
  public :: cylinder_bv_init, set_bv, bndry_src

  type(BoundaryData), public, save :: bdata
  
  !! Parameters for the COIL_H_FIELD function
  !real(kind=rk), private, save :: center(3), radius, length
  !integer,       private, save :: turns
  
  !! Parameters for the CONST_H_FIELD function
  !real(kind=rk), private, save :: const_H(3)
  
  real(kind=rk), parameter :: PI =    3.1415926535897932385_rk
  real(kind=rk), parameter :: TWOPI = 6.2831853071795864769_rk
  
contains

  subroutine cylinder_bv_init (mesh)
  
    use parallel_communication
    use distributed_mesh
    use EM_data_proxy, only: symmetry_axis, get_EM_domain_type
    use EM_utilities
    use GeometricModeler
    use bitfield_type, only: btest
    
    type(dist_mesh), intent(in), target :: mesh
    
    integer :: k, bface(mesh%nface), bx0, by0, bz0, bx1, by1, bz1, br1, b30, b60, b120, b150
    real(kind=rk) :: xmin(3), xmax(3), rmax, xh(3), yh(3), zh(3), tol, rmin, vertex(3), slope
    integer, allocatable :: group(:)
    type(GeometricModel) :: gm
    
    real(kind=rk), parameter :: ORIGIN(3) = (/ 0.0_rk, 0.0_rk, 0.0_rk /)
    
  !!!
  !!! RECOVER BOUNDARY INFORMATION BASED ON A PRIORI KNOWLEDGE OF THE PROBLEM
    
    !! Create a geometric model of the domain
    do k = 1, 3
      xmin(k) = global_minval(mesh%x(k,:)) ! min corner of the tight bounding box
      xmax(k) = global_maxval(mesh%x(k,:)) ! max corner of the tight bounding box
    end do
    
    select case (symmetry_axis())
    case ('X')
      rmax = sqrt(global_maxval(mesh%x(2,:)**2 + mesh%x(3,:)**2))
      xh = YHAT
      yh = ZHAT
      zh = XHAT
    case ('Y')
      rmax = sqrt(global_maxval(mesh%x(1,:)**2 + mesh%x(3,:)**2))
      xh = ZHAT
      yh = XHAT
      zh = YHAT
    case ('Z')
      rmax = sqrt(global_maxval(mesh%x(1,:)**2 + mesh%x(2,:)**2))
      xh = XHAT
      yh = YHAT
      zh = ZHAT
    case default
      INSIST( .false. )
    end select
    
    !! Extract some geometrical info used in the frustum case.
    select case (get_EM_domain_type())
    case ('FRUSTUM')
       select case (symmetry_axis())
       case ('X')
          tol = 1.0d-5 * (xmax(1)-xmin(1))
          !! Radius of frustum bottom
          rmin = maxval(mesh%x(2,:)**2 + mesh%x(3,:)**2, mask=abs(mesh%x(1,:)-xmin(1)) < tol)
          rmin = sqrt(global_maxval(rmin))
          !! Radius of frustum top
          rmax = maxval(mesh%x(2,:)**2 + mesh%x(3,:)**2, mask=abs(mesh%x(1,:)-xmax(1)) < tol)
          rmax = sqrt(global_maxval(rmax))
          slope = (rmax - rmin) / (xmax(1) - xmin(1))
          vertex = (/ xmax(1) - rmax / slope, 0.0_rk, 0.0_rk /)
          slope = abs(slope)
       case ('Y')
          tol = 1.0d-5 * (xmax(2)-xmin(2))
          !! Radius of frustum bottom
          rmin = maxval(mesh%x(1,:)**2 + mesh%x(3,:)**2, mask=abs(mesh%x(2,:)-xmin(2)) < tol)
          rmin = sqrt(global_maxval(rmin))
          !! Radius of frustum top
          rmax = maxval(mesh%x(1,:)**2 + mesh%x(3,:)**2, mask=abs(mesh%x(2,:)-xmax(2)) < tol)
          rmax = sqrt(global_maxval(rmax))
          slope = (rmax - rmin) / (xmax(2) - xmin(2))
          vertex = (/ 0.0_rk, xmax(2) - rmax / slope, 0.0_rk /)
          slope = abs(slope)
       case ('Z')
          tol = 1.0d-5 * (xmax(3)-xmin(3))
          !! Radius of frustum bottom
          rmin = maxval(mesh%x(1,:)**2 + mesh%x(2,:)**2, mask=abs(mesh%x(3,:)-xmin(3)) < tol)
          rmin = sqrt(global_maxval(rmin))
          !! Radius of frustum top
          rmax = maxval(mesh%x(1,:)**2 + mesh%x(2,:)**2, mask=abs(mesh%x(3,:)-xmax(3)) < tol)
          rmax = sqrt(global_maxval(rmax))
          slope = (rmax - rmin) / (xmax(3) - xmin(3))
          vertex = (/ 0.0_rk, 0.0_rk, xmax(3) - rmax / slope /)
          slope = abs(slope)
       case default
          INSIST( .false. )
       end select
    case default
    end select

    call GMTune (gm, tol=1.0e-3_rk*minval(xmax-xmin))
    
    select case (get_EM_domain_type())
      case ('FULL_CYLINDER')
        !! Define the bounding surfaces.
        call AddCylinder (gm, br1, ORIGIN, zh, rmax)  ! z-axial cylindrical surface
        call AddPlane (gm, bz0, xmin, -zh)  ! bottom surface
        call AddPlane (gm, bz1, xmax, zh)   ! top surface

        !! Define surface groups to be treated alike.
        allocate(group(3))
        group(1) = 0                      ! dummy place holder
        group(2) = bit_mask((/br1/))      ! nxH given (source)
        group(3) = bit_mask((/bz0, bz1/)) ! nxH given (source) or nxH=0 (symmetry)

      case ('HALF_CYLINDER')
        !! Define the bounding surfaces.
        call AddCylinder (gm, br1, ORIGIN, zh, rmax)  ! z-axial cylindrical surface
        call AddPlane (gm, bz0, xmin, -zh)    ! bottom surface
        call AddPlane (gm, bz1, xmax, zh)     ! top surface
        call AddPlane (gm, bx0, ORIGIN, -xh)  ! x=0 symmetry plane

        !! Define surface groups to be treated alike.
        allocate(group(3))
        group(1) = bit_mask((/bx0/))      ! nxE = 0 (symmetry)
        group(2) = bit_mask((/br1/))      ! nxH given (source)
        group(3) = bit_mask((/bz0, bz1/)) ! nxH given (source) or nxH=0 (symmetry)

      case ('QUARTER_CYLINDER')
        !! Define the bounding surfaces.
        call AddCylinder (gm, br1, ORIGIN, zh, rmax)  ! z-axial cylindrical surface
        call AddPlane (gm, bz0, xmin, -zh)    ! bottom surface
        call AddPlane (gm, bz1, xmax, zh)     ! top surface
        call AddPlane (gm, bx0, ORIGIN, -xh)  ! x=0 symmetry plane
        call AddPlane (gm, by0, ORIGIN, -yh)  ! y=0 symmetry plane

        !! Define surface groups to be treated alike.
        allocate(group(3))
        group(1) = bit_mask((/bx0, by0/)) ! nxE = 0 (symmetry)
        group(2) = bit_mask((/br1/))      ! nxH given (source)
        group(3) = bit_mask((/bz0, bz1/)) ! nxH given (source) or nxH=0 (symmetry)

      case ('CYLINDER')
        !! Define the bounding surfaces.
        call AddCylinder (gm, br1, ORIGIN, zh, rmax)  ! z-axial cylindrical surface
        call AddPlane (gm, bz0, xmin, -zh)    ! bottom surface
        call AddPlane (gm, bz1, xmax, zh)     ! top surface
        call AddPlane (gm, bx0, ORIGIN, -xh)  ! x=0 symmetry plane
        call AddPlane (gm, by0, ORIGIN, -yh)  ! y=0 symmetry plane
        call AddPlane (gm, b30, ORIGIN, xh - sqrt(3.0_rk)*yh)
        call AddPlane (gm, b60, ORIGIN, sqrt(3.0_rk)*xh - yh)
        call AddPlane (gm, b120, ORIGIN, sqrt(3.0_rk)*xh + yh)
        call AddPlane (gm, b150, ORIGIN, xh + sqrt(3.0_rk)*yh)

        !! Define surface groups to be treated alike.
        allocate(group(3))
        group(1) = bit_mask((/bx0, by0, b60, b120, b150/)) ! nxE = 0 (symmetry)
        group(2) = bit_mask((/br1/))      ! nxH given (source)
        group(3) = bit_mask((/bz0, bz1/)) ! nxH given (source) or nxH=0 (symmetry)

      case ('FRUSTUM')
        !! Define the bounding surfaces.
        call AddCone  (gm, br1, vertex, zh, slope)
        call AddPlane (gm, bz0, xmin, -zh)    ! bottom surface
        call AddPlane (gm, bz1, xmax, zh)     ! top surface
        call AddPlane (gm, bx0, ORIGIN, -xh)  ! x=0 symmetry plane
        call AddPlane (gm, by0, ORIGIN, -yh)  ! y=0 symmetry plane
        call AddPlane (gm, b30, ORIGIN, xh - sqrt(3.0_rk)*yh)   ! other symmetry planes
        call AddPlane (gm, b60, ORIGIN, sqrt(3.0_rk)*xh - yh)
        call AddPlane (gm, b120, ORIGIN, sqrt(3.0_rk)*xh + yh)
        call AddPlane (gm, b150, ORIGIN, xh + sqrt(3.0_rk)*yh)
        
        !! Define surface groups to be treated alike.
        allocate(group(3))
        group(1) = bit_mask((/bx0, by0, b30, b60, b120, b150/)) ! nxE = 0 (symmetry)
        group(2) = bit_mask((/br1/))      ! nxH given (source)
        group(3) = bit_mask((/bz0, bz1/)) ! nxH given (source) or nxH=0 (symmetry)
      
      case ('VERIFICATION1')
        call AddPlane (gm, bx0, xmin, -XHAT)
        call AddPlane (gm, by0, xmin, -YHAT)
        call AddPlane (gm, bz0, xmin, -ZHAT)
        call AddPlane (gm, bx1, xmax, XHAT)
        call AddPlane (gm, by1, xmax, YHAT)
        call AddPlane (gm, bz1, xmax, ZHAT)
        
        !! Define surface groups to be treated alike.
        allocate(group(3))
        group(1) = bit_mask((/by0,by1,bx1/))  ! nxE = 0
        group(2) = bit_mask((/bx0/))          ! nxH given
        group(3) = bit_mask((/bz0,bz1/))      ! nxB = 0 (natural)
        
      case default
        INSIST( .false. )
    end select
    
    call generate_bface (gm, mesh, btest(mesh%face_set_mask,0), group, bface)
    call DestroyGeometricModel (gm)
    deallocate(group)
    
    !! NB: As written, the preceding code is not strict.  A quarter cylinder mesh,
    !! declared QUARTER_CYLINDER, in any of the four quadrants will be properly
    !! handled except that both symmetry plane normals point out of the domain
    !! only in the positive quadrant case (the only officially allowed case).
    !! Similarly, a half cylinder mesh, declared HALF_CYLINDER, in any of the four
    !! coordinate half-spaces will be properly handled with the exception of a
    !! possible inward pointing symmetry plane normal.  Worse, is that a full
    !! cylinder mesh is properly handled if it is declared as either HALF_CYLINDER
    !! or QUARTER_CYLINDER, and a half cylinder mesh properly handled if it is
    !! declared QUARTER_CYLINDER.  In fact, the code could be simplified to
    !! the quarter cylinder case alone (we only examine true boundary faces), and
    !! eliminate the user-declared domain type.  However, it seemed useful to
    !! force the user to declare the domain type in order to highlight to the user
    !! the restrictions on the domain.  And as this is temporary code, it doesn't
    !! seem too useful to go to extra lengths to check that the mesh is strictly
    !! of the type declared and officially documented.  The direction of the normal
    !! on the symmetry planes is irrelevant.
    
    !! The sole result of the preceding calculations is the BFACE array!
    
  !!!
  !!! INITIALIZE THE BOUNDARY DATA STRUCTURE
    
    call create_boundary_data (bdata, mesh, bface, ebgroup=(/1/), hbgroup=(/2,3/))
    call set_Eb_function (bdata, 1, zero_field) ! unnec; defaults to zero.
    call set_Hb_function (bdata, 1, source_H_field)
    call set_Hb_function (bdata, 2, source_H_field)

  end subroutine cylinder_bv_init
  
  !!
  !! These H source fields are used to initialize the boundary data structure.
  !! They use previously set private module variables.
  !!
  
  !function coil_H_field (x) result (H)
  !  use solenoid_fields, only: H_coil
  !  real(kind=rk), intent(in) :: x(:)
  !  real(kind=rk) :: H(3)
  !  H = length * H_coil(x-center, radius, 0.5_rk * length, turns)
  !end function coil_H_field
  
  !function const_H_field (x) result (H)
  !  real(kind=rk), intent(in) :: x(:)
  !  real(kind=rk) :: H(3)
  !  H = const_H
  !end function const_H_field
  
  function source_H_field (x) result (H)
  
    use solenoid_fields, only: H_coil
    use EM_data_proxy, only: uniform_source, solenoid, induction_coils, symmetry_axis
    
    real(kind=rk), intent(in) :: x(:)
    real(kind=rk) :: H(3)
    
    integer :: n
    real(kind=rk) :: y(3)
    character(len=1) :: axis
    type(solenoid), pointer :: coil(:)
    
    y = x
    axis = symmetry_axis()
    select case (axis)
    case ('X')
      y = cshift(y,shift=1)
    case ('Y')
      y = cshift(y,shift=-1)
    end select
    
    H = (/ 0.0_rk, 0.0_rk, uniform_source() /)
    coil => induction_coils()
    do n = 1, size(coil)
      H = H + coil(n)%current * H_coil(y-coil(n)%center, coil(n)%radius, coil(n)%length/2, coil(n)%nturns)
    end do
    
    select case (axis)
    case ('X')
      H = cshift(H,shift=-1)
    case ('Y')
      H = cshift(H,shift=1)
    end select
    
  end function source_H_field
  
  !!
  !! This routine gets passed to the integrator to set the E-field boundary
  !! values when needed.  Not much to do for this problem; the only condition
  !! is nxE = 0.
  !!
  
  subroutine set_bv (t, e)
    real(kind=rk), intent(in) :: t
    real(kind=rk), intent(inout) :: e(:)
    call set_Eb_values (bdata, coef=(/0.0_rk/), e=e)
  end subroutine set_bv
  
  !!
  !! This routine gets passed to the integrator to evaluate the source vector
  !! produced by the integration-by-parts boundary integral in the case of
  !! non-homogeneous natural conditions nxH = f
  !!
  
  subroutine bndry_src (t, s)
    real(kind=rk), intent(in) :: t
    real(kind=rk), intent(out) :: s(:)
    real(kind=rk) :: a
    select case (0)
    case (1)  ! Truncated l2 fit to a square wave
      a = (sin(TWOPI*t)+sin(3*TWOPI*t)/3.0+sin(5*TWOPI*t)/5.0+sin(7*TWOPI*t)/7.0)*(4.0/PI)
    case (2)  ! Non-oscillatory 'square' wave
      a = (1225*sin(TWOPI*t)+245*sin(3*TWOPI*t)+49*sin(5*TWOPI*t)+5*sin(7*TWOPI*t))/1024.0
    case default ! Basic sinusoidal wave form.
      a = sin(TWOPI*t)
    end select
    a = (1.0_rk - exp(-2.0_rk*t**2))*a
    call get_Hb_source (bdata, coef=(/a, a/), bsrc=s)
  end subroutine bndry_src
  
end module EM_boundary_data
