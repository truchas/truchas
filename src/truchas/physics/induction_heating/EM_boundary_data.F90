!!
!!  The EM_BOUNDARY_DATA Module
!!
!!  Neil N. Carlson <nnc@newmexico.com>
!!  Last revised 23 Aug 2003
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module EM_boundary_data

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use MaxwellBoundaryData
  use parallel_communication
  use simpl_mesh_type
  use geometry_model_type
  use bitfield_type
  use parameter_list_type
  implicit none
  private

  public :: cylinder_bv_init, set_bv, bndry_src

  type(BoundaryData), public, save :: bdata

  !! Parameters for the COIL_H_FIELD function
  !real(r8), private, save :: center(3), radius, length
  !integer,       private, save :: turns

  !! Parameters for the CONST_H_FIELD function
  !real(r8), private, save :: const_H(3)

  real(r8), parameter :: PI =    3.1415926535897932385_r8
  real(r8), parameter :: TWOPI = 6.2831853071795864769_r8

contains

  subroutine cylinder_bv_init (mesh, params)

    type(simpl_mesh), intent(inout), target :: mesh
    type(parameter_list), intent(inout) :: params

    integer :: k, bface(mesh%nface), bx0, by0, bz0, bx1, by1, bz1, br1, b30, b60, b120, b150
    real(r8) :: xmin(3), xmax(3), rmax, xh(3), yh(3), zh(3), tol, rmin, vertex(3), slope
    integer, allocatable :: group(:)
    type(geometry_model) :: gm
    character(:), allocatable :: domain_type, symmetry_axis

    real(r8), parameter :: ORIGIN(3) = [ 0.0_r8, 0.0_r8, 0.0_r8 ]

    call params%get('em-domain-type', domain_type)
    call params%get('symmetry-axis', symmetry_axis)

  !!!
  !!! RECOVER BOUNDARY INFORMATION BASED ON A PRIORI KNOWLEDGE OF THE PROBLEM

    !! Create a geometric model of the domain
    do k = 1, 3
      xmin(k) = global_minval(mesh%x(k,:)) ! min corner of the tight bounding box
      xmax(k) = global_maxval(mesh%x(k,:)) ! max corner of the tight bounding box
    end do

    select case (symmetry_axis)
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
    select case (domain_type)
    case ('FRUSTUM')
       select case (symmetry_axis)
       case ('X')
          tol = 1.0d-5 * (xmax(1)-xmin(1))
          !! Radius of frustum bottom
          rmin = maxval(mesh%x(2,:)**2 + mesh%x(3,:)**2, mask=abs(mesh%x(1,:)-xmin(1)) < tol)
          rmin = sqrt(global_maxval(rmin))
          !! Radius of frustum top
          rmax = maxval(mesh%x(2,:)**2 + mesh%x(3,:)**2, mask=abs(mesh%x(1,:)-xmax(1)) < tol)
          rmax = sqrt(global_maxval(rmax))
          slope = (rmax - rmin) / (xmax(1) - xmin(1))
          vertex = [ xmax(1) - rmax / slope, 0.0_r8, 0.0_r8 ]
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
          vertex = [ 0.0_r8, xmax(2) - rmax / slope, 0.0_r8 ]
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
          vertex = [ 0.0_r8, 0.0_r8, xmax(3) - rmax / slope ]
          slope = abs(slope)
       case default
          INSIST( .false. )
       end select
    case default
    end select

    call gm%tune(tol=1.0e-3_r8*minval(xmax-xmin))

    select case (domain_type)
      case ('FULL_CYLINDER')
        !! Define the bounding surfaces.
        call gm%add_cylinder(ORIGIN, zh, rmax, br1)  ! z-axial cylindrical surface
        call gm%add_plane(xmin, -zh, bz0)  ! bottom surface
        call gm%add_plane(xmax,  zh, bz1)  ! top surface

        !! Define surface groups to be treated alike.
        allocate(group(3))
        group(1) = 0                      ! dummy place holder
        group(2) = bit_mask([br1])      ! nxH given (source)
        group(3) = bit_mask([bz0, bz1]) ! nxH given (source) or nxH=0 (symmetry)

      case ('HALF_CYLINDER')
        !! Define the bounding surfaces.
        call gm%add_cylinder(ORIGIN, zh, rmax, br1)  ! z-axial cylindrical surface
        call gm%add_plane(xmin, -zh, bz0)    ! bottom surface
        call gm%add_plane(xmax,  zh, bz1)    ! top surface
        call gm%add_plane(ORIGIN, -xh, bx0)  ! x=0 symmetry plane

        !! Define surface groups to be treated alike.
        allocate(group(3))
        group(1) = bit_mask([bx0])      ! nxE = 0 (symmetry)
        group(2) = bit_mask([br1])      ! nxH given (source)
        group(3) = bit_mask([bz0, bz1]) ! nxH given (source) or nxH=0 (symmetry)

      case ('QUARTER_CYLINDER')
        !! Define the bounding surfaces.
        call gm%add_cylinder(ORIGIN, zh, rmax, br1)  ! z-axial cylindrical surface
        call gm%add_plane(xmin, -zh, bz0)    ! bottom surface
        call gm%add_plane(xmax,  zh, bz1)    ! top surface
        call gm%add_plane(ORIGIN, -xh, bx0)  ! x=0 symmetry plane
        call gm%add_plane(ORIGIN, -yh, by0)  ! y=0 symmetry plane

        !! Define surface groups to be treated alike.
        allocate(group(3))
        group(1) = bit_mask([bx0, by0]) ! nxE = 0 (symmetry)
        group(2) = bit_mask([br1])      ! nxH given (source)
        group(3) = bit_mask([bz0, bz1]) ! nxH given (source) or nxH=0 (symmetry)

      case ('CYLINDER')
        !! Define the bounding surfaces.
        call gm%add_cylinder(ORIGIN, zh, rmax, br1)  ! z-axial cylindrical surface
        call gm%add_plane(xmin, -zh, bz0)    ! bottom surface
        call gm%add_plane(xmax, zh, bz1)     ! top surface
        call gm%add_plane(ORIGIN, -xh, bx0)  ! x=0 symmetry plane
        call gm%add_plane(ORIGIN, -yh, by0)  ! y=0 symmetry plane
        call gm%add_plane(ORIGIN, xh - sqrt(3.0_r8)*yh, b30)
        call gm%add_plane(ORIGIN, sqrt(3.0_r8)*xh - yh, b60)
        call gm%add_plane(ORIGIN, sqrt(3.0_r8)*xh + yh, b120)
        call gm%add_plane(ORIGIN, xh + sqrt(3.0_r8)*yh, b150)

        !! Define surface groups to be treated alike.
        allocate(group(3))
        group(1) = bit_mask([bx0, by0, b60, b120, b150]) ! nxE = 0 (symmetry)
        group(2) = bit_mask([br1])      ! nxH given (source)
        group(3) = bit_mask([bz0, bz1]) ! nxH given (source) or nxH=0 (symmetry)

      case ('FRUSTUM')
        !! Define the bounding surfaces.
        call gm%add_cone(vertex, zh, slope, br1)
        call gm%add_plane(xmin, -zh, bz0)    ! bottom surface
        call gm%add_plane(xmax, zh, bz1)     ! top surface
        call gm%add_plane(ORIGIN, -xh, bx0)  ! x=0 symmetry plane
        call gm%add_plane(ORIGIN, -yh, by0)  ! y=0 symmetry plane
        call gm%add_plane(ORIGIN, xh - sqrt(3.0_r8)*yh, b30)   ! other symmetry planes
        call gm%add_plane(ORIGIN, sqrt(3.0_r8)*xh - yh, b60)
        call gm%add_plane(ORIGIN, sqrt(3.0_r8)*xh + yh, b120)
        call gm%add_plane(ORIGIN, xh + sqrt(3.0_r8)*yh, b150)

        !! Define surface groups to be treated alike.
        allocate(group(3))
        group(1) = bit_mask([bx0, by0, b30, b60, b120, b150]) ! nxE = 0 (symmetry)
        group(2) = bit_mask([br1])      ! nxH given (source)
        group(3) = bit_mask([bz0, bz1]) ! nxH given (source) or nxH=0 (symmetry)

      case ('VERIFICATION1')
        call gm%add_plane(xmin, -XHAT, bx0)
        call gm%add_plane(xmin, -YHAT, by0)
        call gm%add_plane(xmin, -ZHAT, bz0)
        call gm%add_plane(xmax,  XHAT, bx1)
        call gm%add_plane(xmax,  YHAT, by1)
        call gm%add_plane(xmax,  ZHAT, bz1)

        !! Define surface groups to be treated alike.
        allocate(group(3))
        group(1) = bit_mask([by0,by1,bx1])  ! nxE = 0
        group(2) = bit_mask([bx0])          ! nxH given
        group(3) = bit_mask([bz0,bz1])      ! nxB = 0 (natural)

      case default
        INSIST( .false. )
    end select

    call add_face_sets(mesh, gm, group) ! overwrites GROUP
    call generate_bface(mesh, group, bface)

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

    call bdata%init(mesh, bface, ebgroup=[1], hbgroup=[2,3])
    call bdata%set_Eb_function(1, zero_field) ! unnec; defaults to zero.
    call bdata%set_Hb_function(1, source_H_field)
    call bdata%set_Hb_function(2, source_H_field)

  end subroutine cylinder_bv_init

  !! This adds face sets to MESH, one for each surface group in the GROUP
  !! array. The face set ID is the index of the corresponding array element.
  !! This is analogous to the generate_bface subroutine. The elements of the
  !! the GROUP array are overwritten with the corresponding index in the
  !! MESH%FACE_SET_ID array; the index is the bit in the %FACE_SET_MASK array.

  subroutine add_face_sets(mesh, gm, group)

    type(simpl_mesh), intent(inout) :: mesh
    type(geometry_model), intent(in) :: gm
    integer, intent(inout) :: group(:)

    integer :: i, j, fsid(size(group))
    integer, allocatable :: surf(:)

    fsid = [(i + size(mesh%face_set_id), i=1,size(group))]
    mesh%face_set_id = [mesh%face_set_id, fsid]

    do j = 1, mesh%nface
      if (.not.btest(mesh%face_set_mask(j), pos=0)) cycle ! not a boundary face
      associate (face => mesh%x(:,mesh%fnode(:,j)), mask => mesh%face_set_mask(j))
        call gm%get_on_surface_list(face, surf)
        select case (size(surf))
        case (1)
          ASSERT(surf(1) > 0 .and. surf(1) < bit_size(surf(1)))
          do i = 1, size(group)
            if (btest(group(i), surf(1))) mask = ibset(mask, pos=fsid(i))
          end do
        case (2:)
          INSIST(.false.) !TODO: better error handling
        end select
      end associate
    end do

    group = fsid

    block ! check for a complete set of boundary face sets
      use parallel_communication, only: global_any
      use truchas_env, only: output_dir
      use truchas_logging_services
      character(:), allocatable :: vtk_file
      logical :: mask(mesh%nface)
      do j = 1, mesh%nface
        mask(j) = btest(mesh%face_set_mask(j), pos=0) .and. .not.any(btest(mesh%face_set_mask(j), pos=fsid))
      end do
      if (global_any(mask)) then
        vtk_file = trim(output_dir) // 'em_no_bc.vtk'
        call mesh%write_faces_vtk(mask, vtk_file, 'Boundary faces without BC')
        call TLS_fatal('found unexpected EM mesh boundary faces; possible causes:' &
            // new_line('a') &
            // '* invalid EM domain or specification; check EM_domain_type and symmetry_axis' &
            // new_line('a') &
            // '* volumes were not imprinted and merged in Cubit before meshing' &
            // new_line('a') &
            // 'Faces written to the Paraview input file "' // vtk_file // '" for visualization')
      end if
    end block

  end subroutine add_face_sets

  subroutine generate_bface (mesh, group, bface, status)

    use truchas_env, only: output_dir
    use parallel_communication, only: global_any
    use truchas_logging_services

    type(simpl_mesh), intent(in) :: mesh
    integer, intent(in) :: group(:)
    integer, intent(out) :: bface(:)
    integer, intent(out), optional :: status

    integer :: i, j
    logical :: mask(mesh%nface)
    character(:), allocatable :: vtk_file

    ASSERT(size(bface) == mesh%nface)

    bface = 0
    do j = 1, mesh%nface
      if (.not.btest(mesh%face_set_mask(j), pos=0)) cycle ! ignore this face
      do i = 1, size(group)
        if (btest(mesh%face_set_mask(j), pos=group(i))) then
          if (mesh%fcell(2,j) == 0) then ! face normal points out of domain
            bface(j) = i
          else  ! face normal points into the domain
            bface(j) = -i
          end if
          exit  ! NB: the face sets in GROUP are disjoint
        end if
      end do
    end do

    mask = ((bface == 0) .and. btest(mesh%face_set_mask, pos=0))
    if (global_any(mask)) then
      if (present(status)) then
        status = -1
        return
      else
        vtk_file = trim(output_dir) // 'em_no_bc.vtk'
        call mesh%write_faces_vtk(mask, vtk_file, 'Boundary faces without BC')
        call TLS_fatal('found unexpected EM mesh boundary faces; possible causes:' &
            // new_line('a') &
            // '* invalid EM domain or specification; check EM_domain_type and symmetry_axis' &
            // new_line('a') &
            // '* volumes were not imprinted and merged in Cubit before meshing' &
            // new_line('a') &
            // 'Faces written to the Paraview input file "' // vtk_file // '" for visualization')
      end if
    end if

    if (present(status)) status = 0

  end subroutine generate_bface

  !!
  !! These H source fields are used to initialize the boundary data structure.
  !! They use previously set private module variables.
  !!

  !function coil_H_field (x) result (H)
  !  use solenoid_fields, only: H_coil
  !  real(r8), intent(in) :: x(:)
  !  real(r8) :: H(3)
  !  H = length * H_coil(x-center, radius, 0.5_r8 * length, turns)
  !end function coil_H_field

  !function const_H_field (x) result (H)
  !  real(r8), intent(in) :: x(:)
  !  real(r8) :: H(3)
  !  H = const_H
  !end function const_H_field

  function source_H_field (x) result (H)

    use solenoid_fields, only: H_coil
    use EM_data_proxy, only: uniform_source, solenoid, induction_coils, symmetry_axis

    real(r8), intent(in) :: x(:)
    real(r8) :: H(3)

    integer :: n
    real(r8) :: y(3)
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

    H = [ 0.0_r8, 0.0_r8, uniform_source() ]
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
    real(r8), intent(in) :: t
    real(r8), intent(inout) :: e(:)
    call bdata%set_Eb_values(coef=[0.0_r8], e=e)
  end subroutine set_bv

  !!
  !! This routine gets passed to the integrator to evaluate the source vector
  !! produced by the integration-by-parts boundary integral in the case of
  !! non-homogeneous natural conditions nxH = f
  !!

  subroutine bndry_src (t, s)
    real(r8), intent(in) :: t
    real(r8), intent(out) :: s(:)
    real(r8) :: a
    select case (0)
    case (1)  ! Truncated l2 fit to a square wave
      a = (sin(TWOPI*t)+sin(3*TWOPI*t)/3.0+sin(5*TWOPI*t)/5.0+sin(7*TWOPI*t)/7.0)*(4.0/PI)
    case (2)  ! Non-oscillatory 'square' wave
      a = (1225*sin(TWOPI*t)+245*sin(3*TWOPI*t)+49*sin(5*TWOPI*t)+5*sin(7*TWOPI*t))/1024.0
    case default ! Basic sinusoidal wave form.
      a = sin(TWOPI*t)
    end select
    a = (1.0_r8 - exp(-2.0_r8*t**2))*a
    call bdata%get_Hb_source(coef=[a, a], bsrc=s)
  end subroutine bndry_src

end module EM_boundary_data
