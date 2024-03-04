!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module legacy_ih_bc

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use simpl_mesh_type
  use geometry_model_type
  use bitfield_type
  use parameter_list_type
  implicit none
  private

  public :: create_ih_face_sets

  real(r8), parameter :: PI =    3.1415926535897932385_r8
  real(r8), parameter :: TWOPI = 6.2831853071795864769_r8

contains

  subroutine create_ih_face_sets(mesh, params, pec_setid, src_setid)

    type(simpl_mesh), intent(inout), target :: mesh
    type(parameter_list), intent(inout) :: params
    integer, allocatable, intent(out) :: pec_setid(:), src_setid(:)

    integer :: k, bx0, by0, bz0, bx1, by1, bz1, br1, b30, b60, b120, b150
    real(r8) :: xmin(3), xmax(3), rmax, xh(3), yh(3), zh(3), tol, rmin, vertex(3), slope
    integer, allocatable :: group(:), setids(:)
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

    ! GROUP is a list of bit masks, each identifying a set of surfaces
    allocate(setids, mold=group)
    call add_face_sets(mesh, gm, group, setids) ! overwrites GROUP
    ! GROUP is a list of bit positions (in %FACE_SET_MASK), one for eash of the surface sets

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

    !! Side sets for the PEC tangential E and tangential H boundary conditions
    pec_setid = setids(1:1)
    src_setid = setids(2:3)

  end subroutine create_ih_face_sets

  !! This adds face sets to MESH, one for each surface group in the GROUP
  !! array. The face set ID is the index of the corresponding array element.
  !! This is analogous to the generate_bface subroutine. The elements of the
  !! the GROUP array are overwritten with the corresponding index in the
  !! MESH%FACE_SET_ID array; the index is the bit in the %FACE_SET_MASK array.

  subroutine add_face_sets(mesh, gm, group, setids)

    type(simpl_mesh), intent(inout) :: mesh
    type(geometry_model), intent(in) :: gm
    integer, intent(inout) :: group(:)
    integer, intent(out) :: setids(:)

    integer :: i, j, bits(size(group))
    integer, allocatable :: surf(:)

    ASSERT(size(setids) == size(group))

    bits = [(i + size(mesh%face_set_id), i=1,size(group))]
    setids = [(i + max(0, maxval(mesh%face_set_id)), i=1,size(group))]
    mesh%face_set_id = [mesh%face_set_id, setids]

    do j = 1, mesh%nface
      if (.not.btest(mesh%face_set_mask(j), pos=0)) cycle ! not a boundary face
      associate (face => mesh%x(:,mesh%fnode(:,j)), mask => mesh%face_set_mask(j))
        call gm%get_on_surface_list(face, surf)
        select case (size(surf))
        case (1)
          ASSERT(surf(1) > 0 .and. surf(1) < bit_size(surf(1)))
          do i = 1, size(group)
            if (btest(group(i), surf(1))) mask = ibset(mask, pos=bits(i))
          end do
        case (2:)
          INSIST(.false.) !TODO: better error handling
        end select
      end associate
    end do

    group = bits

    block ! check for a complete set of boundary face sets
      use parallel_communication, only: global_any
      use truchas_env, only: output_dir
      use truchas_logging_services
      character(:), allocatable :: vtk_file
      logical :: mask(mesh%nface)
      do j = 1, mesh%nface
        mask(j) = btest(mesh%face_set_mask(j), pos=0) .and. .not.any(btest(mesh%face_set_mask(j), pos=bits))
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

  function bit_mask (list) result (n)
    integer, intent(in) :: list(:)
    integer :: n, j
    n = 0
    do j = 1, size(list)
      n = ibset(n, list(j))
    end do
  end function bit_mask

end module legacy_ih_bc
