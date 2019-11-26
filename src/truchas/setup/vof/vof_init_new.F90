!!
!! VOF_INIT
!!
!! This module provides routines for initializing the VOF scalar based
!! on shapes provided by the user. It uses the divide and conquer method.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! November 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module vof_init_NEW

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use body_identifier_type
  use truchas_logging_services
  implicit none
  private

  public :: vof_initialize_NEW

  ! A divide and conquer cell type to hold the logic for recursively
  ! splitting cells into tets, checking if a cell is mixed, and
  ! computing volume fractions at the lowest level.
  type :: dnc_cell
    private
    real(r8), allocatable :: x(:,:)
    integer, allocatable :: matl_at_node(:)
    integer :: nnode, ntet, recursion_height, cellid
    real(r8) :: volume
    logical :: is_subtet
    integer, pointer :: to_tet(:,:) => null() ! do not own
    type(body_identifier), pointer :: body_id => null() ! do not own
  contains
    procedure :: init => dnc_cell_init
    procedure :: volumes
    procedure, private :: get_subtet
    procedure, private :: prepare_subtet
  end type dnc_cell

  integer, parameter :: tet4_ntet = 8
  integer, parameter :: pyr5_ntet = 4
  integer, parameter :: wed6_ntet = 17
  integer, parameter :: hex8_ntet = 24
  integer, target :: pyr5_to_tet(4,pyr5_ntet), wed6_to_tet(4,wed6_ntet), &
      hex8_to_tet(4,hex8_ntet), tet4_to_tet(4,tet4_ntet)
  data PYR5_TO_TET/1,2,5,6, 2,3,5,6, 3,4,5,6, 1,4,5,6/
  data WED6_TO_TET/1,4,7,9, 2,5,7,8, 3,6,8,9, 7,8,9,10, 7,8,9,11, 1,7,9,10, 2,7,8,10, &
      3,8,9,10, 1,2,7,10, 1,3,9,10, 2,3,8,10, 5,7,8,11, 6,8,9,11, 4,7,9,11, &
      5,6,8,11, 4,5,7,11, 4,6,9,11/
  data HEX8_TO_TET/1,5,9,12, 2,6,9,10, 3,7,10,11, 4,8,11,12, 1,2,9,13, 2,3,10,13, 3,4,11,13, &
      1,4,12,13, 5,6,9,14, 6,7,10,14, 7,8,11,14, 5,8,12,14, 9,10,13,14, &
      10,11,13,14, 11,12,13,14, 9,12,13,14, 5,9,12,14, 6,9,10,14, 7,10,11,14, &
      8,11,12,14, 1,9,12,13, 2,9,10,13, 3,10,11,13, 4,11,12,13/
  data TET4_TO_TET/1,5,6,7, 2,5,8,9, 3,6,8,10, 4,7,9,10, 5,6,7,10, 5,6,8,10, 5,7,9,10, 5,8,9,10/

contains

  ! This subroutine initializes the body volumes in all cells from a
  ! user input function. It calls a divide and conquer algorithm to
  ! calculate the volume fractions given an interface.
  !
  ! TODO: recursion_limit should come from a parameter list.
  subroutine vof_initialize_NEW(mesh, plist, recursion_limit, vof)

    use unstr_mesh_type
    use parameter_list_type
    use timer_tree_type

    type(unstr_mesh), intent(in) :: mesh
    type(parameter_list), intent(in) :: plist
    integer, intent(in) :: recursion_limit
    real(r8), intent(out) :: vof(:,:)

    integer :: i
    type(body_identifier), target :: body_id
    type(dnc_cell) :: cell

    INSIST(recursion_limit > 0)
    ASSERT(size(vof,dim=2) >= mesh%ncell_onP)

    call start_timer("VOF Initialize")

    ! Here we divide by the cell%volume, rather than mesh%volume(i),
    ! to ensure consistency. This routine works by tallying up volumes
    ! from recursively divided tets, and the accumulated volume of
    ! those tets aren't guaranteed to be machine-identical to the
    ! parent volume.
    call body_id%init(plist, mesh) !body_ids)
    ASSERT(size(vof,dim=1) == body_id%nbody)
    do i = 1, mesh%ncell_onP
      call cell%init(i, mesh, body_id, recursion_limit)
      vof(:,i) = cell%volumes() !/ cell%volume
      vof(:,i) = vof(:,i) * mesh%volume(i) / sum(vof(:,i)) ! correct for errors in planar assumption
    end do

    call stop_timer("VOF Initialize")

  end subroutine vof_initialize_NEW


  subroutine dnc_cell_init(this, i, mesh, body_id, recursion_limit)

    use unstr_mesh_type
    use cell_topology

    class(dnc_cell), intent(out) :: this
    class(unstr_mesh), intent(in) :: mesh
    integer, intent(in) :: i, recursion_limit
    type(body_identifier), target, intent(in) :: body_id

    integer :: j, nface

    this%body_id => body_id
    this%recursion_height = recursion_limit
    this%cellid = i
    this%is_subtet = .false.
    this%volume = mesh%volume(i)
    
    associate (cn => mesh%cnode(mesh%xcnode(i):mesh%xcnode(i+1)-1))
      this%nnode = size(cn)
      nface = num_cell_faces(cn)
      
      select case(this%nnode)
      case (4)
        this%to_tet => tet4_to_tet
        allocate(this%x(3,this%nnode+6))
        this%x(:,:this%nnode) = mesh%x(:,cn)
        do j = 1, 6
          this%x(:,4+j) = sum(this%x(:,tet4_edges(:,j)), dim=2) / 2
        end do
      case (5)
        this%to_tet => pyr5_to_tet
        allocate(this%x(3,this%nnode+1))
        this%x(:,:this%nnode) = mesh%x(:,cn)
        ! split the only potentially non-planar face (5)
        this%x(:,this%nnode+1) = sum(mesh%x(:,cn(pyr5_faces(pyr5_xface(5):pyr5_xface(6)-1))),dim=2) / pyr5_fsize(5)
      case (6)
        this%to_tet => wed6_to_tet
        allocate(this%x(3,this%nnode+nface))
        this%x(:,:this%nnode) = mesh%x(:,cn)
        do j = 1, nface
          this%x(:,this%nnode+j) = sum(mesh%x(:,cn(wed6_faces(wed6_xface(j):wed6_xface(j+1)-1))),dim=2) / wed6_fsize(j)
        end do
      case (8)
        this%to_tet => hex8_to_tet
        allocate(this%x(3,this%nnode+nface))
        this%x(:,:this%nnode) = mesh%x(:,cn)
        do j = 1, nface
          this%x(:,this%nnode+j) = sum(mesh%x(:,cn(hex8_faces(hex8_xface(j):hex8_xface(j+1)-1))),dim=2) / hex8_fsize(j)
        end do
      case default
        INSIST(.false.)
      end select

      this%ntet = size(this%to_tet,dim=2)
    end associate

    this%matl_at_node = [(body_id%body_at_point(this%x(:,j), this%cellid), j = 1, size(this%x,dim=2))]
    INSIST(.not.any(this%matl_at_node==0))

  end subroutine dnc_cell_init


  ! Compute the volume fractions of materials in a cell. This will also
  ! update this cell's volume to be the sum of subcell volumes, if
  ! applicable.
  !
  ! Check if all vertices lie within a single material. If so, set the
  ! VOF for that material to 1. If not, divide the cell into tets and
  ! repeat recursively to a given threshold.
  !
  ! TODO: Right now we close by counting how many nodes are inside
  !       the volume, and this is likely what should be done whenever
  !       there are >2 materials in the cell. For 2-material cells,
  !       we could use a plane reconstruction from a signed distance
  !       function to produce a more accurate estimate. This produces
  !       much greater accuracy for very little cost, far less than
  !       adding levels of refinement.
  !
  ! TODO: Interface detection might be improved. By checking only
  !       vertices, we fail to subdivide if an interface intersects a
  !       cell without cutting off any vertices.
  recursive function volumes(this)

    use cell_geometry, only: tet_volume

    class(dnc_cell), intent(inout) :: this
    real(r8) :: volumes(this%body_id%nbody)

    integer :: i
    type(dnc_cell) :: subtet

    volumes = 0

    ! If the cell contains an interface (and therefore has at least
    ! two materials and we haven't yet hit our recursion limit, divide
    ! the cell and repeat.
    if (this%recursion_height > 0 &
        .and. any(this%matl_at_node(:this%nnode) /= this%matl_at_node(1))) then
      ! tally the volumes from refined sub-cells
      call this%prepare_subtet(subtet)
      do i = 1, this%ntet
        call this%get_subtet(subtet, i)
        volumes = volumes + subtet%volumes()
      end do
    else ! if (this%body_id%nbody > 2 .or. .not.this%contains_interface()) then
      ! If we're past the recursion limit or the cell does not contain
      ! an interface, calculate the vof in this cell based on the
      ! materials at its nodes.
      if (this%is_subtet) then
        ! Subtet volume is calculated here for the first time, since
        ! we don't need subtet volumes at every refinement level.
        this%volume = tet_volume(this%x)
        INSIST(this%volume > 0) ! ensure we haven't recursed beyond floating point limits
      end if

      do i = 1, this%body_id%nbody
        volumes(i) = count(this%matl_at_node(:this%nnode) == i)
      end do
      volumes = volumes / this%nnode * this%volume
    ! else
    !   ! For two-material cells, we can reconstruct an interface from
    !   ! 3 interface-points identified by linearly interpolating a
    !   ! signed distance function.
    !   volumes = this%volume_from_level_set()
    end if

  end function volumes


  ! Prepare this structure and an input subtet for subdivision. This
  ! is written so that we can quickly divide into tets, and generate
  ! those tets and work on them one at a time. To this end, this
  ! routine computes edge and face centers so they don't need to be
  ! recalculated for every subdivision, and it initializes the input
  ! subtet with the right array sizes for a tet.
  subroutine prepare_subtet(this, subtet)

    use cell_topology, only: tet4_edges

    class(dnc_cell), intent(inout) :: this
    class(dnc_cell), intent(out) :: subtet

    integer :: j

    ! Compute this cells hanging nodes and find material at that point.
    if (this%is_subtet) then
      do j = 1, 6
        this%x(:,4+j) = sum(this%x(:,tet4_edges(:,j)), dim=2) / 2
        this%matl_at_node(4+j) = this%body_id%body_at_point(this%x(:,4+j), this%cellid)
      end do
      INSIST(.not.any(subtet%matl_at_node==0))
    end if

    ! Subtet structure. The geometry itself is generated by get_subtet.
    subtet%body_id => this%body_id
    subtet%to_tet => tet4_to_tet
    subtet%ntet = tet4_ntet
    subtet%recursion_height = this%recursion_height - 1
    subtet%cellid = this%cellid
    subtet%is_subtet = .true.
    subtet%nnode = 4
    allocate(subtet%x(3,subtet%nnode+6), subtet%matl_at_node(subtet%nnode+6))

  end subroutine prepare_subtet


  ! Return a complete subtet geometry from a given subtet ID.
  subroutine get_subtet(this, subtet, i)
    class(dnc_cell), intent(in) :: this
    class(dnc_cell), intent(inout) :: subtet
    integer, intent(in) :: i
    ASSERT(i > 0 .and. i <= this%ntet)
    subtet%x(:,:4) = this%x(:,this%to_tet(:,i))
    subtet%matl_at_node(:4) = this%matl_at_node(this%to_tet(:,i))
  end subroutine get_subtet


  ! ! For two-material cells, we can reconstruct an interface from
  ! ! three interface-points identified by linearly interpolating a
  ! ! signed distance function. Note that this function ASSUMES it's
  ! ! operating on a cell with only two materials present, and that it
  ! ! is operating on a tet
  ! function volume_from_level_set(this) result(volumes)

  !   use plane_type
  !   use pure_polyhedron_type
  !   use cell_geometry, only: cross_product, normalized

  !   class(dnc_cell), intent(in) :: this
  !   real(r8) :: volumes(this%body_id%nbody)

  !   integer :: m, ierr, i, e, matl(2)
  !   real(r8) :: x(3,3), s1, s2
  !   type(plane) :: P
  !   type(pure_polyhedron) :: poly

  !   volumes = 0
    
  !   ! identify the two materials
  !   matl(1) = this%matl_at_node(1)
  !   do i = 2, this%nnode
  !     if (this%matl_at_node(i) /= matl(1)) then
  !       matl(2) = this%matl_at_node(i)
  !       exit
  !     end if
  !   end do
    
  !   ! get 3 points intersecting the edges of this cell
  !   i = 0
  !   do e = 1, this%nedge
  !     if (this%matl_at_node(this%edges(1,e)) /= this%matl_at_node(this%edges(2,e))) then
  !       s1 = this%body_id%signed_distance(this%x(:,this%edges(1,e)))
  !       s2 = this%body_id%signed_distance(this%x(:,this%edges(2,e)))

  !       i = i + 1
  !       x(:,i) = this%x(:,this%edges(1,e)) + &
  !           s1/(s1-s2) * (this%x(:,this%edges(2,e)) - this%x(:,this%edges(1,e)))
  !       if (i == 3) exit
  !     end if
  !   end do
  !   !ASSERT(e <= this%nedges) ! make sure we found 3 intersection points

  !   if (e <= this%nedges) then
  !     ! Construct a plane from the 3 points.
  !     P%normal = normalized(cross_product(x(:,2) - x(:,1), x(:,3) - x(:,1)))
  !     P%rho = dot_product(x(:,1), P%normal)

  !     ! Orient the interface so the normal points outside relative to matl(1).
  !     do i = 1, this%nnode
  !       s1 = P%signed_distance(this%x(:,i))

  !       if (abs(s1) > 0) then
  !         if ((s1 > 0 .and. this%matl_at_node(i) == matl(1)) .or. &
  !             (s1 < 0 .and. this%matl_at_node(i) /= matl(1))) then
  !           P%normal = -P%normal
  !           P%rho = -P%rho
  !         end if
  !         exit
  !       end if
  !     end do

  !     ! Get the volume fraction behind the plane in this cell.
  !     call poly%init(ierr, this%x, vol=this%volume)
  !     INSIST(ierr == 0)
  !     volumes(matl(1)) = poly%volume_behind_plane(P, ierr)
  !     volumes(matl(2)) = this%volume - volumes(matl(1))
  !     INSIST(ierr == 0)
  !   else
  !     ! Interface doesn't actually intersect this subtet.
  !     ! We can get here if it just barely touches an edge.
  !     volumes(this%matl_at_node(1)) = this%volume
  !   end if
    
  ! end function volume_from_level_set

end module vof_init_NEW
