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
    integer, allocatable :: body_at_node(:)
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
    procedure, private :: contains_interface
    procedure, private :: volumes_from_intersection
    procedure, private :: volumes_from_centroid
  end type dnc_cell

  ! These structures define the cell -> tet subdivision.
  !
  ! Hex -> tet performs nearly the same subdivision as Cubit's HTet structured
  ! conversion, except that it splits every edge. It adds a new node to each
  ! face, edge, and one to the interior, and generates 40 tets per hex.
  !
  ! Note the hex and tet refinement descriptions will generate inverted tets.
  ! This is so we can guarantee symmetric and stable refinement. Tet
  ! refinement Follow's Bey's red refinement algorithm, which specifies a
  ! dominant diagonal. To keep that diagonal consistent, some tets must be
  ! mirrored.
  integer, parameter :: tet4_ntet = 8
  integer, parameter :: pyr5_ntet = 4
  integer, parameter :: wed6_ntet = 17
  !integer, parameter :: hex8_ntet = 28
  !integer, parameter :: hex8_ntet = 40
  !integer, target :: hex8_to_tet(4,hex8_ntet)
  integer, target :: pyr5_to_tet(4,pyr5_ntet), wed6_to_tet(4,wed6_ntet), tet4_to_tet(4,tet4_ntet)
  data PYR5_TO_TET/1,2,5,6, 2,3,5,6, 3,4,5,6, 1,4,5,6/
  data WED6_TO_TET/1,4,7,9, 2,5,7,8, 3,6,8,9, 7,8,9,10, 7,8,9,11, 1,7,9,10, 2,7,8,10, &
      3,8,9,10, 1,2,7,10, 1,3,9,10, 2,3,8,10, 5,7,8,11, 6,8,9,11, 4,7,9,11, &
      5,6,8,11, 4,5,7,11, 4,6,9,11/

  ! ! htet+
  ! data HEX8_TO_TET/&
  !     17,5,9,12, 19,6,9,10, 21,7,10,11, 22,8,11,12, &
  !     17,1,9,12, 19,2,9,10, 21,3,10,11, 22,4,11,12, &

  !     24,5,14,12, 25,6,10,14, 18,2,13,10, 16,1,13,12, &
  !     24,8,14,12, 25,7,10,14, 18,3,13,10, 16,4,13,12, &

  !     23,5,14,9, 26,8,14,11, 20,4,13,11, 15,1,13,9, &
  !     23,6,14,9, 26,7,14,11, 20,3,13,11, 15,2,13,9, &
      
  !     9,12,5,14, 10,9,6,14, 11,10,7,14, 12,11,8,14, &
  !     9,12,1,13, 10,9,2,13, 11,10,3,13, 12,11,4,13, &

  !     9,12,14,27, 10,9,14,27, 11,10,14,27, 12,11,14,27, &
  !     9,12,13,27, 10,9,13,27, 11,10,13,27, 12,11,13,27/

  ! htet
  ! data HEX8_TO_TET/&
  !     1,5,9,12, 2,9,6,10, 3,10,7,11, 4,11,8,12, &
  !     5,12,8,14, 6,7,10,14, 2,10,3,13, 1,4,12,13, &
  !     5,6,9,14, 7,8,11,14, 3,11,4,13, 1,9,2,13, &
      
  !     5,9,12,14, 6,10,9,14, 7,11,10,14, 8,12,11,14, &
  !     1,12,9,13, 2,9,10,13, 3,10,11,13, 4,11,12,13, &
      
  !     14,9,12,15, 14,10,9,15, 14,11,10,15, 14,12,11,15, &
  !     13,12,9,15, 13,9,10,15, 13,10,11,15, 13,11,12,15/

  ! mine
  integer, parameter :: hex8_ntet = 24
  integer, target :: hex8_to_tet(4,hex8_ntet)
  data HEX8_TO_TET/&
      1,2,9,15,  2,6,9,15,  6,5,9,15,  5,1,9,15,&
      2,3,10,15, 3,7,10,15, 7,6,10,15, 6,2,10,15,&
      3,4,11,15, 4,8,11,15, 8,7,11,15, 7,3,11,15,&
      4,1,12,15, 1,5,12,15, 5,8,12,15, 8,4,12,15,&
      2,1,13,15, 1,4,13,15, 4,3,13,15, 3,2,13,15,&
      5,6,14,15, 6,7,14,15, 7,8,14,15, 8,5,14,15/

  ! ! mine+
  ! integer, parameter :: hex8_ntet = 48
  ! integer, target :: hex8_to_tet(4,hex8_ntet)
  ! data HEX8_TO_TET/&
  !     1,2,9,15,  2,6,9,15,  6,5,9,15,  5,1,9,15,&
  !     2,3,10,15, 3,7,10,15, 7,6,10,15, 6,2,10,15,&
  !     3,4,11,15, 4,8,11,15, 8,7,11,15, 7,3,11,15,&
  !     4,1,12,15, 1,5,12,15, 5,8,12,15, 8,4,12,15,&
  !     2,1,13,15, 1,4,13,15, 4,3,13,15, 3,2,13,15,&
  !     5,6,14,15, 6,7,14,15, 7,8,14,15, 8,5,14,15/

  ! 1,5,9,12, 2,9,6,10, 3,10,7,11, 4,11,8,12, 1,9,2,13, 2,10,3,13, 3,11,4,13, &
  ! 1,4,12,13, 5,6,9,14, 6,7,10,14, 7,8,11,14, 5,12,8,14, 9,10,13,14, &
  ! 10,11,13,14, 11,12,13,14, 9,13,12,14, 5,9,12,14, 6,10,9,14, 7,11,10,14, &
  ! 8,12,11,14, 1,12,9,13, 2,9,10,13, 3,10,11,13, 4,11,12,13/
  ! data TET4_TO_TET/1,5,6,7, 2,8,5,9, 3,6,8,10, 7,9,10,4, &
  !     5,6,7,10, 5,8,6,10, 5,7,9,10, 5,9,8,10/
  ! data TET4_TO_TET/1,5,6,7, 5,2,8,9, 6,8,3,10, 7,9,10,4, &
  !     5,6,7,9, 5,6,8,9, 6,7,9,10, 6,8,9,10/
  data TET4_TO_TET/2,8,5,9, 8,3,6,10, 5,6,1,7, 9,10,7,4, &
      8,5,9,10, 8,5,6,10, 5,9,10,7, 5,6,10,7/

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
    if (body_id%nbody > 1) then
      do i = 1, mesh%ncell_onP
        call cell%init(i, mesh, body_id, recursion_limit)
        vof(:,i) = cell%volumes() !/ cell%volume
        ! if (any(i == [323, 324, 325])) then
        !   print *, 'cell: ', i
        !   print '(a,3es14.4)', 'xc:  ', sum(mesh%x(:,mesh%cnode(mesh%xcnode(i):mesh%xcnode(i+1)-1)), dim=2) / 8
        !   print '(a,3es14.4)', 'vof: ', vof(:,i) / mesh%volume(i)
        !   print *
        ! end if
        ! if (all(vof(:,i) > 0)) then
        !   print '(i6,4es13.3)', i, vof(:,i), mesh%volume(i) / sum(vof(:,i))
        !   print *
        ! end if
        vof(:,i) = vof(:,i) * mesh%volume(i) / sum(vof(:,i)) ! correct for errors in planar assumption
      end do
    else
      ! If there's only one body, do something faster.
      vof(1,:) = mesh%volume(:mesh%ncell_onP)
    end if
    !call TLS_Fatal("debugging")
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
        this%x(:,this%nnode+1) = sum(this%x(:,pyr5_faces(pyr5_xface(5):pyr5_xface(6)-1)),dim=2) / pyr5_fsize(5)
      case (6)
        this%to_tet => wed6_to_tet
        allocate(this%x(3,this%nnode+nface))
        this%x(:,:this%nnode) = mesh%x(:,cn)
        do j = 1, nface
          this%x(:,this%nnode+j) = sum(this%x(:,wed6_faces(wed6_xface(j):wed6_xface(j+1)-1)),dim=2) / wed6_fsize(j)
        end do
      case (8)
        this%to_tet => hex8_to_tet

        ! ! htet+
        ! allocate(this%x(3,this%nnode+nface+12+1))
        ! this%x(:,:this%nnode) = mesh%x(:,cn)
        ! do j = 1, nface
        !   this%x(:,this%nnode+j) = sum(this%x(:,hex8_faces(hex8_xface(j):hex8_xface(j+1)-1)),dim=2) / hex8_fsize(j)
        ! end do
        ! do j = 1, 12
        !   this%x(:,this%nnode+nface+j) = sum(this%x(:,hex8_edges(:,j)), dim=2) / 2
        ! end do
        ! this%x(:,this%nnode+nface+13) = sum(this%x(:,:this%nnode),dim=2) / this%nnode

        ! htet & mine
        allocate(this%x(3,this%nnode+nface+1))
        this%x(:,:this%nnode) = mesh%x(:,cn)
        do j = 1, nface
          this%x(:,this%nnode+j) = sum(this%x(:,hex8_faces(hex8_xface(j):hex8_xface(j+1)-1)),dim=2) / hex8_fsize(j)
        end do
        this%x(:,this%nnode+nface+1) = sum(this%x(:,:this%nnode),dim=2) / this%nnode
      case default
        INSIST(.false.)
      end select

      this%ntet = size(this%to_tet,dim=2)
    end associate

    allocate(this%body_at_node(size(this%x,dim=2)))
    do j = 1, this%nnode
      this%body_at_node(j) = body_id%body_at_point(this%x(:,j), this%cellid)
    end do
    !this%body_at_node = [(body_id%body_at_point(this%x(:,j), this%cellid), j = 1, size(this%x,dim=2))]
    INSIST(.not.any(this%body_at_node(:this%nnode)==0))

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

    integer :: i, j, body(3), nbody
    real(r8) :: xc(3), x(3), v(this%body_id%nbody)
    type(dnc_cell) :: subtet

    volumes = 0

    ! if (.not.this%is_subtet .and. any(this%body_at_node(:this%nnode) == 1)) then
    !   do i = 1, this%nnode
    !     print '(i6,3es13.3)', this%body_id%body_at_point(this%x(:,i), this%cellid), this%x(:,i)
    !   end do
    !   print '(8i6)', this%body_at_node(:this%nnode)
    ! end if

    ! If the cell contains an interface (and therefore has at least
    ! two materials and we haven't yet hit our recursion limit, divide
    ! the cell and repeat.
    if (this%recursion_height > 0 .and. this%contains_interface()) then
      ! tally the volumes from refined sub-cells
      call this%prepare_subtet(subtet)
      do i = 1, this%ntet
        ! if (any(this%cellid == [323, 324, 325])) then
        !   print *, i, this%ntet, this%recursion_height
        ! end if
        call this%get_subtet(subtet, i)
        ! v = subtet%volumes()
        !print '(i6,3es13.3)', i, v
        ! if (any(v<0)) then
        !   print *, i, this%ntet, this%recursion_height, subtet%volume
        !   INSIST(.false.)
        ! end if
        volumes = volumes + subtet%volumes()
      end do
    else ! if (this%body_id%nbody > 2 .or. .not.this%contains_interface()) then

      ! If this is a subtet containing two bodies, we can do a plane
      ! reconstruction and intersection. Otherwise, get the volume from
      ! the cell center.
      nbody = 1
      if (this%is_subtet .and. this%contains_interface()) then
        body(1) = this%body_at_node(1)
        do i = 2, this%nnode
          if (all(this%body_at_node(i) /= body(:nbody))) then
            nbody = nbody + 1
            body(nbody) = this%body_at_node(i)
            if (nbody > 2) exit
          end if
        end do
      end if

      if (this%is_subtet) then
        ! Subtet volume is calculated here for the first time, since
        ! we don't need subtet volumes at every refinement level.
        ! We take the absolute value, since some of our refined tets
        ! are inverted (see comment on Bey's red algorithm above).
        this%volume = abs(tet_volume(this%x))
      end if

      if (this%is_subtet .and. nbody == 2) then
        ! For two-material cells, we can reconstruct an interface from
        ! 3 interface-points identified by linearly interpolating a
        ! signed distance function.
        volumes = this%volumes_from_intersection()
      else
        ! If we're past the recursion limit or the cell does not contain
        ! an interface, calculate the vof in this cell based on the
        ! materials at its nodes.
        volumes = this%volumes_from_centroid()
        ! xc = sum(this%x(:,:this%nnode), dim=2) / this%nnode
        ! ! if (any(this%cellid == [323, 324, 325])) then
        ! !   print '(a,3es13.3,i4)', 'vol: ', xc, this%body_id%body_at_point(xc, this%cellid)
        ! !   ! print '(a,3es13.3,i4)', 'x: ', this%x(:,1)
        ! !   ! print '(a,3es13.3,i4)', 'x: ', this%x(:,2)
        ! !   ! print '(a,3es13.3,i4)', 'x: ', this%x(:,3)
        ! !   ! print '(a,3es13.3,i4)', 'x: ', this%x(:,4)
        ! !   ! print *
        ! ! end if
        ! do i = 1, this%nnode
        !   x = (xc + 2*this%x(:,i)) / 3 ! point weighted to just inside the cell
        !   j = this%body_id%body_at_point(x, this%cellid)
        !   volumes(j) = volumes(j) + 1
        ! end do
        ! ! do i = 1, this%body_id%nbody
        ! !   volumes(i) = count(this%body_at_node(:this%nnode) == i)
        ! ! end do
        ! ! if (.not.this%is_subtet .and. all(volumes > 0)) then
        ! !   do i = 1, this%nnode
        ! !     print '(3es13.3)', this%x(:,i)
        ! !   end do
        ! !   print '(2es13.3)', volumes
        ! !   INSIST(.false.)
        ! ! end if
        ! volumes = volumes / this%nnode * this%volume
        ! else
      end if
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

    ! Compute this cell's new nodes and find bodies at those points.
    if (this%is_subtet) then
      do j = 1, 6
        this%x(:,4+j) = sum(this%x(:,tet4_edges(:,j)), dim=2) / 2
        !this%body_at_node(4+j) = this%body_id%body_at_point(this%x(:,4+j), this%cellid)
      end do
      !INSIST(.not.any(this%body_at_node==0))
    end if

    do j = this%nnode+1, size(this%x,dim=2)
      this%body_at_node(j) = this%body_id%body_at_point(this%x(:,j), this%cellid)
    end do
    INSIST(.not.any(this%body_at_node==0))

    ! Subtet structure. The geometry itself is generated by get_subtet.
    subtet%body_id => this%body_id
    subtet%to_tet => tet4_to_tet
    subtet%ntet = tet4_ntet
    subtet%recursion_height = this%recursion_height - 1
    subtet%cellid = this%cellid
    subtet%is_subtet = .true.
    subtet%nnode = 4
    allocate(subtet%x(3,subtet%nnode+6), subtet%body_at_node(subtet%nnode+6))

  end subroutine prepare_subtet


  ! Return a complete subtet geometry from a given subtet ID.
  subroutine get_subtet(this, subtet, i)
    class(dnc_cell), intent(in) :: this
    class(dnc_cell), intent(inout) :: subtet
    integer, intent(in) :: i
    ASSERT(i > 0 .and. i <= this%ntet)
    subtet%x(:,:4) = this%x(:,this%to_tet(:,i))
    subtet%body_at_node(:4) = this%body_at_node(this%to_tet(:,i))
  end subroutine get_subtet


  logical function contains_interface(this)

    class(dnc_cell), intent(in) :: this

    integer :: i, b
    real(r8) :: x(3), xc(3)

    contains_interface = any(this%body_at_node(2:this%nnode) /= this%body_at_node(1))
    ! if (any(this%cellid == [323, 324, 325])) then
    !   print *, 'mixed: ', contains_interface
    ! end if
    ! contains_interface = .false.
    ! xc = sum(this%x(:,:this%nnode), dim=2) / this%nnode
    ! x = (xc + 4*this%x(:,1)) / 5 ! point weighted to just inside the cell
    ! b = this%body_id%body_at_point(x, this%cellid)
    ! do i = 2, this%nnode
    !   x = (xc + 4*this%x(:,i)) / 5 ! point weighted to just inside the cell
    !   if (b /= this%body_id%body_at_point(x, this%cellid)) then
    !     contains_interface = .true.
    !     return
    !   end if
    ! end do
    
  end function contains_interface


  ! For two-material cells, we can reconstruct an interface from
  ! three interface-points identified by linearly interpolating a
  ! signed distance function. Note that this function ASSUMES it's
  ! operating on a cell with only two materials present, and that it
  ! is operating on a tet
  function volumes_from_intersection(this) result(volumes)

    use plane_type
    use pure_polyhedron_type
    use cell_geometry, only: cross_product, normalized, tet_volume
    use cell_topology, only: tet4_edges

    class(dnc_cell), intent(in) :: this
    real(r8) :: volumes(this%body_id%nbody)

    integer, parameter :: nedge = 6
    integer :: ierr, i, j, b1, b2
    real(r8) :: xp(3,3), xn(3,4), signed_distance(4), s1, s2, cell_volume
    type(plane) :: P
    type(pure_polyhedron) :: poly

    ASSERT(this%is_subtet)

    volumes = 0
    
    ! Identify the two bodies present. We will use the
    ! highest-priority body's signed-distance function.
    ! Highest-priority material goes into body(1).
    b1 = this%body_at_node(1)
    do i = 2, this%nnode
      if (this%body_at_node(i) /= b1) then
        b2 = this%body_at_node(i)
        exit
      end if
    end do
    if (b1 > b2) then
      i = b1
      b1 = b2
      b2 = i
    end if
    
    ! Get 3 points intersecting the edges of this cell.
    ! First we look at the nodes, then at intersections along edges.
    ! This prevents inadvertantly adding nodes multiple times if they
    ! are directly intersected, as we can skip intersected nodes while
    ! checking edges.
    j = 0
    do i = 1, this%nnode
      signed_distance(i) = this%body_id%body(b1)%signed_distance(this%x(:,i))
      if (abs(signed_distance(i)) < 1e-10) signed_distance(i) = 0
      if (signed_distance(i) == 0) then
        j = j + 1
        xp(:,j) = this%x(:,i)
        if (j == 3) then
          ! If 3 nodes of a tet are intersected, the cell is not cut by the plane
          volumes = this%volumes_from_centroid()
          return
        end if
      end if
    end do
    if (j < 3) then
      do i = 1, nedge
        ! add edge intersections, skipping node intersections here
        s1 = signed_distance(tet4_edges(1,i))
        s2 = signed_distance(tet4_edges(2,i))
        if (s1*s2 < 0) then
          j = j + 1
          xp(:,j) = this%x(:,tet4_edges(1,i)) + &
              s1/(s1-s2) * (this%x(:,tet4_edges(2,i)) - this%x(:,tet4_edges(1,i)))

          print '(i4,3es13.3,"   ",2es13.3)', i, xp(:,j), s1, s2

          if (j == 3) exit
        end if
      end do
    end if

    if (j == 3) then
      ! Construct a plane from the 3 points.
      P%normal = normalized(cross_product(xp(:,2) - xp(:,1), xp(:,3) - xp(:,1)))
      P%rho = dot_product(xp(:,1), P%normal)

      ! Orient the interface so the normal points outside relative to b.
      ! Signed distance is defined to be positive outside the body.
      do i = 1, this%nnode
        s1 = P%signed_distance(this%x(:,i))

        if (abs(s1) > 0) then
          if ((s1 > 0 .and. this%body_at_node(i) == b1) .or. &
              (s1 < 0 .and. this%body_at_node(i) /= b1)) then
            P%normal = -P%normal
            P%rho = -P%rho
            print *, 'swap'
          end if
          exit
        end if
      end do

      print '(a,4es13.3)', 'normal2: ', P%normal, P%rho
      ASSERT(norm2(P%normal) > 0)

      ! If this volume is inverted, create a normal tet.
      xn = this%x(:,:this%nnode)
      if (tet_volume(xn) < 0) then
        xp(:,1) = xn(:,1)
        xn(:,1) = xn(:,2)
        xn(:,2) = xp(:,1)
      end if

      ! Get the volume fraction behind the plane in this cell.
      call poly%init(ierr, xn, vol=this%volume)
      INSIST(ierr == 0)
      volumes(b1) = poly%volume_behind_plane(P, ierr)
      volumes(b2) = this%volume - volumes(b1)
      INSIST(ierr == 0)
    else
      ! Interface doesn't actually intersect this subtet.
      ! We can get here if it just barely touches an edge.
      volumes = this%volumes_from_centroid()
    end if
    
  end function volumes_from_intersection


  function volumes_from_centroid(this) result(volumes)
    class(dnc_cell), intent(in) :: this
    real(r8) :: volumes(this%body_id%nbody)
    real(r8) :: xc(3)
    volumes = 0
    xc = sum(this%x(:,:this%nnode), dim=2) / this%nnode
    volumes(this%body_id%body_at_point(xc, this%cellid)) = this%volume
  end function volumes_from_centroid

end module vof_init_NEW
