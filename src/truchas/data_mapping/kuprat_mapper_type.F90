!!
!! KUPRAT_MAPPER_TYPE
!!
!! This provides a Truchas-specific class that encapsulates Andrew Kuprat's
!! grid-to-grid mapping package.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module kuprat_mapper_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use data_mapper_class
  use grid_mapping_module
  use parallel_communication, only: is_IOP, collate, distribute, global_sum
  implicit none
  private

  type, extends(data_mapper), public :: kuprat_mapper
    integer :: n1, n2
    type(grid_int_vols) :: gmd
  contains
    procedure :: init
    procedure :: map_field
  end type

  !! Re-export the mapping type options
  public :: LOCALLY_CONSERVATIVE, LOCALLY_BOUNDED, GLOBALLY_CONSERVATIVE

contains

  subroutine init(this, mesh1, mesh2)

    use base_mesh_class

    class(kuprat_mapper), intent(out) :: this
    class(base_mesh), intent(in) :: mesh1, mesh2

    type(gm_mesh) :: gm_mesh1, gm_mesh2

    call base_to_gm_mesh(mesh1, gm_mesh1)
    call base_to_gm_mesh(mesh2, gm_mesh2)

    if (is_IOP) call compute_int_volumes(gm_mesh1, gm_mesh2, this%gmd)

    this%n1 = mesh1%ncell_onP
    this%n2 = mesh2%ncell_onP

  end subroutine init


  subroutine map_field(this, src, dest, defval, map_type, pullback)

    class(kuprat_mapper), intent(in) :: this
    real(r8), intent(in)  :: src(:)
    real(r8), intent(out) :: dest(:)
    real(r8), intent(in)  :: defval
    integer,  intent(in)  :: map_type
    logical,  intent(in), optional :: pullback

    integer :: n
    logical :: reverse_order, preserve_constants, exactly_conservative
    real(r8), allocatable :: col_src(:), col_dest(:)

    select case (map_type)
    case (LOCALLY_CONSERVATIVE)
      preserve_constants   = .false.
      exactly_conservative = .false.
    case (LOCALLY_BOUNDED)
      preserve_constants   = .true.
      exactly_conservative = .false.
    case (GLOBALLY_CONSERVATIVE)
      preserve_constants   = .false.
      exactly_conservative = .true.
    case default
      INSIST(.false.)
    end select

    reverse_order = .false.
    if (present(pullback)) reverse_order = pullback

    if (reverse_order) then
      ASSERT(size(src)  == this%n2)
      ASSERT(size(dest) == this%n1)
    else
      ASSERT(size(src)  == this%n1)
      ASSERT(size(dest) == this%n2)
    end if

    n = global_sum(size(src))
    allocate(col_src(merge(n,0,is_IOP)))
    n = global_sum(size(dest))
    allocate(col_dest(merge(n,0,is_IOP)))

    call collate(col_src, src)

    if (is_IOP) call map_cell_field(col_src, col_dest, this%gmd, defval=defval, &
                                    reverse_order=reverse_order, &
                                    preserve_constants=preserve_constants, &
                                    exactly_conservative=exactly_conservative)

    call distribute(dest, col_dest)

  end subroutine map_field

  !! This auxiliary subroutine converts a distributed Truchas mesh to a
  !! serial GM_MESH object required by Kuprat's grid-to-grid mapping package.

  subroutine base_to_gm_mesh(inmesh, outmesh)

    use base_mesh_class
    use unstr_mesh_type
    use simpl_mesh_type

    class(base_mesh), intent(in) :: inmesh
    type(gm_mesh), intent(out) :: outmesh

    select type (inmesh)
    class is (unstr_mesh)
      call unstr_to_gm_mesh(inmesh, outmesh)
    class is (simpl_mesh)
      call simpl_to_gm_mesh(inmesh, outmesh)
    class default
      INSIST(.false.)
    end select

  end subroutine base_to_gm_mesh

  !! This auxiliary subroutine converts a distributed unstructured mesh stored
  !! as an UNSTR_MESH object to a global GM_MESH object on the IO processor.

  subroutine unstr_to_gm_mesh(inmesh, outmesh)

    use unstr_mesh_type

    type(unstr_mesh), intent(in) :: inmesh
    type(gm_mesh), intent(out) :: outmesh

    integer, parameter :: TET_NODE_MAP(8) = [1,1,2,3,4,4,4,4]
    integer, parameter :: PYR_NODE_MAP(8) = [1,2,3,4,5,5,5,5]
    integer, parameter :: PRI_NODE_MAP(8) = [1,4,5,2,3,6,6,3]

    integer :: j, ncell, nnode
    integer, allocatable :: xcnode(:), cnode(:), blockid(:), tmp(:)

    ncell = inmesh%cell_ip%global_size()
    nnode = inmesh%node_ip%global_size()

    call inmesh%get_global_cnode_array(xcnode, cnode)

    if (is_IOP) then

      outmesh%nnod = nnode
      outmesh%nelt = ncell

      if (all(xcnode(2:)-xcnode(:ncell) == 4)) then
        !! Special case: pure tet mesh
        allocate(outmesh%node_elt(4,ncell))
        do j = 1, ncell
          associate (cnode => cnode(xcnode(j):xcnode(j+1)-1))
            outmesh%node_elt(:,j) = cnode
          end associate
        end do
      else
        !! General case: hexes and degenerate hexes (as needed)
        allocate(outmesh%node_elt(8,ncell))
        do j = 1, ncell
          associate (cnode => cnode(xcnode(j):xcnode(j+1)-1))
            select case (size(cnode))
            case (4)  ! tet
              outmesh%node_elt(:,j) = cnode(TET_NODE_MAP)
            case (5)  ! pyramid
              outmesh%node_elt(:,j) = cnode(PYR_NODE_MAP)
            case (6)  ! prism
              outmesh%node_elt(:,j) = cnode(PRI_NODE_MAP)
            case (8)  ! hex
              outmesh%node_elt(:,j) = cnode
            case default
              INSIST(.false.)
            end select
          end associate
        end do
      end if

    end if

    !! Node coordinate data
    call inmesh%get_global_x_array(outmesh%pos_node)

    !! Element block ID data; need to translate from cell set bitmask
    allocate(blockid(inmesh%ncell_onP))
    do j = 1, inmesh%ncell_onP
      associate (bitmask => inmesh%cell_set_mask(j))
        INSIST(popcnt(bitmask) == 1)
        blockid(j) = inmesh%cell_set_id(trailz(bitmask))
      end associate
    end do
    allocate(tmp(merge(ncell,0,is_IOP)))
    call collate(tmp, blockid)
    call move_alloc(tmp, outmesh%block_elt)

  end subroutine unstr_to_gm_mesh

  !! This auxiliary subroutine converts a distributed tet mesh stored as a
  !! SIMPL_MESH object to a global GM_MESH object on the IO processor.

  subroutine simpl_to_gm_mesh(inmesh, outmesh)

    use simpl_mesh_type

    type(simpl_mesh), intent(in) :: inmesh
    type(gm_mesh), intent(out) :: outmesh

    integer :: j, tmp
    real(r8), allocatable :: volume(:)

    call inmesh%get_global_cnode_array(outmesh%node_elt)
    call inmesh%get_global_cblock_array(outmesh%block_elt)
    call inmesh%get_global_x_array(outmesh%pos_node)
    outmesh%nnod = size(outmesh%pos_node,dim=2)
    outmesh%nelt = size(outmesh%node_elt,dim=2)

    !! Ensure each tet is positively oriented with respect to its volume.
    call inmesh%get_global_volume_array(volume)
    if (is_IOP) then
      do j = 1, size(volume)
        if (volume(j) < 0.0) then
          tmp = outmesh%node_elt(3,j)
          outmesh%node_elt(3,j) = outmesh%node_elt(4,j)
          outmesh%node_elt(4,j) = tmp
        end if
      end do
    end if

  end subroutine simpl_to_gm_mesh

end module kuprat_mapper_type
