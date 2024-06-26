!!
!! PORTAGE_MAPPER_TYPE
!!
!! This provides a Truchas-specific class that encapsulates the unstructured
!! mesh data mapping functionality from the Portage library. Note that a major
!! part of this interface is implemented in a custom C++ interface required by
!! the Portage library. Portage is at https://github.com/laristra/portage
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

module portage_mapper_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_loc, c_int, c_double
  use parallel_communication, only: is_IOP, gather, scatter, global_sum
  use data_mapper_class
  implicit none
  private

  !! This holds an internal serialized copy of a parallel Truchas mesh
  type :: serial_mesh
    integer(c_int) :: num_cell, num_face, num_node
    integer(c_int), allocatable :: xcnode(:), cnode(:)
    integer(c_int), allocatable :: xcface(:), cface(:)
    integer(c_int), allocatable :: xfnode(:), fnode(:)
    integer(c_int), allocatable :: xncell(:), ncell(:)
    real(c_double), allocatable :: coord(:,:)
    integer(c_int), allocatable :: cfdir(:)
    integer(c_int), allocatable :: blockid(:)
  end type serial_mesh

  !! Mesh data needed by Portage
  type, bind(c) :: unstr_mesh_data
    integer(c_int) :: num_cell, num_face, num_node
    type(c_ptr) :: xcnode, cnode   ! int arrays
    type(c_ptr) :: xcface, cface   ! int arrays
    type(c_ptr) :: xfnode, fnode   ! int arrays
    type(c_ptr) :: xncell, ncell   ! int arrays
    type(c_ptr) :: coord           ! double array
    type(c_ptr) :: cfdir           ! int array
    type(c_ptr) :: blockid         ! int array
  end type

  type, extends(data_mapper), public :: portage_mapper
    type(c_ptr) :: mapper = c_null_ptr ! opaque handle to the Portage mapper object
    type(c_ptr) :: mapper_pullback = c_null_ptr ! opaque handle to the Portage mapper object
    type(serial_mesh) :: serial_mesh1, serial_mesh2
    integer :: n1, n2
  contains
    procedure :: init
    procedure :: map_field
    final :: delete
  end type portage_mapper

  interface
    function portage_mapper_new(src_mesh, tgt_mesh) result(mapper) bind(c)
      import unstr_mesh_data, c_ptr
      type(unstr_mesh_data), intent(in) :: src_mesh, tgt_mesh
      type(c_ptr) :: mapper
    end function
    subroutine portage_mapper_delete(mapper) bind(c)
      import c_ptr
      type(c_ptr), value :: mapper
    end subroutine
    subroutine portage_map_field(mapper, src_size, src, dest_size, dest, method) bind(c)
      import c_ptr, c_int, c_double
      type(c_ptr), value :: mapper
      integer(c_int), value :: src_size, dest_size
      real(c_double), intent(in) :: src(*)
      real(c_double), intent(out) :: dest(*)
      integer(c_int), value :: method
    end subroutine
  end interface

  !! Re-export the mapping type options
  public :: LOCALLY_CONSERVATIVE, LOCALLY_BOUNDED, GLOBALLY_CONSERVATIVE

contains

  !! Final subroutine for PORTAGE_MAPPER objects
  subroutine delete(this)
    type(portage_mapper), intent(inout) :: this
    call portage_mapper_delete(this%mapper)
    call portage_mapper_delete(this%mapper_pullback)
  end subroutine


  subroutine init(this, mesh1, mesh2)

    use base_mesh_class

    class(portage_mapper), intent(out) :: this
    class(base_mesh), intent(in) :: mesh1, mesh2

    type(unstr_mesh_data) :: src_mesh, tgt_mesh

    call base_to_serial(mesh1, this%serial_mesh1)
    call base_to_serial(mesh2, this%serial_mesh2)

    if (is_IOP) then
      call serial_to_portage(this%serial_mesh1, src_mesh)
      call serial_to_portage(this%serial_mesh2, tgt_mesh)

      this%mapper = portage_mapper_new(src_mesh, tgt_mesh)
      this%mapper_pullback = portage_mapper_new(tgt_mesh, src_mesh)
    end if

    this%n1 = mesh1%ncell_onP
    this%n2 = mesh2%ncell_onP

  end subroutine init


  subroutine map_field(this, src, dest, defval, map_type, pullback)

    use,intrinsic :: ieee_exceptions

    class(portage_mapper), intent(in) :: this
    real(r8), intent(in)  :: src(:)
    real(r8), intent(out) :: dest(:)
    real(r8), intent(in)  :: defval
    integer,  intent(in)  :: map_type
    logical,  intent(in), optional :: pullback

    integer :: method, n
    logical :: reverse_order, ieee_invalid_mode
    real(r8), allocatable :: col_src(:), col_dest(:)

    select case (map_type)
    case (LOCALLY_CONSERVATIVE)
      method = 0
    case (LOCALLY_BOUNDED)
      method = 1
    case (GLOBALLY_CONSERVATIVE)
      method = 2
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

    call gather(src, col_src)

    if (is_IOP) then
      col_dest = defval
      ! Portage currently performs operations which trip floating
      ! exception handling. Disable while in Portage.
      call ieee_get_halting_mode(ieee_invalid, ieee_invalid_mode)
      call ieee_set_halting_mode(ieee_invalid, .false.)
      if (reverse_order) then
        call portage_map_field(this%mapper_pullback, &
            size(col_src), col_src, size(col_dest), col_dest, method)
      else
        call portage_map_field(this%mapper, &
            size(col_src), col_src, size(col_dest), col_dest, method)
      end if
      call ieee_set_halting_mode(ieee_invalid, ieee_invalid_mode)
    end if

    call scatter(col_dest, dest)

  end subroutine map_field

  !! This auxiliary subroutine converts a distributed Truchas mesh to
  !! a serial UNSTR_MESH_DATA object that will be passed to Portage.

  subroutine base_to_serial(inmesh, outmesh)

    use base_mesh_class
    use unstr_mesh_type
    use simpl_mesh_type

    class(base_mesh), intent(in) :: inmesh
    type(serial_mesh), intent(out) :: outmesh

    select type (inmesh)
    class is (unstr_mesh)
      call unstr_to_serial(inmesh, outmesh)
    class is (simpl_mesh)
      call simpl_to_serial(inmesh, outmesh)
    class default
      INSIST(.false.)
    end select

    if (is_IOP) call add_node_to_cell(outmesh)

  end subroutine base_to_serial

  !! This auxiliary subroutine creates on the IOP a SERIAL_MESH copy of the
  !! given parallel UNSTR_MESH object, which holds just those data components
  !! needed by Portage.

  subroutine unstr_to_serial(inmesh, outmesh)

    use unstr_mesh_type

    type(unstr_mesh), intent(in) :: inmesh
    type(serial_mesh), intent(out) :: outmesh

    if (is_IOP) then
      outmesh%num_cell = inmesh%cell_imap%global_size
      outmesh%num_face = inmesh%face_imap%global_size
      outmesh%num_node = inmesh%node_imap%global_size
    end if

    call inmesh%get_global_cnode_array(outmesh%xcnode, outmesh%cnode)
    call inmesh%get_global_cface_array(outmesh%xcface, outmesh%cface)
    call inmesh%get_global_fnode_array(outmesh%xfnode, outmesh%fnode)
    call inmesh%get_global_x_array(outmesh%coord)

    block !! Convert the packed CFPAR info to the CDIR array.
      integer :: j, k
      integer, allocatable :: cfdir(:)
      allocate(cfdir(inmesh%xcface(inmesh%ncell_onP+1)-1))
      do j = 1, inmesh%ncell_onP
        associate(dir => cfdir(inmesh%xcface(j):inmesh%xcface(j+1)-1))
          do k = 1, size(dir)
            if (btest(inmesh%cfpar(j),pos=k)) then
              dir(k) = -1
            else
              dir(k) = 1
            end if
          end do
        end associate
      end do
      allocate(outmesh%cfdir, mold=outmesh%cface)
      call gather(cfdir, outmesh%cfdir)
    end block

    block !! Convert CELLS_SET_MASK info to the BLOCKID array.
      integer :: j
      integer, allocatable :: blockid(:)
      allocate(blockid(inmesh%ncell_onP))
      do j = 1, inmesh%ncell_onP
        associate (bitmask => inmesh%cell_set_mask(j))
          INSIST(popcnt(bitmask) == 1)
          blockid(j) = inmesh%cell_set_id(trailz(bitmask))
        end associate
      end do
      allocate(outmesh%blockid(merge(outmesh%num_cell,0,is_IOP)))
      call gather(blockid, outmesh%blockid)
    end block

  end subroutine unstr_to_serial

  !! This auxiliary subroutine creates on the IOP a SERIAL_MESH copy of the
  !! given parallel SIMPL_MESH object, which holds just those data components
  !! needed by Portage.

  subroutine simpl_to_serial(inmesh, outmesh)

    use simpl_mesh_type

    type(simpl_mesh), intent(in), target :: inmesh
    type(serial_mesh), intent(out) :: outmesh

    integer :: j
    integer, pointer :: array(:)

    if (is_IOP) then
      outmesh%num_cell = inmesh%cell_imap%global_size
      outmesh%num_face = inmesh%face_imap%global_size
      outmesh%num_node = inmesh%node_imap%global_size
    end if

    !! CNODE -- need to reorient for positive volume
    block
      integer :: tmp
      integer, allocatable, target :: cnode(:,:)
      cnode = inmesh%cnode(:,:inmesh%ncell_onP)
      do j = 1, inmesh%ncell_onP
        if (inmesh%volume(j) < 0) then
          tmp = cnode(3,j)
          cnode(3,j) = cnode(4,j)
          cnode(4,j) = tmp
        end if
      end do
      array(1:size(cnode)) => cnode
      allocate(outmesh%cnode(merge(4*outmesh%num_cell,0,is_IOP)))
      call gather(inmesh%node_imap%global_index(array), outmesh%cnode)
    end block

    !! XCNODE
    allocate(outmesh%xcnode(merge(outmesh%num_cell+1,0,is_IOP)))
    if (is_IOP) then
      do concurrent (j = 0:outmesh%num_cell)
        outmesh%xcnode(j+1) = 1 + 4*j
      end do
    end if

    !! CFACE -- NB: faces no longer opposite corresponding vertex
    allocate(outmesh%cface(merge(4*outmesh%num_cell,0,is_IOP)))
    associate (cface => inmesh%cface(:,:inmesh%ncell_onP))
      array(1:size(cface)) => cface
      call gather(inmesh%face_imap%global_index(array), outmesh%cface)
    end associate

    !! XCFACE
    allocate(outmesh%xcface(merge(outmesh%num_cell+1,0,is_IOP)))
    if (is_IOP) then
      do concurrent (j = 0:outmesh%num_cell)
        outmesh%xcface(j+1) = 1 + 4*j
      end do
    end if

    !! FNODE
    allocate(outmesh%fnode(merge(3*outmesh%num_face,0,is_IOP)))
    associate (fnode => inmesh%fnode(:,:inmesh%nface_onP))
      array(1:size(fnode)) => fnode
      call gather(inmesh%node_imap%global_index(array), outmesh%fnode)
    end associate

    !! XFNODE
    allocate(outmesh%xfnode(merge(outmesh%num_face+1,0,is_IOP)))
    if (is_IOP) then
      do concurrent (j = 0:outmesh%num_face)
        outmesh%xfnode(j+1) = 1 + 3*j
      end do
    end if

    call inmesh%get_global_x_array(outmesh%coord)

    block !! Generate the CDIR array.
      integer :: j
      integer, allocatable, target :: cfdir(:,:)
      allocate(cfdir(4,inmesh%ncell_onP))
      do j = 1, inmesh%ncell_onP
        if (inmesh%volume(j) < 0) then
          cfdir(:,j) = [-1, 1, -1, 1]
        else
          cfdir(:,j) = [1, -1, 1, -1]
        end if
      end do
      allocate(outmesh%cfdir, mold=outmesh%cface)
      array(1:size(cfdir)) => cfdir
      call gather(array, outmesh%cfdir)
    end block

    allocate(outmesh%blockid(merge(outmesh%num_cell,0,is_IOP)))
    call gather(inmesh%cblock(:inmesh%ncell_onP), outmesh%blockid)

  end subroutine simpl_to_serial

  !! This auxiliary subroutine adds the node-to-cell connectivity data
  !! structure to the passed (and initialized) serial mesh object.

  subroutine add_node_to_cell(smesh)
    use graph_type
    type(serial_mesh), intent(inout) :: smesh
    integer :: i, j
    type(graph) :: g
    call g%init(smesh%num_node, smesh%num_cell)
    do j = 1, smesh%num_cell
      associate (cnode => smesh%cnode(smesh%xcnode(j):smesh%xcnode(j+1)-1))
        do i = 1, size(cnode)
          call g%add_edge(cnode(i), j)
        end do
      end associate
    end do
    call g%get_adjacency(smesh%xncell, smesh%ncell)
  end subroutine

  !! This auxiliary subroutine converts a serial unstructured mesh stored as
  !! an UNSTR_MESH object to an UNSTR_MESH_DATA object on the IO processor.

  subroutine serial_to_portage(smesh, pmesh)

    use unstr_mesh_type

    type(serial_mesh), intent(in), target :: smesh
    type(unstr_mesh_data), intent(out) :: pmesh

    pmesh%num_cell = smesh%num_cell
    pmesh%num_face = smesh%num_face
    pmesh%num_node = smesh%num_node
    pmesh%xcnode = c_loc(smesh%xcnode)
    pmesh%cnode  = c_loc(smesh%cnode)
    pmesh%xcface = c_loc(smesh%xcface)
    pmesh%cface  = c_loc(smesh%cface)
    pmesh%xfnode = c_loc(smesh%xfnode)
    pmesh%fnode  = c_loc(smesh%fnode)
    pmesh%xncell = c_loc(smesh%xncell)
    pmesh%ncell  = c_loc(smesh%ncell)
    pmesh%coord  = c_loc(smesh%coord)
    pmesh%cfdir  = c_loc(smesh%cfdir)
    pmesh%blockid = c_loc(smesh%blockid)

  end subroutine serial_to_portage

end module portage_mapper_type
