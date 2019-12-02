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
  use parallel_communication, only: is_IOP, collate, distribute, global_sum
  use data_mapper_class
  implicit none
  private

  !! This holds an internal serialized copy of a parallel Truchas mesh
  type :: serial_mesh
    integer(c_int) :: ncell, nface, nnode
    integer(c_int), allocatable :: xcnode(:), cnode(:)
    integer(c_int), allocatable :: xcface(:), cface(:)
    integer(c_int), allocatable :: xfnode(:), fnode(:)
    real(c_double), allocatable :: coord(:,:)
    integer(c_int), allocatable :: cfdir(:)
    integer(c_int), allocatable :: blockid(:)
  end type serial_mesh

  !! Mesh data needed by Portage
  type, bind(c) :: unstr_mesh_data
    integer(c_int) :: ncell, nface, nnode
    type(c_ptr) :: xcnode, cnode   ! int arrays
    type(c_ptr) :: xcface, cface   ! int arrays
    type(c_ptr) :: xfnode, fnode   ! int arrays
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

    class(portage_mapper), intent(in) :: this
    real(r8), contiguous, intent(in)  :: src(:)
    real(r8), contiguous, intent(out) :: dest(:)
    real(r8), intent(in)  :: defval
    integer,  intent(in)  :: map_type
    logical,  intent(in), optional :: pullback

    integer :: method
    logical :: reverse_order
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

    allocate(col_src(merge(global_sum(size(src)),0,is_IOP)))
    allocate(col_dest(merge(global_sum(size(dest)),0,is_IOP)))

    call collate(col_src, src)

    if (is_IOP) then
      col_dest = defval
      if (reverse_order) then
        call portage_map_field(this%mapper_pullback, &
            size(col_src), col_src, size(col_dest), col_dest, method)
      else
        call portage_map_field(this%mapper, &
            size(col_src), col_src, size(col_dest), col_dest, method)
      end if
    end if

    call distribute(dest, col_dest)

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

  end subroutine base_to_serial

  !! This auxiliary subroutine creates on the IOP a SERIAL_MESH copy of the
  !! given parallel UNSTR_MESH object, which holds just those data components
  !! needed by Portage.

  subroutine unstr_to_serial(inmesh, outmesh)

    use unstr_mesh_type

    type(unstr_mesh), intent(in) :: inmesh
    type(serial_mesh), intent(out) :: outmesh

    if (is_IOP) then
      outmesh%ncell = inmesh%cell_ip%global_size()
      outmesh%nface = inmesh%face_ip%global_size()
      outmesh%nnode = inmesh%node_ip%global_size()
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
      call collate(outmesh%cfdir, cfdir)
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
      allocate(outmesh%blockid(merge(outmesh%ncell,0,is_IOP)))
      call collate(outmesh%blockid, blockid)
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
      outmesh%ncell = inmesh%cell_ip%global_size()
      outmesh%nface = inmesh%face_ip%global_size()
      outmesh%nnode = inmesh%node_ip%global_size()
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
      allocate(outmesh%cnode(merge(4*outmesh%ncell,0,is_IOP)))
      call collate(outmesh%cnode, inmesh%node_ip%global_index(array))
    end block

    !! XCNODE
    allocate(outmesh%xcnode(merge(outmesh%ncell+1,0,is_IOP)))
    if (is_IOP) then
      do concurrent (j = 0:outmesh%ncell)
        outmesh%xcnode(j+1) = 1 + 4*j
      end do
    end if

    !! CFACE -- NB: faces no longer opposite corresponding vertex
    allocate(outmesh%cface(merge(4*outmesh%ncell,0,is_IOP)))
    associate (cface => inmesh%cface(:,:inmesh%ncell_onP))
      array(1:size(cface)) => cface
      call collate(outmesh%cface, inmesh%face_ip%global_index(array))
    end associate

    !! XCFACE
    allocate(outmesh%xcface(merge(outmesh%ncell+1,0,is_IOP)))
    if (is_IOP) then
      do concurrent (j = 0:outmesh%ncell)
        outmesh%xcface(j+1) = 1 + 4*j
      end do
    end if

    !! FNODE
    allocate(outmesh%fnode(merge(3*outmesh%nface,0,is_IOP)))
    associate (fnode => inmesh%fnode(:,:inmesh%nface_onP))
      array(1:size(fnode)) => fnode
      call collate(outmesh%fnode, inmesh%node_ip%global_index(array))
    end associate

    !! XFNODE
    allocate(outmesh%xfnode(merge(outmesh%nface+1,0,is_IOP)))
    if (is_IOP) then
      do concurrent (j = 0:outmesh%nface)
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
      call collate(outmesh%cfdir, array)
    end block

    allocate(outmesh%blockid(merge(outmesh%ncell,0,is_IOP)))
    call collate(outmesh%blockid, inmesh%cblock(:inmesh%ncell_onP))

  end subroutine simpl_to_serial

  !! This auxiliary subroutine converts a serial unstructured mesh stored as
  !! an UNSTR_MESH object to an UNSTR_MESH_DATA object on the IO processor.

  subroutine serial_to_portage(smesh, pmesh)

    use unstr_mesh_type

    type(serial_mesh), intent(in), target :: smesh
    type(unstr_mesh_data), intent(out) :: pmesh

    pmesh%ncell = smesh%ncell
    pmesh%nface = smesh%nface
    pmesh%nnode = smesh%nnode
    pmesh%xcnode = c_loc(smesh%xcnode)
    pmesh%cnode  = c_loc(smesh%cnode)
    pmesh%xcface = c_loc(smesh%xcface)
    pmesh%cface  = c_loc(smesh%cface)
    pmesh%xfnode = c_loc(smesh%xfnode)
    pmesh%fnode  = c_loc(smesh%fnode)
    pmesh%coord  = c_loc(smesh%coord)
    pmesh%cfdir  = c_loc(smesh%cfdir)
    pmesh%blockid = c_loc(smesh%blockid)

  end subroutine serial_to_portage

end module portage_mapper_type
