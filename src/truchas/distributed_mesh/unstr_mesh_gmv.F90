!!
!! UNSTR_MESH_GMV
!!
!! Neil N. Carlson <nnc@lanl.gov> 3/30/2006
!! Adapted for mixed element meshes, July 2015
!!
!! This is a quick-n-dirty high-level layer over GMV's gmvwrite C library
!! that provides some procedures for writing a distributed mesh and cell-
!! based fields over that mesh to a GMV-format graphics file.  Although
!! the output is performed only on the IO processor, the routines can be
!! called in parallel, and when distributed data is passed, must be called
!! in parallel.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! These are the available procedures, presented in the order they should
!! be called:
!!
!! CALL GMV_OPEN (FILE) establishes FILE as the file where the graphics data
!!    will be written.  Any previous contents of the file will be overwritten.
!!
!! CALL GMV_WRITE_UNSTR_MESH (MESH) writes the TYPE(UNSTR_MESH) mesh to the
!!    graphics file.  The cell set information is written as the material
!!    component (see Note 1).  When executing in parallel, mesh partition
!!    information is also written as flag components: cellpart gives the
!!    cell partitioning nodepart the node partitioning.
!!
!! CALL GMV_BEGIN_VARIABLES ([time] [, seq]) prepares the graphics file to
!!    receive data for field variables.  If the optional argument TIME is
!!    present its value is used for the gmv problem time.  If the optional
!!    argument SEQ is present, its value is used for the gmv cycle number.
!!
!! CALL GMV_WRITE_DIST_CELL_VAR (MESH, U, NAME) writes the variable data U to
!!    the graphics file.  U should be a distributed, cell-based field over the
!!    TYPE(UNSTR_MESH) mesh (same as in the earlier call), and NAME is an
!!    arbitrary string used to label the variable; only the first 8 characters
!!    are significant.  This may be called multiple times.
!!
!! CALL GMV_END_VARIABLES () signals that no more variable data will be written.
!!
!! CALL GMV_CLOSE () finalizes the graphics file; nothing more can be written.
!!
!! IMPLEMENTATION NOTES
!!
!! 1. There is no good match in GMV for the mesh cell sets.  A cell may belong
!!    to multiple sets, but GMV requires that each cell belong to exactly one
!!    material; that is, the materials form a disjoint partition of the mesh.
!!    The situation is the same for GMV flags.  In principle a cell set could
!!    be mapped to a flag with two values representing "in" and "not in", but
!!    the number of flags is limited to 10 so this is not a viable solution.
!!    For now we identify materials with cell sets and associate a cell with
!!    the first cell set to which it belongs. This does exactly the right thing
!!    for the current manner in which the mesh is initialized: the disjoint
!!    element blocks are mapped to cell sets, and cell sets are defined in no
!!    other way.
!!

#include "f90_assert.fpp"

module unstr_mesh_gmv

  use kinds, only: r8
  use gmvwrite_c_binding
  use parallel_communication
  implicit none
  private

  public :: gmv_open, gmv_close, gmv_write_unstr_mesh
  public :: gmv_begin_variables, gmv_end_variables
  public :: gmv_write_dist_cell_var, gmv_write_dist_node_var

contains

  subroutine gmv_open (file)
    character(len=*) :: file
    !! 4-byte integer data and 8-byte real data.
    !if (is_IOP) call gmvwrite_openfile_ir (file, 4, 8) ! GMV bug with node ids
    if (is_IOP) call gmvwrite_openfile_ir_ascii_f (file, 4, 8)
  end subroutine gmv_open

  subroutine gmv_close ()
    if (is_IOP) call gmvwrite_closefile_f ()
  end subroutine gmv_close

  subroutine gmv_write_unstr_mesh (mesh)

    use unstr_mesh_type
    use index_partitioning
    use string_utilities, only: i_to_c

    type(unstr_mesh), intent(in) :: mesh

    integer :: j, ncell, nnode
    integer, allocatable :: xcnode(:), cnode(:), map(:), iflag(:)
    integer(kind(mesh%cell_set_mask)), allocatable :: cell_set_mask(:)
    real(r8), allocatable :: x(:,:)

    !! Write the node coordinates
    call mesh%get_global_x_array (x)
    nnode = size(x,dim=2)
    if (is_IOP) call gmvwrite_node_data_f (nnode, x(1,:), x(2,:), x(3,:))
    deallocate(x)

    !! Write the cell connectivity
    call mesh%get_global_cnode_array (xcnode, cnode)
    ncell = size(xcnode) - 1
    if (is_IOP) then
      call gmvwrite_cell_header_f (ncell)
      do j = 1, ncell
        associate (connect => cnode(xcnode(j):xcnode(j+1)-1))
          select case (size(connect))
          case (4)
            call gmvwrite_cell_type_f ('ptet4', size(connect), connect)
          case (5)
            call gmvwrite_cell_type_f ('ppyrmd5', size(connect), connect)
          case (6)
            call gmvwrite_cell_type_f ('pprism6', size(connect), connect)
          case (8)
            call gmvwrite_cell_type_f ('phex8', size(connect), connect)
          case default
          print *, ncell, j, '--', size(connect), '--', connect
            INSIST(.false.)
          end select
        end associate
      end do
    end if
    deallocate(xcnode, cnode)

    !! Write external mesh node numbers as the nodeids -- GMV uses these for display.
    allocate(map(nnode))
    call collate (map, mesh%xnode(:mesh%nnode_onP))
    if (is_IOP) call gmvwrite_nodeids_f (map)
    deallocate (map)

    !! Write external mesh cell numbers as the cellids -- GMV uses these for display.
    allocate(map(ncell))
    call collate (map, mesh%xcell(:mesh%ncell_onP))
    if (is_IOP) call gmvwrite_cellids_f (map)
    deallocate (map)

    !! Write cell materials.  NB: See Note 1
    allocate(cell_set_mask(ncell))
    call collate (cell_set_mask, mesh%cell_set_mask(:mesh%ncell_onP))
    if (is_IOP) then
      call gmvwrite_material_header_f (size(mesh%cell_set_id), CELLDATA)
      do j = 1, size(mesh%cell_set_id)
        call gmvwrite_material_name_f ('block'//i_to_c(mesh%cell_set_id(j)))
      end do
      cell_set_mask = ibclr(cell_set_mask,pos=0)  ! bit 0 may have a special use
      call gmvwrite_material_ids_f (trailz(cell_set_mask), CELLDATA)
    end if
    deallocate(cell_set_mask)

    !! If in parallel write partitioning info as flags.
    if (nPE > 1) then

      if (is_IOP) call gmvwrite_flag_header_f ()

      !! Cell partitioning info ...
      allocate(map(ncell))
      call collate (map, spread(this_PE, dim=1, ncopies=mesh%cell_ip%onP_size()))
      if (is_IOP) then
        call gmvwrite_flag_name_f ('cellpart', nPE, CELLDATA)
        do j = 1, nPE
          call gmvwrite_flag_subname_f('P'//i_to_c(j))
        end do
        call gmvwrite_flag_data_f (CELLDATA, map)
      end if

      !! On/off-process cells for each process
      allocate(iflag(mesh%ncell))
      do j = 1, nPE
        if (this_PE /= j) then
          iflag(:mesh%ncell_onP) = 1
          iflag(mesh%ncell_onP+1:) = 0
        else
          iflag = 2
        end if
        call scatter_boundary_sum(mesh%cell_ip, iflag)
        call collate(map, iflag(:mesh%ncell_onP))
        if (is_IOP) then
          call gmvwrite_flag_name_f('P'//i_to_c(j)//'cells', 3, CELLDATA)
          call gmvwrite_flag_subname_f('other') ! for iflag==1
          call gmvwrite_flag_subname_f('owned') ! for iflag==2
          call gmvwrite_flag_subname_f('ghost') ! for iflag==3
          call gmvwrite_flag_data_f(CELLDATA, map)
        end if
      end do

      !! Node partitioning info ...
      deallocate(map)
      allocate(map(nnode))
      call collate (map, spread(this_PE, dim=1, ncopies=mesh%node_ip%onP_size()))
      if (is_IOP) then
        call gmvwrite_flag_name_f ('nodepart', nPE, NODEDATA)
        do j = 1, nPE
          call gmvwrite_flag_subname_f('P'//i_to_c(j))
        end do
        call gmvwrite_flag_data_f (NODEDATA, map)
      end if
      deallocate(map)

      if (is_IOP) call gmvwrite_flag_endflag_f ()

    end if

  end subroutine gmv_write_unstr_mesh

  subroutine gmv_begin_variables (time, seq)
    real(r8), intent(in), optional :: time
    integer, intent(in), optional :: seq
    if (is_IOP) then
      if (present(time)) call gmvwrite_probtime_f (time)
      if (present(seq))  call gmvwrite_cycleno_f (seq)
      call gmvwrite_variable_header_f()
    end if
  end subroutine gmv_begin_variables

  subroutine gmv_end_variables ()
    if (is_IOP) call gmvwrite_variable_endvars_f()
  end subroutine gmv_end_variables

  subroutine gmv_write_dist_cell_var (mesh, u, name)

    use unstr_mesh_type
    use index_partitioning

    type(unstr_mesh), intent(in) :: mesh
    real(r8), intent(in) :: u(:)
    character(*), intent(in) :: name

    real(r8), pointer :: u_global(:)

    ASSERT(mesh%cell_ip%defined())
    ASSERT(size(u) == mesh%cell_ip%onP_size())

    allocate(u_global(merge(mesh%cell_ip%global_size(),0,is_IOP)))
    call collate (u_global, u)
    if (is_IOP) call gmvwrite_variable_name_data_f (CELLDATA, name, u_global)
    deallocate(u_global)

  end subroutine gmv_write_dist_cell_var

  subroutine gmv_write_dist_node_var (mesh, u, name)

    use unstr_mesh_type
    use index_partitioning

    type(unstr_mesh), intent(in) :: mesh
    real(r8), intent(in) :: u(:)
    character(*), intent(in) :: name

    real(r8), pointer :: u_global(:)

    ASSERT(mesh%node_ip%defined())
    ASSERT(size(u) == mesh%node_ip%onP_size())

    allocate(u_global(merge(mesh%node_ip%global_size(),0,is_IOP)))
    call collate (u_global, u)
    if (is_IOP) call gmvwrite_variable_name_data_f (NODEDATA, name, u_global)
    deallocate(u_global)

  end subroutine gmv_write_dist_node_var

end module unstr_mesh_gmv
