!!
!! SIMPL_MESH_GMV
!!
!! Neil N. Carlson <nnc@lanl.gov> 3/30/2006
!! Last revised 11 Aug 2006.
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
!! CALL GMV_WRITE_SIMPL_MESH (MESH) writes the distributed mesh MESH to the
!!    graphics file.  The cell block ID information is written as the material
!!    component.  When executing in parallel, mesh partition information is
!!    also written as flag components: cellpart gives the cell partitioning
!!    nodepart the node partitioning.
!!
!! CALL GMV_BEGIN_VARIABLES ([time] [, seq]) prepares the graphics file to
!!    receive data for field variables.  If the optional argument TIME is
!!    present its value is used for the gmv problem time.  If the optional
!!    argument SEQ is present, its value is used for the gmv cycle number.
!!
!! CALL GMV_WRITE_DIST_CELL_VAR (MESH, U, NAME) writes the variable data U to
!!    the graphics file.  U should be a distributed, cell-based field over the
!!    distributed mesh MESH (same as in the earlier call), and NAME is an
!!    arbitrary string used to label the variable; only the first 8 characters
!!    are significant.  This may be called multiple times.
!!
!! CALL GMV_END_VARIABLES () signals that no more variable data will be written.
!!
!! CALL GMV_CLOSE () finalizes the graphics file; nothing more can be written.
!!

#include "f90_assert.fpp"

module simpl_mesh_gmv

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use gmvwrite_c_binding
  use parallel_communication
  implicit none
  private

  public :: gmv_open, gmv_close, gmv_write_simpl_mesh
  public :: gmv_begin_variables, gmv_end_variables, gmv_write_dist_cell_var

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

  subroutine gmv_write_simpl_mesh (mesh)

    use simpl_mesh_type
    use string_utilities, only: i_to_c

    type(simpl_mesh), intent(in) :: mesh

    integer :: j
    integer, allocatable :: cnode(:,:), cblock(:)
    integer, pointer :: pdata(:), map(:)
    real(kind=r8), allocatable :: x(:,:)
    character(len=3) :: cell_type

    call mesh%get_global_cnode_array (cnode)
    call mesh%get_global_x_array (x)
    call mesh%get_global_cblock_array (cblock)

    if (is_IOP) then

      !! Infer cell type from the CNODE array (MESH ought to have a type component!)
      select case (size(cnode,1))
      case (4)
        cell_type = 'tet'
      case (8)
        cell_type = 'hex'
      case default
        INSIST( .false. )
      end select

      call gmvwrite_node_data_f (size(x,dim=2), x(1,:), x(2,:), x(3,:))
      call gmvwrite_cell_header_f (size(cnode,dim=2))
      do j = 1, size(cnode,dim=2)
        call gmvwrite_cell_type_f (cell_type, size(cnode,dim=1), cnode(:,j))
      end do

    end if

    !! Write external mesh node numbers as the nodeids -- GMV uses these for display.
    allocate(map(merge(size(x,dim=2),0,is_iop)))
    call gather (mesh%xnode(:mesh%nnode_onP), map)
    if (is_IOP) call gmvwrite_nodeids_f (map)
    deallocate (map)

    !! Write external mesh cell numbers as the cellids -- GMV uses these for display.
    allocate(map(merge(size(cnode,dim=2),0,is_iop)))
    call gather (mesh%xcell(:mesh%ncell_onP), map)
    if (is_IOP) call gmvwrite_cellids_f (map)
    deallocate (map)

    if (is_IOP) then
    
      !! Write the cell block IDs as the cell material.
      call gmvwrite_material_header_f (size(mesh%block_id), CELLDATA)
      do j = 1, size(mesh%block_id)
        call gmvwrite_material_name_f ('block'//i_to_c(mesh%block_id(j)))
      end do
      call gmvwrite_material_ids_f (cblock, CELLDATA)

    end if

    deallocate(cnode, x, cblock)

    !! If in parallel write partitioning info as flags.
    if (nPE > 1) then

      if (is_IOP) call gmvwrite_flag_header_f ()

      !! Cell partitioning info ...
      allocate(pdata(merge(mesh%cell_imap%global_size,0,is_iop)))
      call gather (spread(this_PE, dim=1, ncopies=mesh%cell_imap%onp_size), pdata)
      if (is_IOP) then
        call gmvwrite_flag_name_f ('cellpart', nPE, CELLDATA)
        do j = 1, nPE
          call gmvwrite_flag_subname_f('P'//i_to_c(j))
        end do
        call gmvwrite_flag_data_f (CELLDATA, pdata)
      end if
      deallocate(pdata)

      !! Node partitioning info ...
      allocate(pdata(merge(mesh%node_imap%global_size,0,is_iop)))
      call gather (spread(this_PE, dim=1, ncopies=mesh%node_imap%onp_size), pdata)
      if (is_IOP) then
        call gmvwrite_flag_name_f ('nodepart', nPE, NODEDATA)
        do j = 1, nPE
          call gmvwrite_flag_subname_f('P'//i_to_c(j))
        end do
        call gmvwrite_flag_data_f (NODEDATA, pdata)
      end if
      deallocate(pdata)

      if (is_IOP) call gmvwrite_flag_endflag_f ()

    end if

  end subroutine gmv_write_simpl_mesh

  subroutine gmv_begin_variables (time, seq)
    real(kind=r8), intent(in), optional :: time
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

    use simpl_mesh_type

    type(simpl_mesh), intent(in) :: mesh
    real(kind=r8),   intent(in) :: u(:)
    character(len=*), intent(in) :: name

    real(kind=r8), pointer :: u_global(:)

    ASSERT( mesh%cell_imap%defined() )
    ASSERT( size(u) == mesh%cell_imap%onp_size )

    allocate(u_global(merge(mesh%cell_imap%global_size,0,is_iop)))
    call gather (u, u_global)
    if (is_IOP) call gmvwrite_variable_name_data_f (CELLDATA, name, u_global)
    deallocate(u_global)

  end subroutine gmv_write_dist_cell_var

end module simpl_mesh_gmv
