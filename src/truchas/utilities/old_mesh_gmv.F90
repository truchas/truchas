!!
!! OLD_MESH_GMV
!!
!! This is a quick-n-dirty high-level layer over GMV's gmvwrite C library that
!! provides some procedures for writing the original mesh (from legacy_mesh_api)
!! and cell-based fields on that mesh to a GMV-format graphics file.  It is
!! intended for ad hoc debugging purposes only, enabling the developer to
!! visualize intermediate field data from any place within Truchas.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! April 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! These are the available procedures, presented in the order they should be
!! called.  Although the output is performed only on the IO processor, the
!! procedures should be called by all processes in parallel.
!!
!!  CALL GMV_OPEN (FILE) establishes FILE as the file where the graphics data
!!    will be written.  Any previous contents of the file will be overwritten.
!!
!!  CALL GMV_WRITE_MESH () writes the mesh to the graphics file. When executing
!!    in parallel, mesh partition information is also written as the "cellpart"
!!    flag.
!!
!!  CALL GMV_BEGIN_VARIABLES ([time] [, seq]) prepares the graphics file to
!!    receive data for field variables.  If the optional argument TIME is
!!    present its value is used for the gmv problem time.  If the optional
!!    argument SEQ is present, its value is used for the gmv cycle number.
!!
!!  CALL GMV_WRITE_CELL_VAR (MESH, U, NAME) writes the variable data U to
!!    the graphics file.  U should be a distributed, cell-based field, and
!!    NAME is an arbitrary string used to label the variable; only the first
!!    8 characters are significant.  This may be called multiple times.
!!
!!  CALL GMV_END_VARIABLES () signals that no more variable data will be written.
!!
!!  CALL GMV_CLOSE () finalizes the graphics file; nothing more can be written.
!!

module old_mesh_gmv

  use kinds, only: r8
  use gmvwrite_c_binding
  use parallel_communication
  implicit none
  private

  public :: gmv_open, gmv_close
  public :: gmv_write_mesh
  public :: gmv_begin_variables, gmv_end_variables
  public :: gmv_write_cell_var

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

  subroutine gmv_write_mesh

    use legacy_mesh_api, only: ncells_tot, nnodes_tot, ncells
    use legacy_mesh_api, only: mesh, vertex
    use string_utilities, only: i_to_c

    integer,  allocatable :: cnode(:,:), pdata(:)
    real(r8), allocatable :: x(:,:)
    integer :: j, k

    integer, parameter :: prism_vert_map(6) = [1,4,5,2,3,6]

    allocate(x(3,merge(nnodes_tot,0,is_IOP)))
    do k = 1, size(x,dim=1)
      call collate (x(k,:), vertex%coord(k))
    end do
    if (is_IOP) then
      call gmvwrite_node_data_f (size(x,dim=2), x(1,:), x(2,:), x(3,:))
    end if
    deallocate(x)

    allocate(cnode(8,merge(ncells_tot,0,is_IOP)))
    do k = 1, size(cnode,dim=1)
      call collate (cnode(k,:), mesh%ngbr_vrtx_orig(k))
    end do
    if (is_IOP) then
      call gmvwrite_cell_header_f (size(cnode,dim=2))
      do j = 1, size(cnode,dim=2)
        if (cnode(1,j) == cnode(2,j)) then !  tet
          call gmvwrite_cell_type_f ('ptet4', 4, cnode(2:5,j))
        else if (cnode(5,j) == cnode(6,j)) then ! pyramid
          call gmvwrite_cell_type_f ('ppyrmd5', 5, cnode(:5,j))
        else if (cnode(5,j) == cnode(8,j)) then ! wedge
          call gmvwrite_cell_type_f ('pprism6', 6, cnode(prism_vert_map,j))
        else  ! hex
          call gmvwrite_cell_type_f ('phex8', 8, cnode(:,j))
        end if
      end do
    end if
    deallocate(cnode)

    !! If in parallel write partitioning info as flags.
    if (nPE > 1) then
      if (is_IOP) call gmvwrite_flag_header_f ()
      !! Cell partitioning info ...
      allocate(pdata(merge(ncells_tot,0,is_IOP)))
      call collate (pdata, spread(this_PE, dim=1, ncopies=ncells))
      if (is_IOP) then
        call gmvwrite_flag_name_f ('cellpart', nPE, CELLDATA)
        do j = 1, nPE
          call gmvwrite_flag_subname_f('P'//i_to_c(j))
        end do
        call gmvwrite_flag_data_f (CELLDATA, pdata)
      end if
      deallocate(pdata)
      if (is_IOP) call gmvwrite_flag_endflag_f ()
    end if

  end subroutine gmv_write_mesh

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

  subroutine gmv_write_cell_var (u, name)

    use legacy_mesh_api, only: ncells_tot

    real(r8), intent(in) :: u(:)
    character(*), intent(in) :: name

    real(r8), allocatable :: u_global(:)

    allocate(u_global(merge(ncells_tot,0,is_IOP)))
    call collate (u_global, u)
    if (is_IOP) call gmvwrite_variable_name_data_f (CELLDATA, name, u_global)
    deallocate(u_global)

  end subroutine gmv_write_cell_var

end module old_mesh_gmv
