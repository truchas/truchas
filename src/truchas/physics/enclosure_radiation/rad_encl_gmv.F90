!!
!! RAD_ENCL_GMV
!!
!! This is a quick-n-dirty high level layer over GMV's gmvwrite C library
!! that provides some procedures for writing a radiation enclosure surface
!! mesh and face-based fields over that mesh to a GMV-format graphics file.
!! It is intended for ad hoc use in standalone test programs and debugging
!! situations.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Revised May 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! USAGE
!!
!!  While output is performed only by the IO processor, some subroutines
!!  assume distributed data and must be called in parallel.  For simplicity
!!  all subroutines may be called in parallel, and that is the suggested
!!  practice.  These are the available subroutines, presented in the order
!!  they should be called.  Note that the modules RAD_SOLVER_GMV and
!!  RAD_PROBLEM_GMV extend the generics defined here to additional objects
!!  containing a radiation enclosure.
!!
!!  CALL GMV_OPEN (FILE) establishes FILE as the file where the graphics data
!!    will be written.  Any previous contents of the file will be overwritten.
!!
!!  CALL GMV_WRITE_ENCLOSURE (THIS) writes the distributed radiation enclosure
!!    surface mesh described by the RAD_ENCL type object THIS to the graphics
!!    file.  When executing in parallel, the partitioning of faces is written
!!    as flag data.
!!
!!  CALL GMV_BEGIN_VARIABLES ([TIME] [,SEQ]) prepares the graphics file to
!!    receive data for field variables.  If the optional argument TIME is
!!    present, its value is used for the gmv problem time.  If the optional
!!    argument SEQ is present, its value is used for the gmv cycle number.
!!
!!  CALL GMV_WRITE_VARIABLE (THIS, VAR, NAME) writes the distributed variable
!!    data VAR to the graphics file.  VAR is a face field defined over the
!!    radiation enclosure surface mesh described by the RAD_ENCL type object
!!    THIS.  NAME is an arbitrary string used to label the variable; only the
!!    first 8 characters are significant.  This may be called multiple times.
!!
!!  CALL GMV_END_VARIBLES () signals that no more variable data will be written.
!!
!!  CALL GMV_CLOSE () finalizes and closes the graphics file; nothing more can
!!    be written.
!!

#include "f90_assert.fpp"

module rad_encl_gmv

  use kinds, only: r8
  use gmvwrite_c_binding
  use rad_encl_type
  use index_partitioning
  use parallel_communication
  implicit none
  private

  public :: gmv_open, gmv_close, gmv_write_enclosure
  public :: gmv_begin_variables, gmv_end_variables, gmv_write_variable

  interface gmv_write_enclosure
    procedure gmv_write_encl
  end interface

  interface gmv_write_variable
    procedure gmv_write_var
  end interface

contains

  subroutine gmv_open (file)
    character(len=*) :: file
    !! 4-byte integer data and 8-byte real data.
    if (is_IOP) call gmvwrite_openfile_ir_ascii_f (file, 4, 8)
    !! GMV has a bug with binary nodeids data.
    !if (is_IOP) call gmvwrite_openfile_ir_f (file, 4, 8)
  end subroutine gmv_open

  subroutine gmv_close ()
    if (is_IOP) call gmvwrite_closefile_f ()
  end subroutine gmv_close

  subroutine gmv_begin_variables (time, seq)
    real(r8), intent(in), optional :: time
    integer,  intent(in), optional :: seq
    if (is_IOP) then
      if (present(time)) call gmvwrite_probtime_f (time)
      if (present(seq))  call gmvwrite_cycleno_f (seq)
      call gmvwrite_variable_header_f
    end if
  end subroutine gmv_begin_variables

  subroutine gmv_end_variables ()
    if (is_IOP) call gmvwrite_variable_endvars_f
  end subroutine gmv_end_variables

  subroutine gmv_write_encl (this)

    type(rad_encl), intent(in) :: this

    integer :: j, dimen, nnode, nface, offset
    integer, pointer :: map(:), fsize(:), fnode(:), list(:)
    real(r8), pointer :: x(:), y(:), z(:)
    character(8) :: name

    dimen = size(this%coord,dim=1)
    nnode = this%node_ip%global_size()
    nface = this%face_ip%global_size()

    !! Write the node coordinate data.
    call allocate_collated_array (x, nnode)
    call allocate_collated_array (y, nnode)
    call allocate_collated_array (z, nnode)
    call collate (x, this%coord(1,:this%nnode_onP))
    call collate (y, this%coord(2,:this%nnode_onP))
    call collate (z, this%coord(3,:this%nnode_onP))
    if (is_IOP) call gmvwrite_node_data_f (nnode, x, y, z)
    deallocate(x, y, z)

    !! Write the cell data.
    call allocate_collated_array (fsize, nface)
    call collate (fsize, this%xface(2:this%nface+1)-this%xface(1:this%nface))
    call allocate_collated_array (fnode, sum(fsize))
    call collate (fnode, this%node_ip%global_index(this%fnode))
    if (is_IOP) then
      call gmvwrite_cell_header_f (nface)
      offset = 0
      do j = 1, nface
        list => fnode(offset+1:offset+fsize(j))
        select case (fsize(j))
        case (3)
          call gmvwrite_cell_type_f ('tri', 3, list)
        case (4)
          call gmvwrite_cell_type_f ('quad', 4, list)
        case default
          INSIST(.false.)
        end select
        offset = offset + fsize(j)
      end do
    end if
    deallocate(fsize, fnode)

    !! Write the node map as the nodeids -- GMV uses these for display.
    call allocate_collated_array (map, nnode)
    call collate (map, this%node_map(:this%nnode_onP))
    if (is_IOP) call gmvwrite_nodeids_f (map)
    deallocate (map)

    !! Write the face map as the cellids -- GMV uses these for display.
    call allocate_collated_array (map, nface)
    call collate (map, this%face_map)
    if (is_IOP) call gmvwrite_cellids_f (map)
    deallocate (map)

    !! Write the face block IDs as the cell material.
    call allocate_collated_array (map, nface)
    call collate (map, this%face_block)
    if (is_IOP) then
      call gmvwrite_material_header_f (size(this%face_block_id), CELLDATA)
      do j = 1, size(this%face_block_id)
        write(name,'(a,i0)') 'Block', this%face_block_id(j)
        call gmvwrite_material_name_f (name)
      end do
      call gmvwrite_material_ids_f (map, CELLDATA)
    end if
    deallocate(map)

    if (nPE > 1) then
      !! Write the face partitioning info.
      call allocate_collated_array (map, nface)
      call collate (map, spread(this_PE, dim=1, ncopies=this%nface))
      if (is_IOP) then
        call gmvwrite_flag_header_f ()
        call gmvwrite_flag_name_f ('par-part', nPE, CELLDATA)
        do j = 1, nPE
          write(name,'(a,i0)') 'P', j
          call gmvwrite_flag_subname_f (name)
        end do
        call gmvwrite_flag_data_f (CELLDATA, map)
        call gmvwrite_flag_endflag_f ()
      end if
      deallocate(map)
    end if

  end subroutine gmv_write_encl

  subroutine gmv_write_var (this, u, name)
    type(rad_encl), intent(in) :: this
    real(r8), intent(in) :: u(:)
    character(*), intent(in) :: name
    real(r8), allocatable :: u_global(:)
    ASSERT(size(u) >= this%nface_onP)
    allocate(u_global(merge(this%face_ip%global_size(),0,is_IOP)))
    call collate (u_global, u(:this%nface_onP))
    if (is_IOP) call gmvwrite_variable_name_data_f (CELLDATA, name, u_global)
  end subroutine gmv_write_var

end module rad_encl_gmv
