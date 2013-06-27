#include "f90_assert.fpp"

module ER_solver_gmv

  use kinds, only: r8
  use parallel_communication
  use fgmvwrite
  implicit none
  private

  public :: ERS_gmv_open, ERS_gmv_close, ERS_gmv_write_enclosure
  public :: ERS_gmv_begin_variables, ERS_gmv_end_variables, ERS_gmv_write_var

contains

  subroutine ERS_gmv_open (file)
    character(len=*) :: file
    !! 4-byte integer data and 8-byte real data.
    if (is_IOP) call gmvwrite_openfile_ir_ascii (file, 4, 8)
    !! GMV has a bug with binary nodeids data.
    !if (is_IOP) call gmvwrite_openfile_ir (file, 4, 8)
  end subroutine ERS_gmv_open

  subroutine ERS_gmv_close ()
    if (is_IOP) call gmvwrite_closefile ()
  end subroutine ERS_gmv_close

  subroutine ERS_gmv_write_enclosure (this)

    use ER_solver, only: solver
    use index_partitioning

    type(solver), intent(in) :: this

    integer :: j, dimen, nnode, nface, offset
    integer, pointer :: map(:), fsize(:), fnode(:), list(:)
    real(r8), pointer :: x(:), y(:), z(:)
    character(len=8) :: name

    dimen = size(this%encl%coord,dim=1)
    nnode = global_size(this%encl%node_ip)
    nface = global_size(this%encl%face_ip)

    !! Write the node coordinate data.
    call allocate_collated_array (x, nnode)
    call allocate_collated_array (y, nnode)
    call allocate_collated_array (z, nnode)
    call collate (x, this%encl%coord(1,:this%encl%nnode_onP))
    call collate (y, this%encl%coord(2,:this%encl%nnode_onP))
    call collate (z, this%encl%coord(3,:this%encl%nnode_onP))
    if (is_IOP) call gmvwrite_node_data (nnode, x, y, z)
    deallocate(x, y, z)

    !! Write the cell data.
    call allocate_collated_array (fsize, nface)
    call collate (fsize, this%encl%xface(2:this%encl%nface+1)-this%encl%xface(1:this%encl%nface))
    call allocate_collated_array (fnode, sum(fsize))
    call collate (fnode, global_index(this%encl%node_ip, this%encl%fnode))
    if (is_IOP) then
      call gmvwrite_cell_header (nface)
      offset = 0
      do j = 1, nface
        list => fnode(offset+1:offset+fsize(j))
        select case (fsize(j))
        case (3)
          call gmvwrite_cell_type ('tri', 3, list)
        case (4)
          call gmvwrite_cell_type ('quad', 4, list)
        case default
          INSIST(.false.)
        end select
        offset = offset + fsize(j)
      end do
    end if
    deallocate(fsize, fnode)

    !! Write the node map as the nodeids -- GMV uses these for display.
    call allocate_collated_array (map, nnode)
    call collate (map, this%encl%node_map(:this%encl%nnode_onP))
    if (is_IOP) call gmvwrite_nodeids (map)
    deallocate (map)

    !! Write the face map as the cellids -- GMV uses these for display.
    call allocate_collated_array (map, nface)
    call collate (map, this%encl%face_map)
    if (is_IOP) call gmvwrite_cellids (map)
    deallocate (map)

    !! Write the face block IDs as the cell material.
    call allocate_collated_array (map, nface)
    call collate (map, this%encl%face_block)
    if (is_IOP) then
      call gmvwrite_material_header (size(this%encl%face_block_id), CELLDATA)
      do j = 1, size(this%encl%face_block_id)
        write(name,'(a,i0)') 'Block', this%encl%face_block_id(j)
        call gmvwrite_material_name (name)
      end do
      call gmvwrite_material_ids (map, CELLDATA)
    end if
    deallocate(map)

    if (nPE > 1) then
      !! Write the face partitioning info.
      call allocate_collated_array (map, nface)
      call collate (map, spread(this_PE, dim=1, ncopies=this%nface))
      if (is_IOP) then
        call gmvwrite_flag_header ()
        call gmvwrite_flag_name ('par-part', nPE, CELLDATA)
        do j = 1, nPE
          write(name,'(a,i0)') 'P', j
          call gmvwrite_flag_subname (name)
        end do
        call gmvwrite_flag_data (CELLDATA, map)
        call gmvwrite_flag_endflag ()
      end if
      deallocate(map)
    end if

  end subroutine ERS_gmv_write_enclosure

  subroutine ERS_gmv_begin_variables (time, seq)
    real(r8), intent(in), optional :: time
    integer,  intent(in), optional :: seq
    if (is_IOP) then
      if (present(time)) call gmvwrite_probtime (time)
      if (present(seq))  call gmvwrite_cycleno (seq)
      call gmvwrite_variable_header()
    end if
  end subroutine ERS_gmv_begin_variables

  subroutine ERS_gmv_end_variables ()
    if (is_IOP) call gmvwrite_variable_endvars()
  end subroutine ERS_gmv_end_variables

  subroutine ERS_gmv_write_var (this, u, name)

    use ER_solver, only: solver
    use index_partitioning

    type(solver), intent(in) :: this
    real(r8), intent(in) :: u(:)
    character(len=*), intent(in) :: name

    real(r8), pointer :: u_global(:)

    ASSERT(size(u) == this%nface)

    call allocate_collated_array (u_global, global_size(this%encl%face_ip))
    call collate (u_global, u)
    if (is_IOP) call gmvwrite_variable_name_data (CELLDATA, name, u_global)
    deallocate(u_global)

  end subroutine ERS_gmv_write_var

end module ER_solver_gmv
