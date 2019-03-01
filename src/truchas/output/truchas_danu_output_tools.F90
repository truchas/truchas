!!
!! TRUCHAS_DANU_OUTPUT_TOOLS
!!
!! Neil N. Carlson <nnc@lanl.gov> Apr 2012
!!
!! Prototype module for developing new HDF5 output using the Danu package.
!! This is organized very much like Tbrook_Utilities so that we can more
!! easily drop this in as a replacement -- just the first step in reworking the
!! output.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module truchas_danu_output_tools

  use kinds, only: r8
  use truchas_h5_outfile, only: th5_seq_group
  use parallel_communication
  use,intrinsic :: iso_c_binding, only: c_ptr
  implicit none
  private
  
  public :: write_seq_cell_field, write_seq_node_field
  
  interface write_seq_cell_field
    module procedure write_seq_cell_field_r0, write_seq_cell_field_r1
  end interface
  
  interface write_seq_node_field
    module procedure write_seq_node_field_r0, write_seq_node_field_r1
  end interface

contains
  
  subroutine write_seq_cell_field_r0 (seq, ldata, name, for_viz, viz_name)
  
    use legacy_mesh_api, only: ncells, ncells_tot

    type(th5_seq_group), intent(in) :: seq
    real(r8), intent(in) :: ldata(:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name
    
    INSIST(size(ldata) == ncells)
    
    call seq%write_dist_array(name, ncells_tot, ldata)
    
    if (for_viz) then
      call seq%write_dataset_attr(name, 'FIELDTYPE', 'CELL')
      if (present(viz_name)) then
        call seq%write_dataset_attr(name, 'FIELDNAME', viz_name)
      else
        call seq%write_dataset_attr(name, 'FIELDNAME', name)
      end if
    end if
    
  end subroutine write_seq_cell_field_r0
  
  subroutine write_seq_cell_field_r1 (seq, ldata, name, for_viz, viz_name)
  
    use legacy_mesh_api, only: ncells, ncells_tot
    use string_utilities, only: i_to_c

    type(th5_seq_group), intent(in) :: seq
    real(r8), intent(in) :: ldata(:,:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name(:)

    integer :: n
    
    INSIST(size(ldata,dim=2) == ncells)
    
    call seq%write_dist_array(name, ncells_tot, ldata)
    
    if (for_viz) then
      call seq%write_dataset_attr(name, 'FIELDTYPE', 'CELL')
      INSIST(present(viz_name))
      INSIST(size(viz_name) == size(ldata,1))
      ! Argh, no array attributes!  So we'll do it this way for now.
      !call seq%write_dataset_attr(name, 'FIELDNAME', viz_name)
      do n = 1, size(viz_name)
        call seq%write_dataset_attr(name, 'FIELDNAME'//i_to_c(n), viz_name(n))
      end do
    end if
    
  end subroutine write_seq_cell_field_r1
  
  subroutine write_seq_node_field_r0 (seq, ldata, name, for_viz, viz_name)
  
    use legacy_mesh_api, only: nnodes, nnodes_tot

    type(th5_seq_group), intent(in) :: seq
    real(r8), intent(in) :: ldata(:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name

    INSIST(size(ldata) == nnodes)
    
    call seq%write_dist_array(name, nnodes_tot, ldata)
    
    if (for_viz) then
      call seq%write_dataset_attr(name, 'FIELDTYPE', 'NODE')
      if (present(viz_name)) then
        call seq%write_dataset_attr(name, 'FIELDNAME', viz_name)
      else
        call seq%write_dataset_attr(name, 'FIELDNAME', name)
      end if
    end if
    
  end subroutine write_seq_node_field_r0
  
  subroutine write_seq_node_field_r1 (seq, ldata, name, for_viz, viz_name)
  
    use legacy_mesh_api, only: nnodes, nnodes_tot
    use string_utilities, only: i_to_c

    type(th5_seq_group), intent(in) :: seq
    real(r8), intent(in) :: ldata(:,:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name(:)

    integer :: n
    
    INSIST(size(ldata,dim=2) == nnodes)
    
    call seq%write_dist_array(name, nnodes_tot, ldata)
    
    if (for_viz) then
      call seq%write_dataset_attr(name, 'FIELDTYPE', 'NODE')
      INSIST(present(viz_name))
      INSIST(size(viz_name) == size(ldata,1))
      ! Argh, no array attributes!  So we'll do it this way for now.
      !call seq%write_dataset_attr(name, 'FIELDNAME', viz_name)
      do n = 1, size(viz_name)
        call seq%write_dataset_attr(name, 'FIELDNAME'//i_to_c(n), viz_name(n))
      end do
    end if
    
  end subroutine write_seq_node_field_r1

end module truchas_danu_output_tools
