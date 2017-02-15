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
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module truchas_danu_output_tools

  use kinds, only: r8
  use truchasio, only: DANU_SUCCESS, sequence
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

    class(sequence), intent(in) :: seq
    real(r8), intent(in) :: ldata(:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name

    integer :: stat
    
    INSIST(size(ldata) == ncells)
    
    call seq%data_write (name, ncells_tot, ldata, stat)
    
    INSIST(stat == DANU_SUCCESS)
    
    if (for_viz) then
      INSIST(stat == DANU_SUCCESS)
      call seq%write_dataset_attr(name, 'FIELDTYPE', 'CELL', stat)
      INSIST(stat == DANU_SUCCESS)
      if (present(viz_name)) then
        call seq%write_dataset_attr(name, 'FIELDNAME', viz_name, stat)
      else
        call seq%write_dataset_attr(name, 'FIELDNAME', name, stat)
      end if
      INSIST(stat == DANU_SUCCESS)
    end if
    
  end subroutine write_seq_cell_field_r0
  
  subroutine write_seq_cell_field_r1 (seq, ldata, name, for_viz, viz_name)
  
    use legacy_mesh_api, only: ncells, ncells_tot
    use string_utilities, only: i_to_c

    class(sequence), intent(in) :: seq
    real(r8), intent(in) :: ldata(:,:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name(:)

    integer :: stat, n
    
    INSIST(size(ldata,dim=2) == ncells)
    
    call seq%data_write (name, ncells_tot, ldata, stat)
    
    INSIST(stat == DANU_SUCCESS)
    
    if (for_viz) then
      INSIST(stat == DANU_SUCCESS)
      call seq%write_dataset_attr(name, 'FIELDTYPE', 'CELL', stat)
      INSIST(stat == DANU_SUCCESS)
      INSIST(present(viz_name))
      INSIST(size(viz_name) == size(ldata,1))
      ! Argh, no array attributes!  So we'll do it this way for now.
      !call seq%write_dataset_attr(name, 'FIELDNAME', viz_name, stat)
      !INSIST(stat == DANU_SUCCESS)
      do n = 1, size(viz_name)
        call seq%write_dataset_attr(name, 'FIELDNAME'//i_to_c(n), viz_name(n), stat)
        INSIST(stat == DANU_SUCCESS)
      end do
    end if
    
  end subroutine write_seq_cell_field_r1
  
  subroutine write_seq_node_field_r0 (seq, ldata, name, for_viz, viz_name)
  
    use legacy_mesh_api, only: nnodes, nnodes_tot

    class(sequence), intent(in) :: seq
    real(r8), intent(in) :: ldata(:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name

    integer :: stat
    
    INSIST(size(ldata) == nnodes)
    
    call seq%data_write (name, nnodes_tot, ldata, stat)
    
    INSIST(stat == DANU_SUCCESS)
    
    if (for_viz) then
      INSIST(stat == DANU_SUCCESS)
      call seq%write_dataset_attr(name, 'FIELDTYPE', 'NODE', stat)
      INSIST(stat == DANU_SUCCESS)
      if (present(viz_name)) then
        call seq%write_dataset_attr(name, 'FIELDNAME', viz_name, stat)
      else
        call seq%write_dataset_attr(name, 'FIELDNAME', name, stat)
      end if
      INSIST(stat == DANU_SUCCESS)
    end if
    
  end subroutine write_seq_node_field_r0
  
  subroutine write_seq_node_field_r1 (seq, ldata, name, for_viz, viz_name)
  
    use legacy_mesh_api, only: nnodes, nnodes_tot
    use string_utilities, only: i_to_c

    class(sequence), intent(in) :: seq
    real(r8), intent(in) :: ldata(:,:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name(:)

    integer :: stat, n
    
    INSIST(size(ldata,dim=2) == nnodes)
    
    call seq%data_write (name, nnodes_tot, ldata, stat)
    
    INSIST(stat == DANU_SUCCESS)
    
    if (for_viz) then
      INSIST(stat == DANU_SUCCESS)
      call seq%write_dataset_attr(name, 'FIELDTYPE', 'NODE', stat)
      INSIST(stat == DANU_SUCCESS)
      INSIST(present(viz_name))
      INSIST(size(viz_name) == size(ldata,1))
      ! Argh, no array attributes!  So we'll do it this way for now.
      !call seq%write_dataset_attr(name, 'FIELDNAME', viz_name, stat)
      !INSIST(stat == DANU_SUCCESS)
      do n = 1, size(viz_name)
        call seq%write_dataset_attr(name, 'FIELDNAME'//i_to_c(n), viz_name(n), stat)
        INSIST(stat == DANU_SUCCESS)
      end do
    end if
    
  end subroutine write_seq_node_field_r1

end module truchas_danu_output_tools
