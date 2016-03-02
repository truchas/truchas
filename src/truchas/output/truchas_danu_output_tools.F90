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
  use danu_module
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
  
  subroutine write_seq_cell_field_r0 (seq_id, ldata, name, for_viz, viz_name)
  
    use legacy_mesh_api, only: ncells, ncells_tot

    type(c_ptr), intent(in) :: seq_id
    real(r8), intent(in) :: ldata(:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name

    integer :: stat
    real(r8), allocatable :: gdata(:)
    type(c_ptr) :: dataset_id
    
    INSIST(size(ldata) == ncells)
    
    if (nPE == 1) then
      call simulation_data_write (seq_id, name, ldata, stat)
    else
      if (is_IOP) then
        allocate(gdata(ncells_tot))
      else
        allocate(gdata(0))
      end if
      call collate (gdata, ldata)
      if (is_IOP) call simulation_data_write (seq_id, name, gdata, stat)
      call broadcast (stat)
    end if
    
    INSIST(stat == DANU_SUCCESS)
    
    if (for_viz) then
      if (is_IOP) call simulation_open_data (seq_id, name, dataset_id, stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      if (is_IOP) call attribute_write (dataset_id, 'FIELDTYPE', 'CELL', stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      if (is_IOP) then
        if (present(viz_name)) then
          call attribute_write (dataset_id, 'FIELDNAME', viz_name, stat)
        else
          call attribute_write (dataset_id, 'FIELDNAME', name, stat)
        end if
      end if
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
    end if
    
  end subroutine write_seq_cell_field_r0
  
  subroutine write_seq_cell_field_r1 (seq_id, ldata, name, for_viz, viz_name)
  
    use legacy_mesh_api, only: ncells, ncells_tot
    use string_utilities, only: i_to_c

    type(c_ptr), intent(in) :: seq_id
    real(r8), intent(in) :: ldata(:,:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name(:)

    integer :: stat, n
    real(r8), allocatable :: gdata(:,:)
    type(c_ptr) :: dataset_id
    
    INSIST(size(ldata,dim=2) == ncells)
    
    if (nPE == 1) then
      call simulation_data_write (seq_id, name, ldata, stat)
    else
      if (is_IOP) then
        allocate(gdata(size(ldata,1),ncells_tot))
      else
        allocate(gdata(size(ldata,1),0))
      end if
      call collate (gdata, ldata)
      if (is_IOP) call simulation_data_write (seq_id, name, gdata, stat)
      call broadcast (stat)
    end if
    
    INSIST(stat == DANU_SUCCESS)
    
    if (for_viz) then
      if (is_IOP) call simulation_open_data (seq_id, name, dataset_id, stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      if (is_IOP) call attribute_write (dataset_id, 'FIELDTYPE', 'CELL', stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      INSIST(present(viz_name))
      INSIST(size(viz_name) == size(ldata,1))
      ! Argh, no array attributes!  So we'll do it this way for now.
      !if (is_IOP) call attribute_write (dataset_id, 'FIELDNAME', viz_name, stat)
      !call broadcast (stat)
      !INSIST(stat == DANU_SUCCESS)
      do n = 1, size(viz_name)
        if (is_IOP) call attribute_write (dataset_id, 'FIELDNAME'//i_to_c(n), viz_name(n), stat)
        call broadcast (stat)
        INSIST(stat == DANU_SUCCESS)
      end do
    end if
    
  end subroutine write_seq_cell_field_r1
  
  subroutine write_seq_node_field_r0 (seq_id, ldata, name, for_viz, viz_name)
  
    use legacy_mesh_api, only: nnodes, nnodes_tot

    type(c_ptr), intent(in) :: seq_id
    real(r8), intent(in) :: ldata(:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name

    integer :: stat
    real(r8), allocatable :: gdata(:)
    type(c_ptr) :: dataset_id
    
    INSIST(size(ldata) == nnodes)
    
    if (nPE == 1) then
      call simulation_data_write (seq_id, name, ldata, stat)
    else
      if (is_IOP) then
        allocate(gdata(nnodes_tot))
      else
        allocate(gdata(0))
      end if
      call collate (gdata, ldata)
      if (is_IOP) call simulation_data_write (seq_id, name, gdata, stat)
      call broadcast (stat)
    end if
    
    INSIST(stat == DANU_SUCCESS)
    
    if (for_viz) then
      if (is_IOP) call simulation_open_data (seq_id, name, dataset_id, stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      if (is_IOP) call attribute_write (dataset_id, 'FIELDTYPE', 'NODE', stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      if (is_IOP) then
        if (present(viz_name)) then
          call attribute_write (dataset_id, 'FIELDNAME', viz_name, stat)
        else
          call attribute_write (dataset_id, 'FIELDNAME', name, stat)
        end if
      end if
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
    end if
    
  end subroutine write_seq_node_field_r0
  
  subroutine write_seq_node_field_r1 (seq_id, ldata, name, for_viz, viz_name)
  
    use legacy_mesh_api, only: nnodes, nnodes_tot
    use string_utilities, only: i_to_c

    type(c_ptr), intent(in) :: seq_id
    real(r8), intent(in) :: ldata(:,:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name(:)

    integer :: stat, n
    real(r8), allocatable :: gdata(:,:)
    type(c_ptr) :: dataset_id
    
    INSIST(size(ldata,dim=2) == nnodes)
    
    if (nPE == 1) then
      call simulation_data_write (seq_id, name, ldata, stat)
    else
      if (is_IOP) then
        allocate(gdata(size(ldata,1),nnodes_tot))
      else
        allocate(gdata(size(ldata,1),0))
      end if
      call collate (gdata, ldata)
      if (is_IOP) call simulation_data_write (seq_id, name, gdata, stat)
      call broadcast (stat)
    end if
    
    INSIST(stat == DANU_SUCCESS)
    
    if (for_viz) then
      if (is_IOP) call simulation_open_data (seq_id, name, dataset_id, stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      if (is_IOP) call attribute_write (dataset_id, 'FIELDTYPE', 'NODE', stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      INSIST(present(viz_name))
      INSIST(size(viz_name) == size(ldata,1))
      ! Argh, no array attributes!  So we'll do it this way for now.
      !if (is_IOP) call attribute_write (dataset_id, 'FIELDNAME', viz_name, stat)
      !call broadcast (stat)
      !INSIST(stat == DANU_SUCCESS)
      do n = 1, size(viz_name)
        if (is_IOP) call attribute_write (dataset_id, 'FIELDNAME'//i_to_c(n), viz_name(n), stat)
        call broadcast (stat)
        INSIST(stat == DANU_SUCCESS)
      end do
    end if
    
  end subroutine write_seq_node_field_r1

end module truchas_danu_output_tools
