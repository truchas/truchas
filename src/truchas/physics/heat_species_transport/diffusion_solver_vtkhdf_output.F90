!!
!! DIFFUSION_SOLVER_VTKHDF_OUTPUT
!!
!! VTKHDF graphics output for heat-transfer and species-diffusion fields.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module diffusion_solver_vtkhdf_output

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use diffusion_solver, only: ds_get_phi_view
  use diffusion_solver_data, only: heat_eqn, num_species
  use vtkhdf_mb_file_type, only: vtkhdf_mb_file, vtkhdf_block_handle, &
      vtkhdf_cell_data_handle
  implicit none
  private

  public :: diffusion_solver_vtkhdf_block_data
  public :: diffusion_solver_vtkhdf_register_temporal_data
  public :: diffusion_solver_vtkhdf_output_step

  !! The diffusion-solver VTKHDF data registered for one mesh block.
  type :: diffusion_solver_vtkhdf_block_data
    private
    type(vtkhdf_cell_data_handle) :: temperature
    type(vtkhdf_cell_data_handle) :: enthalpy
    type(vtkhdf_cell_data_handle), allocatable :: species(:)
  end type

contains

  subroutine diffusion_solver_vtkhdf_register_temporal_data(file, block, block_cells, &
      moving, data)

    use string_utilities, only: i_to_c
    use zone_module, only: zone

    type(vtkhdf_mb_file), intent(inout) :: file
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:)
    logical, intent(in) :: moving
    type(diffusion_solver_vtkhdf_block_data), intent(inout) :: data

    integer :: m

    if (heat_eqn .or. moving) then
      data%temperature = file%register_temporal_cell_data(block, 'T', 0.0_r8)
      data%enthalpy = file%register_temporal_cell_data(block, 'Enthalpy', 0.0_r8)
    else
      call file%write_cell_data(block, 'T', zone(block_cells)%temp)
      call file%write_cell_data(block, 'Enthalpy', zone(block_cells)%enthalpy)
    end if

    if (num_species > 0) then
      allocate(data%species(num_species))
      do m = 1, num_species
        data%species(m) = file%register_temporal_cell_data(block, &
            'phi'//i_to_c(m), 0.0_r8)
      end do
    end if

  end subroutine diffusion_solver_vtkhdf_register_temporal_data

  subroutine diffusion_solver_vtkhdf_output_step(file, block, block_cells, moving, &
      block_data)

    use zone_module, only: zone

    type(vtkhdf_mb_file), intent(inout) :: file
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:)
    logical, intent(in) :: moving
    type(diffusion_solver_vtkhdf_block_data), intent(inout) :: block_data

    integer :: m
    real(r8), pointer :: phi(:)

    if (heat_eqn .or. moving) then
      call file%write_cell_data(block, block_data%temperature, zone(block_cells)%temp)
      call file%write_cell_data(block, block_data%enthalpy, zone(block_cells)%enthalpy)
    end if

    do m = 1, num_species
      call ds_get_phi_view(m, phi)
      call file%write_cell_data(block, block_data%species(m), phi(block_cells))
    end do

  end subroutine diffusion_solver_vtkhdf_output_step

end module diffusion_solver_vtkhdf_output
