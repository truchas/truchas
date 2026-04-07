!!
!! FLOW_VTKHDF_OUTPUT
!!
!! VTKHDF graphics output for flow fields.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module flow_vtkhdf_output

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use flow_driver, only: flow_enabled, flow_vel_cc_view, flow_P_cc_view
  use vtkhdf_mb_file_type, only: vtkhdf_mb_file, vtkhdf_block_handle, &
      vtkhdf_cell_data_handle
  implicit none
  private

  public :: flow_vtkhdf_block_data
  public :: flow_vtkhdf_register_temporal_data
  public :: flow_vtkhdf_output_step

  !! The flow VTKHDF data registered for one mesh block.
  type :: flow_vtkhdf_block_data
    private
    type(vtkhdf_cell_data_handle) :: velocity
    type(vtkhdf_cell_data_handle) :: pressure
  end type

contains

  subroutine flow_vtkhdf_register_temporal_data(file, block, data)

    type(vtkhdf_mb_file), intent(inout) :: file
    type(vtkhdf_block_handle), intent(in) :: block
    type(flow_vtkhdf_block_data), intent(inout) :: data

    real(r8) :: vector(3)

    if (.not.flow_enabled()) return

    data%velocity = file%register_temporal_cell_data(block, 'Velocity', vector)
    data%pressure = file%register_temporal_cell_data(block, 'P', 0.0_r8)

  end subroutine flow_vtkhdf_register_temporal_data

  subroutine flow_vtkhdf_output_step(file, block, block_cells, data)

    type(vtkhdf_mb_file), intent(inout) :: file
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:)
    type(flow_vtkhdf_block_data), intent(inout) :: data

    real(r8), pointer :: p(:), v(:,:)

    if (.not.flow_enabled()) return

    v => flow_vel_cc_view()
    p => flow_P_cc_view()
    call file%write_cell_data(block, data%velocity, v(:,block_cells))
    call file%write_cell_data(block, data%pressure, p(block_cells))

  end subroutine flow_vtkhdf_output_step

end module flow_vtkhdf_output
