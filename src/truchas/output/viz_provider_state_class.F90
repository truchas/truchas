!!
!! VIZ_PROVIDER_STATE_CLASS
!!
!! Abstract stream-local provider state for VTKHDF output. A provider state
!! is the provider-specific object owned by one viz_stream_state. It stores
!! VTKHDF dataset handles and any stream-local setup or selected-field state
!! needed by that provider. The stream calls it at the three output lifecycle
!! points: registering primary mesh-block temporal data, writing primary
!! mesh-block timestep values, and optionally defining/writing provider-defined
!! blocks such as solid-mechanics interface blocks.
!!

module viz_provider_state_class

  use vtkhdf_mb_file_type
  implicit none
  private

  type, abstract, public :: viz_provider_state
  contains
    procedure :: configure_blocks
    procedure(register_mesh_block_temporal_data_ifc), deferred :: register_mesh_block_temporal_data
    procedure(write_mesh_block_timestep_ifc), deferred :: write_mesh_block_timestep
    procedure :: define_provider_blocks
    procedure :: write_provider_blocks_timestep
  end type

  type, public :: viz_provider_state_box
    class(viz_provider_state), allocatable :: p
  end type

  abstract interface

    subroutine register_mesh_block_temporal_data_ifc(this, file, iblock, block, &
        block_cells, block_nodes)
      import :: vtkhdf_mb_file, vtkhdf_block_handle, viz_provider_state
      class(viz_provider_state), intent(inout) :: this
      type(vtkhdf_mb_file), intent(inout) :: file
      integer, intent(in) :: iblock
      type(vtkhdf_block_handle), intent(in) :: block
      integer, intent(in) :: block_cells(:), block_nodes(:)
    end subroutine

    subroutine write_mesh_block_timestep_ifc(this, file, iblock, block, &
        block_cells, block_nodes)
      import :: vtkhdf_mb_file, vtkhdf_block_handle, viz_provider_state
      class(viz_provider_state), intent(inout) :: this
      type(vtkhdf_mb_file), intent(inout) :: file
      integer, intent(in) :: iblock
      type(vtkhdf_block_handle), intent(in) :: block
      integer, intent(in) :: block_cells(:), block_nodes(:)
    end subroutine

  end interface

contains

  ! Optional pre-registration hook for block-level output policy that is owned
  ! by the stream rather than intrinsic to the provider's fields.
  subroutine configure_blocks(this, block_requires_temporal_output)
    class(viz_provider_state), intent(inout) :: this
    logical, intent(in) :: block_requires_temporal_output(:)
  end subroutine

  ! Optional stream-level setup hook for provider-defined blocks that are not
  ! naturally handled by the primary mesh-block registration path.
  subroutine define_provider_blocks(this, file)
    class(viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
  end subroutine

  ! Optional stream-level timestep write hook paired with
  ! define_provider_blocks() for provider-defined block output.
  subroutine write_provider_blocks_timestep(this, file)
    class(viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
  end subroutine

end module viz_provider_state_class
