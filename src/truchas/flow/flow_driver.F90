!!
!! FLOW_DRIVER
!!

#include "f90_assert.fpp"

module flow_driver

  use kinds, only: r8
  use unstr_mesh_type
  use flow_model
  use vof_model
  use parameter_list_type
  use truchas_logging_services
  use timing_tree
  implicit none
  private

  public :: read_flow_namelist
  public :: flow_init, flow_update_vof, flow_update_velocity
  public :: flow_set_initial_vof
  public :: flow_destroy
  !public :: flow_read_checkpoint, flow_skip_checkpoint
  public :: flow_enabled

  !! Bundle up all the driver state data as a singleton THIS of private
  !! derived type.  All procedures use/modify this object.
  type :: flow_state
    type(unstr_mesh), pointer :: mesh => null()  ! reference only -- do not own
    !integer, allocatable :: sol_matid(:), liq_matid(:)
    type(flow_model_t) :: flow_t
    type(vof_model_t) :: vof_t
 end type flow_state
  type(flow_state), allocatable, save :: state

  !! Input data cached in a private parameter list.
  type(parameter_list), save :: params

contains

  subroutine flow_destroy
    if (allocated(state)) deallocate(state)
  end subroutine flow_destroy

  logical function flow_enabled ()
    flow_enabled = allocated(state)
  end function flow_enabled

  !! Current Truchas design requires that parameter input and object
  !! initialization be separated and occur at distinct execution stages.
  !! The following procedure reads the MICROSTRUCTURE namelist and stores
  !! the data in private module data.

  subroutine read_flow_namelist (lun)
    integer, intent(in) :: lun

    print*, 'read_flow_namelist(lun): ', lun
    allocate(state)
  end subroutine read_flow_namelist

  !! This initializes the driver.  It should only be called if heat transfer is
  !! enabled, and after its initialization (mesh_interop data).

  subroutine flow_init (t)
    real(r8), intent(in) :: t
    print *, 'flow_init(t) :', t
  end subroutine flow_init

  subroutine flow_update_vof (t)
    real(r8), intent(in) :: t
    print *, 'flow_update_vof(t): ', t
  end subroutine flow_update_vof

  subroutine flow_update_velocity (t)
    real(r8), intent(in) :: t
    print *, 'flow_update_velocity(t): ', t
  end subroutine flow_update_velocity

  subroutine flow_set_initial_vof (mesh, nbody, vof)
    type(unstr_mesh), intent(in) :: mesh
    integer, intent(in) :: nbody
    real(r8), intent(out) :: vof(:,:)
    print *, 'flow_set_initial_vof'
    call state%vof_t%initial_state_all(mesh, nbody, vof)
  end subroutine flow_set_initial_vof
end module flow_driver
