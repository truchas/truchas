#include "f90_assert.fpp"

module flow_driver_type

  use kinds, only: r8
  use flow_mesh_type
  use flow_type
  use parameter_list_type
  use truchas_logging_services
  use truchas_timers
  use index_partitioning
  implicit none
  private

  public :: read_flow_namelist
  public :: flow_driver_init, flow_step, flow_final, flow_enabled

  type :: flow_driver
    type(flow_mesh), pointer :: mesh => null()
    type(flow) :: flow
    type(flow_props) :: props
  end type flow_driver
  type(flow_driver), allocatable :: this
  type(parameter_list) :: params

contains

  subroutine flow_final
    if (allocated(this)) deallocate(this)
  end subroutine flow_final


  logical function flow_enabled()
    flow_enabled = allocated(this)
  end function flow_enabled


  subroutine read_flow_namelist(lun)

    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use input_utilities, only: seek_to_namelist

    integer, intent(in) :: lun

    .......................
  end subroutine read_flow_namelist


  subroutine flow_driver_init(t, mesh, ...)
    real(r8), intent(in) :: t
    type(unstr_mesh), pointer, intent(in) :: mesh



    ..............
  end subroutine flow_driver_init

  subroutine flow_step(t, dt)
    real(r8), intent(in) :: t, dt
    !-

    ! need to duplicate the matl query code here.  When the vof and flow drivers are subsumed
    ! into a more inteligent driver structure, this should be reworked.

    call props%update(vof, temperature_cc, initial=.false.)
    call flow%zero_out_solid_velocities(props)

  end subroutine flow_step

end module flow_driver_type
