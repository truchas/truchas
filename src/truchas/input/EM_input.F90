!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module EM_input

  use mesh_manager, only: enable_mesh
  use electromagnetics_namelist, only: read_electromagnetics_namelist, params
  use induction_coil_namelist, only: read_induction_coil_namelists
  use electromagnetic_bc_namelist, only: read_electromagnetic_bc_namelists
  use truchas_logging_services
  implicit none
  private

  public :: read_EM_input

contains

  subroutine read_EM_input(lun)
    use parameter_list_type
    integer, intent(in) :: lun
    type(parameter_list), pointer :: plist
    logical :: exists
    call read_electromagnetics_namelist(lun)
    plist => params%sublist('coils')
    call read_induction_coil_namelists(lun, plist)
    plist => params%sublist('bc')
    call read_electromagnetic_bc_namelists(lun, plist)
    call enable_mesh('em', exists)
    if (.not.exists) call TLS_fatal('EM_MESH namelist was not specified')
  end subroutine

end module EM_input
