!!
!! PROBES_DRIVER
!!
!! Driver subroutines for solution probes and holds the probes object.
!! Consider this part of the Truchas multiphysics driver.
!!
!! Michael Hall <hall@lanl.gov>
!! July 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module probes_driver

  use probes_type
  implicit none
  private

  public :: probes_init, probes_write

  type(probes) :: p

contains

  subroutine probes_init

    use unstr_mesh_type
    use truchas_probe_field_factory_type
    use probe_namelist, only: params
    use mesh_manager, only: unstr_mesh_ptr
    use truchas_env, only: output_dir
    use truchas_logging_services

    type(unstr_mesh), pointer :: mesh
    type(truchas_probe_field_factory) :: pf_fac
    integer :: stat
    character(:), allocatable :: errmsg

    mesh => unstr_mesh_ptr('MAIN')

    call TLS_info('')
    call TLS_info('Initializing solution probes ...')
    call p%init(mesh, pf_fac, trim(output_dir), params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('error initializing probes: ' // errmsg)

  end subroutine probes_init

  subroutine probes_write(time)
    use,intrinsic :: iso_fortran_env, only: r8 => real64
    real(r8), intent(in) :: time
    call p%write(time)
  end subroutine probes_write

end module probes_driver
