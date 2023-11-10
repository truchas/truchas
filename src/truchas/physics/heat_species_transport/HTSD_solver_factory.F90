!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module HTSD_solver_factory

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use HTSD_model_type
  use HTSD_solver_type
  use matl_mesh_func_type
  use parameter_list_type
  implicit none
  private

  public :: create_HTSD_solver

contains

  function create_HTSD_solver(mmf, model, params) result(solver)

    use enclosure_radiation_namelist, only: er_params => params
    use parallel_communication
    use truchas_env, only: output_file_name

    type(matl_mesh_func), intent(in), target :: mmf
    type(HTSD_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    type(HTSD_solver), pointer :: solver

    integer :: j, n, lun
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    character(:), allocatable :: string
    character(16), allocatable :: vfr_precon_coupling(:)
    real(r8), allocatable :: rad_tol(:)
    logical :: verbose_stepping

    call params%get('verbose-stepping', verbose_stepping)
    if (verbose_stepping) then
      lun = -1
      if (is_IOP) open(newunit=lun,file=output_file_name('bdf2.out'),position='rewind',action='write')
      call params%set('output-unit', lun)
      plist => params%sublist('norm')
      call plist%set('unit', lun)
    end if

    !TODO: read enclosure radiation namelists directly into sublist
    if (associated(model%ht)) then
      if (associated(model%ht%vf_rad_prob)) then
        !TODO! Assumes model%ht%vf_rad_prob array was initialized under the same
        !TODO! loop so that that array and vfr_precon_coupling are in correspondence.
        !TODO! This is fragile and needs to be fixed.
        piter = parameter_list_iterator(er_params)
        n = piter%count()
        allocate(vfr_precon_coupling(n), rad_tol(n))
        do j = 1, n
          plist => piter%sublist()
          call plist%get('precon-coupling-method', string, default='BACKWARD GS')
          vfr_precon_coupling(j) = string
          call piter%next
        end do
        plist => params%sublist('precon')
        call plist%set('vfr-precon-coupling', vfr_precon_coupling)
      end if
    end if

    allocate(solver)
    call HTSD_solver_init (solver, mmf, model, params)

  end function create_HTSD_solver

end module HTSD_solver_factory
