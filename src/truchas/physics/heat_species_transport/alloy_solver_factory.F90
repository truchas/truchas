!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module alloy_solver_factory

  use alloy_model_type
  use alloy_solver_type
  use matl_mesh_func_type
  use parameter_list_type
  implicit none
  private

  public :: create_alloy_solver

contains

  function create_alloy_solver(model, params, stat, errmsg) result(solver)

    use parallel_communication
    use truchas_env, only: output_file_name

    type(alloy_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    type(alloy_solver), pointer :: solver
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: lun
    type(parameter_list), pointer :: plist
    logical :: verbose_stepping

    call params%get('verbose-stepping', verbose_stepping)
    if (verbose_stepping) then
      lun = -1
      if (is_IOP) open(newunit=lun,file=output_file_name('bdf2.out'),position='rewind',action='write')
      call params%set('output-unit', lun)
      plist => params%sublist('norm')
      call plist%set('unit', lun)
    end if

    allocate(solver)
    call solver%init(model, params, stat, errmsg)

  end function create_alloy_solver

end module alloy_solver_factory
