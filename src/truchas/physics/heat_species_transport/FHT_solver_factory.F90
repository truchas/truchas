!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module FHT_solver_factory

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use FHT_model_type
  use FHT_solver_type
  use matl_mesh_func_type
  use parameter_list_type
  implicit none
  private

  public :: create_FHT_solver

contains

  function create_FHT_solver(mmf, model, params) result(solver)

    use enclosure_radiation_namelist, only: er_params => params
    use cutoffs_module, only: cutvof
    use parallel_communication
    use truchas_env, only: output_file_name

    type(matl_mesh_func), intent(in), target :: mmf
    type(FHT_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    type(FHT_solver), pointer :: solver

    integer :: j, n, lun
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    character(:), allocatable :: string
    character(16), allocatable :: vfr_precon_coupling(:)
    real(r8), allocatable :: rad_tol(:)
    logical :: verbose_stepping

    !FIXME! we need some user-settable parameters here, especially DELTA
    call params%set('tofh-tol', 0.0d0)
    call params%set('tofh-delta', 1.0d-3)
    call params%set('tofh-max-try', 50)

    call params%get('verbose-stepping', verbose_stepping)
    if (verbose_stepping) then
      lun = -1
      if (is_IOP) open(newunit=lun,file=output_file_name('bdf2.out'),position='rewind',action='write')
      call params%set('unit', lun)
      plist => params%sublist('norm')
      call plist%set('unit', lun)
    end if

    !TODO: read enclosure radiation namelists directly into sublist
    !! Nonlinear solver and preconditioner parameters connected to enclosure radiation.
    if (associated(model%vf_rad_prob)) then
      !TODO! Assumes model%vf_rad_prob array was initialized under the same
      !TODO! loop so that that array and vfr_precon_coupling are in correspondence.
      !TODO! This is fragile and needs to be fixed.
      piter = parameter_list_iterator(er_params, sublists_only=.true.)
      n = piter%count()
      allocate(vfr_precon_coupling(n), rad_tol(n))
      do j = 1, n
        plist => piter%sublist()
        call plist%get('precon-coupling-method', string, default='BACKWARD GS')
        vfr_precon_coupling(j) = string
        call plist%get('error-tol', rad_tol(j), default=1d-3)
        call piter%next
      end do
      plist => params%sublist('precon')
      call plist%set('vfr-precon-coupling', vfr_precon_coupling)
      plist => params%sublist('norm')
      call plist%set('rad-tol', rad_tol)
    end if

    allocate(solver)
    call FHT_solver_init(solver, mmf, model, params)

  end function create_FHT_solver

end module FHT_solver_factory
