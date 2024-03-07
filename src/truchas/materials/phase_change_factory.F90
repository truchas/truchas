!!
!! PHASE_CHANGE_FACTORY
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module phase_change_factory

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use phase_change_class
  use parameter_list_type
  implicit none
  private

  public :: alloc_phase_change

contains

  subroutine alloc_phase_change(pc, params, stat, errmsg)

    use smooth_phase_change_type
    use tabular_phase_change_type

    class(phase_change), allocatable, intent(out) :: pc
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: latent_heat, tlo, thi
    real(r8), allocatable :: table(:,:)

    if (params%is_parameter('solid-frac-table')) then

      call params%get('solid-frac-table', table, stat, errmsg)
      if (stat /= 0) return
      if (size(table,dim=1) /= 2) then
        stat = 1
        errmsg = 'malformed table'
        return
      end if
      call alloc_tabular_phase_change(pc, table(1,:), table(2,:), errmsg)
      if (.not.allocated(pc)) then
        stat = 1
        errmsg = 'tabular phase change: ' // errmsg
        return
      end if

    else if (params%is_parameter('solidus-temp') .and. &
             params%is_parameter('liquidus-temp')) then

      call params%get('solidus-temp', tlo, stat, errmsg)
      if (stat /= 0) return
      call params%get('liquidus-temp', thi, stat, errmsg)
      if (stat /= 0) return
      if (tlo >= thi) then
        errmsg = 'solidus-temp must be < liquidus-temp'
        return
      end if
      call alloc_smooth_phase_change(pc, tlo, thi)
      INSIST(allocated(pc))

    else

      stat = 1
      errmsg = 'unrecognized phase change specification'
      return

    end if

    if (params%is_parameter('latent-heat')) then
      call params%get('latent-heat', latent_heat, stat, errmsg)
      if (stat /= 0) return
      if (latent_heat < 0) then
        stat = 1
        errmsg = 'latent-heat must be > 0'
        return
      end if
      pc%latent_heat = latent_heat
    end if

  end subroutine alloc_phase_change

end module phase_change_factory
