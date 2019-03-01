!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module porous_drag_data

  use kinds, only: r8
  implicit none
  private

  ! PHYSICS namelist variables -
  ! Flag for enabling/disabling the porous drag model.
  logical, public, save :: porous_flow

  ! NUMERICS namelist input
  ! Flag for weighting time-integrator
  real(r8), public, save :: porous_implicitness

end module porous_drag_data
