!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
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
