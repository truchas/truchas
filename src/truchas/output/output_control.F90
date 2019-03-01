!!
!! OUTPUT_CONTROL
!!
!! This module was created to hold the variables from OUTPUT_DATA_MODULE and
!! OUTPUT_MODULE that control when the solution was output.  Everything else
!! in those modules is tied to TBrook and they are slated for deletion when
!! TBrook is ultimately deleted from the code.  I expect this module to exist
!! only temporarily until the control of output is refactored/reimplemented.
!!
!! Neil Carlson <nnc@lanl.gov>
!! October 2012
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module output_control

  use kinds, only: r8
  use parameter_module, only: mops
  use toolpath_type
  implicit none
  private

  real(r8), save, public :: output_dt(mops)
  real(r8), save, public :: output_t(mops+1)
  integer,  save, public :: output_dt_multiplier(mops) = 1
  integer,  save, public :: nops, next_op
  logical,  save, public :: precise_output
  logical,  save, public :: retain_last_step = .false.

  integer, allocatable, public :: part(:)
  type(toolpath), pointer, public :: part_path => null()

end module output_control
