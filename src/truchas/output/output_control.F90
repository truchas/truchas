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
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
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
  real(r8), save, public :: face_dump_time = -1e10
  logical,  save, public :: face_dumped = .false.
  real(r8), save, public :: face_dump_bbox(6)
  integer,  save, public :: temp_dump_freq
  integer,  save, public :: nops, next_op
  logical,  save, public :: precise_output
  logical,  save, public :: retain_last_step = .false.

  integer, allocatable, public :: part(:)
  type(toolpath), pointer, public :: part_path => null()

end module output_control
