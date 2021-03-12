!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This used to hold the DIFFUSION_SOLVER namelist input variables and a
! a subroutine for reading it. That has been extracted into a new procedure
! that repackages the input into a parameter list. What remains here is the
! leftover stuff that didn't belong with the new procedure. This ultimately
! needs to go elsewhere as well.

#include "f90_assert.fpp"

module diffusion_solver_data

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  public :: ds_data_init

  !! These variables are defined manually prior to calling READ_DS_NAMELIST
  logical, public :: ds_enabled = .false.
  character(32), public :: system_type = ''
  integer,  public :: num_species = 0

  !! These variables are defined by DS_DATA_INIT.
  logical, public :: heat_eqn = .false.
  integer, public :: ds_sys_type = 0
  integer, parameter, public :: DS_SPEC_SYS = 1
  integer, parameter, public :: DS_TEMP_SYS = 2
  integer, parameter, public :: DS_TEMP_SPEC_SYS = 3

  !! This will be set later using input data from the DIFFUSION_SOLVER namelist
  real(r8), public :: void_temperature

  character(4), parameter, public :: mesh_name = 'MAIN'

contains

  !! This (temporary!) stub is what remains after extracting the namelist
  !! reading portion of read_ds_namelist to a separate procedure.

  subroutine ds_data_init

    use string_utilities

    INSIST(ds_enabled)
    system_type = raise_case(adjustl(system_type))
    select case (system_type)
    case ('SPECIES')
      ds_sys_type = DS_SPEC_SYS
      heat_eqn = .false.
      INSIST(num_species > 0)
    case ('THERMAL')
      ds_sys_type = DS_TEMP_SYS
      heat_eqn = .true.
      num_species = 0
    case ('THERMAL+SPECIES')
      ds_sys_type = DS_TEMP_SPEC_SYS
      heat_eqn = .true.
      INSIST(num_species > 0)
    case default
      INSIST(.false.)
    end select

  end subroutine

end module diffusion_solver_data
