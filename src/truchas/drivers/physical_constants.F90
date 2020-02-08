!!
!! PHYSICAL_CONSTANTS
!!
!! This module provides physical constants used by various physics modules.
!! These are fixed constants whose value depend on the particular system of
!! units being used.  The following constants are provided and, by default,
!! have been assigned their SI unit values.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  STEFAN_BOLTZMANN -- the Stefan-Boltzmann constant that appears in thermal
!!    radiation models.
!!  ABSOLUTE_ZERO -- the value of absolute zero in the temperature scale; used
!!    in thermal radiation models.
!!  VACUUM_PERMITTIVITY, VACUUM_PERMEABILITY -- the electric permittivity and
!!    magnetic permeability of free space; used by the electromagnetic solver.
!!
!! The following subroutine is provided to overwrite the default values with
!! user-specified values read from an input file.
!!
!!  CALL READ_PHYSICAL_CONSTANTS (LUN) reads the above variables from the
!!    first occurrence of the physical_constants namelist in the file opened
!!    on logical unit LUN.  If the namelist is not found no action is taken.
!!    Any variable not specified in the namelist retains its default value.
!!
!! Note to maintainers:  These constants should be considered read-only, and
!! their values should only be defined by methods provided by this module.
!! Please adhere to this policy when adding new constants.
!!

module physical_constants

  use kinds, only: r8
  implicit none
  private

  public :: read_physical_constants

  !! Stefan-Boltzmann constant; default in SI units.
  real(r8), public, save :: stefan_boltzmann = 5.67e-8_r8

  !! Absolute zero; default for Kelvin temperature scale.
  real(r8), public, save :: absolute_zero = 0.0_r8

  !! Permittivity of free space; default in SI units.
  real(r8), public, save :: vacuum_permittivity = 8.854188e-12_r8

  !! Permeability of free space; default in SI units.
  real(r8), public, save :: vacuum_permeability = 1.256637e-6_r8

contains

  subroutine read_physical_constants (lun)

    use input_utilities, only: seek_to_namelist
    use parallel_communication
    use truchas_logging_services

    integer, intent(in) :: lun

    integer :: ios
    logical :: found

    namelist /physical_constants/ stefan_boltzmann, absolute_zero, &
                                  vacuum_permittivity, vacuum_permeability

    !! Position the input file at the first occurrence of the namelist.
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist (lun, 'PHYSICAL_CONSTANTS', found, iostat=ios)
    end if
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading input file')

    !! Return if we didn't find it -- it's an optional namelist.
    call broadcast (found)
    if (.not.found) return

    !! Read the namelist; variables not read retain their default value.
    call TLS_info ('')
    call TLS_info (' Reading PHYSICAL_CONSTANTS namelist ...')
    if (is_IOP) read(lun,nml=physical_constants,iostat=ios)
    call broadcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading PHYSICAL_CONSTANTS namelist')

    !! Replicate the values on all processes.
    call broadcast (Stefan_Boltzmann)
    call broadcast (Absolute_Zero)
    call broadcast (Vacuum_Permittivity)
    call broadcast (Vacuum_Permeability)

    !! Check the values.
    if (Stefan_Boltzmann <= 0) call TLS_fatal ('STEFAN_BOLTZMANN must be > 0.0')
    if (vacuum_permeability <= 0) call TLS_fatal ('VACUUM_PERMEABILITY must be > 0.0')
    if (vacuum_permittivity <= 0) call TLS_fatal ('VACUUM_PERMITTIVITY must be > 0.0')

  end subroutine read_physical_constants

end module physical_constants
