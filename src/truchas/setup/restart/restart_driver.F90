!!
!! RESTART_DRIVER
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 26 Apr 2005; revised 28 July 2006.
!!
!! This module provides procedures for coordinating the reading of Truchas
!! restart file.  The restart file is read in parts from scattered places in
!! Truchas, so the primary purposes of the module is to serve as the focal
!! point of the process and document what must be read, and in what order.
!! To avoid circular module dependencies, the general restart variables have
!! been placed in the separate module RESTART_VARIABLES.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  o Several segments of the restart file are optional.  The following
!!    logical variables indicate whether the data segments are present:
!!
!!      HAVE_ALLOY_DATA, HAVE_FLUID_FLOW_DATA, HAVE_JOULE_HEAT_DATA,
!!      HAVE_SOLID_MECHANICS_DATA.
!!
!!    These variables are defined by OPEN_RESTART_MESH.
!!
!!  o OPEN_RESTART_MESH also defines the variables read from the header of the
!!    restart file:
!!
!!      RESTART_T, RESTART_DT (real);
!!      RESTART_NCELLS, RESTART_NNODES, RESTART_CYCLE_NUMBER (integer).
!!
!!  o The variables RESTART and RESTART_FILE must be set prior to calling any
!!    any procedure.  (These are set by PROCESS_COMMANDLINE.)
!!
!!  o The following public procedures must be called in order:
!!
!!    CALL OPEN_RESTART_FILE () opens the restart file and reads the header.
!!      Called from READ_INPUT prior to calling MESH_SIZES which requires the
!!      number of nodes and cells for the restart mesh read by this procedure.
!!
!!    The remaining routines should be called after the input file is read.
!!
!!    CALL RESTART_MESH () reads the mesh connectivity and vertex coordinates.
!!      Called from MESH_GEN.  The mesh and vertex partitions and associated
!!      permutations are determined after this procedure is called.  The
!!      remaining procedures permute the data they read accordingly.
!!
!!    CALL RESTART_SIDE_SETS () reads the side set data, if any.
!!      Called from MESH_GEN.
!!
!!    CALL RESTART_MATLZONE () reads zone structure data, and volume fraction
!!      data from which the matl structure is defined.  Called from INITIAL.
!!
!!    CALL RESTART_SOLID_MECHANICS () reads (or skips) the solid mechanics
!!      data, if any.  Called from INITIAL.
!!
!!    CALL RESTART_SPECIES ([PHI]) reads the species data and returns the
!!      data in PHI if PHI is present.  Otherwise if PHI is not present it
!!      skips the species data, if any.  PHI is a rank-2 real array whose
!!      extent in the first dimension is the number of cells, and extent in
!!      the second dimension is the number of species.  It is an error for
!!      PHI to be present if the restart file does not contain species data,
!!      or if the extents of PHI don't match the restart data.
!!
!!    CALL RESTART_JOULE_HEAT (DEFINED) reads (or skips) the joule heat data,
!!      if any.  DEFINED returns true if the Joule heat data was read;
!!      otherwise it returns false.
!!
!!    CALL RESTART_USTRUC () reads (or skips) the microstructure analysis
!!      checkpoint data, if any.
!!
!!    CALL CLOSE_RESTART_FILE () closes the restart file.  Called from SETUP.
!!

#include "f90_assert.fpp"

module restart_driver

  use restart_variables

  implicit none
  private

  public :: open_restart_file, skip_restart_mesh, skip_restart_side_sets, restart_matlzone, &
            restart_solid_mechanics, restart_species, restart_joule_heat, restart_ustruc, &
            close_restart_file

  !! Private data; defined by OPEN_RESTART_FILE.
  integer, save :: unit, version

contains

  subroutine open_restart_file ()

    use parallel_info_module, only: p_info
    use pgslib_module, only: pgslib_bcast
    use restart_utilities, only: read_var, info, halt
    use string_utilities, only: i_to_c

    integer :: j, n, ios
    character(len=8) :: fileid
    character(len=32) :: feature
    character(len=1024) :: string

    call info (' Reading restart file ' // trim(restart_file) // ' ...')

    !! Open the restart file.
    unit = -1
    if (p_info%IOP) open(newunit=unit,file=restart_file,form='unformatted',action='read',status='old',position='rewind',iostat=ios)
    call pgslib_bcast (ios)
    if (ios /= 0) call halt ('Unable to open restart file ' // trim(restart_file) // &
                             ' for unformatted reading: iostat=' // i_to_c(ios))

    !! Read and check the initial FILEID record.
    call read_var (unit, fileid, 'OPEN_RESTART_FILE: error reading the initial record!')
    read(fileid(5:8),fmt=*,iostat=ios) version
    if (ios /= 0 .or. fileid(1:4) /= 'TRF-') call halt ('OPEN_RESTART_FILE: Not a valid restart file.')

    !! Check the version number.
    if (all(version /= (/2,3/))) call halt ('OPEN_RESTART_FILE: unsupported restart version: ' // i_to_c(version))

    !! Read the feature strings.
    call read_var (unit, n, 'OPEN_RESTART_FILE: error reading NFEAT record')
    do j = 1, n
      call read_var (unit, feature, 'OPEN_RESTART_FILE: error reading FEATURE record')
      select case (feature)
      case ('fluid_flow')
        have_fluid_flow_data = .true.
      case ('joule_heat')
        have_joule_heat_data = .true.
      case ('solid_mechanics')
        have_solid_mechanics_data = .true.
      case ('temperature')
        !! NNC, 18 Jan 2006.
        !have_temperature_data = .true.
        !! The restart file always has both temperature and enthalpy field data.
        !! The issue is that they may not be compatible because of a mapped restart.
        !! The presence (absence) of the 'temperature' feature simply indicates a
        !! normal (mapped) restart file.  Truchas now ignores this flag; one of the
        !! fields is always computed from the other.
        !! FIXME: this feature should be dropped from the restart file format or
        !! renamed if it is of interest to distinguish a mapped restart from a normal.
      case ('species')
        have_species_data = .true.
      case ('microstructure')
        have_microstructure_data = .true.
      case default
        call halt ('OPEN_RESTART_FILE: unknown feature: ' // trim(feature))
      end select
    end do

    !! Read and echo the simulation specification strings.
    call info ('  Restart file simulation specification strings:')
    call read_var (unit, n, 'OPEN_RESTART_FILE: error reading NSS record')
    do j = 1, n
      call read_var (unit, string, 'OPEN_RESTART_FILE: error reading SS record')
      call info ('   ' // trim(string))
    end do

    !! Read the remaining global values.
    call read_var (unit, restart_t, 'OPEN_RESTART_FILE: error reading T record')
    call read_var (unit, restart_dt, 'OPEN_RESTART_FILE: error reading DT record')
    call read_var (unit, restart_cycle_number, 'OPEN_RESTART_FILE: error reading CYCLE record')
    call read_var (unit, restart_ncells, 'OPEN_RESTART_FILE: error reading NCELLS record')
    call read_var (unit, restart_nnodes, 'OPEN_RESTART_FILE: error reading NNODES record')

  end subroutine open_restart_file

  !!
  !! Read the mesh data segment: initializes certain components of the module
  !! structures MESH and VERTEX defined in MESH_MODULE, and the module array
  !! MESH_FACE_SET defined in BC_DATA_MODULE.
  !!

  subroutine skip_restart_mesh ()
    use restart_utilities, only: read_var, skip_records
    integer :: n
    call skip_records (unit, 8, 'SKIP_RESTART_MESH: error skipping VERTEX records')
    call read_var (unit, n, 'SKIP_RESTART_MESH: error reading NCBLOCK record')
    if (n /= 0) call skip_records (unit, 1, 'SKIP_RESTART_MESH: error skipping CBLOCKID record')
    call skip_records (unit, 3, 'SKIP_RESTART_MESH: error skipping COORD records')
  end subroutine skip_restart_mesh

  subroutine skip_restart_side_sets ()
    use bc_data_module, only: skip_side_set_data
    !! Starting with version 3, this info is no longer present.
    !! call read_side_set_data (unit, version)
    !! Skip this data in older version files.
    call skip_side_set_data (unit, version)
  end subroutine skip_restart_side_sets

  !!
  !! Read the core data segment: initializes the module structures ZONE and
  !! MATL defined in ZONE_MODULE and MATL_MODULE, respectively.
  !!

  subroutine restart_matlzone ()
    use zone_module, only: read_zone_data
    use matl_utilities, only: read_matl_data
    use fluid_data_module, only: read_flow_data
    call read_zone_data (unit, version)
    call read_flow_data (unit, version)
    call read_matl_data (unit, version)
  end subroutine restart_matlzone

  !!
  !! Read the optional solid mechanics data segment: initializes the module
  !! structures SMECH_IP, SMECH_CELL, RHS, THERMAL_STRAIN, PC_STRAIN, and
  !! DISPLACEMENT defined in SOLID_MECHANICS_DATA.
  !!

  subroutine restart_solid_mechanics ()
    use solid_mechanics_output, only: read_SM_data, skip_SM_data
    use solid_mechanics_input, only: solid_mechanics
    if (have_solid_mechanics_data) then
      if (solid_mechanics .and. .not.ignore_solid_mechanics) then
        call read_SM_data (unit, version)
      else
        call skip_SM_data (unit, version)
      end if
    end if
  end subroutine restart_solid_mechanics

  !!
  !! Read (or skip) the species concentration data.
  !!

  subroutine restart_species (phi, found)

    use kinds, only: r8
    use restart_utilities, only: read_var, read_dist_array, skip_records, halt
    use legacy_mesh_api, only: pcell => unpermute_mesh_vector
    
    real(r8), intent(out), optional :: phi(:,:)
    logical,  intent(out), optional :: found

    integer :: i, n, m
    character(len=80) :: errmsg

    if (present(found)) found = have_species_data

    if (have_species_data) then
      call read_var (unit, n, 'RESTART_SPECIES: error reading NSPECIES record')
      if (present(phi)) then
        m = size(phi,dim=2)
        if (n /= m) then
          write(errmsg,'(a,i0,a,i0)') 'RESTART_SPECIES: expecting ', m, ' species, found ', n
          call halt (errmsg)
        end if
        INSIST(size(phi,dim=1) == restart_ncells)
        do i = 1, n
          call read_dist_array (unit, phi(:,i), pcell, 'RESTART_SPECIES: error reading PHI records')
        end do
      else
        call skip_records (unit, n, 'RESTART_SPECIES: error skipping the species data')
      end if
    else if (present(phi) .and. .not.present(found)) then
      call halt ('RESTART_SPECIES: restart file does not contain the requested species data')
    end if

  end subroutine restart_species

  !!
  !! Read the optional microstructure data segment
  !!

  subroutine restart_ustruc
    use restart_utilities, only: info
    use ustruc_driver, only: ustruc_enabled, ustruc_read_checkpoint, ustruc_skip_checkpoint
    if (have_microstructure_data) then
      if (ustruc_enabled())  then
        call ustruc_read_checkpoint (unit)
      else
        call ustruc_skip_checkpoint (unit)
      end if
    else if (ustruc_enabled()) then
      call info ('  WARNING: restart file contains no microstructure state data')
    end if
  end subroutine restart_ustruc

  !!
  !! Read the optional Joule heat data segment.
  !!

  subroutine restart_joule_heat (defined)
    use EM_data_proxy, only: EM_is_on, read_joule_data, skip_joule_data
    logical, intent(out) :: defined
    defined = .false.
    if (have_joule_heat_data) then
      if (EM_is_on() .and. .not.ignore_joule_heat) then
        call read_joule_data (unit, version)
        defined = .true.
      else
        call skip_joule_data (unit, version)
      end if
    end if
  end subroutine restart_joule_heat

  subroutine close_restart_file ()
    use restart_utilities, only: info
    use parallel_info_module, only: p_info
    if (p_info%IOP) close(unit)
    call info ('')
    call info ('Closing restart file ' // trim(restart_file))
  end subroutine close_restart_file

end module restart_driver
