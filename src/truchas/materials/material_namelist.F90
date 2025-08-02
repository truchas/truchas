!!
!! MATERIAL_NAMELIST
!!
!! This provides a subroutine that reads the MATERIAL and PHASE namelists from
!! the Truchas input file. The data read from the namelists is assembled into
!! a parameter list that provides the actual input to the code. This parameter
!! list is accessible as a public module variable for later consumption.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2019
!!
!! NB: The long-term plan is to replace namelist-based input with JSON-format
!! parameter list input. This module provides a temporary bridge until that
!! is realized.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module material_namelist

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use input_utilities
  use string_utilities, only: i_to_c
  use parameter_list_type
  use scalar_func_table
  use truchas_logging_services
  implicit none
  private

  public :: read_material_namelists

  type(parameter_list), public :: params  ! the output of read_material_namelists

  integer, parameter :: MAXPHASES=16
  integer, parameter :: MAX_SPECIES = 9

  type(parameter_list) :: phase_list, pc_list  ! temporaries

contains

  subroutine read_material_namelists(lun)

    integer, intent(in) :: lun

    !! MATERIAL namelist specific variables
    character(32) :: name, phases(MAXPHASES)
    namelist /material/ name, phases

    !! Additional property variables (almost the same as for PHASE)
    logical :: is_fluid
    real(r8) :: density, specific_heat, ref_temp, ref_enthalpy, conductivity, viscosity, &
        thermal_expan_coef, expan_ref_temp, electrical_conductivity, electric_susceptibility, &
        magnetic_susceptibility, diffusivity(MAX_SPECIES), soret_coef(MAX_SPECIES), &
        tm_ref_density, tm_ref_temp, tm_linear_cte, tm_lame1, tm_lame2
    character(32) :: specific_enthalpy_func, specific_heat_func, conductivity_func, &
        density_delta_func, viscosity_func, electrical_conductivity_func, &
        electric_susceptibility_func, magnetic_susceptibility_func, &
        diffusivity_func(MAX_SPECIES), soret_coef_func(MAX_SPECIES), &
        tm_linear_cte_func, tm_lame1_func, tm_lame2_func
    namelist /material/ is_fluid, density, ref_temp, ref_enthalpy, specific_heat, &
        specific_heat_func, specific_enthalpy_func, conductivity, conductivity_func, &
        density_delta_func, thermal_expan_coef, expan_ref_temp, viscosity, viscosity_func, &
        electrical_conductivity, electrical_conductivity_func, &
        electric_susceptibility, electric_susceptibility_func, &
        magnetic_susceptibility, magnetic_susceptibility_func, &
        diffusivity, diffusivity_func, soret_coef, soret_coef_func, &
        tm_ref_density, tm_ref_temp, tm_linear_cte, tm_linear_cte_func, &
        tm_lame1, tm_lame1_func, tm_lame2, tm_lame2_func

    logical :: found
    integer :: i, n, ios, num_phases
    character(128) :: iom
    character(:), allocatable :: label
    type(parameter_list), pointer :: matl_plist, matl_phase_plist, plist

    call TLS_info('Reading MATERIAL namelists ...')

    if (is_IOP) rewind(lun)

    n = 0
    do  ! until all MATERIAL namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'MATERIAL', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'MATERIAL[' // i_to_c(n) // ']'

      !! Assign default values before reading the namelist
      name = NULL_C
      phases = NULL_C

      is_fluid = .false.
      density = NULL_R
      specific_heat = NULL_R
      specific_heat_func = NULL_C
      ref_temp = NULL_R
      ref_enthalpy = NULL_R
      specific_enthalpy_func = NULL_C
      conductivity = NULL_R
      conductivity_func = NULL_C
      density_delta_func = NULL_C
      thermal_expan_coef = NULL_R
      expan_ref_temp = NULL_R
      viscosity = NULL_R
      viscosity_func = NULL_C
      electrical_conductivity = NULL_R
      electrical_conductivity_func = NULL_C
      electric_susceptibility = NULL_R
      electric_susceptibility_func = NULL_C
      magnetic_susceptibility = NULL_R
      magnetic_susceptibility_func = NULL_C
      diffusivity = NULL_R
      diffusivity_func = NULL_C
      soret_coef = NULL_R
      soret_coef_func = NULL_C
      tm_ref_density = NULL_R
      tm_ref_temp = NULL_R
      tm_linear_cte = NULL_R
      tm_linear_cte_func = NULL_C
      tm_lame1 = NULL_R
      tm_lame1_func = NULL_C
      tm_lame2 = NULL_R
      tm_lame2_func = NULL_C

      if (is_IOP) read(lun,nml=material,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      !! Broadcast the namelist variables.
      call broadcast(name)
      call broadcast(phases)

      call broadcast(is_fluid)
      call broadcast(density)
      call broadcast(specific_heat)
      call broadcast(specific_heat_func)
      call broadcast(ref_temp)
      call broadcast(ref_enthalpy)
      call broadcast(specific_enthalpy_func)
      call broadcast(conductivity)
      call broadcast(conductivity_func)
      call broadcast(density_delta_func)
      call broadcast(thermal_expan_coef)
      call broadcast(expan_ref_temp)
      call broadcast(viscosity)
      call broadcast(viscosity_func)
      call broadcast(electrical_conductivity)
      call broadcast(electrical_conductivity_func)
      call broadcast(electric_susceptibility)
      call broadcast(electric_susceptibility_func)
      call broadcast(magnetic_susceptibility)
      call broadcast(magnetic_susceptibility_func)
      call broadcast(diffusivity)
      call broadcast(diffusivity_func)
      call broadcast(soret_coef)
      call broadcast(soret_coef_func)
      call broadcast(tm_ref_density)
      call broadcast(tm_ref_temp)
      call broadcast(tm_linear_cte)
      call broadcast(tm_linear_cte_func)
      call broadcast(tm_lame1)
      call broadcast(tm_lame1_func)
      call broadcast(tm_lame2)
      call broadcast(tm_lame2_func)

      !! Check the user-supplied name for the material.
      select case (name)
      case (NULL_C)
        call TLS_fatal(label // ': NAME not specified')
      case ('')
        call TLS_fatal(label // ': invalid 0-length NAME')
      case ('SOLID', 'VOID')
        call TLS_fatal(label // ': invalid NAME; ' // trim(name) // ' is a reserved keyword')
      end select
      if (params%is_sublist(name)) then
        call TLS_fatal(label // ': another MATERIAL namelist has this name: ' // trim(name))
      else
        matl_plist => params%sublist(trim(name))
      end if

      !! Check PHASES for basic correctness.
      num_phases = count(phases /= NULL_C)
      select case (num_phases)
      case (1)
        call TLS_fatal(label // ': PHASES requires two or more values')
      case (2:)
        if (any(phases(:num_phases) == NULL_C)) then
          call TLS_fatal(label // ': phase names must be packed into the PHASES array')
        end if
        !! Check PHASES names for uniqueness.
        matl_phase_plist => matl_plist%sublist('phases')
        do i = 1, num_phases
          if (matl_phase_plist%is_sublist(trim(phases(i)))) then
            call TLS_fatal(label // ': repeated PHASES name: ' // trim(phases(i)))
          else if (phase_list%is_sublist(phases(i))) then
            call TLS_fatal(label // ': PHASES name belongs to another MATERIAL: ' // trim(phases(i)))
          else
            plist => matl_phase_plist%sublist(trim(phases(i))) ! empty placeholder for properties
            plist => phase_list%sublist(phases(i))  ! record the phase name and
            call plist%set('material', name)  ! record the material that defined the phase
          end if
        end do
        !! Record the phase changes we will need data for.
        do i = 1, num_phases-1
          plist => pc_list%sublist(phases(i))
          call plist%set('index', i)
          call plist%set('other-phase', phases(i+1))
          call plist%set('material', name)
        end do
      end select

      !! Process any property variables that were defined.
      plist => matl_plist%sublist('properties')

      !! Boolean properties
      if (is_fluid) call plist%set('fluid', is_fluid)  ! default is F if omitted

      !! Process thermal properties
      call process2(plist, density, NULL_C, 'density', 'density', label)
      call process2(plist, ref_temp, NULL_C, 'ref_temp', 'ref-temp', label)
      call process2(plist, ref_enthalpy, NULL_C, 'ref_enthalpy', 'ref-enthalpy', label)
      call process2(plist, specific_heat, specific_heat_func, 'SPECIFIC_HEAT', 'specific-heat', label)
      call process2(plist, NULL_R, specific_enthalpy_func, 'SPECIFIC_ENTHALPY', 'specific-enthalpy', label)
      call process2(plist, conductivity, conductivity_func, 'CONDUCTIVITY', 'conductivity', label)

      !! Process fluid flow properties
      call process2(plist, NULL_R, density_delta_func, 'DENSITY_DELTA', 'density-delta', label)
      call process2(plist, thermal_expan_coef, NULL_C, 'THERMAL_EXPAN_COEF', 'thermal-expan-coef', label)
      call process2(plist, expan_ref_temp, NULL_C, 'EXPAN_REF_TEMP', 'expan-ref-temp', label)
      call process2(plist, viscosity, viscosity_func, 'VISCOSITY', 'viscosity', label)

      !! Process electromagnetic properties
      call process2(plist, electrical_conductivity, electrical_conductivity_func, 'ELECTRICAL_CONDUCTIVITY', 'electrical-conductivity', label)
      call process2(plist, electric_susceptibility, electric_susceptibility_func, 'ELECTRIC_SUSCEPTIBILITY', 'electric-susceptibility', label)
      call process2(plist, magnetic_susceptibility, magnetic_susceptibility_func, 'MAGNETIC_SUSCEPTIBILITY', 'magnetic-susceptibility', label)

      !! Process solid mechanics variables
      call process2(plist, tm_ref_temp, NULL_C, 'TM_REF_TEMP', 'tm-ref-temp', label)
      call process2(plist, tm_ref_density, NULL_C, 'TM_REF_DENSITY', 'tm-ref-density', label)
      call process2(plist, tm_linear_cte, tm_linear_cte_func, 'TM_LINEAR_CTE', 'tm-linear-cte', label)
      call process2(plist, tm_lame1, tm_lame1_func, 'TM_LAME1', 'tm-lame1', label)
      call process2(plist, tm_lame2, tm_lame2_func, 'TM_LAME2', 'tm-lame2', label)

      !! Process species diffusion variables
      do i = 1, MAX_SPECIES
        call process2i(plist, diffusivity(i), diffusivity_func(i), 'DIFFUSIVITY', 'diffusivity', i, label)
        call process2i(plist, soret_coef(i), soret_coef_func(i), 'SORET_COEF', 'soret-coef', i, label)
      end do

      call TLS_info('  read namelist "' // trim(name) // '"')
    end do

    call read_phase_namelists(lun)
    call read_phase_change_namelists(lun)

  end subroutine read_material_namelists


  subroutine read_phase_namelists(lun)

    integer, intent(in) :: lun

    !! Namelist variables
    logical :: is_fluid
    real(r8) :: specific_heat, conductivity, viscosity, thermal_expan_coef, expan_ref_temp, &
        electrical_conductivity, electric_susceptibility, magnetic_susceptibility, &
        diffusivity(MAX_SPECIES), soret_coef(MAX_SPECIES), tm_ref_density, &
        tm_ref_temp, tm_linear_cte, tm_lame1, tm_lame2
    character(32) :: name, specific_heat_func, specific_enthalpy_func, conductivity_func, &
        density_delta_func, viscosity_func, electrical_conductivity_func, &
        electric_susceptibility_func, magnetic_susceptibility_func, &
        diffusivity_func(MAX_SPECIES), soret_coef_func(MAX_SPECIES), &
        tm_linear_cte_func, tm_lame1_func, tm_lame2_func
    namelist /phase/ name, is_fluid, &
        specific_heat, specific_heat_func, specific_enthalpy_func, conductivity, &
        conductivity_func, density_delta_func, thermal_expan_coef, expan_ref_temp, &
        viscosity, viscosity_func, &
        electrical_conductivity, electrical_conductivity_func, &
        electric_susceptibility, electric_susceptibility_func, &
        magnetic_susceptibility, magnetic_susceptibility_func, &
        diffusivity, diffusivity_func, soret_coef, soret_coef_func, &
        tm_ref_density, tm_ref_temp, tm_linear_cte, tm_linear_cte_func, &
        tm_lame1, tm_lame1_func, tm_lame2, tm_lame2_func

    integer :: n, i, ios
    logical :: found
    character(:), allocatable :: label, matl
    character(128) :: iom
    type(parameter_list), pointer :: matl_plist, plist

    call TLS_info('Reading PHASE namelists ...')

    if (is_IOP) rewind(lun)

    n = 0
    do  ! until all PHASE namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'PHASE', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'PHASE[' // i_to_c(n) // ']'

      !! Assign default values
      name = NULL_C
      is_fluid = .false.
      specific_heat = NULL_R
      specific_heat_func = NULL_C
      specific_enthalpy_func = NULL_C
      conductivity = NULL_R
      conductivity_func = NULL_C
      density_delta_func = NULL_C
      thermal_expan_coef = NULL_R
      expan_ref_temp = NULL_R
      viscosity = NULL_R
      viscosity_func = NULL_C
      electrical_conductivity = NULL_R
      electrical_conductivity_func = NULL_C
      electric_susceptibility = NULL_R
      electric_susceptibility_func = NULL_C
      magnetic_susceptibility = NULL_R
      magnetic_susceptibility_func = NULL_C
      diffusivity = NULL_R
      diffusivity_func = NULL_C
      soret_coef = NULL_R
      soret_coef_func = NULL_C
      tm_ref_density = NULL_R
      tm_ref_temp = NULL_R
      tm_linear_cte = NULL_R
      tm_linear_cte_func = NULL_C
      tm_lame1 = NULL_R
      tm_lame1_func = NULL_C
      tm_lame2 = NULL_R
      tm_lame2_func = NULL_C

      if (is_IOP) read(lun,nml=phase,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(name)
      call broadcast(is_fluid)
      call broadcast(specific_heat)
      call broadcast(specific_heat_func)
      call broadcast(specific_enthalpy_func)
      call broadcast(conductivity)
      call broadcast(conductivity_func)
      call broadcast(density_delta_func)
      call broadcast(thermal_expan_coef)
      call broadcast(expan_ref_temp)
      call broadcast(viscosity)
      call broadcast(viscosity_func)
      call broadcast(electrical_conductivity)
      call broadcast(electrical_conductivity_func)
      call broadcast(electric_susceptibility)
      call broadcast(electric_susceptibility_func)
      call broadcast(magnetic_susceptibility)
      call broadcast(magnetic_susceptibility_func)
      call broadcast(diffusivity)
      call broadcast(diffusivity_func)
      call broadcast(soret_coef)
      call broadcast(soret_coef_func)
      call broadcast(tm_ref_density)
      call broadcast(tm_ref_temp)
      call broadcast(tm_linear_cte)
      call broadcast(tm_linear_cte_func)
      call broadcast(tm_lame1)
      call broadcast(tm_lame1_func)
      call broadcast(tm_lame2)
      call broadcast(tm_lame2_func)

      !! Check the NAME
      if (name == NULL_C) then
        call TLS_fatal(label // ': NAME not specified')
      else if (name == '') then
        call TLS_fatal(label // ': invalid 0-length string for NAME')
      else if (.not.phase_list%is_sublist(name)) then
        call TLS_fatal(label // ': unknown PHASE name: ' // trim(name))
      end if

      !! Get the material parameter list that references this phase
      plist => phase_list%sublist(name)
      call plist%get('material', matl)
      call plist%set('defined', .true.)
      matl_plist => params%sublist(matl)

      !! Get the sublist to populate with the phase properties.
      plist => matl_plist%sublist('phases')
      plist => plist%sublist(name)

      !! Boolean properties
      if (is_fluid) call plist%set('fluid', is_fluid)  ! default is F if omitted

      !! Process thermal properties
      call process2(plist, specific_heat, specific_heat_func, 'SPECIFIC_HEAT', 'specific-heat', label)
      call process2(plist, NULL_R, specific_enthalpy_func, 'SPECIFIC_ENTHALPY', 'specific-enthalpy', label)
      call process2(plist, conductivity, conductivity_func, 'CONDUCTIVITY', 'conductivity', label)

      !! Process fluid flow properties
      call process2(plist, NULL_R, density_delta_func, 'DENSITY_DELTA', 'density-delta', label)
      call process2(plist, thermal_expan_coef, NULL_C, 'THERMAL_EXPAN_COEF', 'thermal-expan-coef', label)
      call process2(plist, expan_ref_temp, NULL_C, 'EXPAN_REF_TEMP', 'expan-ref-temp', label)
      call process2(plist, viscosity, viscosity_func, 'VISCOSITY', 'viscosity', label)

      !! Process electromagnetic properties
      call process2(plist, electrical_conductivity, electrical_conductivity_func, 'ELECTRICAL_CONDUCTIVITY', 'electrical-conductivity', label)
      call process2(plist, electric_susceptibility, electric_susceptibility_func, 'ELECTRIC_SUSCEPTIBILITY', 'electric-susceptibility', label)
      call process2(plist, magnetic_susceptibility, magnetic_susceptibility_func, 'MAGNETIC_SUSCEPTIBILITY', 'magnetic-susceptibility', label)

      !! Process solid mechanics variables
      call process2(plist, tm_ref_temp, NULL_C, 'TM_REF_TEMP', 'tm-ref-temp', label)
      call process2(plist, tm_ref_density, NULL_C, 'TM_REF_DENSITY', 'tm-ref-density', label)
      call process2(plist, tm_linear_cte, tm_linear_cte_func, 'TM_LINEAR_CTE', 'tm-linear-cte', label)
      call process2(plist, tm_lame1, tm_lame1_func, 'TM_LAME1', 'tm-lame1', label)
      call process2(plist, tm_lame2, tm_lame2_func, 'TM_LAME2', 'tm-lame2', label)

      !! Process species diffusion variables
      do i = 1, MAX_SPECIES
        call process2i(plist, diffusivity(i), diffusivity_func(i), 'DIFFUSIVITY', 'diffusivity', i, label)
        call process2i(plist, soret_coef(i), soret_coef_func(i), 'SORET_COEF', 'soret-coef', i, label)
      end do

      call TLS_info('  read namelist "' // trim(name) // '"')
    end do

    if (n == 0) call TLS_info('  none found')

  end subroutine read_phase_namelists

  subroutine read_phase_change_namelists(lun)

    integer, intent(in) :: lun

    !! Namelist variables
    character(32) :: low_temp_phase, high_temp_phase
    real(r8) :: latent_heat, solidus_temp, liquidus_temp, solid_frac_table(2,101)
    namelist /phase_change/ low_temp_phase, high_temp_phase, latent_heat, &
        solidus_temp, liquidus_temp, solid_frac_table

    integer :: n, ios, j
    logical :: found
    character(:), allocatable :: label, name, matl
    character(128) :: iom
    type(parameter_list), pointer :: matl_plist, plist

    call TLS_info('Reading PHASE_CHANGE namelists ...')

    if (is_IOP) rewind(lun)

    n = 0
    do  ! until all PHASE_CHANGE namelists have been read or an error occurs

      if (is_IOP) call seek_to_namelist(lun, 'PHASE_CHANGE', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'PHASE_CHANGE[' // i_to_c(n) // ']'

      !! Assign default values
      low_temp_phase = NULL_C
      high_temp_phase = NULL_C
      latent_heat = NULL_R
      solidus_temp = NULL_R
      liquidus_temp = NULL_R
      solid_frac_table = NULL_R

      if (is_IOP) read(lun,nml=phase_change,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(low_temp_phase)
      call broadcast(high_temp_phase)
      call broadcast(latent_heat)
      call broadcast(solidus_temp)
      call broadcast(liquidus_temp)
      call broadcast(solid_frac_table)

      !! Check LOW_TEMP_PHASE
      if (low_temp_phase == NULL_C) then
        call TLS_fatal(label // ': LOW_TEMP_PHASE not specified')
      else if (low_temp_phase == '') then
        call TLS_fatal(label // ': invalid 0-length string for LOW_TEMP_PHASE')
      else if (.not.phase_list%is_sublist(low_temp_phase)) then
        call TLS_fatal(label // ': LOW_TEMP_PHASE: unknown phase name: ' // trim(low_temp_phase))
      end if

      !! Check HIGH_TEMP_PHASE
      if (high_temp_phase == NULL_C) then
        call TLS_fatal(label // ': HIGH_TEMP_PHASE not specified')
      else if (high_temp_phase == '') then
        call TLS_fatal(label // ': invalid 0-length string for HIGH_TEMP_PHASE')
      else if (.not.phase_list%is_sublist(high_temp_phase)) then
        call TLS_fatal(label // ': HIGH_TEMP_PHASE: unknown PHASE name: ' // trim(high_temp_phase))
      end if

      !! Get the material parameter list that references this phase
      plist => pc_list%sublist(low_temp_phase)
      call plist%get('other-phase', name, default='')
      if (name /= high_temp_phase) call TLS_fatal(label // ': invalid phase change')
      call plist%get('material', matl)
      call plist%set('defined', .true.)
      matl_plist => params%sublist(matl)

      plist => matl_plist%sublist('phase-changes')
      plist => plist%sublist(trim(low_temp_phase) // ':' // trim(high_temp_phase))

      if (any(solid_frac_table /= NULL_R)) then

        do j = size(solid_frac_table,dim=2), 1, -1
          if (solid_frac_table(1,j) /= NULL_R) exit
        end do

        if (any(solid_frac_table(:,:j) == NULL_R) .or. &
            any(solid_frac_table(:,j+1:) /= NULL_R) .or. j < 2) &
            call TLS_fatal('invalid table format')

        call plist%set('solid-frac-table', solid_frac_table(:,:j))

      else

        if (solidus_temp == NULL_R) then
          call TLS_fatal(label // ': SOLIDUS_TEMP not specified')
        else
          call plist%set('solidus-temp', solidus_temp)
        end if

        if (liquidus_temp == NULL_R) then
          call TLS_fatal(label // ': LIQUIDUS_TEMP not specified')
	else if (liquidus_temp <= solidus_temp) then
	  call TLS_fatal(label // ': LIQUIDUS_TEMP must be > SOLIDUS_TEMP')
        else
          call plist%set('liquidus-temp', liquidus_temp)
        end if

      end if

      if (latent_heat /= NULL_R) then
        if (latent_heat < 0) call TLS_fatal(label // ': LATENT_HEAT must be >= 0.0')
        call plist%set('latent-heat', latent_heat)
      end if

    end do

    select case (n)
    case (0)
      call TLS_info('  none found')
    case (1)
      call TLS_info('  read 1 PHASE_CHANGE namelist')
    case default
      call TLS_info('  read ' // i_to_c(n) // ' PHASE_CHANGE namelists')
    end select

    block
#ifdef NAG_BUG20220103
      use parameter_list_type ! to get the parameter_list_iterator constructor
#endif
      type(parameter_list_iterator) :: piter
      character(:), allocatable :: name
      piter = parameter_list_iterator(pc_list, sublists_only=.true.)
      do while (.not.piter%at_end())
        plist => piter%sublist()
        call plist%get('defined', found, default=.false.)
        if (.not.found) then
          call plist%get('other-phase', name)
          call TLS_fatal('missing PHASE_CHANGE namelist for ' // piter%name() &
              // ' / ' // name)
        end if
        call piter%next
      end do
    end block

  end subroutine read_phase_change_namelists

  !! This auxiliary subroutine checks that at most one of CONST and FNAME has
  !! been assigned a value, and if one has, sets the PLIST parameter PNAME to
  !! that value. If any errors are encountered a fatal error message is written.
  !! BNAME is the base name of the corresponding namelist variables, and LABEL
  !! is a prefix to use with the error message.

  subroutine process2(plist, const, fname, bname, pname, label)
    type(parameter_list), intent(inout) :: plist
    real(r8),     intent(in) :: const ! possible constant value
    character(*), intent(in) :: fname ! possible function name value
    character(*), intent(in) :: bname ! namelist variable base name
    character(*), intent(in) :: pname ! parameter list variable name
    character(*), intent(in) :: label
    if (const /= NULL_R .and. fname /= NULL_C) then
      call TLS_fatal(label // ': cannot specify both' // bname // ' and ' // bname // '_FUNC')
    else if (const /= NULL_R) then
      call plist%set(pname, const)
    else if (fname /= NULL_C) then
      if (known_func(fname)) then
        call plist%set(pname, trim(fname))
      else
        call TLS_fatal(label // ': unknown function for ' // bname // '_FUNC: ' // trim(fname))
      end if
    end if
  end subroutine

  !! This auxiliary subroutine checks that at most one of CONST and FNAME has
  !! been assigned a value, and if one has, sets a PLIST parameter to that
  !! value. In this version the actual arguments for CONST and FNAME are
  !! supposed to be the INDEX element of arrays. The PLIST parameter name is
  !! PNAME appended with INDEX as a string. If any errors are encountered a
  !! fatal error message is written. BNAME is the base name of the
  !! corresponding namelist variables, and LABEL is a prefix to use with the
  !! error message.

  subroutine process2i(plist, const, fname, bname, pname, index, label)
    type(parameter_list), intent(inout) :: plist
    real(r8),     intent(in) :: const ! possible constant value
    character(*), intent(in) :: fname ! possible function name value
    character(*), intent(in) :: bname ! namelist variable base name
    character(*), intent(in) :: pname ! parameter list variable name
    integer,      intent(in) :: index ! array index
    character(*), intent(in) :: label
    if (const /= NULL_R .and. fname /= NULL_C) then
      call TLS_fatal(label // ': cannot specify both ' // bname //  '(' // &
          i_to_c(index) // ') and ' // bname // '_FUNC(' // i_to_c(index) // ')')
    else if (const /= NULL_R) then
      call plist%set(pname // i_to_c(index), const)
    else if (fname /= NULL_C) then
      if (known_func(fname)) then
        call plist%set(pname // i_to_c(index), trim(fname))
      else
        call TLS_fatal(label // ': unknown function for ' // bname // &
            '_FUNC(' // i_to_c(index) // '): ' // trim(fname))
      end if
    end if
  end subroutine

end module material_namelist
