!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!
!!  The EM_input Module
!!
!!    Robert Ferrell <ferrell@diablotech.com>, original version
!!    Neil N. Carlson <nnc@newmexico.com>
!!
!!  Input parameters for electromagnetic simulation.
!!
!!  The module provides the following procedure:
!!
!!    * call read_EM_input ()
!!
!!      Reads the electromagnetics namelist, checks the validity of the
!!      input, and broadcasts the namelist variables.
!!
!!  NB: Only EM_data_proxy is allowed direct access to the data stored in
!!  this module.  All other direct access is strictly forbidden!  Access
!!  must go through EM_data_proxy.
!!

#include "f90_assert.fpp"

module EM_input

  use kinds, only: r8
  use parameter_module, only: string_len, MAXSV
  use string_utilities, only: i_to_c
  use truchas_logging_services
  implicit none
  private

  public :: read_EM_input

  !! Domain type: 'FULL_CYLINDER', 'HALF_CYLINDER', or 'QUARTER_CYLINDER'
  character(string_len), public :: EM_Domain_Type
  character, public :: Symmetry_Axis

  !! Container for the INDUCTION_COIL namelist data.
  type, public :: coil_data
    real(r8) :: center(3)
    real(r8) :: radius
    real(r8) :: length
    integer  :: nturns
    real(r8), allocatable :: current(:)
    !real(r8), pointer :: frequency(:) => null()
    !real(r8), pointer :: times(:)     => null()
  end type coil_data

  !! Magnetic source field parameters

  real(r8), allocatable, public :: src_time(:)
  real(r8), allocatable, public :: src_freq(:)
  real(r8), allocatable, public :: unif_src(:)
  type(coil_data), pointer, public, save :: coil_array(:) => null()

  !! EM solver control parameters
  integer,  public :: Steps_Per_Cycle
  integer,  public :: Maximum_Source_Cycles
  real(r8), public :: SS_Stopping_Tolerance
  integer,  public :: Maximum_CG_Iterations
  real(r8), public :: CG_Stopping_Tolerance
  real(r8), public :: Num_Etasq
  real(r8), public :: Material_Change_Threshold

  !! EM output control parameters
  integer,  public :: Output_Level
  logical,  public :: Graphics_Output
  real(r8), public :: Probe_Points(3,10)
  integer,  public :: num_probes

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ THE ELECTROMAGNETIC NAMELIST, CHECK INPUT, AND BROADCAST
 !!

  subroutine read_EM_input (lun)

    use mesh_manager, only: enable_mesh
    integer, intent(in) :: lun
    logical :: exists

    call read_electromagnetics_namelist(lun)
    call read_induction_coil_namelists(lun)

    call enable_mesh('alt', exists)
    if (.not.exists) call TLS_fatal('ALTMESH namelist was not specified')

  end subroutine read_EM_input

  subroutine read_electromagnetics_namelist(lun)

    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use parallel_communication, only: is_IOP, broadcast
    use string_utilities, only: raise_case

    integer, intent(in) :: lun

    integer :: n, ios
    logical :: found, exists
    character(128) :: iom
    real(r8) :: source_times(MAXSV-1)
    real(r8) :: source_frequency(MAXSV)
    real(r8) :: uniform_source(MAXSV)

    namelist /electromagnetics/ EM_Domain_Type, Symmetry_Axis, &
      Source_Times, Source_Frequency, Uniform_Source, &
      Steps_Per_Cycle, Maximum_Source_Cycles, SS_Stopping_Tolerance, &
      Maximum_CG_Iterations, CG_Stopping_Tolerance, Material_Change_Threshold, &
      Num_Etasq, Output_Level, Graphics_Output, Probe_Points

    call TLS_info('Reading ELECTROMAGNETICS namelist ...')

    if (is_IOP) rewind(lun)

    if (is_IOP) call seek_to_namelist(lun, 'electromagnetics', found, iostat=ios)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
    call broadcast(found)
    if (.not.found) call TLS_fatal('ELECTROMAGNETICS namelist not found')

    em_domain_type = NULL_C
    symmetry_axis  = NULL_C

    source_times = NULL_R
    source_frequency = NULL_R
    uniform_source = NULL_R

    steps_per_cycle = NULL_I
    maximum_source_cycles = NULL_I
    ss_stopping_tolerance = NULL_R
    maximum_cg_iterations = NULL_I
    cg_stopping_tolerance = NULL_R
    material_change_threshold = NULL_R
    num_etasq = NULL_R

    output_level = NULL_I
    graphics_output = .false.
    probe_points = NULL_R

    if (is_IOP) read(lun,nml=electromagnetics,iostat=ios,iomsg=iom)
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('error reading ELECTROMAGNETICS namelist: ' // trim(iom))

    call broadcast(em_domain_type)
    call broadcast(symmetry_axis)

    call broadcast(source_times)
    call broadcast(source_frequency)
    call broadcast(uniform_source)

    call broadcast(steps_per_cycle)
    call broadcast(maximum_source_cycles)
    call broadcast(ss_stopping_tolerance)
    call broadcast(maximum_cg_iterations)
    call broadcast(cg_stopping_tolerance)
    call broadcast(material_change_threshold)
    call broadcast(num_etasq)

    call broadcast(output_level)
    call broadcast(graphics_output)
    call broadcast(probe_points)

    em_domain_type = raise_case(trim(adjustl(em_domain_type)))
    select case (em_domain_type)
    case ('FULL_CYLINDER')
    case ('HALF_CYLINDER')
    case ('QUARTER_CYLINDER')
    case ('CYLINDER')
    case ('FRUSTUM')
    case ('VERIFICATION1')
    case (NULL_C)
      call TLS_fatal('EM_DOMAIN_TYPE must be assigned a value')
    case default
      call TLS_fatal('invalid EM_DOMAIN_TYPE: ' // trim(em_domain_type))
    end select

    symmetry_axis = raise_case(trim(adjustl(symmetry_axis)))
    select case (symmetry_axis)
    case ('X', 'Y', 'Z')
    case (NULL_C)
      call TLS_info('  using default value "Z" for SYMMETRY_AXIS')
      symmetry_axis = 'Z'
    case default
      call TLS_fatal('invalid SYMMETRY_AXIS: ' // trim(symmetry_axis))
    end select

    src_time = pack(source_times, mask=(source_times /= NULL_R))
    if (size(src_time) > 1) then
      n = size(src_time)
      if (any(src_time(2:n) <= src_time(:n-1))) &
          call TLS_fatal('SOURCE_TIMES values must be strictly increasing')
    end if

    src_freq = pack(source_frequency, mask=(source_frequency /= NULL_R))
    if (size(src_freq) == 0) then
      call TLS_fatal('SOURCE_FREQUENCY must be assigned a value')
    else if (any(src_freq <= 0.0_r8)) then
      call TLS_fatal('SOURCE_FREQUENCY values must be > 0.0')
    else if (size(src_freq) /= size(src_time) + 1) then
      call TLS_fatal('wrong number of values provided for SOURCE_FREQUENCY')
    end if

    unif_src = pack(uniform_source, mask=(uniform_source /= NULL_R))
    if (size(unif_src) > 0) then
      if (size(unif_src) /= size(src_freq)) then
        call TLS_fatal('wrong number of values provided for UNIFORM_SOURCE')
      end if
    end if

    if (steps_per_cycle == NULL_I) then
      steps_per_cycle = 20
      call TLS_info('  using default value 20 for STEPS_PER_CYCLE')
    else if (steps_per_cycle < 1) then
      call TLS_fatal('STEPS_PER_CYCLE must be > 0')
    else if (steps_per_cycle < 10) then
      call TLS_warn('for decent accuracy, STEPS_PER_CYCLE should be at least 10')
    end if

    if (maximum_source_cycles == NULL_I) then
      maximum_source_cycles = 10
      call TLS_info('  using default value 10 for MAXIMUM_SOURCE_CYCLES')
    else if (maximum_source_cycles < 1) then
      call TLS_fatal('MAXIMUM_SOURCE_CYCLES must be > 0')
    end if

    if (ss_stopping_tolerance == NULL_R) then
      ss_stopping_tolerance = 1.0e-2_r8
      call TLS_info('  using default value 0.01 for SS_STOPPING_TOLERANCE')
    else if (ss_stopping_tolerance <= 0.0_r8) then
      call TLS_fatal('SS_STOPPING_TOLERANCE must be > 0.0')
    else if (ss_stopping_tolerance > 0.1_r8) then
      call TLS_fatal('SS_STOPPING_TOLERANCE is very loose; consider decreasing the value')
    end if

    if (maximum_cg_iterations == NULL_I) then
      call TLS_info('  using default value 500 for MAXIMUM_CG_ITERATIONS')
      maximum_cg_iterations = 500
    else if (maximum_cg_iterations < 1) then
      call TLS_fatal('MAXIMUM_CG_ITERATIONS must be > 0')
    end if

    if (cg_stopping_tolerance == NULL_R) then
      cg_stopping_tolerance = 1.0e-5_r8
      call TLS_info('  using default value 1.0e-5 for CG_STOPPING_TOLERANCE')
    else if (cg_stopping_tolerance <= 0.0_r8 .or. cg_stopping_tolerance >= 0.1_r8) then
      call TLS_fatal('CG_STOPPING_TOLERANCE must be > 0.0 and < 0.1')
    else if (cg_stopping_tolerance > 1.0e-4_r8) then
      call TLS_warn('CG_STOPPING_TOLERANCE is very loose and may lead to the build up of errors.')
    else if (cg_stopping_tolerance < 1.e-3_r8 * epsilon(1.0_r8)) then
      call TLS_warn('CG_STOPPING_TOLERANCE is too tight; CG iterations are unlikely to converge.')
    end if

    if (material_change_threshold == NULL_R) then
      material_change_threshold = 0.3_r8
      call TLS_info('  using default value 0.3 for MATERIAL_CHANGE_THRESHOLD')
    else if (material_change_threshold <= 0.0_r8) then
      call TLS_fatal('MATERIAL_CHANGE_THRESHOLD must be > 0.0')
    end if

    if (num_etasq == NULL_R) then
      num_etasq = 0.0_r8
    else if (num_etasq < 0.0_r8) then
      call TLS_fatal('NUM_ETASQ must be >= 0.0')
    end if

    if (output_level == NULL_I) then
      output_level = 1
    else if (output_level < 1 .or. output_level > 4) then
      call TLS_fatal('OUTPUT_LEVEL must be >= 1 and <= 4')
    end if

    num_probes = 0
    do n = 1, size(probe_points,dim=2)
      if (any(probe_points(:,n) == NULL_R)) exit
      num_probes = n
    end do

  end subroutine read_electromagnetics_namelist

  subroutine read_induction_coil_namelists(lun)

    use input_utilities, only: seek_to_namelist, NULL_I, NULL_R
    use parallel_communication, only: is_IOP, broadcast

    integer, intent(in)  :: lun

    integer :: n, ios, j
    logical :: found
    character(128) :: iom
    character(:), allocatable :: label

    !! Namelist variables
    integer :: nturns
    real(r8) :: center(3), radius, length, current(MAXSV)
    namelist /induction_coil/ center, radius, length, nturns, current

    type :: list_node
      type(coil_data) :: coil
      type(list_node), pointer :: next => null()
    end type list_node
    type(list_node), pointer :: list, old_list

    call TLS_info('Reading INDUCTION_COIL namelists ...')

    if (is_IOP) rewind(lun)

    list => null()
    n = 0 ! namelist counter
    do ! until all INDUCTION_COIL namelists have been read

      if (is_IOP) call seek_to_namelist(lun, 'induction_coil', found, iostat=ios)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading input file: iostat=' // i_to_c(ios))
      call broadcast(found)
      if (.not.found) exit

      n = n + 1
      label = 'INDUCTION_COIL[' // i_to_c(n) // ']'

      center = NULL_R
      radius = NULL_R
      length = NULL_R
      nturns = NULL_I
      current = NULL_R

      if (is_IOP) read(lun,nml=induction_coil,iostat=ios,iomsg=iom)
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading ' // label // ' namelist: ' // trim(iom))

      call broadcast(center)
      call broadcast(radius)
      call broadcast(length)
      call broadcast(nturns)
      call broadcast(current)

      !! Prepend a new instance of the coil data structure to the list.
      old_list => list
      allocate(list)
      list%next => old_list

      if (any(center == NULL_R)) then
        if (all(center == NULL_R)) then
          call TLS_info('  using default CENTER value (0,0,0) for ' // label)
          center = 0.0_r8
        else
          call TLS_fatal(label // ': CENTER requires 3 values')
        end if
      end if
      list%coil%center = center

      if (nturns == NULL_I) then
        call TLS_fatal(label // ': NTURNS must be assigned a value')
      else if (nturns <= 0) then
        call TLS_fatal(label // ': NTURNS must be > 0')
      end if
      list%coil%nturns = nturns

      if (nturns == 1) then
        if (length /= NULL_R) then
          call TLS_info('  ignoring LENGTH value for ' // label)
        end if
        length = 0.0_r8
      end if

      if (length == NULL_R) then
        call TLS_fatal(label // ': LENGTH must be assigned a value')
      else if (length < 0.0_r8) then
        call TLS_fatal(label // ': LENGTH must be >= 0.0')
      end if
      list%coil%length = length

      if (radius == NULL_R) then
        call TLS_fatal(label // ': RADIUS must be assigned a value')
      else if (radius <= 0.0_r8) then
        call TLS_fatal(label // ': RADIUS must be > 0.0')
      end if
      list%coil%radius = radius

      list%coil%current = pack(current, mask=(current/=NULL_R))
      if (size(list%coil%current) /= size(src_freq)) then
        call TLS_fatal(label // ': inconsistent number of CURRENT values')
      end if
    end do

    select case (n)
    case (0)
      call TLS_info('  none found')
    case (1)
      call TLS_info('  read 1 INDUCTION_COIL namelist')
    case default
      call TLS_info('  read ' // i_to_c(n) // ' INDUCTION_COIL namelists')
    end select

    !! Copy the coil data into the output array and deallocate the list.
    allocate(coil_array(n))
    do j = n, 1, -1
      coil_array(j) = list%coil
      old_list => list%next
      deallocate(list)
      list => old_list
    end do

  end subroutine read_induction_coil_namelists

end module EM_input
