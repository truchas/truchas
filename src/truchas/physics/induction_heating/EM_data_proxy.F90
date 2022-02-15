!!
!! EM_DATA_PROXY
!!
!! This module provides procedures and persistent data for exchanging data
!! between the electromagnetic solver and the core of Truchas.  It also
!! provides accessors to the electromagnetics and induction_coil namelist
!! input data.
!!
!! Neil N. Carlson <nnc@newmexico.com>
!! Last revised 26 May 2004
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! I. Data Proxy Functions.  These procedures are involved in exchanging data
!!  arrays between the hex-mesh-based core of Truchas and the tet-mesh-based
!!  EM solver.  These are all global procedures and none should be called
!!  unless the EM simulation has been enabled.
!!
!!  CALL INIT_EM_DATA_PROXY () initializes the module, enabling the data store/
!!    retrieval procedures that follow.  Specifically it allocates private
!!    module arrays used to cache the data, and it invokes the computation of
!!    the grid-to-grid mappings used in the data exchange.  This must not be
!!    called until the meshes have been completely initialized.
!!
!!  CALL SET_PERMITTIVITY (VALUES)
!!  CALL SET_PERMEABILITY (VALUES)
!!  CALL SET_CONDUCTIVITY (VALUES)
!!
!!    These routines store the so-named EM material parameter data arrays for
!!    later retrieval by the EM solver.  VALUES is a rank-1 INTENT(IN) real
!!    array whose size equals the number of hex-mesh cells on this processor;
!!    VALUES(j) is the value of the parameter on local hex j.  Internally,
!!    these values are mapped to cell-based values on the tet-mesh, and it is
!!    these values that are actually cached.
!!
!!  PERMITTIVITY()
!!  PERMEABILITY()
!!  CONDUCTIVITY()
!!
!!    These functions return a pointer to a rank-1 real array containing the
!!    values of the so-named EM parameter on the cells of the tet mesh used
!!    by the EM solver.  The size of the result equals the number of tet cells
!!    on this processor.  The result can be used in an expression or as the
!!    target of a pointer assignment.  However, the target of the result should
!!    never be deallocated.  The corresponding 'set' routine must have been
!!    previously called.
!!  
!!  CALL SET_JOULE_POWER_DENSITY (VALUES) stores the values of the Joule power
!!    density on the tet-mesh cells, which is computed by the EM solver, for
!!    later retrieval by the hex-mesh side of Truchas.  VALUES is a rank-1
!!    INTENT(IN) real array whose size equals the number of tet-mesh cells on
!!    this processor; VALUES(j) is the average Joule power density on cell j.
!!    Internally, these values are mapped to cell-based values on the hex-mesh,
!!    and it is these values that are actually cached.
!!
!!  JOULE_POWER_DENSITY() returns a pointer to a rank-1 real array containing the
!!    values of the Joule power density on the cells of the hex mesh.  The size
!!    of the result equals the number of hex cells on this processor.  The result
!!    can be used in an expression or as the target of a pointer assignment.
!!    However, the target of the result should never be deallocated.  If called
!!    before the corresponding 'set' routine, the function returns an array of
!!    zero values.
!!
!! II. Input Data Accessors.  The following functions return the value of the
!!  so-named input parameter.  See the electromagnetics namelist section of the
!!  Reference Manual for a description of the parameters.  These are local
!!  procedures, however each processor will have the same value for parameters.
!!  None of these procedures should be used unless the EM simulation has been
!!  enabled.
!!   
!!  Real-valued functions:
!!    GET_EPSILON_0(), GET_MU_0(), GET_SS_STOPPING_TOLERANCE(),
!!    GET_CG_STOPPING_TOLERANCE(), GET_NUM_ETASQ()
!!
!!  Integer-valued functions:
!!    GET_MAXIMUM_CG_ITERATIONS(), GET_MAXIMUM_SOURCE_CYCLES(),
!!    GET_NUM_PROBE_POINTS(), GET_OUTPUT_LEVEL(), GET_STEPS_PER_CYCLE()
!!
!!  Character-valued functions:
!!    GET_EM_DOMAIN_TYPE()
!!
!!  Logical-valued functions:
!!    GRAPHICS() returns the value of Graphics_Output.
!!
!!  Real array-valued function:
!!    GET_PROBE_POINTS() returns a (3,GET_NUM_PROBE_POINTS()) array
!!
!! III. Source Field Data Procedures.
!!
!!  INDUCTION_COILS() returns a pointer to an array of TYPE(SOLENOID) which is
!!    a container for the input values describing the physical characteristics
!!    of the driving coil, and the current amplitude in the coil.
!!
!!    TYPE(SOLENOID) has the following components:
!!      o CENTER(3) is the position of the center of the coil;
!!      o RADIUS is the radius of the coil;
!!      o LENGTH is the length of the coil;
!!      o NTURNS is the number of turns of the coil.
!!      o CURRENT is the amplitude of the (total) current flowing in the coil.
!!    A physical coil is represented numerically as an equi-spaced array of
!!    circular current loops, and not as a helical current line.
!!
!!  SOURCE_FREQUENCY() returns the (linear) frequency of the sinusoidally
!!    varying current.  All coils operate at the same frequency (and phase).
!!    The value returned is set by a call to SET_SOURCE_PROPERTIES.
!!
!!  UNIFORM_SOURCE() returns the magnitude of the uniform magnetic
!!    source field. The value returned is set by a call to SET_SOURCE_PROPERTIES.
!!
!!  CALL SET_SOURCE_PROPERTIES(T) sets the values returned by SOURCE_FREQUENCY
!!    and UNIFORM_SOURCE, and the CURRENT components of the array returned by
!!    INDUCTION_COILS, to the appropriate values for time T.  These values --
!!    individual coil currents, current frequency, and strength of an optional
!!    uniform source field -- may be described as time-dependent functions in
!!    the input file.  It is expected that after calling this procedure and
!!    setting these values, that the Joule heat will be recalculated (using
!!    these values) and stored via a call to SET_JOULE_POWER_DENSITY, bringing
!!    the result of JOULE_POWER_DENSITY current
!!
!!  SOURCE_HAS_CHANGED() returns the value true if any of the magnetic source
!!    field properties returned by SOURCE_FREQUENCY, UNIFORM_SOURCE, and the
!!    CURRENT components of the array returned by INDUNCTION_COIL differ from
!!    those values used to compute the Joule heat currently being returned by
!!    JOULE_POWER_DENSITY.  Its purpose is to determine whether the Joule heat
!!    ought to be recomputed.
!!
!!  MATERIAL_HAS_CHANGED() returns the value true if the material parameter
!!    arrays returned by CONDUCTIVITY and PERMEABILITY differ significantly
!!    from the values used to compute the Joule heat currently being returned
!!    by JOULE_POWER_DENSITY.  The maximum relative change is taken as the
!!    difference measure, and when this difference exceeds the value of the
!!    input variable MATERIAL_CHANGE_THRESHOLD, the difference is considered
!!    significant.  For conductivity, only the conducting region (where the
!!    value is positive) is considered when computing the difference, and an
!!    underlying assumption is that this region remains fixed throughout the
!!    simulation.
!!
!!  NO_SOURCE_FIELD() returns the value true if the source field parameters as
!!    returned by UNIFORM_SOURCE and INDUCTION_COILS indicate that there is no
!!    magnetic source field present.  Otherwise it returns the value false.
!!    This would be normally called after a call to SET_SOURCE_PROPERTIES to
!!    determine if it were necessary to actually calculate the Joule heat. If
!!    there is no source field, there is no Joule heat, and the following
!!    procedure can be called to set it accordingly.
!!
!!  CALL ZERO_JOULE_POWER_DENSITY() sets the values returned by the function
!!    JOULE_POWER_DENSITY to zero.  Typically this would be called if
!!    NO_SOURCE_FIELD returns the value true, instead of invoking a procedure
!!    to actually compute the Joule heat.
!!
!! IV. Miscellaneous.
!!
!!  SET_EM_SIMULATION_ON_OR_OFF(ON) enables the electromagnetic simulation if
!!  the INTENT(IN) logical argument ON has the value true; otherwise the EM
!!  simulation is disabled.  The enabled/disabled state of the EM simulation
!!  is described by a private module variable.  The default state is disabled.
!!  This is typically called when the physics namelist is read.  This is a
!!  global procedure; only the value of ON on the IO processor is relevant
!!
!!  EM_IS_ON() returns the value true if the EM simulation is enabled; it
!!    returns the value false otherwise.  This may be called at any time.
!!    All processors have the same state.
!!
!!  EM_IS_OFF() returns the value false if the EM simulation is enabled; it
!!    returns the value true otherwise.  This may be called at any time.
!!    All processors have the same state.
!!
!!  IN_TIME_INTERVAL(T, N) returns the value true if the time T lies within
!!    the Nth source time interval (1-based numbering); otherwise it returns
!!    the value false.  This is provided for the use of sensitivity analysis.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! IMPLEMENTATION NOTES
!!
!! o  The services in section II and III are completely independent of the data
!!    proxy services in section I, and probably ought to be split off as a
!!    separate module or incorporated into a new 'EM physics' driver module.
!!    (Something like the module EM, but with a refactored COMPUTE_JOULE_HEAT.)
!!    The services of section IV belong there as well.
!!
!! o  The reason the data retrieval functions of section I return a pointer to
!!    internal private module data instead of simply an array, is to avoid the
!!    undesirable copying of data, although this will permit the malicous
!!    corruption of this private data.
!!
  
#include "f90_assert.fpp"

module EM_data_proxy

  use kinds, only : rk => r8
  use data_mapper_class
  use truchas_logging_services
  use truchas_timers
  
  implicit none
  private
  
  !! EM simulation state procedures
  public :: set_EM_simulation_on_or_off, EM_is_on, EM_is_off
  
  !! Initialization procedure
  public :: init_EM_data_proxy
  
  !! Data store/retrieval procedures 
  public :: set_permittivity, permittivity
  public :: set_permeability, permeability
  public :: set_conductivity, conductivity
  public :: set_joule_power_density, joule_power_density
  public :: EM_mesh
  
  !! Accessor functions for the electromagnetic namelist input parameters.
  public :: get_epsilon_0, get_mu_0
  public :: get_EM_domain_type, symmetry_axis
  public :: get_steps_per_cycle, get_maximum_source_cycles, get_ss_stopping_tolerance
  public :: get_maximum_cg_iterations, get_cg_stopping_tolerance
  public :: get_num_etasq, get_output_level, graphics
  public :: get_num_probe_points, get_probe_points
  
  !! Magnetic source field procedures
  public :: source_has_changed, set_source_properties
  public :: source_frequency, uniform_source, induction_coils
  public :: no_source_field, zero_joule_power_density
  public :: source_is_scaled, scale_joule_power_density
  public :: in_time_interval
  public :: time_interval
  public :: material_has_changed
  
  !! Restart procedures
  public :: read_joule_data, skip_joule_data, danu_write_joule
  
  !! Container for the description of a driving coil
  type, public :: solenoid
     real(kind=rk) :: current
     !real(kind=rk) :: frequency
     !real(kind=rk) :: axis(3)
     real(kind=rk) :: center(3)
     real(kind=rk) :: length
     real(kind=rk) :: radius
     integer       :: nturns
  end type solenoid
  
  !! Hex-tet grid mapping data.
  class(data_mapper), allocatable :: ht2em
  
  !! Distributed tet-mesh cell-based vectors: EM material parameters.
  real(kind=rk), allocatable, target, save :: eps(:), mu(:), sigma(:)
  real(kind=rk), allocatable, target, save :: eps_q(:), mu_q(:), sigma_q(:)
  
  !! Distributed hex-mesh cell-based vector: Joule power density.
  real(kind=rk), allocatable, target, save :: joule(:)

  !! Source field parameters.
  real(kind=rk), save :: uhfs   ! Uniform H-field strength returned by UNIFORM_SOURCE, and
  real(kind=rk), save :: uhfs_q !   that value used to compute JOULE.
  real(kind=rk), save :: freq   ! Source frequency returned by SOURCE_FREQUENCY, and
  real(kind=rk), save :: freq_q !   that value used to compute JOULE.
  type(solenoid), pointer, save :: coil(:) => null()  ! Coil properties returned by INDUCTION_COILS,
  type(solenoid), pointer, save :: coil_q(:) => null()!   and those values used to compute JOULE.
  
  !! Source field parameters corresponding to JOULE
  
  !! EM simulation state variable
  logical, save :: EM_enabled = .false.
  
CONTAINS
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SET_EM_SIMULATION_ON_OR_OFF, EM_IS_ON, EM_IS_OFF
 !!
 !! Set and query the state of the EM simulation.  When setting the state,
 !! only the value of the argument on the IO processor is relevant.
 !!
 
  subroutine set_EM_simulation_on_or_off (on)
    use parallel_communication, only : broadcast
    logical, intent(in) :: on
    EM_enabled = on
    call broadcast (EM_enabled)
  end subroutine set_EM_simulation_on_or_off

  logical function EM_is_on ()
    EM_is_on = EM_enabled
  end function EM_is_on

  logical function EM_is_off ()
    EM_is_off = .not.EM_enabled
  end function EM_is_off

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! INIT_EM_DATA_PROXY
 !!
 !! Initializes the module for the data store/retrieval procedures.
 !! Specifically it allocates module arrays used to cache the data and invokes
 !! the computation of the grid-to-grid mappings used in the exchange of data.
 !!
 !! NOTES
 !!
 !! o Ideally the main mesh would be handled exactly the same as the alt mesh:
 !!   the (distributed) mesh would be accessed through the mesh broker, and
 !!   collated components formed as needed by the grid-to-grid mesh data
 !!   structure.  For now we kludge it by 'importing' the collated main mesh
 !!   directly and accessing a core Truchas module to get the number of hex
 !!   cells on this processor.
 !!

  subroutine init_EM_data_proxy ()
  
    use parallel_communication
    use base_mesh_class
    use mesh_manager, only: named_mesh_ptr
    use EM_input, only: coil_array
    use kuprat_mapper_type
#ifdef USE_PORTAGE
    use portage_mapper_type
#endif
    use altmesh_namelist, only: data_mapper_kind
    
    integer :: j
    class(base_mesh), pointer :: ht_mesh, em_mesh
    
    !! Generate the mapping between the HT and EM meshes
    ht_mesh => named_mesh_ptr('main')
    em_mesh => named_mesh_ptr('alt')
    select case (data_mapper_kind)
    case ('portage')
#ifdef USE_PORTAGE
      call TLS_info('  Creating Portage mesh-to-mesh mapper ...')
      allocate(portage_mapper :: ht2em)
#else
      INSIST(.false.)
#endif
    case default
      call TLS_info('  Creating Kuprat mesh-to-mesh mapper ...')
      allocate(kuprat_mapper :: ht2em)
    end select
    call ht2em%init(ht_mesh, em_mesh)
    call TLS_info('    mesh-to-mesh mapper created')
    
    !! Allocate the module's tet-mesh data store arrays.
    if (allocated(eps)) deallocate(eps)
    if (allocated(mu)) deallocate(mu)
    if (allocated(sigma)) deallocate(sigma)
    allocate(eps(em_mesh%ncell), mu(em_mesh%ncell), sigma(em_mesh%ncell))
    
    if (allocated(eps_q)) deallocate(eps_q)
    if (allocated(mu_q)) deallocate(mu_q)
    if (allocated(sigma_q)) deallocate(sigma_q)
    allocate(eps_q(em_mesh%ncell), mu_q(em_mesh%ncell), sigma_q(em_mesh%ncell))

    !! Allocate the module's hex-mesh data store array.
    if (allocated(joule)) deallocate(joule)
    allocate(joule(ht_mesh%ncell_onP))
    
    !! Copy the input coil physical properties to our own data structure.
    allocate(coil(size(coil_array)), coil_q(size(coil_array)))
    do j = 1, size(coil)
      coil(j)%center = coil_array(j)%center
      coil(j)%radius = coil_array(j)%radius
      coil(j)%length = coil_array(j)%length
      coil(j)%nturns = coil_array(j)%nturns
      coil(j)%current = 0.0_rk
    end do
    coil_q = coil

    !! The initial joule heat is set later in INITIALIZE_EM.
    
  end subroutine init_EM_data_proxy

  function EM_mesh () result (ptr)
    use mesh_manager, only: simpl_mesh_ptr
    use simpl_mesh_type
    type(simpl_mesh), pointer :: ptr
    ptr => simpl_mesh_ptr('alt')
  end function EM_mesh
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Truchas-side data store procedures
 !!
 !! Procedures map cell-based, hex-mesh input data to cell-based data on the
 !! tet mesh, storing the results as module data for later retrieval.
 !! 
  
  subroutine set_permittivity (values)
    use simpl_mesh_type
    real(kind=rk), intent(in) :: values(:)
    type(simpl_mesh), pointer :: mesh
    ASSERT( allocated(eps) )
    mesh => EM_mesh()
    call start_timer('mesh-to-mesh mapping')
    call ht2em%map_field(values, eps(:mesh%ncell_onP), defval=1.0_rk, map_type=LOCALLY_BOUNDED)
    call mesh%cell_ip%gather_offp(eps)
    call stop_timer('mesh-to-mesh mapping')
  end subroutine set_permittivity

  subroutine set_permeability (values)
    use simpl_mesh_type
    real(kind=rk), intent(in) :: values(:)
    type(simpl_mesh), pointer :: mesh
    ASSERT( allocated(mu) )
    mesh => EM_mesh()
    call start_timer('mesh-to-mesh mapping')
    call ht2em%map_field(values, mu(:mesh%ncell_onP), defval=1.0_rk, map_type=LOCALLY_BOUNDED)
    call mesh%cell_ip%gather_offp(mu)
    call stop_timer('mesh-to-mesh mapping')
  end subroutine set_permeability

  subroutine set_conductivity (values)
    use simpl_mesh_type
    real(kind=rk), intent(in) :: values(:)
    type(simpl_mesh), pointer :: mesh
    ASSERT( allocated(sigma) )
    mesh => EM_mesh()
    call start_timer('mesh-to-mesh mapping')
    call ht2em%map_field(values, sigma(:mesh%ncell_onP), defval=0.0_rk, map_type=LOCALLY_BOUNDED)
    call mesh%cell_ip%gather_offp(sigma)
    call stop_timer('mesh-to-mesh mapping')
  end subroutine set_conductivity
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Solver-side data retrieval procedures
 !!
 !! These functions return pointers to previously stored module data, which are
 !! cell-based arrays on the tet mesh.  Functions may be used in expressions,
 !! or as the target of a pointer assignment, however the target of the result
 !! should never be deallocated.
 !! 
 
  function permittivity () result (ptr)
    real(kind=rk), pointer :: ptr(:)
    ASSERT( allocated(eps) )
    ptr => eps
  end function permittivity
  
  function permeability () result (ptr)
    real(kind=rk), pointer :: ptr(:)
    ASSERT( allocated(mu) )
    ptr => mu
  end function permeability

  function conductivity () result (ptr)
    real(kind=rk), pointer :: ptr(:)
    ASSERT( allocated(sigma) )
    ptr => sigma
  end function conductivity
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Solver-side data store procedures
 !!
 !! Procedures map cell-based, tet-mesh input data to cell-based data on the
 !! hex mesh, storing the results as module data for later retrieval.
 !!
 
  subroutine set_joule_power_density (values)
  
    use simpl_mesh_type

    real(kind=rk), intent(in) :: values(:)
    type(simpl_mesh), pointer :: mesh

    ASSERT( allocated(joule) )
    
    mesh => EM_mesh()
    call start_timer('mesh-to-mesh mapping')
    call ht2em%map_field(values(:mesh%ncell_onP), joule, defval=0.0_rk, &
                         map_type=GLOBALLY_CONSERVATIVE, pullback=.true.)
    call stop_timer('mesh-to-mesh mapping')
    
    !! Record the data that gave rise to this Joule heat field.
    uhfs_q  = uhfs
    freq_q  = freq
    coil_q  = coil
    eps_q   = eps
    mu_q    = mu
    sigma_q = sigma

  end subroutine set_joule_power_density 

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Truchas-side data retrieval procedures
 !!
 !! These functions return pointers to previously stored module data, which are
 !! cell-based arrays on the hex mesh.  Functions may be used in expressions,
 !! or as the target of a pointer assignment, however the target of the result
 !! should never be deallocated.
 !! 
 
  function joule_power_density () result (ptr)
    real(kind=rk), pointer :: ptr(:)
    ASSERT( allocated(joule) )
    ptr => joule
  end function joule_power_density
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! EM namelist input accessor functions
 !!
  
  function get_epsilon_0 ()
    use physical_constants, only: vacuum_permittivity
    real(kind=rk) :: get_epsilon_0
    get_epsilon_0 = vacuum_permittivity
  end function get_epsilon_0
  
  function get_mu_0 ()
    use physical_constants, only: vacuum_permeability
    real(kind=rk) :: get_mu_0
    get_mu_0 = vacuum_permeability
  end function get_mu_0
  
  function get_EM_domain_type () result (string)
    use EM_input, only: EM_Domain_Type
    character(len=len_trim(EM_Domain_Type)) :: string
    string = EM_Domain_Type
  end function get_EM_domain_type
  
  function symmetry_axis ()
    use EM_input, only: sa => Symmetry_Axis
    character(len=1) :: symmetry_axis
    symmetry_axis = sa
  end function symmetry_axis
  
  integer function get_steps_per_cycle ()
    use EM_input, only: steps_per_cycle
    get_steps_per_cycle = steps_per_cycle
  end function get_steps_per_cycle
  
  integer function get_maximum_source_cycles ()
    use EM_input, only: Maximum_Source_Cycles
    get_maximum_source_cycles = Maximum_Source_Cycles
  end function get_maximum_source_cycles
  
  function get_ss_stopping_tolerance ()
    use EM_input, only: SS_Stopping_Tolerance
    real(kind=rk) :: get_ss_stopping_tolerance
    get_ss_stopping_tolerance = SS_Stopping_Tolerance
  end function get_ss_stopping_tolerance
  
  integer function get_maximum_cg_iterations ()
    use EM_input, only: Maximum_CG_Iterations
    get_maximum_cg_iterations = Maximum_CG_Iterations
  end function get_maximum_cg_iterations
  
  function get_cg_stopping_tolerance ()
    use EM_input, only: CG_Stopping_Tolerance
    real(kind=rk) :: get_cg_stopping_tolerance
    get_cg_stopping_tolerance = CG_Stopping_Tolerance
  end function get_cg_stopping_tolerance
  
  function get_num_etasq ()
    use EM_input, only: Num_Etasq
    real(kind=rk) :: get_num_etasq
    get_num_etasq = Num_Etasq
  end function get_num_etasq
  
  integer function get_output_level ()
    use EM_input, only: Output_Level
    get_output_level = Output_Level
  end function get_output_level
  
  logical function graphics ()
    use EM_input, only: Graphics_Output
    graphics = Graphics_Output
  end function graphics
  
  integer function get_num_probe_points ()
    use EM_input, only: num_probes
    get_num_probe_points = num_probes
  end function get_num_probe_points
  
  function get_probe_points () result (array)
    use EM_input, only: Probe_Points, num_probes
    real(kind=rk) :: array(3,num_probes)
    array = Probe_Points(:,:num_probes)
  end function get_probe_points
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SOURCE_HAS_CHANGED
 !!
 !! This procedure returns the value true if any of the values of the magnetic
 !! source field parameters (as returned by SOURCE_FREQUENCY, UNIFORM_SOURCE,
 !! and the CURRENT components of the array returned by INDUCTION_COILS) differ
 !! from those used to compute the Joule heat being returned by
 !! JOULE_POWER_DENSITY; otherwise the procedure returns the value false.
 !! It is used to determine whether the Joule heat needs to be recomputed.
 !! At present, the values of the coil currents, current frequency, and uniform
 !! field strength are used in the determination.
 !!
 !! At initialization, the Joule heat is set to zero, and the source field
 !! parameters set accordingly.  If the input specifies a nonzero source
 !! field (at some time), then this procedure will detect a change from the
 !! default pre-existing state as desired.
 !!

  logical function source_has_changed ()

    source_has_changed = .false.
    if (any(coil%current /= coil_q%current)) then
      source_has_changed = .true.   ! at least one coil current changed
    else if (uhfs /= uhfs_q) then
      source_has_changed = .true.   ! the uniform field strength changed
    else if (any(coil%current /= 0.0_rk) .or. (uhfs /= 0.0_rk)) then
      if (freq /= freq_q) then
        source_has_changed = .true. ! the frequency changed with a non-zero source field
      end if
    end if

  end function source_has_changed
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! MATERIAL_HAS_CHANGED
 !!
 !! This procedure returns the value true if the EM material parameters
 !! returned by CONDUCTIVITY and PERMEABILITY differ significantly from the
 !! values used to compute the Joule heat returned by JOULE_POWER_DENSITY;
 !! otherwise the procedure returns the value false.
 !!
 !! The maximum relative change is taken as the difference measure. Only
 !! when this difference exceeds MATERIAL_CHANGE_THRESHOLD is the difference
 !! considered significant.
 !!
 !! For conductivity, only the conducting region (sigma>0) is considered.
 !! An underlying assumption is that the conducting region remains fixed
 !! throughout the simulation.
 !!
 !! The permittivity only enters the equations through the displacement
 !! current term, which is exceedingly small in this magnetostatic regime
 !! and ought to be dropped entirely.  Thus the solution is essentially
 !! independent of the permittivity and so ignore any changes in its value.
 !!

  logical function material_has_changed ()

    use EM_input, only: Material_Change_Threshold
    use parallel_communication, only: global_maxval

    integer :: j
    real(kind=rk) :: dmu, dsigma
    character(len=80) :: string

    material_has_changed = .false.

    dmu = global_maxval(abs(mu-mu_q)/mu)

    dsigma = -huge(1.0_rk)
    do j = 1, size(sigma)
      if (sigma(j) > 0.0_rk) dsigma = max(dsigma, abs(sigma(j)-sigma_q(j))/sigma(j))
    end do
    dsigma = global_maxval(dsigma)

    if (max(dmu, dsigma) > Material_Change_Threshold) material_has_changed = .true.

    write(string,fmt='(3x,2(a,es10.3))') 'maximum relative change: sigma=', dsigma, ', mu=', dmu
    call TLS_info(string)

  end function material_has_changed

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SET_SOURCE_PROPERTIES
 !!
 !! This procedure sets the values returned by the functions SOURCE_FREQUENCY
 !! and UNIFORM_SOURCE, and the current-component of the value returned
 !! by INDUCTION_COILS, to the appropriate values for time T.
 !!
  
  subroutine set_source_properties (t)
  
    real(kind=rk), intent(in) :: t
    
    uhfs = uhfs_at_time(t)
    freq = freq_at_time(t)
    coil%current = curr_at_time(t)
    
  end subroutine set_source_properties
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! NO_SOURCE_FIELD, ZERO_JOULE_POWER_DENSITY
 !!
 !! The first procedure returns the value true if the values returned by
 !! UNIFORM_SOURCE and INDUCTION_COILS indicate that no magnetic source field
 !! is present.  In this case the second procedure can be called to simply set
 !! the resulting Joule power density field to zero (on the hex mesh),
 !! eliminating the need to 'compute' the Joule power density.
 !!
 
  logical function no_source_field ()
    no_source_field = (uhfs == 0.0_rk) .and. all(coil%current == 0.0_rk)
  end function no_source_field
  
  subroutine zero_joule_power_density ()
    ASSERT( allocated(joule) )
    joule = 0.0_rk
    !! Data from which the preceding joule heat derives.
    uhfs_q = 0.0_rk
    freq_q = 0.0_rk
    coil_q%current = 0.0_rk
    eps_q   = eps
    mu_q    = mu
    sigma_q = sigma
  end subroutine zero_joule_power_density 

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SOURCE_IS_SCALED, SCALE_JOULE_POWER_DENSITY
 !!
 !! The first function returns the value true if the *only* change to the
 !! magnetic source field parameters (see SOURCE_HAS_CHANGED) is a common
 !! scaling of the uniform source and current values; the observed scale
 !! factor is returned in the argument S.  In this case, the new Joule heat
 !! is simply S**2 times the previously computed Joule heat, and the second
 !! procedure can be called, passing S, to perform this scaling, avoiding a
 !! costly Joule heat computation.
 !!
 !! In practice this scaling is most likely to occur when there is only a
 !! uniform source or only a single coil.  However, the determination of
 !! the scale factor applies as well for any combination of coils and uniform
 !! source, though in this case the tolerance (1.0e-6) comes into play and
 !! its value is relatively arbitrary.
 !!

  logical function source_is_scaled (s)

    real(kind=rk), intent(out) :: s

    real(kind=rk) :: a, err

    !! Best scale factor in least-squares sense.
    a = sqrt(uhfs_q**2 + sum(coil_q%current**2))
    if (a > 0.0_rk) then
      s = (uhfs*uhfs_q + sum(coil%current*coil_q%current)) / a**2
    else
      s = 0.0_rk
    end if

    !! L2 error in best scaling.
    err = sqrt((uhfs - s*uhfs_q)**2 + sum((coil%current - s*coil_q%current)**2))

    !! Only change is a common scaling of the currents and uniform field strength.
    source_is_scaled = (freq == freq_q) .and. (err <= a*1.0e-6)

  end function source_is_scaled

  subroutine scale_joule_power_density (s)

    real(kind=rk), intent(in) :: s

    ASSERT( allocated(joule) )

    !! Scale the source strength: in the absence of roundoff
    !! UHFS == S*UHFS_Q and COIL%CURRENT == S * COIL_Q%CURRENT
    !! but we must ensure absolute equality, so instead ...
    uhfs_q = uhfs
    coil_q%current = coil%current

    !! Scale the joule heat accordingly.
    joule = s**2 * joule

 end subroutine scale_joule_power_density

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !!
 !! Solver-side source field parameter accessor functions
 !!
 !! These functions return data previously set in the module.  The physical
 !! property components of the coil data are static and are set during the 
 !! initialization of the module.  The electrical current component of the
 !! of the coil data, and the remaining data are time-dependent and are set
 !! by SET_SOURCE_PROPERTIES for a specific time.
 !!

  function induction_coils ()
    type(solenoid), pointer :: induction_coils(:)
    induction_coils => coil
  end function induction_coils
  
  function source_frequency ()
    real(kind=rk) :: source_frequency
    source_frequency = freq
  end function source_frequency
  
  function uniform_source ()
    real(kind=rk) :: uniform_source
    uniform_source = uhfs
  end function uniform_source

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! FREQ_AT_TIME, CURR_AT_TIME, UHFS_AT_TIME
 !!
 !! These auxillary procedures scan the time-dependent frequency, coil current,
 !! and uniform-source data that were read from the input file and returns the
 !! desired quantity at the specified time T.  If no uniform-source data was
 !! specified in the input file, UHFS_AT_TIME returns 0.
 !!
 !! If the time dependent quantity f is defined by the values f(1), ..., f(n+1)
 !! and times t(1), ..., t(n), then these functions return  f(1) for t < t(1),
 !! f(n+1) for t > t(n), and f(j+1) for t(j) <= t < t(j+1).  In the constant
 !! case (n=0), the functions simply return f(1).
 !!

  function freq_at_time (t) result (freq)

    use EM_input, only: src_freq

    real(kind=rk), intent(in) :: t
    real(kind=rk) :: freq

    freq = src_freq(time_interval(t))

  end function freq_at_time


  function curr_at_time (t) result (curr)

    use EM_input, only: coil_array

    real(kind=rk), intent(in) :: t
    real(kind=rk) :: curr(size(coil_array))

    integer :: n, j

    n = time_interval(t)
    do j = 1, size(curr)
      curr(j) = coil_array(j)%current(n)
    end do

  end function curr_at_time


  function uhfs_at_time (t) result (uhfs)

    use EM_input, only: unif_src

    real(kind=rk), intent(in) :: t
    real(kind=rk) :: uhfs

    if (size(unif_src) > 0) then
      uhfs = unif_src(time_interval(t))
    else
      uhfs = 0.0_rk
    end if

  end function uhfs_at_time


  function time_interval (t) result (n)

    use EM_input, only: src_time

    real(kind=rk), intent(in) :: t
    integer :: n

    n = size(src_time)
    do while (n > 0)
      if (t >= src_time(n)) exit
      n = n - 1
    end do
    n = n + 1

  end function time_interval

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! IN_TIME_INTERVAL
 !!
 !! This utility function returns the value true if the time T lies within
 !! the Nth source time interval (1-based numbering); otherwise it returns
 !! the value false.  This is provided for the use of sensitivity analysis.
 !!

  logical function in_time_interval (t, n)

    real(kind=rk), intent(in) :: t
    integer, intent(in) ::n

    in_time_interval = (n == time_interval(t))

  end function in_time_interval

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_JOULE_DATA
 !!
 !! This subroutine reads the Joule heat segment of a restart file opened (and
 !! prepositioned) on UNIT, and initializes the module variables FREQ_Q, UHFS_Q,
 !! COIL_Q, MU_Q, SIGMA_Q, and JOULE with this data (properly distributed and
 !! permuted).  It is an error if that data is not compatible with the mesh
 !! sizes and number of coils specified through the input file.  VERSION is the
 !! version number of the restart file format.
 !!

  subroutine read_joule_data (unit, version)

    use mesh_manager, only: simpl_mesh, simpl_mesh_ptr
    use legacy_mesh_api, only: pcell => unpermute_mesh_vector
    use restart_utilities, only: read_var, read_dist_array, halt
    use string_utilities, only: i_to_c
    use parallel_communication, only: global_sum

    integer, intent(in) :: unit,  version
    type(simpl_mesh), pointer :: mesh => null()
    
    integer :: n

    call TLS_info ('  reading the Joule heat data from the restart file')

    call read_var (unit, freq_q, 'READ_JOULE_DATA: error reading FREQ record')
    call read_var (unit, uhfs_q, 'READ_JOULE_DATA: error reading UHFS record')

    call read_var (unit, n, 'READ_JOULE_DATA: error reading NCOIL record')
    if (n /= size(coil_q)) call halt ('READ_JOULE_DATA: incompatible NCOIL value: ' // i_to_c(n))
    do n = 1, size(coil_q)
      call read_var (unit, coil_q(n)%current, 'READ_JOULE_DATA: error reading CURRENT record')
      call read_var (unit, coil_q(n)%center,  'READ_JOULE_DATA: error reading CENTER record')
      call read_var (unit, coil_q(n)%length,  'READ_JOULE_DATA: error reading LENGTH record')
      call read_var (unit, coil_q(n)%radius,  'READ_JOULE_DATA: error reading RADIUS record')
      call read_var (unit, coil_q(n)%nturns,  'READ_JOULE_DATA: error reading NTURNS record')
    end do

    mesh => simpl_mesh_ptr('alt')

    call read_var (unit, n, 'READ_JOULE_DATA: error reading NMU record')
    if (n /= mesh%cell_ip%global_size) call halt ('READ_JOULE_DATA: incompatible NMU value: ' // i_to_c(n))
    call read_dist_array (unit, mu_q(:mesh%ncell_onP), mesh%xcell(:mesh%ncell_onP), 'READ_JOULE_DATA: error reading MU record')
    call mesh%cell_ip%gather_offp(mu_q)

    call read_var (unit, n, 'READ_JOULE_DATA: error reading NSIGMA record')
    if (n /= mesh%cell_ip%global_size) call halt ('READ_JOULE_DATA: incompatible NSIGMA value: ' // i_to_c(n))
    call read_dist_array (unit, sigma_q(:mesh%ncell_onP), mesh%xcell(:mesh%ncell_onP), 'READ_JOULE_DATA: error reading SIGMA record')
    call mesh%cell_ip%gather_offp(sigma_q)

    call read_var (unit, n, 'READ_JOULE_DATA: error reading NJOULE record')
    if (n /= global_sum(size(joule))) call halt ('READ_JOULE_DATA: incompatible NJOULE value: ' // i_to_c(n))
    call read_dist_array (unit, joule, pcell, 'READ_JOULE_DATA: error reading JOULE record')
      
  end subroutine read_joule_data


  subroutine skip_joule_data (unit, version)
    use restart_utilities, only: skip_records, read_var
    integer, intent(in) :: unit, version
    integer :: n
    call skip_records (unit, 2, 'SKIP_JOULE_DATA: error skipping the Joule heat data')
    call read_var (unit, n, 'SKIP_JOULE_DATA: error skipping the Joule heat data')
    call skip_records (unit, 5*n + 6, 'SKIP_JOULE_DATA: error skipping the Joule heat data')
  end subroutine skip_joule_data

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DANU_WRITE_JOULE
 !!
 !! This subroutine writes the Joule heat array, and associated data used to
 !! compute it, to the HDF output file.  This is a global procedure.
 !! This is the direct analog of XML_WRITE_JOULE.
 !!
 !! Any I/O error that occurs is considered fatal, and this subroutine will
 !! abort execution of the code in such an event.
 !!
 !! NB: At this time the tet mesh is not written to the xml file, thus the
 !! fields over the tet mesh must be permuted to the external ordering prior
 !! to being written.  Conventionally, the xml parser would handle the
 !! permutation to external ordering using data written with the mesh.
 !!

  subroutine danu_write_joule (t)

    use parallel_communication
    use permutations
    use mesh_manager, only: simpl_mesh, simpl_mesh_ptr
    use truchas_h5_outfile, only: th5_sim_group
    use truchas_danu_output_data, only: outfile, io_group_size
    use truchas_env, only: output_file_name
    
    real(rk), intent(in) :: t

    integer :: n
    type(simpl_mesh), pointer :: mesh => null()
    integer, pointer :: cell_perm(:) => null()
    real(kind=rk), pointer :: col_mu(:)    => null()
    real(kind=rk), pointer :: col_sigma(:) => null()

    integer, save :: sim_num = 0
    character(32) :: sim_name
    type(th5_sim_group) :: sim

    ASSERT( allocated(mu_q) )
    ASSERT( allocated(sigma_q) )
    ASSERT( size(mu_q) == size(sigma_q) )
    ASSERT( allocated(joule) )
    ASSERT( associated(coil_q) )

    call outfile%reopen (output_file_name('h5'), io_group_size, is_IOP)
    
    mesh => simpl_mesh_ptr('alt')
    n = global_sum(mesh%ncell_onP)

    !! NNC, 2/2017. Not entirely sure why the reordering for MU and SIGMA is
    !! necessary.  I believe if we wrote the cell map for the tet mesh, then
    !! post-processing tools could do it when needed, as is done for the main
    !! mesh, and we could dispense with the collation here and truly write in
    !! parallel.  FIXME

    !! Collate the cell permutation array.
    allocate(cell_perm(merge(n,0,is_IOP)))
    call collate (cell_perm, mesh%xcell(:mesh%ncell_onP))

    !! Collate the cell-based MU array on the tet mesh, and restore it to the external order.
    allocate(col_mu(merge(n,0,is_IOP)))
    call collate (col_mu, mu_q(:mesh%ncell_onP))
    if (is_IOP) call reorder (col_mu, cell_perm, forward=.true.)

    !! Collate the cell-based SIGMA array on the tet mesh, and restore it to the external order.
    allocate(col_sigma(merge(n,0,is_IOP)))
    call collate (col_sigma, sigma_q(:mesh%ncell_onP))
    if (is_IOP) call reorder (col_sigma, cell_perm, forward=.true.)

    !! Write the data.
    sim_num = sim_num + 1
    write(sim_name,'(a,i3.3)') 'EM', sim_num
    call outfile%add_sim_group(trim(sim_name), sim)
    call sim%write_attr ('TIME', t)
    call TLS_info ('  writing EM restart data for ' // trim(sim_name))
    call sim%write_repl_data('FREQ', freq_q)
    call sim%write_repl_data('UHFS', uhfs_q)
    call sim%write_repl_data('COILS', solenoid_serialize(coil_q))
    call sim%write_repl_data('MU', col_mu)
    call sim%write_repl_data('SIGMA', col_sigma)
    call sim%write_dist_array('JOULE', joule, global_sum(size(joule)))

    call outfile%close()

    deallocate(cell_perm, col_mu, col_sigma)

  end subroutine danu_write_joule
  
  function solenoid_serialize (coil) result (array)
    type(solenoid), intent(in) :: coil(:)
    real(rk) :: array(7,size(coil))
    integer :: j
    do j = 1, size(coil)
      array(1,j) = coil(1)%current
      array(2,j) = coil(j)%center(1)
      array(3,j) = coil(j)%center(2)
      array(4,j) = coil(j)%center(3)
      array(5,j) = coil(j)%length
      array(6,j) = coil(j)%radius
      array(7,j) = coil(j)%nturns ! THIS IS FRAGILE
    end do
  end function solenoid_serialize

end module EM_data_proxy
