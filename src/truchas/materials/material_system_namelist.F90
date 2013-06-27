
#include "f90_assert.fpp"

module material_system_namelist
  
  use kinds
  use parallel_communication
  use material_system

  implicit none
  private

  public :: read_material_system_namelists

  integer, parameter :: MAXPHASES=16

  integer,   parameter :: NULL_I = HUGE(1)
  real(r8),  parameter :: NULL_R = HUGE(1.0_r8)
  character, parameter :: NULL_C = char(0)

contains

  subroutine read_material_system_namelists (lun)

    use material_table
    use phase_property_table
    use input_utilities, only: seek_to_namelist
    use string_utilities, only: i_to_c
    use truchas_logging_services
  
    integer, intent(in) :: lun

    !! the namelist variables
    character(len=MT_MAX_NAME_LEN)  :: name
    character(len=PPT_MAX_NAME_LEN) :: phases(MAXPHASES)
    integer                         :: number_of_components
    logical                         :: temperature_dependent
    real(r8)                        :: reference_temp
    real(r8)                        :: reference_enthalpy
    !character(len=MAX_NAME_LENGTH)  :: phase_diagram_type
    real(r8)                        :: transition_temps_low(MAXPHASES-1)
    real(r8)                        :: transition_temps_high(MAXPHASES-1)
    real(r8)                        :: latent_heat(MAXPHASES-1)
    real(r8)                        :: smoothing_radius(MAXPHASES-1)
    
    namelist /material_system/ name, phases, number_of_components, &
         temperature_dependent, reference_temp, reference_enthalpy, &!phase_diagram_type, &
         transition_temps_low, transition_temps_high, smoothing_radius, latent_heat

    !! local variables
    logical :: found
    integer :: i, j, n, stat, num_phases, id, num_nml
    integer, allocatable :: phase_id(:), mid(:)
    integer, pointer :: pid(:)
    character(len=10) :: string
    type(mat_system) :: ms

    call TLS_info ('')
    call TLS_info ('Reading MATERIAL_SYSTEM namelists ...')

    if (is_IOP) rewind(lun)
    num_nml = 0
    
    do  ! until all MATERIAL_SYSTEM namelists have been read or an error occurs
    
      if (is_IOP) call seek_to_namelist (lun, 'MATERIAL_SYSTEM', found, iostat=stat)
      call broadcast (stat)
      if (stat /= 0) call TLS_fatal ('error reading input file')

      call broadcast (found)
      if (.not.found) return  ! no further MATERIAL_SYSTEM namelists found

      num_nml = num_nml + 1
      call TLS_info ('  Reading MATERIAL_SYSTEM namelist #' // i_to_c(num_nml))

      !! Read the namelist variables, assigning default values first.
      if (is_IOP) then
        name =                 NULL_C
        phases =               NULL_C
        number_of_components = NULL_I
        reference_temp =       NULL_R
        reference_enthalpy =   NULL_R
        !phase_diagram_type =   NULL_C
        transition_temps_low = NULL_R
        transition_temps_high= NULL_R
        latent_heat =          NULL_R
        smoothing_radius =     NULL_R
        temperature_dependent = .true.
        read(lun, nml=material_system, iostat=stat)
      endif
      call broadcast (stat)
      if (stat /= 0) call TLS_fatal ('error reading MATERIAL_SYSTEM namelist')

      !! Broadcast the namelist variables.
      call broadcast (name)
      call broadcast (phases)
      call broadcast (number_of_components)
      call broadcast (reference_temp)
      call broadcast (reference_enthalpy)
      !call broadcast (phase_diagram_type)
      call broadcast (transition_temps_low)
      call broadcast (transition_temps_high)
      call broadcast (latent_heat)
      call broadcast (smoothing_radius)
      call broadcast (temperature_dependent)
      
      !! Check the user-supplied name for the material system.
      if (name == NULL_C .or. name == '') then
        call TLS_fatal ('NAME must be assigned a nonempty value')
      else if (mt_has_material(name)) then
        call TLS_fatal ('already read a MATERIAL_SYSTEM namelist with this name: ' // trim(name))
      end if

      !! Check PHASES for basic correctness.
      num_phases = count(phases /= NULL_C)
      if (num_phases < 1) then
        call TLS_fatal ('at least one phase name must be assigned to PHASES')
      else if (any(phases(:num_phases) == NULL_C)) then
        call TLS_fatal ('phase names must be packed into the PHASES array')
      end if

      !! Check that these are known phases and retrieve their IDs.
      allocate(phase_id(num_phases))
      do i = 1, num_phases
        if (ppt_has_phase(phases(i))) then
          phase_id(i) = ppt_phase_id(phases(i))
        else
          call TLS_fatal ('unknown phase: "' // trim(phases(i)) // '"')
        end if
      end do
      
      !! Check that the phases are distinct.
      do i = 1, num_phases-1
        do j = i+1, num_phases
          if (phases(i) == phases(j)) then
            call TLS_fatal ('phase name appears more than once: "' // trim(phases(i)) // '"')
          end if
        end do
      end do
      
      !! Check that none of the phases belong to a previous material system.
      allocate(mid(mt_num_material()))
      call mt_get_material_ids (mid)
      do i = 1, size(mid)
        call ms_get_phase_id (mt_get_material(mid(i)), pid)
        do j = 1, size(phase_id)
          if (any(phase_id(j) == pid)) then
            call TLS_fatal ('phase "' // trim(phases(j)) // '" already belongs to material system "' &
                                      // trim(mt_material_name(mid(i))) // '"')
          end if
        end do
        deallocate(pid)
      end do
      deallocate(mid)

      !! Check the number of components, supplying the default of 1 if necessary.
      if (number_of_components == NULL_I) then
        number_of_components = 1
      else if (number_of_components < 1) then
        call TLS_fatal ('NUMBER_OF_COMPONENTS must be positive')
      end if

      if (num_phases > 1) then  ! have a phase diagram
      
        if (.not.temperature_dependent) then
          call TLS_fatal ('multi-phase material systems must be temperature-dependent')
        end if
      
        if (number_of_components == 1) then ! unary phase diagram
      
          n = num_phases - 1  ! the number of transition intervals

          !! Check TRANSITION_TEMPS_LOW and TRANSITION_TEMPS_HIGH for correctness.
          if (any(transition_temps_low(n+1:) /= NULL_R)) then
            call TLS_fatal ('found extraneous TRANSITION_TEMPS_LOW values')
          end if
          if (any(transition_temps_low(:n) == NULL_R)) then
            call TLS_fatal ('too few TRANSITION_TEMPS_LOW values; expecting ' // i_to_c(n))
          end if

          if (any(transition_temps_high(n+1:) /= NULL_R)) then
            call TLS_fatal ('found extraneous TRANSITION_TEMPS_HIGH values')
          end if
          if (any(transition_temps_high(:n) == NULL_R)) then
            call TLS_fatal ('too few TRANSITION_TEMPS_HIGH values; expecting ' // i_to_c(n))
          end if

          if (any(transition_temps_high(:n) <= transition_temps_low(:n))) then
            call TLS_fatal ('all transition interval widths must be positive')
          end if

          !! Check for properly ordered transition intervals.
          if (any(transition_temps_high(1:n-1) >= transition_temps_low(2:n))) then
            call TLS_fatal('the sequence of transition intervals are not monotonically increasing')
          end if

          !! Assign default values for SMOOTHING_RADIUS where needed.
          do i = 1, n
            if (smoothing_radius(i) == NULL_R) then
              smoothing_radius(i) = 0.25_r8*(transition_temps_high(i) - transition_temps_low(i))
              if (i > 1) then
                smoothing_radius(i) = min(smoothing_radius(i), &
                     0.25_r8*(transition_temps_low(i)-transition_temps_high(i-1)))
              end if
              if (i < n) then
                smoothing_radius(i) = min(smoothing_radius(i), &
                     0.25_r8*(transition_temps_low(i+1)-transition_temps_high(i)))
              end if
              write(string,'(es10.4)') smoothing_radius(i)
              call TLS_info ('    using default SMOOTHING_RADIUS(' // i_to_c(i) // ') = ' // string)
            end if
          end do

          !! Check SMOOTHING_RADIUS for correctness.
          if (any(smoothing_radius(n+1:) /= NULL_R)) then
            call TLS_fatal ('found extraneous SMOOTHING_RADIUS values')
          end if
          do i = 1, n
            if (smoothing_radius(i) < 0.0_r8) then
              call TLS_fatal ('all SMOOTHING_RADIUS values must be nonnegative')
            end if
            if (smoothing_radius(i) > (transition_temps_high(i) - transition_temps_low(i))/2) then
              call TLS_fatal ('SMOOTHING_RADIUS exceeds half the transition width')
            end if
          end do

          !! Check that the smoothed transition intervals remain separated.
          do i = 1, n-1
            if (transition_temps_high(i) + smoothing_radius(i) > &
                transition_temps_low(i+1) - smoothing_radius(i+1)) then
              call TLS_fatal ('SMOOTHING_RADIUS is too large; adjacent transition intervals overlap')
            end if
          end do

          !! Check LATENT_HEAT for correctness.
          if (any(latent_heat(n+1:) /= NULL_R)) then
            call TLS_fatal ('found extraneous LATENT_HEAT values')
          end if
          if (any(latent_heat(:n) == NULL_R)) then
            call TLS_fatal ('too few LATENT_HEAT values; expecting ' // i_to_c(n))
          end if
          if (any(latent_heat(:n) <= 0.0_r8)) then
            call TLS_fatal('all LATENT_HEAT values must be positive')
          end if

          !if (phase_diagram_type == NULL_C) then
          !   call halt('a phase_diagram_type must be specified')
          !else if (phase_diagram_type /= 'unary') then
          !   call halt('unary is currently the only supported phase diagram type')
          !endif
        
        else  ! multi-component phase diagram

          call TLS_fatal ('multi-phase, multi-component material systems are not yet implemented')
           
        end if
        
      else  ! this is a single-phase material system -- no phase diagram

        if (number_of_components == 1 .and. .not.temperature_dependent) then
          call TLS_fatal ('single-component material systems must be temperature-dependent')
        end if
      
        if (any(transition_temps_low /= NULL_R)) then
          call TLS_fatal ('use of TRANSITION_TEMPS_LOW is not relevant to single-phase systems')
        endif
        if (any(transition_temps_high /= NULL_R)) then
          call TLS_fatal ('use of TRANSITION_TEMPS_HIGH is not relevant to single-phase systems')
        endif
        if (any(smoothing_radius /= NULL_R)) then
          call TLS_fatal ('use of SMOOTHING_RADIUS is not relevant to single-phase systems')
        endif
        if (any(latent_heat /= NULL_R)) then
          call TLS_fatal ('use of LATENT_HEAT is not relevant to single-phase systems')
        endif
        !if (phase_diagram_type /= NULL_C) then
        !  call TLS_fatal ('use of PHASE_DIAGRAM_TYPE is not relevant to single-phase systems')
        !endif
        
      endif

      if (reference_temp == NULL_R) then
        reference_temp = 0.0_r8
        write(string,'(es10.3)') reference_temp
        call TLS_info ('    using default REFERENCE_TEMP = ' // string)
      endif

      if (reference_enthalpy == NULL_R) then
        reference_enthalpy = 0.0_r8
        write(string,'(es10.3)') reference_enthalpy
        call TLS_info ('    using default REFERENCE_ENTHALPY = ' // string)
      endif

      !! now initialize the current material_system object

      call ms_create(ms,number_of_components,temperature_dependent,phase_id,&
           transition_temps_low(1:num_phases-1), transition_temps_high(1:num_phases-1),&
           latent_heat(1:num_phases-1),smoothing_radius(1:num_phases-1),&
           reference_temp,reference_enthalpy)
      call mt_add_material (name, ms, id)
      call destroy (ms)

      deallocate(phase_id)
    
    end do
 
  end subroutine read_material_system_namelists

end module material_system_namelist
