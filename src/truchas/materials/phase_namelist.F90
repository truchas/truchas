
#include "f90_assert.fpp"

module phase_namelist

  use kinds
  use parallel_communication

  implicit none
  private
  
  public :: read_phase_namelists

  integer, parameter :: MAX_PHASE_PROPS = 10
  
  integer,   parameter :: NULL_I = HUGE(1)
  real(r8),  parameter :: NULL_R = HUGE(1.0_r8)
  character, parameter :: NULL_C = char(0)

contains

  subroutine read_phase_namelists (lun)

    use scalar_functions
    use function_table
    use phase_property_table
    use input_utilities, only: seek_to_namelist
    use string_utilities, only: i_to_c
    use ds_utilities
    
    integer, intent(in) :: lun
    
    !! namelist variables
    character(len=PPT_MAX_NAME_LEN) :: name
    character(len=PPT_MAX_NAME_LEN) :: property_name(MAX_PHASE_PROPS)
    character(len=FT_MAX_NAME_LEN)  :: property_function(MAX_PHASE_PROPS)
    real(r8)                        :: property_constant(MAX_PHASE_PROPS)

    namelist /phase/ name, property_name, property_function, property_constant

    integer :: n, stat, phase_id, prop_id, i
    logical :: found
    type(scafun) :: f

    call ds_info ('')
    call ds_info ('Reading PHASE namelists ...')

    if (is_IOP) rewind(lun)
    n = 0
    
    do  ! until all PHASE namelists have been read or an error occurs
    
      if (is_IOP) call seek_to_namelist (lun, 'PHASE', found, iostat=stat)
      call broadcast(stat)
      if (stat /= 0) call ds_halt ('error reading input file')

      call broadcast (found)
      if (.not.found) return  ! no further PHASE namelists found

      n = n + 1
      call ds_info ('  Reading PHASE namelist #' // i_to_c(n))

      !! Read the namelist variables, assigning default values first.
      if (is_IOP) then
        name              = NULL_C
        property_name     = NULL_C
        property_function = NULL_C
        property_constant = NULL_R
        read(lun, nml=phase, iostat=stat)
      end if

      call broadcast(stat)
      if (stat /= 0) call ds_halt('error reading PHASE namelist')

      !! Broadcast the namelist variables.
      call broadcast (name)
      call broadcast (property_name)
      call broadcast (property_function)
      call broadcast (property_constant)

      !! Check the phase name and add it to the phase property table.
      if (name == NULL_C .or. name == '') call ds_halt ('NAME must be assigned a nonempty value')
      if (ppt_has_phase(name)) then
        call ds_halt ('already read a PHASE namelist with this name: ' // trim(name))
      end if
      call ppt_add_phase (name, phase_id)

      !! Warn of any extraneous PROPERTY_CONSTANT or PROPERTY_FUNCTION values.
      if (any((property_name == NULL_C) .and. (property_constant /= NULL_R))) then
        call ds_warn ('found PROPERTY_CONSTANT values without a matching PROPERTY_NAME value')
      end if
      if (any((property_name == NULL_C) .and. (property_function /= NULL_C))) then
        call ds_warn ('found PROPERTY_FUNCTION values without a matching PROPERTY_NAME value')
      end if

      do i = 1, size(property_name)
      
        if (property_name(i) == NULL_C) cycle
        
        !! Check the property name and add it to the phase property table.
        if (property_name(i) == '') then
          call ds_halt ('empty string assigned to PROPERTY_NAME(' // i_to_c(i) // ')')
        end if
        call ppt_add_property (property_name(i), prop_id)
        if (ppt_has_phase_property(phase_id, prop_id)) then
          call ds_halt ('attempting to redefine phase property: ' // trim(property_name(i)))
        end if

        !! Check that only one of PROPERTY_FUNCTION and PROPERTY_CONSTANT was specified.
        if (property_constant(i) == NULL_R .eqv. property_function(i) == NULL_C) then
          call ds_halt ('exactly one of PROPERTY_CONSTANT(' // i_to_c(i) // &
                        ') and PROPERTY_FUNCTION(' // i_to_c(i) // ') must be assigned a value')
        end if

        !! Assign the specified phase property.
        if (property_constant(i) /= NULL_R) then
          call create_scafun_const (f, c=property_constant(i))
          call ppt_assign_phase_property (phase_id, prop_id, f)
          call destroy (f)
        else
          if (ft_has_function(property_function(i))) then
            call ppt_assign_phase_property (phase_id, prop_id, ft_get_function(property_function(i)))
          else
            call ds_halt ('unknown function name: ' // trim(property_function(i)))
          end if
        end if

      end do
    end do
    
  end subroutine read_phase_namelists

end module phase_namelist
