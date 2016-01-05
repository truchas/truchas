!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module material_utilities

  use phase_property_table
  use material_table
  use material_system
  implicit none
  private
  
  public :: required_property_check, optional_property_check, constant_property_check
  public :: define_enthalpy_density_property, request_property
  
contains

  subroutine required_property_check (matids, prop, stat, errmsg)

    integer, intent(in)  :: matids(:)
    character(*), intent(in) :: prop
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg

    integer :: i, n, propid
    integer, pointer :: phaseid(:)

    if (ppt_has_property(prop)) then

      propid = ppt_property_id(prop)
      do i = 1, size(matids)  ! for each material system
        ASSERT(mt_valid_material(matids(i)))
        call ms_get_phase_id (mt_get_material(matids(i)), phaseid)
        do n = 1, size(phaseid)  ! loop through the phases of the material
          ASSERT(ppt_valid_phase(phaseid(n)))
          if (ppt_has_phase_property(phaseid(n), propid)) cycle  ! we're good
          stat = -1
          errmsg = 'missing property "' // trim(prop) // &
                   '" for phase "' // trim(ppt_phase_name(phaseid(n))) // '"'
          deallocate(phaseid)
          return
        end do
        deallocate(phaseid)
      end do

      stat = 0
      errmsg = ''

    else

      stat = -2
      errmsg = 'required property "' // trim(prop) // '" is unknown'

    end if

  end subroutine required_property_check

  subroutine optional_property_check (matids, prop, stat, errmsg)

    integer, intent(in)  :: matids(:)
    character(*), intent(in) :: prop
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg

    integer :: i, n, propid
    integer, pointer :: phaseid(:)
    logical :: some_have, some_lack

    if (ppt_has_property(prop)) then

      propid = ppt_property_id(prop)
      some_have = .false.; some_lack = .false.
      do i = 1, size(matids)  ! for each material system
        ASSERT(mt_valid_material(matids(i)))
        call ms_get_phase_id (mt_get_material(matids(i)), phaseid)
        do n = 1, size(phaseid)  ! loop through the phases of the material
          ASSERT(ppt_valid_phase(phaseid(n)))
          if (ppt_has_phase_property(phaseid(n), propid)) then
            some_have = .true.
          else
            some_lack = .true.
          end if
          if (some_have .and. some_lack) then
            stat = -1
            errmsg = 'incomplete specification of property "' // trim(prop) // '"'
            deallocate(phaseid)
            return
          end if
        end do
        deallocate(phaseid)
      end do

      if (some_have) then ! all phases have the property
        stat = 0
        errmsg = ''
      else  ! none have have the property
        stat = 1
        errmsg = ''
      end if

    else  ! none of the phases have the property

      stat = 1
      errmsg = ''

    end if

  end subroutine optional_property_check

 !!
 !! CONSTANT_PROPERTY_CHECK
 !!
 !! This subroutine verifies that the specified property is constant for each
 !! phase of the specified material systems and that within a material system
 !! the phase properties have the same constant values.
 !!
 !! Precondition: REQUIRED_PROPERTY_CHECK returns STAT=0 for the same arguments.
 !!

  subroutine constant_property_check (matids, prop, stat, errmsg)

    use kinds, only: r8
    use scalar_func_class
    use scalar_func_tools, only: is_const

    integer, intent(in)  :: matids(:)
    character(*), intent(in) :: prop
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: i, n, propid
    real(r8) :: value, state(0)
    integer, pointer :: phaseid(:)
    class(scalar_func), allocatable :: propfun

    ASSERT(ppt_has_property(prop))
    propid = ppt_property_id(prop)

    do i = 1, size(matids)  ! for each material system
      ASSERT(mt_valid_material(matids(i)))
      call ms_get_phase_id (mt_get_material(matids(i)), phaseid)
      do n = 1, size(phaseid)  ! loop through the phases of the material
        ASSERT(ppt_valid_phase(phaseid(n)))
        call ppt_get_phase_property (phaseid(n), propid, propfun)
        ASSERT(allocated(propfun))
        if (.not.is_const(propfun)) then
          stat = -2
          errmsg = 'property "' // trim(prop) // '" for phase "' // &
              trim(ppt_phase_name(phaseid(n))) // '" is required to be constant'
          deallocate(phaseid)
          return
        end if
        if (n == 1) then  ! save its value for comparison with remaining phases
          value = propfun%eval(state)
        else if (value /= propfun%eval(state)) then
          stat = -1
          errmsg = 'property "' // trim(prop) // '" for material system "' // &
              trim(mt_material_name(i)) // '" is required to be uniformly constant'
          deallocate(phaseid)
          return
        end if
      end do
      deallocate(phaseid)
    end do

    stat = 0
    errmsg = ''

  end subroutine constant_property_check
  
  !! Ensure that every material phase has the specified property
  !! defined, adding it if necessary with the given default value.

  subroutine request_property (matids, prop, default)

    use kinds, only: r8
    use scalar_func_factories
    use truchas_logging_services

    integer, intent(in)  :: matids(:)
    character(*), intent(in) :: prop
    real(r8), intent(in), optional :: default

    integer :: i, n, propid
    integer, pointer :: phaseid(:)
    class(scalar_func), allocatable :: f_default
    character(128) :: message
    
    if (ppt_has_property(prop)) then
      propid = ppt_property_id(prop)
    else
      if (present(default)) then
        call ppt_add_property (prop, propid)
      else
        call TLS_fatal ('required property "' // trim(prop) // '" is unknown')
      end if
    end if

    do i = 1, size(matids)  ! for each material system
      ASSERT(mt_valid_material(matids(i)))
      call ms_get_phase_id (mt_get_material(matids(i)), phaseid)
      do n = 1, size(phaseid)  ! loop through the phases of the material
        ASSERT(ppt_valid_phase(phaseid(n)))
        if (ppt_has_phase_property(phaseid(n),propid)) cycle  ! we're good
        if (present(default)) then  ! go ahead and define the property
          call alloc_const_scalar_func (f_default, default)
          call ppt_assign_phase_property (phaseid(n), propid, f_default)
          write(message, '(2x,3a,es10.3,3a)') 'Using default value "', trim(prop), &
              '" =', default, ' for phase "', trim(ppt_phase_name(phaseid(n))), '"'
          call TLS_info (message)
        else
          message = 'missing property "' // trim(prop) // '" for phase "' // &
                    trim(ppt_phase_name(phaseid(n))) // '"'
          deallocate(phaseid)
          call TLS_fatal (message)
        end if
      end do
      deallocate(phaseid)
    end do
      
  end subroutine request_property

  subroutine define_enthalpy_density_property (matid, stat, errmsg)

    use kinds, only: r8
    use scalar_func_class
    use scalar_func_tools

    integer, intent(in)  :: matid(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, n, cpid, hid, rhoid, hdid
    integer, pointer :: phaseid(:)
    real(r8) :: t0, h0
    type(mat_system), pointer :: ms
    class(scalar_func), allocatable :: h, hd, cp, rho

    ASSERT(ppt_has_property('specific heat'))
    cpid = ppt_property_id('specific heat')
    ASSERT(ppt_has_property('density'))
    rhoid = ppt_property_id('density')
    call ppt_add_property ('enthalpy', hid)
    call ppt_add_property ('enthalpy density', hdid)

    do i = 1, size(matid)  ! for each material system
      ASSERT(mt_valid_material(matid(i)))
      ms => mt_get_material(matid(i))
      call ms_get_phase_id (ms, phaseid)
      t0 = ms_ref_temp(ms)
      h0 = ms_ref_enthalpy(ms)
      do n = 1, size(phaseid)  ! loop through the phases of the material
        ASSERT(ppt_valid_phase(phaseid(n)))
        call ppt_get_phase_property (phaseid(n), cpid, cp)
        ASSERT(allocated(cp))
        call alloc_scalar_func_antideriv (cp, t0, h0, h, stat, errmsg)
        if (stat /= 0) then
          stat = -1
          errmsg = 'unable to integrate the "specific heat" property for phase "' // &
              trim(ppt_phase_name(phaseid(n))) // '"'
          deallocate(phaseid)
          return
        end if
        if (n < size(phaseid)) then
          t0 = ms_temp_hi(ms, n)
          h0 = ms_latent_heat(ms, n) + h%eval([t0])
        end if
        call ppt_get_phase_property (phaseid(n), rhoid, rho)
        ASSERT(allocated(rho))
        call alloc_scalar_func_product (rho, h, hd, stat, errmsg)
        ASSERT(stat == 0)
        if (ppt_has_phase_property(phaseid(n), hid)) then
          stat = -2
          errmsg = 'found conflicting "enthalpy" property for phase "' // &
              trim(ppt_phase_name(phaseid(n))) // '"'
          deallocate(phaseid)
          return
        end if
        call ppt_assign_phase_property (phaseid(n), hid, h)
        if (ppt_has_phase_property(phaseid(n), hdid)) then
          stat = -3
          errmsg = 'found conflicting "enthalpy density" property for phase "' // &
              trim(ppt_phase_name(phaseid(n))) // '"'
          deallocate(phaseid)
          return
        end if
        call ppt_assign_phase_property (phaseid(n), hdid, hd)
      end do
      deallocate(phaseid)
    end do

    stat = 0
    errmsg = ''

  end subroutine define_enthalpy_density_property

end module material_utilities
