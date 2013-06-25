!!
!! MATERIAL_INTEROP
!!
!! This module provides some data and procedures that facilitate the exchange
!! of information between the new material software model used by the diffusion
!! solver and the aging material/property implementation still used by the rest
!! of Truchas.
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL GENERATE_MATERIAL_MAPPINGS (STAT, ERRMSG) generates various mappings
!!    that link Truchas material numbers with the material system and phase IDs
!!    of the new material software model employed by the diffusion solver.
!!    The mappings are provided as the following integer module arrays:
!!
!!    VOID_MATERIAL_INDEX
!!        The Truchas material number of the "void" material (density==0);
!!        this has value 0 if no void material exists.
!!
!!    MATERIAL_TO_PHASE(:)
!!        Maps Truchas material numbers to phase IDs.
!!        The map value is 0 for the void material if it exists.
!!
!!    MATERIAL_TO_SYSTEM(:)
!!        Maps Truchas material numbers to material system IDs.
!!        The map value is 0 for the void material if it exists.
!!
!!    PHASE_TO_MATERIAL(:)
!!        Maps phase IDs to Truchas material numbers.  A negative value marks
!!        an invalid phase ID, and a zero value marks a phase with no
!!        corresponding truchas material.
!!
!!    PHASE_TO_SYSTEM(:)
!!        Maps phase IDs to the ID of the material system to which it belongs.
!!        A negative value marks an invalid phase ID, and a zero value marks
!!        a phase not associated with a material system.
!!
!!    A 1-1 correspondence between materials and used phases is required that
!!    is based on the user-assigned names of the corresponding MATERIAL and
!!    PHASE namelists.  The only exception is for the void material (at most
!!    one is allowed) which is identified by a zero density; it must not have
!!    a matching phase.  A phase is "used" if it belongs to one of the material
!!    systems.  In addition, every phase that corresponds to a material must
!!    belong to some material system.  STAT returns a nonzero value if one of
!!    these constraints is violated, and an explanatory error message is
!!    assigned to the character string ERRMSG.
!!
!!    Other constraints, which will trip assertions if violated, are that the
!!    Truchas material names must be distinct from each other and that a
!!    phase can belong to at most one material system.  These constraints
!!    should be checked elsewhere (when processing the input, for example).
!!
!! The following routines evaluate selected phase properties used by the
!! diffusion solver given the Truchas material number and state variable
!! values.  These are intended for use solely by the advection solver that
!! computes the increment of heat advected into/out of a cell.  That solver
!! operates with Truchas material numbers but needs to use the same material
!! properties used by the diffusion solver.
!!
!!  DS_ENTHALPY(T_MATN, TEMP) returns the value of the enthalpy, as used by the
!!    diffusion solver, of Truchas material number T_MATN at temperature TEMP.
!!
!!  DS_DENSITY(T_MATN, TEMP) returns the value of the density, as used by the
!!    diffusion solver, of Truchas material number T_MATN at temperature TEMP.
!!
!!  DS_ENTHALPY_DENSITY(T_MATN, TEMP) returns the value of the enthalpy density,
!!    as used by the diffusion solver, of Truchas material number T_MATN at
!!    temperature TEMP.
!!

#include "f90_assert.fpp"

module material_interop

  use kinds, only: r8
  use phase_property_table
  use parameter_module, only: nmat
  use property_data_module, only: material_name, density
  use scalar_functions

  implicit none
  private
  save

  public :: generate_material_mappings
  public :: ds_density, ds_enthalpy, ds_enthalpy_density

  integer, allocatable, public :: material_to_phase(:), material_to_system(:)
  integer, allocatable, public :: phase_to_material(:), phase_to_system(:)
  integer, public :: void_material_index

  !! Private module data used by the enthalpy and density routines.
  type :: box
    type(scafun), pointer :: fn => null()
  end type box
  type(box), allocatable :: enthalpy_mat(:), density_mat(:), enthalpy_density_mat(:)

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! GENERATE_MATERIAL_MAPPINGS
 !!

  subroutine generate_material_mappings (stat, errmsg)

    use material_system
    use material_table

    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: i, n
    integer, pointer :: phase_id(:)
    integer, allocatable :: material_id(:)
    type(mat_system), pointer :: ms

    !! Identify the void material if any.
    void_material_index = 0
    do i = 1, nmat
      if (density(i) == 0.0_r8) then
        if (void_material_index > 0) then
          stat = -1
          errmsg = "found multiple void materials: " // &
              trim(material_name(void_material_index)) // ", " // trim(material_name(i))
          return
        else
          void_material_index = i
        end if
      end if
    end do

    !! Mapping MATERIAL_TO_PHASE from Truchas material numbers to phase IDs.
    allocate(material_to_phase(nmat))
    do i = 1, nmat
      if (i == void_material_index) then
        material_to_phase(i) = 0  ! marks the void material; no corresponding phase
      else if (ppt_has_phase(material_name(i))) then
        material_to_phase(i) = ppt_phase_id(material_name(i))
      else
        stat = -1
        errmsg = 'no matching phase for material "' // trim(material_name(i)) // '"'
        return
      end if
    end do

    !! Allocate mapping arrays from phase IDs.  A value of -1 marks invalid
    !! phase IDs and a value of 0 marks a valid phase ID for which no mapping
    !! has been assigned.  (associative arrays are what's needed here :(
    call ppt_get_phase_ids (phase_id)
    n = max(0, maxval(phase_id))
    allocate(phase_to_material(n), phase_to_system(n))
    phase_to_material  = -1; phase_to_material(phase_id)  = 0
    phase_to_system    = -1; phase_to_system(phase_id) = 0
    deallocate(phase_id)

    !! Mapping PHASE_TO_MATERIAL from phase IDs to Truchas material numbers.
    do i = 1, nmat
      if (material_to_phase(i) > 0) phase_to_material(material_to_phase(i)) = i
    end do

    !! Mapping PHASE_TO_SYSTEM from phase IDs to material system IDs.
    allocate(material_id(mt_num_material()))
    call mt_get_material_ids (material_id)
    do i = 1, size(material_id)
      ms => mt_get_material(material_id(i))
      call ms_get_phase_id (ms, phase_id)
      INSIST(all(phase_to_system(phase_id) == 0))
      phase_to_system(phase_id) = material_id(i)
    end do
    deallocate(material_id)

    !! Every used phase must map back to a Truchas material.
    do i = 1, size(phase_to_system)
      if (phase_to_system(i) > 0 .and. phase_to_material(i) <= 0) then
        stat = -1
        errmsg = 'no matching material for phase "' // trim(ppt_phase_name(i)) // '"'
        return
      end if
    end do

    !! Mapping MATERIAL_TO_SYSTEM from Truchas material numbers to material system IDs.
    allocate(material_to_system(nmat))
    do i = 1, nmat
      if (i == void_material_index) then
        material_to_system(i) = 0
      else
        material_to_system(i) = phase_to_system(material_to_phase(i))
      end if
    end do

    !! Every (non-void) material must map to a material system.
    do i = 1, nmat
      if (i /= void_material_index .and. material_to_system(i) == 0) then
        stat = -1
        errmsg = 'no matching material system for material "' // trim(material_name(i)) // '"'
        return
      end if
    end do

    stat = 0
    errmsg = ''

  end subroutine generate_material_mappings

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DS_ENTHALPY
 !!

  function ds_enthalpy (t_matn, state) result (value)

    integer,  intent(in) :: t_matn    ! Truchas material number
    real(r8), intent(in) :: state(:)  ! state variable values
    real(r8) :: value

    integer :: i, property_id, phase_id

    ASSERT(t_matn > 0 .and. t_matn <= nmat)

    if (.not.allocated(enthalpy_mat)) then
      allocate(enthalpy_mat(nmat))
      property_id = ppt_property_id('enthalpy')
      do i = 1, nmat
        if (i == void_material_index) cycle
        phase_id = material_to_phase(i)
        call ppt_get_phase_property (phase_id, property_id, enthalpy_mat(i)%fn)
      end do
    end if

    if (t_matn == void_material_index) then
      value = 0.0_r8
    else
      ASSERT(associated(enthalpy_mat(t_matn)%fn))
      value = eval(enthalpy_mat(t_matn)%fn, state)
    end if

  end function ds_enthalpy

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DS_DENSITY
 !!

  function ds_density (t_matn, state) result (value)

    integer,  intent(in) :: t_matn    ! Truchas material number
    real(r8), intent(in) :: state(:)  ! state variable values
    real(r8) :: value

    integer :: i, property_id, phase_id

    ASSERT(t_matn > 0 .and. t_matn <= nmat)

    if (.not.allocated(density_mat)) then
      allocate(density_mat(nmat))
      property_id = ppt_property_id('density')
      do i = 1, nmat
        if (i == void_material_index) cycle
        phase_id = material_to_phase(i)
        call ppt_get_phase_property (phase_id, property_id, density_mat(i)%fn)
      end do
    end if

    if (t_matn == void_material_index) then
      value = 0.0_r8
    else
      ASSERT(associated(density_mat(t_matn)%fn))
      value = eval(density_mat(t_matn)%fn, state)
    end if

  end function ds_density

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DS_ENTHALPY_DENSITY
 !!

  function ds_enthalpy_density (t_matn, state) result (value)

    integer,  intent(in) :: t_matn    ! Truchas material number
    real(r8), intent(in) :: state(:)  ! state variable values
    real(r8) :: value

    integer :: i, property_id, phase_id

    ASSERT(t_matn > 0 .and. t_matn <= nmat)

    if (.not.allocated(enthalpy_density_mat)) then
      allocate(enthalpy_density_mat(nmat))
      property_id = ppt_property_id('enthalpy density')
      do i = 1, nmat
        if (i == void_material_index) cycle
        phase_id = material_to_phase(i)
        call ppt_get_phase_property (phase_id, property_id, enthalpy_density_mat(i)%fn)
      end do
    end if

    if (t_matn == void_material_index) then
      value = 0.0_r8
    else
      ASSERT(associated(enthalpy_density_mat(t_matn)%fn))
      value = eval(enthalpy_density_mat(t_matn)%fn, state)
    end if

  end function ds_enthalpy_density

end module material_interop
