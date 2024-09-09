!!
!! EM_PROPERTIES
!!
!! This module provides procedures that deal with the EM properties, including
!! subroutines that evaluate the temperature-dependent properties on the heat
!! transfer mesh.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NNC, Dec 2019. This has been extracted from the original properties_module.
!!

#include "f90_assert.fpp"

module em_properties

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use material_model_driver, only:  matl_model
  implicit none
  private

  public :: define_default_em_properties
  public :: get_permittivity, get_permittivity_im, get_permeability, get_conductivity
  public :: permittivity_is_const, permittivity_im_is_const, permeability_is_const, &
      conductivity_is_const

contains

  !! Define default values for the EM properties read from input. EM properties
  !! only need to be specified in the input where they differ from the default.

  subroutine define_default_em_properties
    use material_utilities, only: define_property_default
    call define_property_default(matl_model, 'electrical-conductivity', 0.0_r8)
    call define_property_default(matl_model, 'electric-susceptibility', 0.0_r8)
    call define_property_default(matl_model, 'electric-susceptibility-im', 0.0_r8)
    call define_property_default(matl_model, 'magnetic-susceptibility', 0.0_r8)
  end subroutine

  !! Returns the real part of the electric permittivity at the current temperature.
  subroutine get_permittivity(value)
    use physical_constants, only: vacuum_permittivity
    use zone_module, only: zone
    real(r8), intent(out) :: value(:)
    call compute_cell_property('electric-susceptibility', zone%temp, value)
    value = vacuum_permittivity*(1.0_r8 + value)
  end subroutine

  !! Returns the imaginary part of the electric permittivity at the current temperature.
  subroutine get_permittivity_im(value)
    use physical_constants, only: vacuum_permittivity
    use zone_module, only: zone
    real(r8), intent(out) :: value(:)
    call compute_cell_property('electric-susceptibility-im', zone%temp, value)
    value = vacuum_permittivity*value
  end subroutine

  !! Returns true if the real part of the electric permittivity is constant
  !! with respect to time/temperature; otherwise it returns false.
  logical function permittivity_is_const()
    use material_utilities, only: constant_property_check
    integer :: stat
    character(:), allocatable :: errmsg
    call constant_property_check(matl_model, 'electric-susceptibility', stat, errmsg)
    permittivity_is_const = (stat == 0)
  end function

  !! Returns true if the imaginary part of the electric permittivity is
  !! constant with respect to time/temperature; otherwise it returns false.
  logical function permittivity_im_is_const()
    use material_utilities, only: constant_property_check
    integer :: stat
    character(:), allocatable :: errmsg
    call constant_property_check(matl_model, 'electric-susceptibility-im', stat, errmsg)
    permittivity_im_is_const = (stat == 0)
  end function

  !! Returns the magnetic permeability at the current temperature.
  subroutine get_permeability(value)
    use physical_constants, only: vacuum_permeability
    use zone_module, only: zone
    real(r8) :: value(:)
    call compute_cell_property('magnetic-susceptibility', zone%temp, value)
    value = vacuum_permeability*(1.0_r8 + value)
  end subroutine


  !! Returns true of the magnetic permeability is constant with respect to
  !! time/temperature; otherwise it returns false.
  logical function permeability_is_const()
    use material_utilities, only: constant_property_check
    integer :: stat
    character(:), allocatable :: errmsg
    call constant_property_check(matl_model, 'magnetic-susceptibility', stat, errmsg)
    permeability_is_const = (stat == 0)
  end function

  !! Returns the electrical conductivity at the current temperature.
  subroutine get_conductivity(value)
    use zone_module, only: zone
    real(r8) :: value(:)
    call compute_cell_property('electrical-conductivity', zone%temp, value)
  end subroutine

  !! Returns true if the electrical conductivity is constant with respect to
  !! time/temperature; otherwise it returns false.
  logical function conductivity_is_const()
    use material_utilities, only: constant_property_check
    integer :: stat
    character(:), allocatable :: errmsg
    call constant_property_check(matl_model, 'electrical-conductivity', stat, errmsg)
    conductivity_is_const = (stat == 0)
  end function

  !! This computes the specified property for the cells on the heat transfer
  !! mesh. The routine handles the material averaging of the property over a
  !! cell using the volume fraction data.
  !!
  !! NB: VOID material is ignored in this computation, which is suitable for
  !! the electrical conductivity and electric and magnetic susceptibilities
  !! which are all zero for void/vacuum. This would not be the case for the
  !! electric permittivity, for example.

  subroutine compute_cell_property(prop, temp, value)

    use scalar_func_class
    use scalar_func_tools, only: is_const
    use legacy_matl_api, only: gather_vof

    character(*), intent(in) :: prop
    real(r8), intent(in) :: temp(:)
    real(r8), intent(out) :: value(:)

    integer :: m, j
    real(r8) :: vofm (size(temp)), state(1), mval
    class(scalar_func), allocatable :: prop_fun

    ASSERT(size(value) == size(temp))

    value = 0.0_r8
    do m = 1, matl_model%nphase_real
      call matl_model%get_phase_prop(m, prop, prop_fun)
      ASSERT(allocated(prop_fun))
      call gather_vof (m, vofm)
      if (is_const(prop_fun)) then
        mval = prop_fun%eval(state)  ! state is ignored, but required
        value = value + mval*vofm
      else
        do j = 1, size(value)
          if (vofm(j) > 0.0_r8) then
            state(1) = temp(j)
            mval = prop_fun%eval(state)
            value(j) = value(j) + mval*vofm(j)
          end if
        end do
      end if
    end do

  end subroutine compute_cell_property

end module em_properties
