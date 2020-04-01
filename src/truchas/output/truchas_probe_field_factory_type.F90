!!
!! TRUCHAS_PROBE_FIELD_FACTORY_TYPE
!!
!! This defines the derived type TRUCHAS_PROBE_FIELD_FACTORY that implements
!! a factory method that creates concrete PROBE_FIELD class objects for the
!! types of probe solution fields supported by Truchas. Consider this part of
!! the Truchas multiphysics driver.
!!
!! Michael Hall <hall@lanl.gov>
!! July 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module truchas_probe_field_factory_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use probes_type, only: probe_field, probe_field_factory
  implicit none
  private

  type, extends(probe_field_factory), public :: truchas_probe_field_factory
  contains
    procedure :: alloc_field
  end type truchas_probe_field_factory

  !! For scalar cell/node-based field arrays that can be a pointer target.
  type, extends(probe_field) :: probe_field_vector
    real(r8), pointer :: vector(:) => null()
  contains
    procedure :: value => value_vector
  end type probe_field_vector

  !! For vector cell/node-based field arrays that can be a pointer target.
  type, extends(probe_field) :: probe_field_array
    real(r8), pointer :: array(:,:) => null()
  contains
    procedure :: value => value_array
  end type probe_field_array

  !! Special for the zone cell-centered velocity field.
  type, extends(probe_field) :: probe_field_zone_vc
  contains
    procedure :: value => value_zone_vc
  end type probe_field_zone_vc

contains

  function value_vector(this, index)
    class(probe_field_vector), intent(in) :: this
    integer, intent(in) :: index
    real(r8), allocatable :: value_vector(:)
    value_vector = [this%vector(index)]
  end function

  function value_array(this, index)
    class(probe_field_array), intent(in) :: this
    integer, intent(in) :: index
    real(r8), allocatable :: value_array(:)
    value_array = this%array(:,index)
  end function

  function value_zone_vc(this, index)
    use zone_module, only: zone
    class(probe_field_zone_vc), intent(in) :: this
    integer, intent(in) :: index
    real(r8), allocatable :: value_zone_vc(:)
    value_zone_vc = zone(index)%vc
  end function

  subroutine alloc_field(this, field, data)

    use string_utilities, only: lower_case
    use flow_driver, only: flow_enabled, flow_vel_cc_view, flow_p_cc_view
    use zone_module, only: zone

    class(truchas_probe_field_factory), intent(in) :: this
    class(probe_field), allocatable, intent(out) :: field
    character(*), intent(in) :: data

    select case (lower_case(data))
    case ('temperature')

      block
        type(probe_field_vector), allocatable :: pf
        allocate(pf)
        pf%label  = ', temperature'
        pf%kind   = 'cell'
        pf%vector => zone%temp
        call move_alloc(pf, field)
      end block

    case ('enthalpy')

      block
        type(probe_field_vector), allocatable :: pf
        allocate(pf)
        pf%label  = ', enthalpy'
        pf%kind   = 'cell'
        pf%vector => zone%enthalpy
        call move_alloc(pf, field)
      end block

    case ('pressure')

      if (flow_enabled()) then
        block
          type(probe_field_vector), allocatable :: pf
          allocate(pf)
          pf%label  = ', pressure'
          pf%kind   = 'cell'
          pf%vector => flow_p_cc_view()
          call move_alloc(pf, field)
        end block
      else  ! legacy flow
        block
          type(probe_field_vector), allocatable :: pf
          allocate(pf)
          pf%label  = ', pressure'
          pf%kind   = 'cell'
          pf%vector => zone%p
          call move_alloc(pf, field)
        end block
      end if

    case ('velocity')

      if (flow_enabled()) then
        block
          type(probe_field_array), allocatable :: pf
          allocate(pf)
          pf%label = ', vx, vy, vz'
          pf%kind  = 'cell'
          pf%array => flow_vel_cc_view()
          call move_alloc(pf, field)
        end block
      else  ! legacy flow
        block
          type(probe_field_zone_vc), allocatable :: pf
          allocate(pf)
          pf%label = ', vx, vy, vz'
          pf%kind  = 'cell'
          call move_alloc(pf, field)
        end block
      end if

    end select

  end subroutine

end module truchas_probe_field_factory_type
