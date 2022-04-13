!! LASER_IRRAD_FACTORY
!!
!! Provides a procedure for instantiating a LASER_IRRAD class object.
!! The dynamic type of the object and the parameters that define it are
!! specified through a parameter list.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! March 2017
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module laser_irrad_factory

  use laser_irrad_class
  use gauss_laser_irrad_type
  use beam_laser_irrad_type
  use parameter_list_type
  implicit none
  private

  public :: alloc_laser_irrad

contains

  subroutine alloc_laser_irrad(this, params)

    class(laser_irrad), allocatable, intent(out) :: this
    type(parameter_list) :: params

    integer :: stat
    character(:), allocatable :: irrad_type

    call params%get('type', irrad_type)
    select case (irrad_type)
    case ('gaussian')
      allocate(gauss_laser_irrad :: this)
    case ('gaussian beam')
      allocate(beam_laser_irrad :: this)
    case default
      INSIST(.false.)
    end select
    call this%init(params, stat)

  end subroutine alloc_laser_irrad

end module laser_irrad_factory
