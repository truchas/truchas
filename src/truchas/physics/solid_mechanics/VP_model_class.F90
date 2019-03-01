!!
!! VP_model_class
!!
!! This module defines the abstract base class VP_model that defines
!! the viscoplastic model interface expected by the solid mechanics kernel.
!! The power law and MTS models will extend this base class.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module VP_model_class

  use kinds, only: r8
  implicit none
  private
  
  type, abstract, public :: VP_model
  contains
    procedure(strain_rate), deferred :: strain_rate
  end type
  
  abstract interface
    function strain_rate (this, stress, temp, dt)
      import :: VP_model, r8
      class(VP_model), intent(in) :: this
      real(r8), intent(in) :: stress, temp, dt
      real(r8) :: strain_rate
    end function
  end interface

end module VP_model_class
