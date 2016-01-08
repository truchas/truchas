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
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
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
