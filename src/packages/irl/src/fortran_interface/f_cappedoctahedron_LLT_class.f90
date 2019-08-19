!  This file is part of the Interface Reconstruction Library (IRL),
!  a library for interface reconstruction and computational geometry operations.
!
!  Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
!
!  This Source Code Form is subject to the terms of the Mozilla Public
!  License, v. 2.0. If a copy of the MPL was not distributed with this
!  file, You can obtain one at https://mozilla.org/MPL/2.0/.

!> \file f_CappedOctahedron_LLT_class.f90
!!
!! This file contains the Fortran interface for the
!! CappedDodecahedron class.

!> \brief A fortran type class that allows the creation of
!! IRL's CappedOctahedron class along with enabling
!! some of its methods.
module f_CapOcta_LLT_class
  use, intrinsic :: iso_c_binding
  use f_DefinedTypes
  implicit none

  type, public, bind(C) :: c_CapOcta_LLT
    type(C_PTR), private :: object = C_NULL_PTR
  end type c_CapOcta_LLT

  type, public :: CapOcta_LLT_type
    type(c_CapOcta_LLT) :: c_object
  contains
    final :: CapOcta_LLT_class_delete
  end type CapOcta_LLT_type

  interface new
    module procedure CapOcta_LLT_class_new
  end interface

  interface construct
    module procedure CapOcta_LLT_class_construct
  end interface

  interface adjustCapToMatchVolume
    module procedure CapOcta_LLT_class_adjustCapToMatchVolume
  end interface

  interface calculateVolume
    module procedure CapOcta_LLT_class_calculateVolume
  end interface
  
  interface getBoundingPts
    module procedure CapOcta_LLT_class_getBoundingPts
  end interface

  interface printToScreen
    module procedure CapOcta_LLT_class_printToScreen
  end interface 

  interface getPt
    module procedure CapOcta_LLT_class_getPt
  end interface

  interface

    subroutine F_CapOcta_LLT_new(this) &
      bind(C, name="c_CapOcta_LLT_new")
      import
      implicit none
      type(c_CapOcta_LLT) :: this
    end subroutine F_CapOcta_LLT_new

    subroutine F_CapOcta_LLT_delete(this) &
      bind(C, name="c_CapOcta_LLT_delete")
      import
      implicit none
      type(c_CapOcta_LLT) :: this
    end subroutine F_CapOcta_LLT_delete

    subroutine F_CapOcta_LLT_construct(this, a_octahedron) &
      bind(C, name="c_CapOcta_LLT_construct")
      import
      implicit none
      type(c_CapOcta_LLT) :: this
      real(C_DOUBLE), dimension(*), intent(in) :: a_octahedron ! dimension(1:3,1:7)
    end subroutine F_CapOcta_LLT_construct

    subroutine F_CapOcta_LLT_adjustCapToMatchVolume(this, a_correct_signed_volume) &
      bind(C, name="c_CapOcta_LLT_adjustCapToMatchVolume")
      import
      implicit none
      type(c_CapOcta_LLT) :: this
      real(C_DOUBLE), intent(in) :: a_correct_signed_volume ! scalar
    end subroutine F_CapOcta_LLT_adjustCapToMatchVolume

    function F_CapOcta_LLT_calculateVolume(this) result(a_octahedron_volume) &
      bind(C, name="c_CapOcta_LLT_calculateVolume")
      import
      implicit none
      type(c_CapOcta_LLT) :: this
      real(C_DOUBLE) :: a_octahedron_volume
    end function F_CapOcta_LLT_calculateVolume
    
    subroutine F_CapOcta_LLT_printToScreen(this) &
      bind(C, name="c_CapOcta_LLT_printToScreen")
      import
      implicit none
      type(c_CapOcta_LLT) :: this
    end subroutine F_CapOcta_LLT_printToScreen
    
    subroutine F_CapOcta_LLT_getBoundingPts(this, a_lower_pt, a_upper_pt) &
      bind(C, name="c_CapOcta_LLT_getBoundingPts")
      import
      implicit none
      type(c_CapOcta_LLT) :: this
      real(C_DOUBLE), dimension(*), intent(out) :: a_lower_pt ! dimension(1:3)
      real(C_DOUBLE), dimension(*), intent(out) :: a_upper_pt ! dimension(1:3)
    end subroutine F_CapOcta_LLT_getBoundingPts

    subroutine F_CapOcta_LLT_getPt(this, a_index, a_pt) &
      bind(C, name="c_CapOcta_LLT_getPt")
      import
      implicit none
      type(c_CapOcta_LLT) :: this
      integer(C_INT) :: a_index
      real(C_DOUBLE), dimension(*), intent(out) :: a_pt ! dimension(1:3)
    end subroutine F_CapOcta_LLT_getPt

  end interface


  contains

    impure elemental subroutine CapOcta_LLT_class_delete(this)
      implicit none
      type(CapOcta_LLT_type), intent(in) :: this
      call F_CapOcta_LLT_delete(this%c_object)
    end subroutine CapOcta_LLT_class_delete

    subroutine CapOcta_LLT_class_new(this)
      implicit none
      type(CapOcta_LLT_type), intent(inout) :: this
      call F_CapOcta_LLT_new(this%c_object)
    end subroutine CapOcta_LLT_class_new

    subroutine CapOcta_LLT_class_construct(this, a_octahedron)
      implicit none
      type(CapOcta_LLT_type), intent(inout) :: this
      real(IRL_double), dimension(1:3,1:7), intent(in) :: a_octahedron
      call F_CapOcta_LLT_construct(this%c_object, a_octahedron)
    end subroutine CapOcta_LLT_class_construct

    subroutine CapOcta_LLT_class_adjustCapToMatchVolume(this, a_correct_signed_volume)
      implicit none
      type(CapOcta_LLT_type), intent(inout) :: this
      real(IRL_double), intent(in) :: a_correct_signed_volume
      call F_CapOcta_LLT_adjustCapToMatchVolume(this%c_object, a_correct_signed_volume)
    end subroutine CapOcta_LLT_class_adjustCapToMatchVolume

    function CapOcta_LLT_class_calculateVolume(this) result(a_octahedron_volume)
      implicit none
      type(CapOcta_LLT_type), intent(inout) :: this
      real(IRL_double) :: a_octahedron_volume
      a_octahedron_volume = F_CapOcta_LLT_calculateVolume(this%c_object)
    end function CapOcta_LLT_class_calculateVolume
    
    subroutine CapOcta_LLT_class_printToScreen(this)
      implicit none
      type(CapOcta_LLT_type), intent(inout) :: this
      call F_CapOcta_LLT_printToScreen(this%c_object)
    end subroutine CapOcta_LLT_class_printToScreen
    
    subroutine CapOcta_LLT_class_getBoundingPts(this, a_lower_pt, a_upper_pt)
      implicit none
      type(CapOcta_LLT_type), intent(inout) :: this
      real(IRL_double), dimension(1:3), intent(out) :: a_lower_pt
      real(IRL_double), dimension(1:3), intent(out) :: a_upper_pt
      call F_CapOcta_LLT_getBoundingPts(this%c_object, a_lower_pt, a_upper_pt)
    end subroutine CapOcta_LLT_class_getBoundingPts

    function CapOcta_LLT_class_getPt(this, a_index) result(a_pt)
      implicit none
      type(CapOcta_LLT_type), intent(in) :: this
      integer(IRL_UnsignedIndex_t), intent(in) :: a_index
      real(IRL_double), dimension(3) :: a_pt
      call F_CapOcta_LLT_getPt(this%c_object, a_index, a_pt)
      return
    end function CapOcta_LLT_class_getPt

end module f_CapOcta_LLT_class
