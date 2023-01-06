!!
!! USTRUC_COMP_CLASS
!!
!! This module provides an extension USTRUC_COMP of the abstract base class
!! USTRUC_ANALYSIS from which optional analysis components will be derived.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!! 1) While the USTRUC_COMP class should be regarded as an abstract type that
!! cannot/should not be instantiated, it is not actually abstract.  This is
!! because extensions of this type override the type bound procedures with
!! their specific functionality and then call back to the original procedures
!! defined here to pass control off to the next analysis component in the
!! chain.  Fortran requires that this type not be abstract in order to call
!! these procedures (even though they aren't deferred).
!!
!! 2) The canonical method for creating the ultimate USTRUC_ANALYSIS object is
!! to allocate a USTRUC_ANALYSIS pointer via the NEW_USTRUC_CORE function, pass
!! this to a NEW_USTRUC_* call to wrap it with an optional analysis component,
!! and take the resulting USTRUC_ANALYSIS pointer and continue to wrap with
!! additional optional analysis components. As long as no reference is retained
!! to any of the intermediate pointers, the final pointer can be considered to
!! own the chain it begins, and when it is deleted the finalization defined
!! here will do the right thing.
!!

#include "f90_assert.fpp"

module ustruc_comp_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ustruc_analysis_class
  use ustruc_core_type
  implicit none
  private

  public :: ustruc_analysis, ustruc_core  ! re-export

  type, extends(ustruc_analysis), public :: ustruc_comp ! see Note 1
    private
    class(ustruc_analysis), pointer :: next => null() ! the next component in the chain
    type(ustruc_core), pointer, public :: core => null() ! reference only -- do not own
  contains
    procedure :: init
    procedure :: set_state
    procedure :: update_state
    procedure :: get_comp_list
    procedure :: has
    procedure :: getl1
    procedure :: geti1
    procedure :: getr1
    procedure :: getr2
    procedure :: serialize
    procedure :: deserialize
    final :: delete_chain
  end type

  !! Unique IDs for the existing analysis components
  integer, parameter, public :: USTRUC_CORE_ID = 1
  integer, parameter, public :: USTRUC_GL_ID   = 2
  integer, parameter, public :: USTRUC_LDRD_ID = 3

contains

  !! Finalizer for objects of class USTRUC_COMP.  This effectively
  !! walks down the chain of CLASS(USTRUC_ANALYSIS) and deallocates each
  !! starting with the one at the end.  The assumption here is that the
  !! chain owns each of its objects.  See NOTE 2.

  recursive subroutine delete_chain(this)
    type(ustruc_comp) :: this
    if (associated(this%next)) deallocate(this%next)
  end subroutine

  subroutine init(this, comp)
    class(ustruc_comp), intent(out) :: this
    class(ustruc_analysis), pointer, intent(in) :: comp
    ASSERT(associated(comp))
    this%next => comp
    select type (comp)
    type is (ustruc_core)
      this%core => comp
    class is (ustruc_comp)
      ASSERT(associated(comp%core))
      this%core => comp%core
    end select
    this%n = this%core%n
  end subroutine

  recursive subroutine set_state(this, t, temp, temp_grad, frac, invalid)
    class(ustruc_comp), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)
    ASSERT(associated(this%next))
    call this%next%set_state(t, temp, temp_grad, frac, invalid)
  end subroutine

  recursive subroutine update_state(this, t, temp, temp_grad, frac, invalid)
    class(ustruc_comp), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)
    ASSERT(associated(this%next))
    call this%next%update_state(t, temp, temp_grad, frac, invalid)
  end subroutine

  recursive subroutine get_comp_list(this, list)
    class(ustruc_comp), intent(in) :: this
    integer, allocatable, intent(out) :: list(:)
    call this%next%get_comp_list(list)
  end subroutine

  recursive logical function has(this, name)
    class(ustruc_comp), intent(in) :: this
    character(*), intent(in) :: name
    ASSERT(associated(this%next))
    has = this%next%has(name)
  end function

  recursive subroutine getl1(this, name, array)
    class(ustruc_comp), intent(in) :: this
    character(*), intent(in) :: name
    logical, intent(out) :: array(:)
    ASSERT(associated(this%next))
    call this%next%get(name, array)
  end subroutine

  recursive subroutine geti1(this, name, array, invalid)
    class(ustruc_comp), intent(in) :: this
    character(*), intent(in) :: name
    integer, intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    ASSERT(associated(this%next))
    call this%next%get(name, array, invalid)
  end subroutine

  recursive subroutine getr1(this, name, array, invalid)
    class(ustruc_comp), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    ASSERT(associated(this%next))
    call this%next%get(name, array, invalid)
  end subroutine

  recursive subroutine getr2(this, name, array, invalid)
    class(ustruc_comp), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:,:)
    logical, intent(out), optional :: invalid(:)
    ASSERT(associated(this%next))
    call this%next%get(name, array, invalid)
  end subroutine

  recursive subroutine serialize(this, cid, array)
    use,intrinsic :: iso_fortran_env, only: int8
    class(ustruc_comp), intent(in) :: this
    integer, intent(in) :: cid
    integer(int8), allocatable, intent(out) :: array(:,:)
    ASSERT(associated(this%next))
    call this%next%serialize(cid, array)
  end subroutine

  recursive subroutine deserialize(this, cid, array)
    use,intrinsic :: iso_fortran_env, only: int8
    class(ustruc_comp), intent(inout) :: this
    integer, intent(in) :: cid
    integer(int8), intent(in) :: array(:,:)
    call this%next%deserialize(cid, array)
  end subroutine

end module ustruc_comp_class
