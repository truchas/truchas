!!
!! USTRUC_PLUGIN_CLASS
!!
!! This module provides an extension USTRUC_PLUGIN of the abstract base class
!! USTRUC_COMP from which optional analysis components will be derived.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2014
!!
!! PROGRAMMING INTERFACE
!!
!! For flexibility, the microstructure model is decomposed into individual
!! analysis components that can be combined dynamically at run time.  The
!! decorator programming pattern is used to implement this design.  This
!! type corresponds to the "decorator" class in that pattern.  There is a
!! non-optional analysis component implemented by USTRUC_CORE, which is
!! wrapped by one or more optional "decorator" analysis components, forming
!! a chain of components (all of class USTRUC_COMP) with the core component
!! at the end.  USTRUC_PLUGIN type encapsulates this chain.
!!
!! NOTES
!!
!! 1) While the USTRUC_PLUGIN class should be regarded as an abstract type that
!! cannot/should not be instantiated, it is not actually abstract.  This is
!! because extensions of this type override the type bound procedures with
!! their specific functionality and then call back to the original procedures
!! defined here to pass control off to the next analysis component in the
!! chain.  Fortran requires that this type not be abstract in order to call
!! these procedures (even though they aren't deferred).
!!
!! 2) The canonical method for creating the ultimate USTRUC_COMP object is
!! to allocate a USTRUC_COMP pointer via the NEW_USTRUC_CORE function, pass
!! this to a NEW_USTRUC_* call to wrap it with an optional analysis component,
!! and take the resulting USTRUC_COMP pointer and continue to wrap with
!! additional optional analysis components.  As long as no reference is
!! retained to any of the intermediate pointers, the final pointer can be
!! considered to own the chain it begins, and when it is deleted the finalization
!! defined here will do the right thing.
!!

#include "f90_assert.fpp"

module ustruc_plugin_class

  use kinds, only: r8
  use ustruc_comp_class
  use ustruc_core_type
  implicit none
  private

  public :: ustruc_comp, ustruc_core  ! re-export

  type, extends(ustruc_comp), public :: ustruc_plugin
    private
    class(ustruc_comp), pointer :: comp => null() ! the next component in the chain
    type(ustruc_core),  pointer, public :: core => null() ! reference only -- do not own
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
  integer, parameter, public :: USTRUC_GV0_ID  = 2
  integer, parameter, public :: USTRUC_GV1_ID  = 3

contains

  !! Finalizer for objects of class USTRUC_PLUGIN.  This effectively
  !! walks down the chain of CLASS(USTRUC_OBJECTS) and deallocates each
  !! starting with the one at the end.  The assumption here is that the
  !! chain owns each of its objects.  See NOTE 2.

  recursive subroutine delete_chain (this)
    type(ustruc_plugin) :: this
    if (associated(this%comp)) deallocate(this%comp)
  end subroutine delete_chain

  subroutine init (this, comp)
    class(ustruc_plugin), intent(out) :: this
    class(ustruc_comp), pointer, intent(in) :: comp
    ASSERT(associated(comp))
    this%comp => comp
    select type (comp)
    type is (ustruc_core)
      this%core => comp
    class is (ustruc_plugin)
      ASSERT(associated(comp%core))
      this%core => comp%core
    end select
    this%n = this%core%n
  end subroutine init

  recursive subroutine set_state (this, t, temp, temp_grad, frac, invalid)
    class(ustruc_plugin), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)
    ASSERT(associated(this%comp))
    call this%comp%set_state (t, temp, temp_grad, frac, invalid)
  end subroutine set_state

  recursive subroutine update_state (this, t, temp, temp_grad, frac, invalid)
    class(ustruc_plugin), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)
    ASSERT(associated(this%comp))
    call this%comp%update_state (t, temp, temp_grad, frac, invalid)
  end subroutine update_state

  recursive subroutine get_comp_list (this, list)
    class(ustruc_plugin), intent(in) :: this
    integer, allocatable, intent(out) :: list(:)
    call this%comp%get_comp_list (list)
  end subroutine get_comp_list

  recursive logical function has (this, name)
    class(ustruc_plugin), intent(in) :: this
    character(*), intent(in) :: name
    ASSERT(associated(this%comp))
    has = this%comp%has (name)
  end function has

  recursive subroutine getl1 (this, name, array)
    class(ustruc_plugin), intent(in) :: this
    character(*), intent(in) :: name
    logical, intent(out) :: array(:)
    ASSERT(associated(this%comp))
    call this%comp%get (name, array)
  end subroutine getl1

  recursive subroutine geti1 (this, name, array, invalid)
    class(ustruc_plugin), intent(in) :: this
    character(*), intent(in) :: name
    integer, intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    ASSERT(associated(this%comp))
    call this%comp%get (name, array, invalid)
  end subroutine geti1

  recursive subroutine getr1 (this, name, array, invalid)
    class(ustruc_plugin), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    ASSERT(associated(this%comp))
    call this%comp%get (name, array, invalid)
  end subroutine getr1

  recursive subroutine getr2 (this, name, array, invalid)
    class(ustruc_plugin), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:,:)
    logical, intent(out), optional :: invalid(:)
    ASSERT(associated(this%comp))
    call this%comp%get (name, array, invalid)
  end subroutine getr2

  recursive subroutine serialize (this, cid, array)
    use,intrinsic :: iso_fortran_env, only: int8
    class(ustruc_plugin), intent(in) :: this
    integer, intent(in) :: cid
    integer(int8), allocatable, intent(out) :: array(:,:)
    ASSERT(associated(this%comp))
    call this%comp%serialize (cid, array)
  end subroutine serialize

  recursive subroutine deserialize (this, cid, array)
    use,intrinsic :: iso_fortran_env, only: int8
    class(ustruc_plugin), intent(inout) :: this
    integer, intent(in) :: cid
    integer(int8), intent(in) :: array(:,:)
    call this%comp%deserialize (cid, array)
  end subroutine deserialize

end module ustruc_plugin_class
