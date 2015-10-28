!!
!! USTRUC_COMP_CLASS
!!
!! This module provides the abstract base class USTRUC_COMP that defines the
!! interface to the low-level microstructure analysis component used by the
!! microstructure modeling kernel as implemented by the USTRUC_MODEL type.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2014
!!
!! PROGRAMMING INTERFACE
!!
!! This interface provides for point microstructure models that depend on
!! local quantities like temperature, solid fraction, and their gradients
!! at a point.  The models are applied independently at each point in a
!! collection of points; at this level there is no reference to any spatial
!! mesh.  The models are expected to be time dependent however, requiring
!! concrete implementations to maintain state through time as needed.
!!
!! For flexibility, the microstructure model is decomposed into individual
!! analysis components that can be combined dynamically at run time.  The
!! decorator programming pattern is used to implement this design.  The
!! core component which holds the state used by other components is defined
!! by the concrete implementation USTRUC_CORE of this abstract class.  Other
!! optional analysis components which add specific functionality will be
!! derived from the USTRUC_PLUGIN class which also derives from this abstract
!! class.
!!
!! The abstract base class specifies the following type bound procedures.
!! The data component N is the number of points the analysis component is
!! being applied to.
!!
!!  SET_STATE (T, TEMP, TEMP_GRAD, FRAC, FRAC_GRAD, INVALID) sets the initial
!!    state, removing and resetting any existing state back to its starting
!!    condition.  All arguments are intent(in).  The scalar T is time, TEMP
!!    and FRAC are the rank-1 arrays of length N giving the temperature and
!!    solid fraction, and TEMP_GRAD and FRAC_GRAD are corresponding rank-2
!!    arrays with first dimension size 3, giving the spatial gradients.  The
!!    corresponding rank-1 logical array INVALID marks those points for which
!!    the preceding arrays do not have valid data.
!!
!!  UPDATE_STATE (T, TEMP, TEMP_GRAD, FRAC, FRAC_GRAD, INVALID) updates the
!!    state.  The interface is identical to SET_STATE.  The difference here
!!    is that the analysis component is expected to be advancing its internal
!!    state using the previous state and this new passed state.  SET_STATE
!!    must be called before this subroutine is called.
!!
!!  GET (NAME, ARRAY) returns the named data in the provided array.  NAME
!!    is a intent(in) character variable with an implementation-defined
!!    value that identifies the requested data.  The data is returned in
!!    the intent(out) ARRAY. It is a rank-1 array of length N, or a rank-2
!!    array of shape (3,N), and of real, integer, or logical type.  It is
!!    an error if the name is not known.  Concrete implementations must
!!    provide implementations of the specific subroutines GETL1, GETI1, GETR1,
!!    and GETR2.
!!
!! The abstract class also defines the type bound function VECTOR_MAGNITUDE
!! which computes the magnitude of the passed 2 or 3-vector using a numerically
!! robust procedure.  Implementations should use this function when computing
!! the norm of a gradient.
!!

#include "f90_assert.fpp"

module ustruc_comp_class

  use kinds, only: r8
  use,intrinsic :: iso_fortran_env, only: int8
  implicit none
  private

  type, abstract, public :: ustruc_comp
    integer :: n  ! number of points
  contains
    procedure(update), deferred :: set_state
    procedure(update), deferred :: update_state
    procedure(has),    deferred :: has
    procedure(comp_list), deferred :: get_comp_list
    generic :: get => getl1, geti1, getr1, getr2
    procedure(getl1), deferred :: getl1
    procedure(geti1), deferred :: geti1
    procedure(getr1), deferred :: getr1
    procedure(getr2), deferred :: getr2
    procedure(serialize), deferred :: serialize
    procedure(deserialize), deferred :: deserialize
    procedure, nopass :: vector_magnitude
  end type

  abstract interface
    subroutine update (this, t, temp, temp_grad, frac, frac_grad, invalid)
      import ustruc_comp, r8
      class(ustruc_comp), intent(inout) :: this
      real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:), frac_grad(:,:)
      logical,  intent(in) :: invalid(:)
    end subroutine
    logical function has (this, name)
      import ustruc_comp
      class(ustruc_comp), intent(in) :: this
      character(*), intent(in) :: name
    end function has
    subroutine comp_list (this, list)
      import ustruc_comp
      class(ustruc_comp), intent(in) :: this
      integer, allocatable, intent(out) :: list(:)
    end subroutine
    subroutine getl1 (this, name, array)
      import ustruc_comp
      class(ustruc_comp), intent(in) :: this
      character(*), intent(in) :: name
      logical, intent(out) :: array(:)
    end subroutine
    subroutine geti1 (this, name, array, invalid)
      import ustruc_comp, r8
      class(ustruc_comp), intent(in) :: this
      character(*), intent(in) :: name
      integer, intent(out) :: array(:)
      logical, intent(out), optional :: invalid(:)
    end subroutine
    subroutine getr1 (this, name, array, invalid)
      import ustruc_comp, r8
      class(ustruc_comp), intent(in) :: this
      character(*), intent(in) :: name
      real(r8), intent(out) :: array(:)
      logical, intent(out), optional :: invalid(:)
    end subroutine
    subroutine getr2 (this, name, array, invalid)
      import ustruc_comp, r8
      class(ustruc_comp), intent(in) :: this
      character(*), intent(in) :: name
      real(r8), intent(out) :: array(:,:)
      logical, intent(out), optional :: invalid(:)
    end subroutine
    subroutine serialize (this, cid, array)
      import ustruc_comp, int8
      class(ustruc_comp), intent(in) :: this
      integer, intent(in) :: cid
      integer(int8), allocatable, intent(out) :: array(:,:)
    end subroutine
    subroutine deserialize (this, cid, array)
      import ustruc_comp, int8
      class(ustruc_comp), intent(inout) :: this
      integer, intent(in) :: cid
      integer(int8), intent(in) :: array(:,:)
    end subroutine
  end interface

contains

  !! Numerically robust procedure to compute the magnitude of a 2 or 3-vector.

  pure function vector_magnitude (v) result (vmag)

    real(r8), intent(in) :: v(:)  ! length 2 or 3
    real(r8) :: vmag

    real(r8) :: a, b, c, t

    select case (size(v))
    case (2)
      a = abs(v(1))
      b = abs(v(2))
      !! Swap largest value to A.
      if (b > a) then
        t = a
        a = b
        b = t
      end if
      vmag = 0.0_r8
      if (a > 0.0_r8) vmag = a * sqrt(1.0_r8 + (b/a)**2)
    case (3)
      a = abs(v(1))
      b = abs(v(2))
      c = abs(v(3))
      !! Swap largest value to A.
      if (b > a) then
        if (c > b) then
          t = a
          a = c
          c = t
        else
          t = a
          a = b
          b = t
        end if
      else if (c > a) then
        t = a
        a = c
        c = t
      end if
      vmag = 0.0_r8
      if (a > 0.0_r8) vmag = a * sqrt(1.0_r8 + ((b/a)**2 + (c/a)**2))
    end select

  end function vector_magnitude

end module ustruc_comp_class
