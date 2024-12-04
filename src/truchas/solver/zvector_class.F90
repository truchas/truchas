!!
!! ZVECTOR_CLASS
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! March 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! CLASS METHODS
!!
!!  CLASS(ZVECTOR) :: V
!!
!!  CALL V%CLONE(C [,N]) creates a clone C of V. The allocatable CLASS(ZVECTOR)
!!  variable C is allocated with the same dynamic type and internal structure
!!  as V, but the values of its vector elements are not defined (see COPY).
!!  C may be a rank-1 array, in which case N must also be specified, and its
!!  value is the size to which C is allocated.
!!
!!    This is a generic interface with deferred specific procedures CLONE1 and
!!    CLONE2 that subclasses must implement.
!!
!!  CALL V%COPY(SRC) copies the vector element values from the CLASS(ZVECTOR)
!!  variable SRC to V. V and SRC must have the same dynamic type and be
!!  compatibly defined, as is the case if V were a clone of SRC.
!!
!!    This is implemented using the non-virtual interface (NVI) pattern.
!!    Subclasses implement the private deferred procedure COPY_ and are
!!    assured that SRC and V will have the same dynamic type.
!!
!!  CALL V%SETVAL(VAL) sets the vector element values to the complex scalar VAL.
!!
!!    This is a deferred procedure implemented by subclasses.
!!
!!  CALL V%SCALE(A) multiplies the vector elements of V by the complex scalar A.
!!
!!    This is a deferred procedure implemented by subclasses.
!!
!!  CALL V%UPDATE(A, X [,B [,Y [,C]]]) performs SAXPY-like vector updates:
!!    V <-- A*X + V
!!    V <-- A*X + B*V
!!    V <-- A*X + B*Y + V
!!    V <-- A*X + B*Y + C*V
!!  V, X, and Y must all have the same dynamic type and be compatibly defined.
!!  A, B, and C are complex scalars.
!!
!!    This is a generic procedure with specific procedures implemented using
!!    the NVI pattern. Subclasses implement the private deferred procedures
!!    UPDATE1_, UPDATE2_, UPDATE3_, UPDATE4_ corresponding to the different
!!    versions of the update, and are assured that V, X, and Y will all have
!!    the same dynamic type.
!!
!!  V%DOTC(Y) is the dot product of the conjugate of V with the CLASS(ZVECTOR)
!!  variable Y. V and Y must have the same dynamic type and be compatibly
!!  defined. Note that this form matches the behavior of DOT_PRODUCT in the
!!  case of complex arguments, and the behavior of the BLAS CDOTC function.
!!
!!    This is a deferred procedure implemented by subclasses.
!!
!!  V%NORM1() is the l1 norm of V; i.e., SUM(ABS(V)) where the sum runs over
!!  the elements that would be included in a reduction operation like DOTC.
!!
!!    This is a deferred procedure implemented by subclasses.
!!
!!  V%NORM2() is the l2 norm of V; i.e., SQRT(V%DOTC(V))
!!
!!    This is a deferred procedure implemented by subclasses.
!!
!!  V%NORM_MAX() is the max norm of V; i.e., MAXVAL(V) where the max runs over
!!  the elements that would be included in a reduction operation.
!!
!!    This is a deferred procedure implemented by subclasses.
!!
!!  V%CHECKSUM([FULL]) is a string that is a checksum of the vector elements
!!  of V. The default behavior is to include only the significant elements in
!!  the checksum. These are the elements that would be included in a reduction
!!  operation like NORM2, for example. Otherwise, when the optional argument
!!  FULL is present with value TRUE, the checksum may include additional data,
!!  like duplicated ghost values, at the discretion of the implementation.
!!

module zvector_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  type, abstract, public :: zvector
  contains
    generic :: clone => clone1, clone2
    procedure(clone1), deferred :: clone1 ! private except to subclass
    procedure(clone2), deferred :: clone2 ! private except to subclass
    procedure :: copy
    procedure(setval), deferred :: setval
    procedure(setzero), deferred :: setzero
    procedure :: conjg
    procedure(scale), deferred :: scale
    generic :: update => update1, update2, update3, update4
    procedure, private :: update1, update2, update3, update4
    procedure :: dotc
    procedure(norm), deferred :: norm1
    procedure(norm), deferred :: norm2
    procedure(norm), deferred :: norm_max
    procedure(checksum), deferred :: checksum
    ! remaining are all private except to subclass
#if defined(INTEL_SR04329123) || defined(__GFORTRAN__)
    procedure(cpy), deferred :: copy_
    procedure(upd1), deferred :: update1_
    procedure(upd2), deferred :: update2_
    procedure(upd3), deferred :: update3_
    procedure(upd4), deferred :: update4_
    procedure(dt), deferred :: dotc_
#else
    procedure(copy), deferred :: copy_
    procedure(update1), deferred :: update1_
    procedure(update2), deferred :: update2_
    procedure(update3), deferred :: update3_
    procedure(update4), deferred :: update4_
    procedure(dotc), deferred :: dotc_
#endif
    procedure(conjg1), deferred :: conjg1_
    procedure(conjg2), deferred :: conjg2_
  end type zvector

  abstract interface
    subroutine clone1(this, clone)
      import zvector
      class(zvector), intent(in)  :: this
      class(zvector), allocatable, intent(out) :: clone
    end subroutine
    subroutine clone2(this, clone, n)
      import zvector
      class(zvector), intent(in)  :: this
      class(zvector), allocatable, intent(out) :: clone(:)
      integer, intent(in) :: n
    end subroutine
  end interface

  abstract interface
    subroutine setval(this, val)
      import zvector, r8
      class(zvector), intent(inout) :: this
      complex(r8), intent(in) :: val
    end subroutine
    subroutine setzero(this)
      import zvector
      class(zvector), intent(inout) :: this
    end subroutine
  end interface

  abstract interface
    subroutine conjg1(this)
      import zvector
      class(zvector), intent(inout) :: this
    end subroutine
    subroutine conjg2(this, src)
      import zvector
      class(zvector), intent(inout) :: this
      class(zvector), intent(in) :: src
    end subroutine
  end interface

  abstract interface
    subroutine scale(this, a)
      import zvector, r8
      class(zvector), intent(inout) :: this
      complex(r8), intent(in) :: a
    end subroutine
  end interface

  abstract interface
    function norm(this)
      import zvector, r8
      class(zvector), intent(in) :: this
      real(r8) :: norm
    end function
  end interface

  abstract interface
    function checksum(this, full) result(string)
      import zvector
      class(zvector), intent(in) :: this
      logical, intent(in), optional :: full ! default is FALSE
      character(:), allocatable :: string
    end function
  end interface

#if defined(INTEL_SR04329123) || defined(__GFORTRAN__)
  abstract interface
    subroutine cpy(dest, src)
      import zvector
      class(zvector), intent(inout) :: dest
      class(zvector), intent(in) :: src
    end subroutine cpy
    subroutine upd1(this, a, x)
      import zvector, r8
      class(zvector), intent(inout) :: this
      class(zvector), intent(in) :: x
      complex(r8), intent(in) :: a
    end subroutine
    subroutine upd2(this, a, x, b)
      import zvector, r8
      class(zvector), intent(inout) :: this
      class(zvector), intent(in) :: x
      complex(r8), intent(in) :: a, b
    end subroutine
    subroutine upd3(this, a, x, b, y)
      import zvector, r8
      class(zvector), intent(inout) :: this
      class(zvector), intent(in) :: x, y
      complex(r8), intent(in) :: a, b
    end subroutine
    subroutine upd4(this, a, x, b, y, c)
      import zvector, r8
      class(zvector), intent(inout) :: this
      class(zvector), intent(in) :: x, y
      complex(r8), intent(in) :: a, b, c
    end subroutine
    function dt(x, y)
      import zvector, r8
      class(zvector), intent(in) :: x, y
      complex(r8) :: dt
    end function
  end interface
#endif

contains

  recursive subroutine conjg(this, src)
    class(zvector), intent(inout) :: this
    class(zvector), intent(in), optional :: src
    if (present(src)) then
      if (same_type_as(this, src)) then
        call this%conjg2_(src)
      else
        error stop 'incompatible arguments to ZVECTOR%CONJG'
      end if
    else
      call this%conjg1_()
    end if
  end subroutine

  recursive subroutine copy(dest, src)
    class(zvector), intent(inout) :: dest
    class(zvector), intent(in) :: src
    if (same_type_as(dest, src)) then
      call dest%copy_(src)
    else
      error stop 'incompatible arguments to ZVECTOR%COPY'
    end if
  end subroutine

  recursive function dotc(x, y)
    class(zvector), intent(in) :: x, y
    complex(r8) :: dotc
    if (same_type_as(x, y)) then
      dotc = x%dotc_(y)
    else
      error stop 'incompatible arguments to ZVECTOR%DOTC'
    end if
  end function

  !! Conventional SAXPY procedure: y <-- y + a*x
  recursive subroutine update1(this, a, x)
    class(zvector), intent(inout) :: this
    class(zvector), intent(in) :: x
    complex(r8), intent(in) :: a
    if (a == 0.0_r8) return
    if (same_type_as(this, x)) then
      call this%update1_(a, x)
    else
      error stop 'incompatible arguments to ZVECTOR%UPDATE'
    end if
  end subroutine

  !! SAXPY-like procedure: y <-- a*x + b*y
  recursive subroutine update2(this, a, x, b)
    class(zvector), intent(inout) :: this
    class(zvector), intent(in) :: x
    complex(r8), intent(in) :: a, b
    if (a == 0.0_r8) then
      call this%scale(b)
    else if (same_type_as(this, x)) then
      call this%update2_(a, x, b)
    else
      error stop 'incompatible arguments to ZVECTOR%UPDATE'
    end if
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + z
  recursive subroutine update3(this, a, x, b, y)
    class(zvector), intent(inout) :: this
    class(zvector), intent(in) :: x, y
    complex(r8), intent(in) :: a, b
    if (a == 0.0_r8) then
      call update1(this, b, y)
    else if (b == 0.0_r8) then
      call update1(this, a, x)
    else if (same_type_as(this, x) .and. same_type_as(this, y)) then
      call this%update3_(a, x, b, y)
    else
      error stop 'incompatible arguments to ZVECTOR%UPDATE'
    end if
  end subroutine

  !! SAXPY-like procedure: z <-- a*x + b*y + c*z
  recursive subroutine update4(this, a, x, b, y, c)
    class(zvector), intent(inout) :: this
    class(zvector), intent(in) :: x, y
    complex(r8), intent(in) :: a, b, c
    if (a == 0.0_r8) then
      call update2(this, b, y, c)
    else if (b == 0.0_r8) then
      call update2(this, a, x, c)
    else if (same_type_as(this, x) .and. same_type_as(this, y)) then
      call this%update4_(a, x, b, y, c)
    else
      error stop 'incompatible arguments to ZVECTOR%UPDATE'
    end if
  end subroutine

end module zvector_class
