!!
!! USTRUC_CORE_TYPE
!!
!! A concrete implementation the abstract base class USTRUC_CORE. This
!! implementation defines the non-optional core microstructure analysis
!! component that is referenced by the optional analysis components.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module ustruc_core_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ustruc_analysis_class
  use truchas_logging_services
  implicit none
  private

  public :: new_ustruc_core

  type, extends(ustruc_analysis), public :: ustruc_core
    real(r8) :: t = 0.0_r8
    real(r8), allocatable :: temp(:)
    real(r8), allocatable :: temp_grad(:,:)
    real(r8), allocatable :: frac(:)
    logical,  allocatable :: invalid(:)
  contains
    procedure, private :: init
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
  end type ustruc_core

#ifndef INTEL_BUG20200721
  !! Number of bytes (per cell) of internal state for serialization/deserialization
  type(ustruc_core), allocatable :: dummy  ! only use is in the following parameter declaration
  integer, parameter :: NBYTES = 0
#endif

contains

  function new_ustruc_core(n, params) result(this)
    use parameter_list_type
    integer, intent(in) :: n
    type(parameter_list) :: params
    type(ustruc_core), pointer :: this
    allocate(this)
    call this%init(n, params)
  end function

  subroutine init(this, n, params)
    use parameter_list_type
    class(ustruc_core), intent(out) :: this
    integer, intent(in) :: n
    type(parameter_list) :: params
    ASSERT(n >= 0)
    this%n = n
    allocate(this%temp(n), this%temp_grad(3,n), this%frac(n), this%invalid(n))
  end subroutine

  subroutine set_state(this, t, temp, temp_grad, frac, invalid)

    class(ustruc_core), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)

    ASSERT(size(temp) == this%n)
    ASSERT(size(temp_grad,1) == 3)
    ASSERT(size(temp_grad,2) == this%n)
    ASSERT(size(frac) == this%n)
    ASSERT(size(invalid) == this%n)

    this%t = t
    this%temp = temp
    this%temp_grad = temp_grad
    this%frac = frac
    this%invalid = invalid

  end subroutine

  subroutine update_state(this, t, temp, temp_grad, frac, invalid)

    class(ustruc_core), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)

    ASSERT(size(temp) == this%n)
    ASSERT(size(temp_grad,1) == 3)
    ASSERT(size(temp_grad,2) == this%n)
    ASSERT(size(frac) == this%n)
    ASSERT(size(invalid) == this%n)

    this%t = t
    this%temp = temp
    this%temp_grad = temp_grad
    this%frac = frac
    this%invalid = invalid

  end subroutine

  subroutine get_comp_list(this, list)
    class(ustruc_core), intent(in) :: this
    integer, allocatable, intent(out) :: list(:)
    list = [1]
  end subroutine

  !! USTRUC_CORE is provides no public data and is end-of-the line for data
  !! requests. If execution reaches any of these GET procedures the specified
  !! data name was not recognized, and a fatal error message is written and
  !! execution is halted.

  logical function has(this, name)
    class(ustruc_core), intent(in) :: this
    character(*), intent(in) :: name
    has = .false.
  end function

  subroutine getl1(this, name, array)
    class(ustruc_core), intent(in) :: this
    character(*), intent(in) :: name
    logical, intent(out) :: array(:)
    call TLS_fatal('USTRUC_ANALYSIS%GET: unknown name: ' // name)
  end subroutine

  subroutine geti1(this, name, array, invalid)
    class(ustruc_core), intent(in) :: this
    character(*), intent(in) :: name
    integer, intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    call TLS_fatal('USTRUC_ANALYSIS%GET: unknown name: ' // name)
  end subroutine

  subroutine getr1(this, name, array, invalid)
    class(ustruc_core), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    call TLS_fatal('USTRUC_ANALYSIS%GET: unknown name: ' // name)
  end subroutine

  subroutine getr2(this, name, array, invalid)
    class(ustruc_core), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:,:)
    logical, intent(out), optional :: invalid(:)
    call TLS_fatal('USTRUC_ANALYSIS%GET: unknown name: ' // name)
  end subroutine


  subroutine serialize(this, cid, array)

    use,intrinsic :: iso_fortran_env, only: int8
    use serialization_tools, only: copy_to_bytes

    class(ustruc_core), intent(in) :: this
    integer, intent(in) :: cid
    integer(int8), allocatable, intent(out) :: array(:,:)

    integer :: j, offset
#ifdef INTEL_BUG20200721
    integer :: NBYTES
    NBYTES = 0
#endif

    if (cid == 1) then
      allocate(array(NBYTES,this%n))
    end if

  end subroutine

  subroutine deserialize(this, cid, array)

    use,intrinsic :: iso_fortran_env, only: int8
    use serialization_tools, only: copy_from_bytes

    class(ustruc_core), intent(inout) :: this
    integer, intent(in) :: cid
    integer(int8), intent(in) :: array(:,:)

    integer :: j, offset
#ifdef INTEL_BUG20200721
    integer :: NBYTES
    NBYTES = 0
#endif

    if (cid == 1) then
      INSIST(size(array,1) == NBYTES)
      INSIST(size(array,2) == this%n)
    end if

  end subroutine

end module ustruc_core_type
