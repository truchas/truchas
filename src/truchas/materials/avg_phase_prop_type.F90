!!
!! AVG_PHASE_PROP_TYPE
!!
!! A freqently encountered pattern is the evaluation of a linear combination
!! function of the form $f(x, w) = \sum_{i=1}^m w_i f_i(x)$ for a fixed set
!  of functions $f_i$. This module defines the type AVG_PHASE_PROP that
!! encapsulates this linear combination for the case of SCALAR_FUNC class
!! phase property functions
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module avg_phase_prop_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scalar_func_class
  implicit none
  private

  type :: scalar_func_box
    class(scalar_func), allocatable :: func
  end type

  type, public :: avg_phase_prop
    private
    type(scalar_func_box), allocatable :: phase(:)
  contains
    generic :: init => init_list, init_all
    generic :: compute_value => compute_value_1
    procedure, private :: init_list, init_all
    procedure, private :: compute_value_1
  end type

contains

  subroutine init_list(this, name, pids, model, stat, errmsg, void_value)
    use material_model_type
    use scalar_func_factories, only: alloc_const_scalar_func
    class(avg_phase_prop), intent(out) :: this
    character(*), intent(in) :: name
    integer, intent(in) :: pids(:)
    type(material_model), intent(in) :: model
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), intent(in), optional :: void_value
    integer :: n
    stat = 0
    allocate(this%phase(size(pids)))
    do n = 1, size(pids)
      if (pids(n) == model%void_index) then
        if (present(void_value)) then
          call alloc_const_scalar_func(this%phase(n)%func, void_value)
        else
          call alloc_const_scalar_func(this%phase(n)%func, 0.0_r8)
        end if
      else
        call model%get_phase_prop(pids(n), name, this%phase(n)%func)
        if (.not.allocated(this%phase(n)%func)) then
          stat = 1
          errmsg = name // ' undefined for phase ' // model%phase_name(pids(n))
          return
        end if
      end if
    end do
  end subroutine

  subroutine init_all(this, name, model, stat, errmsg, void_value)
    use material_model_type
    use scalar_func_factories, only: alloc_const_scalar_func
    class(avg_phase_prop), intent(out) :: this
    character(*), intent(in) :: name
    type(material_model), intent(in) :: model
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), intent(in), optional :: void_value
    integer :: n
    call init_list(this, name, [(n, n=1,model%nphase)], model, stat, errmsg, void_value)
  end subroutine

  subroutine compute_value_1(this, w, state, value)
    class(avg_phase_prop), intent(in) :: this
    real(r8), intent(in) :: w(:), state(:)
    real(r8), intent(out) :: value
    integer :: n
    real(r8) :: tmp2
    ASSERT(size(w) >= size(this%phase))
    tmp2 = 0.0_r8
    do n = 1, size(this%phase)
      if (w(n) /= 0) tmp2 = tmp2 + w(n)*this%phase(n)%func%eval(state)
    end do
    value = tmp2
  end subroutine

end module avg_phase_prop_type
