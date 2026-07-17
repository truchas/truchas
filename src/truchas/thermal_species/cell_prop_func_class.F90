!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cell_prop_func_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  implicit none
  private

  public :: cell_prop_func

  type, abstract :: cell_prop_func
  contains
    procedure(compute_value_t_iface), deferred :: compute_value_t
    procedure(compute_value_c_iface), deferred :: compute_value_c
    procedure(compute_value_tc_iface), deferred :: compute_value_tc
    procedure(compute_value_cell_iface), deferred :: compute_value_cell
    procedure(compute_value_cell_t_iface), deferred :: compute_value_cell_t
    procedure(compute_value_cell_tc_iface), deferred :: compute_value_cell_tc
    procedure(compute_deriv_t_iface), deferred :: compute_deriv_t
    procedure(compute_deriv_c_iface), deferred :: compute_deriv_c
    procedure(compute_deriv_tc_iface), deferred :: compute_deriv_tc
    generic :: compute_value => compute_value_t, compute_value_c, compute_value_tc, &
                                compute_value_cell, compute_value_cell_t, &
                                compute_value_cell_tc
    generic :: compute_deriv => compute_deriv_t, compute_deriv_c, compute_deriv_tc
  end type cell_prop_func

  abstract interface

    subroutine compute_value_t_iface(this, temp, value)
      import :: cell_prop_func, r8
      class(cell_prop_func), intent(in) :: this
      real(r8), intent(in)  :: temp(:)
      real(r8), intent(out) :: value(:)
    end subroutine compute_value_t_iface

    subroutine compute_value_c_iface(this, state, value)
      import :: cell_prop_func, r8
      class(cell_prop_func), intent(in) :: this
      real(r8), intent(in)  :: state(:,:)
      real(r8), intent(out) :: value(:)
    end subroutine compute_value_c_iface

    subroutine compute_value_tc_iface(this, temp, conc, value)
      import :: cell_prop_func, r8
      class(cell_prop_func), intent(in) :: this
      real(r8), intent(in)  :: temp(:)
      real(r8), intent(in)  :: conc(:,:)
      real(r8), intent(out) :: value(:)
    end subroutine compute_value_tc_iface

    subroutine compute_value_cell_iface(this, n, state, value)
      import :: cell_prop_func, r8
      class(cell_prop_func), intent(in) :: this
      integer, intent(in) :: n
      real(r8), intent(in)  :: state(:)
      real(r8), intent(out) :: value
    end subroutine compute_value_cell_iface

    subroutine compute_value_cell_t_iface(this, n, temp, value)
      import :: cell_prop_func, r8
      class(cell_prop_func), intent(in) :: this
      integer, intent(in) :: n
      real(r8), intent(in)  :: temp
      real(r8), intent(out) :: value
    end subroutine compute_value_cell_t_iface

    subroutine compute_value_cell_tc_iface(this, n, temp, conc, value)
      import :: cell_prop_func, r8
      class(cell_prop_func), intent(in) :: this
      integer, intent(in) :: n
      real(r8), intent(in)  :: temp
      real(r8), intent(in)  :: conc(:)
      real(r8), intent(out) :: value
    end subroutine compute_value_cell_tc_iface

    subroutine compute_deriv_t_iface(this, temp, index, value)
      import :: cell_prop_func, r8
      class(cell_prop_func), intent(in) :: this
      real(r8), intent(in)  :: temp(:)
      integer,  intent(in)  :: index
      real(r8), intent(out) :: value(:)
    end subroutine compute_deriv_t_iface

    subroutine compute_deriv_c_iface(this, state, index, value)
      import :: cell_prop_func, r8
      class(cell_prop_func), intent(in) :: this
      real(r8), intent(in)  :: state(:,:)
      integer,  intent(in)  :: index
      real(r8), intent(out) :: value(:)
    end subroutine compute_deriv_c_iface

    subroutine compute_deriv_tc_iface(this, temp, conc, index, value)
      import :: cell_prop_func, r8
      class(cell_prop_func), intent(in) :: this
      real(r8), intent(in)  :: temp(:)
      real(r8), intent(in)  :: conc(:,:)
      integer,  intent(in)  :: index
      real(r8), intent(out) :: value(:)
    end subroutine compute_deriv_tc_iface

  end interface

end module cell_prop_func_class
