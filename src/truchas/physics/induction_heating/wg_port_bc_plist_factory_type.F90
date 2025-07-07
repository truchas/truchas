!!
!! WG_PORT_BC_PLIST_FACTORY_TYPE
!!
!! This module defines a derived type that encapsulates the time-dependent
!! power of waveguide port feed BCs with a method for setting the power in
!! the BC parameter list at a given time. This is intended for use in
!! microwave heating simulations.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! June 2025
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module wg_port_bc_plist_factory_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  implicit none
  private

  type, public :: wg_port_bc_plist_factory
    private
    real(r8), allocatable, public :: times(:) ! protected
    real(r8), allocatable :: powers(:,:)
    character(:), allocatable :: wg_port_bc(:)
    type(parameter_list), pointer :: bc_params => null() ! unowned reference
  contains
    procedure :: init
    procedure :: set_plist_power
    procedure :: power_differs
    procedure :: power_is_scaled
    procedure :: power_is_zero
    procedure :: power_data
    procedure, private :: data_index
  end type

contains

  subroutine init(this, params, bc_params, stat, errmsg)

    class(wg_port_bc_plist_factory), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    type(parameter_list), intent(inout), target :: bc_params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    character(:), allocatable :: bc_type
    type(parameter_list), pointer :: plist

    this%bc_params => bc_params

    call params%get('wg-port-bc', this%wg_port_bc, stat, errmsg)
    if (stat /= 0) return

    !! Check that the specified wg-port BC are actually defined
    do n = 1, size(this%wg_port_bc)
      if (bc_params%is_sublist(this%wg_port_bc(n))) then
        plist => bc_params%sublist(this%wg_port_bc(n))
        call plist%get('type', bc_type, stat, errmsg)
        if (stat /= 0) return
        if (bc_type /= 'wg-port') then
          stat = 1
          errmsg = 'not a "wg-port" type bc: ' // trim(this%wg_port_bc(n))
          return
        end if
      else
        stat = 1
        errmsg = 'unknown electromagnetic BC: ' // trim(this%wg_port_bc(n))
        return
      end if
    end do

    !TODO? check that no other wg-port BC are defined?

    if (params%is_parameter('times')) then
      call params%get('times', this%times, stat, errmsg)
      if (stat /= 0) return
      n = size(this%times)
      if (n > 1) then
        if (any(this%times(2:n) <= this%times(1:n-1))) then
          stat = 1
          errmsg = 'times values not strictly increasing'
          return
        end if
      end if
    else
      allocate(this%times(0))
      n = 0
    end if

    call params%get('powers', this%powers, stat, errmsg)
    if (stat /= 0) return
    if (any(shape(this%powers) /= [size(this%wg_port_bc),n+1])) then
      stat = 1
      errmsg = 'shape of powers is incompatible with wg-port-bc and times'
      return
    end if

  end subroutine

  subroutine set_plist_power(this, t)
    class(wg_port_bc_plist_factory), intent(in) :: this
    real(r8), intent(in) :: t
    integer :: n, i
    type(parameter_list), pointer :: plist
    n = this%data_index(t)
    do i = 1, size(this%wg_port_bc)
      plist => this%bc_params%sublist(trim(this%wg_port_bc(i)))
      call plist%set('power', this%powers(i,n))
    end do
  end subroutine

  function power_data(this, t) result(power)
    class(wg_port_bc_plist_factory), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), allocatable :: power(:)
    power = this%powers(:,this%data_index(t))
  end function

  !! Return true if any of the powers at time T2 differ from P1.

  pure logical function power_differs(this, p1, t2) result(differs)
    class(wg_port_bc_plist_factory), intent(in) :: this
    real(r8), intent(in) :: p1(:), t2
    integer :: n2
    n2 = this%data_index(t2)
    associate (p2 => this%powers(:,n2))
      differs = any(p1 /= p2)
    end associate
  end function

  !! Return true if the powers at time T2 are a multiple of P1
  !! and return the scale factor in the argument SCF.

  logical function power_is_scaled(this, p1, t2, scf) result(scaled)
    class(wg_port_bc_plist_factory), intent(in) :: this
    real(r8), intent(in) :: p1(:), t2
    real(r8), intent(out) :: scf
    integer :: n2
    real(r8) :: a, err
    n2 = this%data_index(t2)
    associate (p2 => this%powers(:,n2))
      ! Best scale factor in least square sense
      a = norm2(p1)
      if (a > 0) then
        scf = dot_product(p1, p2) / a**2
      else
        scf = 0.0_r8
      end if
      err = norm2(p2 - scf*p1)  ! l2 error in best scaling
      scaled = (err <= a*1.0d-6)
    end associate
  end function

  !! Return true if the powers are 0 at the given time T.

  pure logical function power_is_zero(this, t)
    class(wg_port_bc_plist_factory), intent(in) :: this
    real(r8), intent(in) :: t
    integer :: n
    n = this%data_index(t)
    power_is_zero = all(this%powers(:,n) == 0.0_r8)
  end function

  !! Auxiliary function that locates the time interval containing time T.

  pure integer function data_index(this, t) result(n)
    class(wg_port_bc_plist_factory), intent(in) :: this
    real(r8), intent(in) :: t
    n = size(this%times)
    do while (n > 0)
      if (t >= this%times(n)) exit
      n = n - 1
    end do
    n = n + 1
  end function

end module wg_port_bc_plist_factory_type
