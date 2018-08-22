!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_scalar_func_copy
  use kinds
  use scalar_func_containers
  use scalar_func_class
  use scalar_func_factories
  use parallel_util_module, only: parallel_init
  use parallel_communication
  use pgslib_module
  use parameter_list_type
  implicit none

  character(PGSLib_CL_MAX_TOKEN_LENGTH), pointer :: argv(:) => null()
  type(scalar_func_box) :: sf_box(4)
  type(parameter_list) :: sf_const, sf_poly, sf_tab, sf_ramp
  real(r8) :: results(4), in(4)
  integer :: i

  call parallel_init(argv)
  !call init_parallel_communication()

  call sf_const%set("type", "constant")
  call sf_const%set("value", -1.0_r8)

  call sf_poly%set("type", "polynomial")
  call sf_poly%set("polynomial coefficients", [1.0_r8, 2.0_r8, 3.0_r8])
  call sf_poly%set("polynomial exponents", [1, 2, 3])

  call sf_tab%set("type", "tabular")
  call sf_tab%set("tabular x", [1.0_r8, 2.0_r8, 3.0_r8])
  call sf_tab%set("tabular y", [10.0_r8, 11.0_r8, 12.0_r8])

  call sf_ramp%set("type", "smooth ramp")
  call sf_ramp%set("begin time", 0.0_r8)
  call sf_ramp%set("begin value", 1000.0_r8)
  call sf_ramp%set("end time", 7.0_r8)
  call sf_ramp%set("end value", 4.0_r8)

  call alloc_scalar_func(sf_box(1)%f, sf_const)
  call alloc_scalar_func(sf_box(2)%f, sf_poly)
  call alloc_scalar_func(sf_box(3)%f, sf_tab)
  call alloc_scalar_func(sf_box(4)%f, sf_ramp)

  in = [0.0_r8, 4.0_r8, 1.2_r8, 1.2_r8]
  do i = 1, size(in)
    results(i) = sf_box(i)%eval(in(i:i))
  end do

  !  call sf_copy(sf_box, in, results)
  call sf_move(sf_box, in, results)

  call halt_parallel_communication()
contains

!!$  subroutine sf_copy(sf_box, in, results)
!!$    type(scalar_func_box), intent(inout) :: sf_box(:)
!!$    real(r8), intent(in) :: in(:), results(:)
!!$    !-
!!$    real(r8) :: copy_results(size(in))
!!$    type(scalar_func_box) :: sf_box_copy(size(in))
!!$    integer :: i
!!$
!!$    do i = 1, size(in)
  ! this line fails with:
  ! In an intrinsic assignment statement, variable shall not be polymorphic.
!!$      sf_box_copy(i)%f = sf_box(i)%f
!!$      copy_results(i) = sf_box_copy(i)%f%eval(in(i:i))
!!$    end do
!!$
!!$    do i = 1, size(in)
!!$      if (abs(copy_results(i)-results(i)) > epsilon(1.0_r8)) then
!!$        write(*,*) "ERROR: function ", i
!!$        write(*,*) "ERROR: baseline results: ", results(i)
!!$        write(*,*) "ERROR: copied resuts   : ", copy_results(i)
!!$      end if
!!$    end do
!!$  end subroutine sf_copy


  subroutine sf_move(sf_box, in, results)
    type(scalar_func_box), intent(inout) :: sf_box(:)
    real(r8), intent(in) :: in(:), results(:)
    !-
    real(r8) :: move_results(size(in))
    type(scalar_func_box) :: sf_box_move(size(in))
    integer :: i

    do i = 1, size(in)
      call move_alloc(sf_box(i)%f, sf_box_move(i)%f)
      move_results(i) = sf_box_move(i)%f%eval(in(i:i))
    end do

    do i = 1, size(in)
      if (abs(move_results(i)-results(i)) > epsilon(1.0_r8)) then
        write(*,*) "ERROR: function ", i
        write(*,*) "ERROR: baseline results: ", results(i)
        write(*,*) "ERROR: moved resuts   : ", move_results(i)
      end if
    end do
  end subroutine sf_move

end program test_scalar_func_copy
