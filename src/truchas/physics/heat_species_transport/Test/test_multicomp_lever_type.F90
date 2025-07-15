program test_alloy_lever_rule_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use parameter_list_json
  use material_class
  use alloy_lever_rule_type
  
  type(parameter_list), pointer :: matl_params, foo_params
  character(:), allocatable :: errmsg
  class(material), allocatable :: matl
  type(alloy_lever_rule) :: foo
  integer :: stat, j
  real(r8) :: T, g
  
  
  character(*), parameter :: matl_input = &
      '{"phases":{"solid":{"enthalpy":{"type":"polynomial","poly-coef":[1.0],"poly-powers":[1]}},&
      &"liquid":{"enthalpy":{"type":"polynomial","poly-coef":[2.0,1.0],"poly-powers":[0,1]}}},&
      &"phase-changes":{"solid:liquid":{"solidus-temp":1.0,"liquidus-temp":2.0}}}'

  call parameter_list_from_json_string(matl_input, matl_params, errmsg)
  if (.not.associated(matl_params)) then
    write(*,*) errmsg
    stop 1
  end if

  block
    use material_factory, only: alloc_material
    call alloc_material(matl, 'stuff', matl_params, stat, errmsg)
    if (stat /= 0) then
      write(*,*) errmsg
      stop 2
    end if
  end block

  call ternary

contains

  subroutine binary

    character(*), parameter :: foo_input = &
        '{"num-comp":1,"temp-fusion":2.0,"temp-eutectic":1.0,"liq-slope":[-1.0],"part-coef":[0.5]}'

    call parameter_list_from_json_string(foo_input, foo_params, errmsg)
    if (.not.associated(foo_params)) then
      write(*,*) errmsg
      stop 3
    end if
    call foo%init(matl, foo_params, stat, errmsg)
    if (stat /= 0) then
      write(*,*) errmsg
      stop 4
    end if

    call foo%solve(4.0_r8, [0.5_r8], 0.0_r8, 1.0_r8, T, g)
    write(*,*) T, g

    call foo%solve(1.5_r8, [0.1_r8], 0.0_r8, 5.0_r8, T, g)
    write(*,*) T, g

    call foo%solve(2.5_r8, [0.375_r8], 0.0_r8, 5.0_r8, T, g)
    write(*,*) T, g

    call foo%solve(2.0_r8, [0.75_r8], 0.0_r8, 5.0_r8, T, g)
    write(*,*) T, g

    call foo%solve(1.5_r8, [0.75_r8], 0.0_r8, 5.0_r8, T, g)
    write(*,*) T, g

    call foo%solve(0.5_r8, [0.75_r8], 0.0_r8, 5.0_r8, T, g)
    write(*,*) T, g
  
  !  do j = 0, 32
  !    call foo%solve(0.125_r8*j, [1.0_r8], 0.0_r8, 5.0_r8, T, g)
  !    write(*,*) T, g
  !  end do

  end subroutine

  subroutine ternary

    character(*), parameter :: foo_input = &
        '{"num-comp":2,"temp-fusion":2.0,"temp-eutectic":1.0,&
        &"liq-slope":[-1.0,-1.0],"part-coef":[0.5,0.5]}'

    call parameter_list_from_json_string(foo_input, foo_params, errmsg)
    if (.not.associated(foo_params)) then
      write(*,*) errmsg
      stop 3
    end if
    call foo%init(matl, foo_params, stat, errmsg)
    if (stat /= 0) then
      write(*,*) errmsg
      stop 4
    end if

    call foo%solve(4.0_r8, [0.25_r8, 0.25_r8], 0.0_r8, 1.0_r8, T, g)
    write(*,*) T, g

    call foo%solve(2.5_r8, [0.375_r8,0.0_r8], 0.0_r8, 5.0_r8, T, g)
    write(*,*) T, g

    call foo%solve(2.5_r8, [0.0_r8, 0.375_r8], 0.0_r8, 5.0_r8, T, g)
    write(*,*) T, g

    call foo%solve(2.5_r8, [0.1875_r8, 0.1875_r8], 0.0_r8, 5.0_r8, T, g)
    write(*,*) T, g

  end subroutine
  
end program
