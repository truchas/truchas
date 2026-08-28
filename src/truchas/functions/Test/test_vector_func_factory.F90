program test_vector_func_factory

  use,intrinsic :: iso_fortran_env, only: r8 => real64, error_unit
  use parameter_list_type
  use parameter_list_json
  use vector_func_class, only: vector_func
  use vector_func_factories, only: alloc_vector_func
  implicit none

  character(*), parameter :: library_path = TEST_VECTOR_FUNC_LIBRARY
  type(parameter_list), pointer :: params
  class(vector_func), allocatable :: f
  character(:), allocatable :: errmsg, input
  real(r8) :: value(2)

  input = '{"type":"library","library-path":"' // library_path // &
      '","library-symbol":"test_vector_func","dimension":2,"parameters":[0.5,-2.0]}'
  call parameter_list_from_json_string(input, params, errmsg)
  if (.not.associated(params)) call fail(errmsg)

  call alloc_vector_func(f, params)
  value = f%eval([1.0_r8, 2.0_r8, 3.0_r8])
  if (any(abs(value - [2.5_r8, -6.0_r8]) > 1.0e-14_r8)) call fail('incorrect library vector function value')

  input = '{"type":"polynomial","poly-coef":[[-0.4330127018922193,-0.25],[0.75,0.4330127018922193]],' // &
      '"poly-powers":[[1,0],[0,1]]}'
  call parameter_list_from_json_string(input, params, errmsg)
  if (.not.associated(params)) call fail(errmsg)

  call alloc_vector_func(f, params)
  value = f%eval([2.0_r8, 3.0_r8, 4.0_r8])
  if (any(abs(value - [1.3839745962155614_r8, 0.7990381056766580_r8]) > 1.0e-14_r8)) &
      call fail('incorrect polynomial vector function value')
  if (abs(f%eval_comp(2, [2.0_r8, 3.0_r8, 4.0_r8]) - 0.7990381056766580_r8) > 1.0e-14_r8) &
      call fail('incorrect polynomial vector function component value')

contains

  subroutine fail(message)
    character(*), intent(in) :: message
    write(error_unit,'(a)') message
    error stop 1
  end subroutine fail

end program test_vector_func_factory
