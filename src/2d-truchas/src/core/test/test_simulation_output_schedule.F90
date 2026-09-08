program test_simulation_output_schedule

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use simulation_output_schedule, only: get_output_times
  implicit none

  type(parameter_list) :: params
  type(parameter_list) :: default_params
  type(parameter_list), pointer :: schedule
  type(parameter_list), pointer :: default_schedule
  real(r8), allocatable :: output_times(:)
  character(:), allocatable :: errmsg
  integer :: stat

  call params%set_path('sim-control')
  schedule => params%sublist('output-times')
  call schedule%set('times', [0.0_r8, 1.0_r8, 2.0_r8])
  call schedule%set('subintervals', [2, 4])
  call get_output_times(params, 0.0_r8, output_times, stat, errmsg)
  if (stat /= 0) error stop errmsg
  call require(size(output_times) == 6, 'unexpected number of expanded output times')
  call require(all(abs(output_times - [0.5_r8, 1.0_r8, 1.25_r8, 1.5_r8, 1.75_r8, 2.0_r8]) < 1.0e-14_r8), &
      'incorrect expanded output times')

  call params%set('final-time', 3.0_r8)
  call get_output_times(params, 1.0_r8, output_times, stat, errmsg)
  if (stat /= 0) error stop errmsg
  call require(size(output_times) == 5, 'unexpected number of filtered output times')
  call require(all(abs(output_times - [1.25_r8, 1.5_r8, 1.75_r8, 2.0_r8, 3.0_r8]) < 1.0e-14_r8), &
      'incorrect filtered output times')

  call default_params%set_path('sim-control')
  default_schedule => default_params%sublist('output-times')
  call default_schedule%set('times', [0.0_r8, 1.0_r8, 2.0_r8])
  call get_output_times(default_params, 0.0_r8, output_times, stat, errmsg)
  if (stat /= 0) error stop errmsg
  call require(size(output_times) == 2, 'unexpected number of default output times')
  call require(all(abs(output_times - [1.0_r8, 2.0_r8]) < 1.0e-14_r8), &
      'incorrect default output times')

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not.condition) error stop message
  end subroutine require

end program test_simulation_output_schedule
