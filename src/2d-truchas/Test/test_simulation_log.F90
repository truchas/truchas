program test_simulation_log

  use simulation_log_type
  use parallel_communication
  implicit none

  type(simulation_log) :: simlog
  character(:), allocatable :: errmsg
  integer :: stat

  call init_parallel_communication

  call simlog%init('simulation_log_test.log', stat, errmsg, LOG_NORMAL, .false.)
  call require(stat == 0, 'could not initialize normal log')
  call simlog%info('normal message')
  call simlog%info('detail message', LOG_DETAIL)
  call simlog%warn('warning message')
  call simlog%error('error message')
  call simlog%close()

  call require(log_contains('simulation_log_test.log', 'normal message'), 'normal message missing')
  call require(.not.log_contains('simulation_log_test.log', 'detail message'), 'detail message was not filtered')
  call require(log_contains('simulation_log_test.log', 'Warning: warning message'), 'warning message missing')
  call require(log_contains('simulation_log_test.log', 'ERROR: error message'), 'error message missing')

  call simlog%init('simulation_log_test.log', stat, errmsg, LOG_DETAIL, .false.)
  call require(stat == 0, 'could not initialize detail log')
  call simlog%info('detail message', LOG_DETAIL, terminal=.false.)
  call simlog%close()
  call require(log_contains('simulation_log_test.log', 'detail message'), 'detail message missing')

  call halt_parallel_communication

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message

    if (global_any(.not.condition)) then
      if (is_IOP) write(*, '(2a)') 'ERROR: ', message
      error stop 1
    end if
  end subroutine


  logical function log_contains(filename, text) result(found)
    character(*), intent(in) :: filename, text

    character(256) :: line
    integer :: ios, unit

    if (is_IOP) then
      found = .false.
      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios == 0) then
        do
          read(unit, '(a)', iostat=ios) line
          if (ios /= 0) exit
          if (index(line, text) > 0) then
            found = .true.
            exit
          end if
        end do
        close(unit)
      end if
    end if
    call broadcast(found)
  end function

end program test_simulation_log
