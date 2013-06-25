program test_truchas_logging_services
  use truchas_logging_services
  implicit none
  call test1
  call test_warn
  call test_error
  call test_info_advance
contains
  subroutine test1
    character(80) :: message(3)
    call TLS_initialize
    call TLS_info ('a normal message')
    call TLS_info ('another normal message', TLS_VERB_NOISY)
    TLS_verbosity = TLS_VERB_NOISY
    call TLS_info ('another normal message', TLS_VERB_NOISY)
    message(1) = 'line 1'
    message(2) = ' line 2'
    message(3) = '  line 3'
    call TLS_info (message)
    call TLS_finalize
  end subroutine test1
  subroutine test_warn
    character(80) :: message(2)
    call TLS_initialize
    call TLS_warn ('you did something questionable')
    message(1) = 'you did something else questionable'
    message(2) = 'and here is some more information'
    call TLS_warn (message)
    call TLS_finalize
  end subroutine test_warn
  subroutine test_error
    character(80) :: message(2)
    call TLS_initialize
    call TLS_error ('you did something bad')
    message(1) = 'you did something else bad'
    message(2) = 'and here is some more information'
    call TLS_error (message)
    call TLS_finalize
  end subroutine test_error
  subroutine test_info_advance
    call TLS_initialize
    TLS_verbosity = TLS_VERB_NORMAL
    call TLS_info ('starting some process...         ', advance=.false., verbosity=TLS_VERB_NOISY)
    call TLS_info ('done.')
    call TLS_finalize
  end subroutine
end program
