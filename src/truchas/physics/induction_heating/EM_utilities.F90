module EM_utilities
  
  use truchas_logging_services, only: EM_info => TLS_info
  implicit none
  private
  
  public :: EM_error, EM_warning, EM_info
  
contains
  
  subroutine EM_error (name, message)
    use truchas_logging_services, only: TLS_fatal
    character(len=*), intent(in) :: name, message
    call TLS_fatal (message)
  end subroutine EM_error
  
  subroutine EM_warning (name, message)
    use truchas_logging_services, only: TLS_warn
    character(*), intent(in) :: name, message
    call TLS_warn (trim(name) // ': ' // message)
  end subroutine EM_warning
  
end module EM_utilities
    
