module ds_utilities

  use truchas_logging_services, only: ds_info => TLS_info, ds_warn => TLS_warn
  implicit none
  private

  public :: ds_info, ds_warn, ds_halt, ds_check_stat

contains

  subroutine ds_halt (errmsg)
    use truchas_logging_services, only: TLS_fatal
    character(*), intent(in) :: errmsg
    call TLS_fatal (errmsg)
  end subroutine ds_halt

  subroutine ds_check_stat (stat, errmsg)
    use parallel_communication, only: global_any
    integer, intent(in) :: stat
    character(len=*), intent(in) :: errmsg
    if (global_any(stat /= 0)) call ds_halt (errmsg)
  end subroutine ds_check_stat

end module ds_utilities
