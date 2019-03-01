!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module pcsr_precon_factory

  use pcsr_precon_class
  use pcsr_matrix_type
  use parameter_list_type
  implicit none
  private

  public :: alloc_pcsr_precon

contains

  subroutine alloc_pcsr_precon (this, A, params)

    use pcsr_precon_ssor_type
    use pcsr_precon_boomer_type
    use truchas_logging_services
    use string_utilities, only: raise_case

    class(pcsr_precon), allocatable, intent(out) :: this
    type(pcsr_matrix), intent(in), target :: A
    type(parameter_list) :: params

    integer :: stat
    character(:), allocatable :: method, errmsg
    type(parameter_list), pointer :: plist

    call params%get ('method', method, stat=stat, errmsg=errmsg)
    if (stat /= 0) call TLS_fatal ('ALLOC_PCSR_PRECON: ' // errmsg)

    select case (raise_case(method))
    case ('SSOR')
      allocate (pcsr_precon_ssor :: this)
    case ('BOOMERAMG')
      allocate (pcsr_precon_boomer :: this)
    case default
      call TLS_fatal ('ALLOC_PCSR_PRECON: unknown "method": '// method)
    end select

    if (params%is_sublist('params')) then
      plist => params%sublist('params')
      call this%init (A, plist)
    else
      call TLS_fatal ('ALLOC_PCSR_PRECON: missing "params" sublist parameter')
    end if

  end subroutine alloc_pcsr_precon

end module pcsr_precon_factory
