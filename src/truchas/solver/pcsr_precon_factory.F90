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

  subroutine alloc_pcsr_precon(this, A, params, stat, errmsg)

    use pcsr_precon_ssor_type
    use pcsr_precon_boomer_type
    use string_utilities, only: raise_case

    class(pcsr_precon), allocatable, intent(out) :: this
    type(pcsr_matrix), intent(in), target :: A
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: context, method
    type(parameter_list), pointer :: plist

    context = 'processing ' // params%path() // ': '

    call params%get('method', method, stat, errmsg, default='BOOMERAMG')
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if

    select case (raise_case(method))
    case ('SSOR')
      allocate (pcsr_precon_ssor :: this)
    case ('BOOMERAMG')
      allocate (pcsr_precon_boomer :: this)
    case default
      stat = 1
      errmsg = context // 'unknown "method": ' // method
      return
    end select

    if (params%is_sublist('params')) then
      plist => params%sublist('params')
      call this%init(A, plist, stat, errmsg)
    else
      stat = 1
      errmsg = context // 'missing "params" sublist parameter'
    end if

  end subroutine alloc_pcsr_precon

end module pcsr_precon_factory
