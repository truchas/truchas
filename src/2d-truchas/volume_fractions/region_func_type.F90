module region_func_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use region_class
  implicit none
  private

  type, public :: region_func
    private
    type(region_box), allocatable :: reg(:)
  contains
    procedure :: init
    procedure :: num_region
    procedure :: region_index
  end type

contains

  subroutine init(this, mesh, params, stat, errmsg)

    use base_mesh_class
    use parameter_list_type
    use region_factory

    class(region_func), intent(out) :: this
    class(base_mesh), intent(in) :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    
    integer :: n
    type(parameter_list), pointer :: plist
    type(parameter_list_iterator) :: piter

    piter = parameter_list_iterator(params, sublists_only=.true.)
    allocate(this%reg(piter%count()))
    do n = 1, size(this%reg)
      plist => piter%sublist()
      call alloc_region(this%reg(n)%reg, mesh, plist, stat, errmsg)
      if (stat /= 0) return
      !select type (reg => this%reg(n))
      !type is (background_region)
      !  if (n /= size(this%reg)) then
      !    ! Warn that remaining regions will be unused?
      !  end if
      !end select
      call piter%next
    end do

  end subroutine

  pure integer function num_region(this)
    class(region_func), intent(in) :: this
    num_region = size(this%reg)
  end function

  pure integer function region_index(this, x, bitmask) result(i)
    class(region_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: bitmask
    do i = 1, size(this%reg)
      if (this%reg(i)%encloses(x, bitmask)) return
    end do
    i = 0
  end function

end module region_func_type
