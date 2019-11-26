!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module body_identifier_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use body_class
  use unstr_mesh_type
  use truchas_logging_services
  implicit none
  private

  type, public :: body_identifier
    private
    integer, public :: nbody
    class(body_box), allocatable :: body(:)
    type(unstr_mesh), pointer :: mesh => null() ! do not own

    logical :: legacy = .true. ! This should go away
  contains
    procedure :: init
    procedure :: body_at_point
  end type body_identifier

  integer, parameter :: POINT = 1
  integer, parameter :: CELL = 2

contains

  subroutine init(this, plist, mesh)

    use parameter_list_type
    use body_factories
    use interfaces_module, only: nbody ! legacy, need to remove

    class(body_identifier), intent(out) :: this
    type(parameter_list), intent(in) :: plist
    type(unstr_mesh), target, intent(in) :: mesh

    integer :: i, matl_id, matl_index, stat
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: bparams
    character(:), allocatable :: context, errmsg

    ! if legacy, don't initialize the body box
    if (this%legacy) then
      this%nbody = nbody
      return
    end if
    
    this%mesh => mesh
    this%nbody = plist%count()
    allocate(this%body(this%nbody))
    piter = parameter_list_iterator(plist)
    do i = 1, this%nbody
      context = 'processing ' // piter%name() // ': '
      bparams => piter%sublist()
      call alloc_body(mesh, bparams, this%body(i)%f)
      call piter%next()
    end do
    
  end subroutine init


  ! Return the ID of the body which contains the given point x.
  integer function body_at_point(this, x, cellid)

    use vof_init, only: body_id_from_vertex ! legacy, need to remove

    class(body_identifier), intent(in) :: this
    real(r8), intent(in) :: x(3)
    integer, intent(in) :: cellid

    if (this%legacy) then
      body_at_point = body_id_from_vertex(x)
      return
    end if

    do body_at_point = 1, this%nbody
      if (this%body(body_at_point)%eval(x, cellid)) return
    end do

    ! TODO: this should return an error status to be handled at a higher level.
    ASSERT(.false.)

  end function body_at_point

end module body_identifier_type
