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
    class(body_box), allocatable, public :: body(:)
    type(unstr_mesh), pointer :: mesh => null() ! do not own
  contains
    procedure :: init
    procedure :: body_at_point
  end type body_identifier

contains

  subroutine init(this, plist, mesh)

    use parameter_list_type
    use body_factories
    use background_body_type

    class(body_identifier), intent(out) :: this
    type(parameter_list), intent(in) :: plist
    type(unstr_mesh), target, intent(in) :: mesh

    integer :: i, j, matl_id, matl_index, stat, nback
    class(body), allocatable :: tmp
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: bparams
    character(:), allocatable :: context, errmsg

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

    ! If there is a background body which isn't lowest priority, move it to the
    ! lowest priority. In principle, we can move the background body to the end
    ! to get the right priority ordering for the body_volume_initialize function,
    ! but the order of volumes in the output of that function corresponds to the
    ! namelist ordering of bodies, and matl_init relies on that ordering. So for
    ! now, leave the body ordering as is and expect the user to put background
    ! at the end. Return a warning message if that is not the case.
    nback = 0
    do i = 1, this%nbody - 1
      associate(f => this%body(i)%f)
        select type (f)
        type is (background_body)
          nback = nback + 1
          ! call move_alloc(this%body(i)%f, tmp)
          ! do j = i+1, this%nbody
          !   call move_alloc(this%body(j)%f, this%body(j-1)%f)
          ! end do
          ! call move_alloc(tmp, this%body(this%nbody)%f)
        end select
      end associate
    end do
    ! if (nback > 1) call TLS_fatal("Too many background bodies.")
    if (nback > 0) call TLS_warn("Background body is not lowest priority. This may produce unexpected results.")
    
  end subroutine init


  ! Return the ID of the body which contains the given point x.
  integer function body_at_point(this, x, cellid)

    class(body_identifier), intent(in) :: this
    real(r8), intent(in) :: x(3)
    integer, intent(in) :: cellid

    do body_at_point = 1, this%nbody
      if (this%body(body_at_point)%eval(x, cellid)) return
    end do

    ! TODO: this should return an error status to be handled at a higher level.
    INSIST(.false.)

  end function body_at_point

end module body_identifier_type
