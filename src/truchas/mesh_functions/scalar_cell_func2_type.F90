!!
!! SCALAR_CELL_FUNC2_TYPE
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module scalar_cell_func2_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use scalar_mesh_func2_class
  use scalar_func_containers
  use cell_group_builder_type
  implicit none
  private

  type, extends(scalar_mesh_func2), public :: scalar_cell_func2
    private
    class(unstr_mesh), pointer :: mesh => null()  ! reference only -- not owned
    logical :: evaluated = .false.
    integer :: ngroup
    integer, allocatable :: xgroup(:), index(:)
    type(scalar_func_box), allocatable :: farray(:)
    integer, allocatable :: hint(:)
    ! construction phase temporaries
    type(cell_group_builder), allocatable :: builder
    type(scalar_func_list) :: flist
  contains
    procedure :: init
    procedure :: add
    procedure :: assemble
    procedure :: compute
  end type scalar_cell_func2

  !! Optimization hint values.
  integer, parameter :: HINT_NONE    = 0
  integer, parameter :: HINT_CONST   = 1
  integer, parameter :: HINT_X_INDEP = 3

contains

  subroutine init(this, mesh)
    class(scalar_cell_func2), intent(out) :: this
    class(unstr_mesh), intent(in), target :: mesh
    this%mesh => mesh
    call this%mesh%init_cell_centroid
    allocate(this%builder)
    call this%builder%init(mesh)
  end subroutine

  subroutine add(this, f, setids, stat, errmsg)
    class(scalar_cell_func2), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: f
    integer, intent(in)  :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_cell_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%flist%append(f)
  end subroutine add

  subroutine assemble(this)
    use const_scalar_func_type
    class(scalar_cell_func2), intent(inout) :: this
    integer :: n
    ASSERT(allocated(this%builder))
    call this%builder%get_cell_groups(this%ngroup, this%xgroup, this%index)
    deallocate(this%builder)
    call scalar_func_list_to_box_array(this%flist, this%farray)
    allocate(this%hint(this%ngroup))
    do n = 1, this%ngroup
      select type (f => this%farray(n)%f)
      type is (const_scalar_func)
        this%hint(n) = HINT_CONST
      class default
        this%hint(n) = HINT_NONE
      end select
    end do
  end subroutine assemble


  subroutine compute(this, t, v)

    class(scalar_cell_func2), intent(inout) :: this
    real(r8), intent(in) :: t, v(:)

    integer :: j, n
    real(r8) :: args(0:size(this%mesh%x,dim=1))

    args(0) = t
    do n = 1, this%ngroup
      associate (index => this%index(this%xgroup(n):this%xgroup(n+1)-1))
        select case (this%hint(n))
        case (HINT_CONST)
          if (.not.this%evaluated) this%value(index) = this%farray(n)%eval(args)
        case (HINT_X_INDEP)
          this%value(index) = this%farray(n)%eval(args)
        case default
          do j = 1, size(index)
            args(1:3) = this%mesh%cell_centroid(:,index(j))
            args(4) = v(index(j))
            this%value(index(j)) = this%farray(n)%eval(args)
          end do
        end select
      end associate
    end do
    this%evaluated = .true.

  end subroutine compute

end module scalar_cell_func2_type
