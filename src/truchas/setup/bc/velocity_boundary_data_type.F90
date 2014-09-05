!!
!! VELOCITY_BOUNDARY_DATA_TYPE
!!
!! A time-dependent Dirichlet velocity BC data hack.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 16 Jan 2014
!!
!! Until now (version 2.8) Truchas has only allowed constant velocity Dirichlet
!! conditions to be specified.  That data is stored in the array BC_Vel(:,:,:)
!! (3 by 6 by number of cells), together with some additional mask arrays. This
!! module is an implementation of time-dependent velocity boundary conditions
!! that is designed to slip into the existing code with *minimal* structural
!! changes to both the BC initialization code and the flow physics kernel.
!! A much more efficient implementation would be possible with a comprehensive
!! refactoring of the flow physics kernel.
!!
!! USAGE
!!
!! No data is stored here; only the VELOCITY_BOUNDARY_DATA derived type and
!! its type bound procedures are defined.  Other Truchas code will declare a
!! variable of this type and initialize/use it, by executing the following
!! in order:
!!
!! (1) Call the INIT method.  Must be done after the final mesh distribution
!!     has been performed.
!! (2) Call the SET_VEL_FUNC method for each boundary surface where a Dirichlet
!!     velocity condition was specified and pass the corresponding function.
!! (3) Call SET_NO_SLIP and SET_DIRICHLET methods as needed.  The API was
!!     chosen to slip directly into the existing BC_INIT code.
!! (4) Now calls to the GET method can be made.  The API was chosen so that
!!     it matches closely access to the original BC_Vel array, including
!!     return values for array elements where no valid data exists.
!!

#include "f90_assert.fpp"

module velocity_boundary_data_type

  use kinds, only: r8
  use vector_func_class
  implicit none
  private

  type :: vfbox
    class(vector_func), pointer :: vf => null()
  contains
    final :: vfbox_delete
  end type vfbox

  type, public :: velocity_boundary_data
    private
    integer, allocatable :: index(:,:)
    type(vfbox), allocatable :: vel_func_array(:)
  contains
    procedure :: init
    procedure :: set_vel_func
    procedure :: set_no_slip
    procedure :: set_dirichlet
    procedure, private :: get_one
    procedure, private :: get_comp_one
    generic :: get => get_one, get_comp_one
    procedure :: get_all => get_comp_all
  end type velocity_boundary_data

contains

  !! Final subroutine for the private vfbox type.
  subroutine vfbox_delete (this)
    type(vfbox), intent(inout) :: this
    if (associated(this%vf)) deallocate(this%vf)
  end subroutine vfbox_delete

  !! Allocates and initializes the internal data structures;
  !! NSURF is the number of BC namelists.
  subroutine init (this, ncell, nsurf)
    class(velocity_boundary_data), intent(out) :: this
    integer, intent(in) :: ncell  ! number of local cells
    integer, intent(in) :: nsurf  ! number of boundary surfaces (globally)
    ASSERT(ncell >= 0)
    ASSERT(nsurf >= 0)
    allocate(this%index(6,ncell)) ! N.B. hardwiring for hex cells here.
    this%index = -1 ! default is no data -- get returns NULL_R
    allocate(this%vel_func_array(nsurf))
  end subroutine

  !! Sets the velocity function to use for the specified boundary surface.
  !! The object assumes ownership of the function; pointer is returned null.
  !! Not all boundary surfaces will have a function (nor necessarily any).
  subroutine set_vel_func (this, p, vf)
    class(velocity_boundary_data), intent(inout) :: this
    integer, intent(in) :: p  ! the boundary surface index
    class(vector_func), pointer :: vf ! pointer to the velocity function
    ASSERT(p > 0 .and. p <= size(this%vel_func_array))
    ASSERT(vf%dim == 3)
    if (associated(this%vel_func_array(p)%vf)) deallocate(this%vel_func_array(p)%vf)
    this%vel_func_array(p)%vf => vf
    nullify(vf) ! the object assumes ownership
  end subroutine

  !! Sets the faces specified by f and mask to return 0 velocity;
  !! index value 0 represents this special case.
  subroutine set_no_slip (this, f, mask)
    class(velocity_boundary_data), intent(inout) :: this
    integer, intent(in) :: f  ! cell face index
    logical, intent(in) :: mask(:)  ! cell-based mask array
    ASSERT(size(mask) == size(this%index,dim=2))
    where (mask) this%index(f,:) = 0
  end subroutine

  !! Sets the faces specified by f and mask to return the velocity
  !! computed by the velocity function associated with boundary surface p.
  subroutine set_dirichlet (this, p, f, mask)
    class(velocity_boundary_data), intent(inout) :: this
    integer, intent(in) :: p        ! the boundary surface index
    integer, intent(in) :: f        ! cell face index
    logical, intent(in) :: mask(:)  ! cell-based mask array
    ASSERT(p > 0 .and. p <= size(this%vel_func_array))
    ASSERT(associated(this%vel_func_array(p)%vf))
    ASSERT(f > 0 .and. f <= size(this%index,dim=1))
    ASSERT(size(mask) == size(this%index,dim=2))
    where (mask) this%index(f,:) = p
  end subroutine

  !! Returns the velocity 3-vector for the face specified by f and n
  !! at the given time t.
  function get_one (this, f, n, t) result (v)
    use input_utilities, only: NULL_R
    class(velocity_boundary_data), intent(in) :: this
    integer,  intent(in) :: f ! cell face index
    integer,  intent(in) :: n ! cell number
    real(r8), intent(in) :: t ! time
    real(r8) :: v(3)
    integer p
    ASSERT(f > 0 .and. f <= size(this%index,dim=1))
    ASSERT(n > 0 .and. n <= size(this%index,dim=2))
    p = this%index(f,n)
    select case (p)
    case (0)  ! no slip
      v = 0.0_r8
    case (1:) ! general Dirichlet using pth function
      v = this%vel_func_array(p)%vf%eval([t])
    case default
      v = NULL_R
    end select
  end function get_one

  !! Returns the specified velocity component i for the face specified
  !! by f and n and the specified time t.
  function get_comp_one (this, i, f, n, t) result (vi)
    use input_utilities, only: NULL_R
    class(velocity_boundary_data), intent(in) :: this
    integer,  intent(in) :: i ! component index
    integer,  intent(in) :: f ! cell face index
    integer,  intent(in) :: n ! cell number
    real(r8), intent(in) :: t ! time
    real(r8) :: vi
    integer p
    ASSERT(f > 0 .and. f <= size(this%index,dim=1))
    ASSERT(n > 0 .and. n <= size(this%index,dim=2))
    p = this%index(f,n)
    select case (p)
    case (0)  ! no slip
      vi = 0.0_r8
    case (1:) ! general Dirichlet using pth function
      vi = this%vel_func_array(p)%vf%eval_comp(i,[t])
    case default
      vi = NULL_R
    end select
  end function get_comp_one

  !! Returns the specified velocity component for face f on all cells.
  !! I would rather not do this, but there is at least one place that would
  !! require significant structural rearrangement in order to use the other get.
  function get_comp_all (this, i, f, t) result (v)
    class(velocity_boundary_data), intent(in) :: this
    integer,  intent(in) :: i ! component index
    integer,  intent(in) :: f ! cell face index
    real(r8), intent(in) :: t ! time
    real(r8) :: v(size(this%index,dim=2))
    integer :: n
    do n = 1, size(v)
      v(n) = get_comp_one(this, i, f, n, t)
    end do
  end function get_comp_all

end module velocity_boundary_data_type
