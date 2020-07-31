!!
!! SOURCE_MESH_FUNCTION
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Last revised 3 Apr 2009
!!
!! This module implements the derived type SOURCE_MF and associated evaluation
!! method SMF_EVAL that are required by the heat transfer/species diffusion
!! (HTSD) solver.  An instance of this type is intended to describe a time and
!! space-dependent external source term that has been discretized over the
!! cells of a distributed mesh.  Given the time value, the evaluation method
!! returns a cell-based array with the cell-averaged source values at that time.
!! This is the only method expected by the HTSD solver, which also treats
!! instances as opaque objects.
!!
!! This particular implementation supports sources that are the superposition
!! of an intrinsic heat source, arising from the explicit coupling of the HTSD
!! solver to other physics solvers, with an extrinsic user-specifed functional
!! source.  Examples of the former type of source are advected heat and Joule
!! heat; this type can be specfied as a cell-based array of source values.
!! The latter functional source can be specified in a piecewise manner: SCAFUN-
!! type functions associated to disjoint sets of cells described via mesh cell
!! set IDs.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! The derived type SOURCE_MF is an opaque structure with private components.
!! No defined assignment is provided; do not use instances of this type in an
!! assignment statement unless you really know what the default assignment is
!! doing (it's almost certainly not what you want).  The extrinsic functional
!! source term is defined by an initial call to SMF_PREP, zero or more calls
!! to SMF_ADD_FUNCTION and an optional call to SMF_SET_DEFAULT, and a final
!! call to SMF_DONE.  This definition is mandatory.  At this stage the object
!! may be evaluated with SMF_EVAL.  Also at any time after this stage the
!! instrinsic source may be set using SMF_SET_EXTRA_SOURCE as often as needed;
!! the specified value will be incorporated into the result returned by
!! subsequent calls to SMF_EVAL.  This definition is optional.
!!
!!  CALL SMF_PREP (THIS, MESH) prepares the SOURCE_MF object THIS to begin
!!    receiving the data associated with one or more source functions.  MESH
!!    is a pointer to the distributed mesh object over which the source is to
!!    be discretized.  The implementation assumes that the mesh cell data will
!!    never be modified.
!!
!!  CALL SMF_ADD_FUNCTION (THIS, F, SETID, STAT, ERRMSG) associates the SCAFUN
!!    function object F to the cell sets with IDs given in the rank-1 integer
!!    array SETID.  The possible cell sets are those defined in MESH.  This
!!    function object should expect (/ t, x, y, z /) as the array argument.
!!    An independent copy of F is stored so that the actual argument may be
!!    destroyed afterwards if desired.  A nonzero value is assigned to the
!!    integer STAT if an error condition occurs, and an explanatory message
!!    is assigned to the character string ERRMSG. Possible errors are referring
!!    to an invalid cell set ID or attempting to associate a cell to more than
!!    one function.  This routine may be called multiple times (or even not at
!!    all), as long as the specified cell sets do not overlap from one call to
!!    the next.
!!
!!  CALL SMF_SET_DEFAULT (THIS, DEFAULT) sets the default source value to
!!    DEFAULT.  This value will be used for any cell that is not associated
!!    to a function.  Calling this routine is optional.  If no default value
!!    is set, then every cell must be associated to a function.
!!
!!  CALL SMF_DONE (THIS, STAT, ERRMSG) performs the final configuration of the
!!    SOURCE_MF object THIS after all the desired calls to SMF_ADD_FUNCTION
!!    and SMF_ADD_DEFAULT have been made.  A nonzero value is assigned to the
!!    integer STAT if an error condition occurs, and an explanatory error
!!    message is assigned to the character string ERRMSG.  It is an error for
!!    a cell to be without a function if no default value was set.  After this
!!    call the object THIS may be evaluated using SMF_EVAL.
!!
!!  CALL SMF_SET_EXTRA_SOURCE (THIS, Q) sets the optional intrinsic source term
!!    to the value of the cell-based array Q.  The value will be incorportated
!!    into the result returned by subsequent calls to SMF_EVAL. This subroutine
!!    may be called as needed.
!!
!!  CALL SMF_EVAL (THIS, T, Q) evaluates the SOURCE_MF object THIS at time T
!!    and returns the computed source values in the cell-based array Q.  The
!!    returned value is the sum of the source specified by SMF_SET_EXTRA_SOURCE,
!!    if any, and the cell-averaged source functions.  A simple 1-point (cell-
!!    centroid) quadrature is currently used.
!!
!!  CALL SMF_DESTROY (THIS) frees the resources associated with the SOURCE_MF
!!    object THIS, returning it to its default initialization state.
!!
!! IMPLEMENTATION NOTE
!!
!! The evaluation procedure is written to short-circuit some potentially costly
!! calculations as directed by some optimization "hinting".  However, this
!! hinting is not currently exposed in the API for the user to control.  But
!! since we can directly interrogate SCAFUN objects to see if they are of the
!! constant flavor, we do set the optimization hinting for those functions at
!! least.
!!

#include "f90_assert.fpp"

module source_mesh_function

  use kinds
  use scalar_func_containers
  use unstr_base_mesh_class
  implicit none
  private

  public :: source_mf, smf_eval   ! all that the HTSD solver needs
  public :: smf_prep, smf_add_function, smf_set_default, smf_done
  public :: smf_set_extra_source, smf_destroy

  type :: source_mf
    private
    class(unstr_base_mesh), pointer :: mesh => null()
    real(r8) :: tlast = -huge(1.0_r8)
    real(r8) :: default = 0.0_r8
    real(r8), pointer :: fvalue(:) => null()
    real(r8), pointer :: xvalue(:) => null()
    integer :: ngroup = -1
    integer, pointer :: xgroup(:) => null()
    integer, pointer :: index(:) => null()
    type(scalar_func_box), allocatable :: farray(:)
    integer, pointer :: hint(:) => null()
    !! temporary components used during initialization
    logical :: no_default = .true.
    integer, pointer :: tag(:) => null()
    type(scalar_func_list) :: flist
  end type source_mf

  integer, parameter, public :: SMF_HINT_NONE    = 0
  integer, parameter, public :: SMF_HINT_CONST   = 1
  integer, parameter, public :: SMF_HINT_T_INDEP = 2
  integer, parameter, public :: SMF_HINT_X_INDEP = 3

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SMF_PREP
 !!

  subroutine smf_prep (this, mesh)

    type(source_mf), intent(out) :: this
    class(unstr_base_mesh), intent(in), target :: mesh

    this%mesh => mesh
    this%ngroup = 0
    allocate(this%tag(mesh%ncell))
    this%tag = 0

  end subroutine smf_prep

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SMF_ADD_FUNCTION
 !!

  subroutine smf_add_function (this, f, setID, stat, errmsg)

    type(source_mf), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: f
    integer, intent(in) :: setID(:)
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    !! Verify that THIS is in the correct state.
    INSIST(this%ngroup>=0 .and. associated(this%tag))

    !! Tag the cells associated with this function.
    call set_tag_array (this, setID, stat, errmsg)
    if (stat /= 0) return

    !! Append this function to the function list.
    call this%flist%append (f)

  end subroutine smf_add_function

  !!
  !! This auxillary subroutine tags the cells identified by the list of
  !! cell set IDs, checking for error conditions in the process.
  !!

  subroutine set_tag_array (this, setID, stat, errmsg)

    type(source_mf), intent(inout) :: this
    integer, intent(in) :: setID(:)
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: i, j
    logical :: mask(this%mesh%ncell)
    integer(kind(this%mesh%cell_set_mask)) :: bitmask

    !! Create the bitmask corresponding to setID.
    bitmask = 0
    do i = 1, size(setID)
      do j = size(this%mesh%cell_set_ID), 1, -1
        if (setID(i) == this%mesh%cell_set_ID(j)) exit
      end do
      if (j == 0) then
        stat = 1
        write(errmsg,'(a,i0)') 'unknown cell set ID: ', setID(i)
        return
      end if
      bitmask = ibset(bitmask, j)
    end do

    !! Identify the cells specified by setID.
    mask = (iand(bitmask, this%mesh%cell_set_mask) /= 0)

    !! Check that these cells don't overlap those from preceding calls.
    if (any(mask .and. this%tag /= 0)) then
      stat = 1
      errmsg = 'cell already associated with a function'
      return
    end if

    !! Set the tag array.
    this%ngroup = 1 + this%ngroup
    this%tag = merge(this%ngroup, this%tag, mask)

    stat = 0
    errmsg = ''

  end subroutine set_tag_array

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SMF_SET_DEFAULT
 !!

  subroutine smf_set_default (this, default)

    type(source_mf), intent(inout) :: this
    real(r8), intent(in) :: default

    !! Verify that THIS is in the correct state.
    INSIST(this%ngroup>=0 .and. associated(this%tag))

    this%no_default = .false.
    this%default = default

  end subroutine smf_set_default

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SMF_DONE
 !!

  subroutine smf_done (this, stat, errmsg)
  
    use const_scalar_func_type

    type(source_mf), intent(inout) :: this
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: n, j

    !! Verify that THIS is in the correct state.
    INSIST(this%ngroup>=0 .and. associated(this%tag))

    ASSERT(minval(this%tag)>= 0 .and. maxval(this%tag)<=this%ngroup)

    !! Verify that every cell has been associated with a function.
    n = count(this%tag > 0)
    if (n /= this%mesh%ncell .and. this%no_default) then
      stat = 1
      errmsg = 'no default and some cells not associated with a function'
      return
    end if

    allocate(this%index(n), this%xgroup(this%ngroup+1))

    !! Prepare XGROUP: cells of group N will be INDEX(XGROUP(N):XGROUP(N+1)-1).
    this%xgroup(1) = 1
    do n = 1, this%ngroup
      this%xgroup(n+1) = this%xgroup(n) + count(this%tag == n)
    end do

    !! Fill the INDEX array; XGROUP(N) stores the next free location for group N.
    do j = 1, size(this%tag)
      n = this%tag(j)
      if (n == 0) cycle
      this%index(this%xgroup(n)) = j
      this%xgroup(n) = 1 + this%xgroup(n)
    end do
    deallocate(this%tag)

    !! Restore XGROUP: XGROUP(N) is now the start of group N+1 instead of N.
    do n = this%ngroup, 1, -1
      this%xgroup(n+1) = this%xgroup(n)
    end do
    this%xgroup(1) = 1

    !! Convert the function list into the final function array.
    call scalar_func_list_to_box_array (this%flist, this%farray)

    !! For now we don't expose optimization hinting to the user,
    !! but we can determine directly which functions are constant.
    allocate(this%hint(this%ngroup))
    do n = 1, this%ngroup
      select type (f => this%farray(n)%f)
      type is (const_scalar_func)
        this%hint(n) = SMF_HINT_CONST
      class default
        this%hint(n) = SMF_HINT_NONE
      end select
    end do

    stat = 0
    errmsg = ''

  end subroutine smf_done

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SMF_SET_EXTRA_SOURCE
 !!

  subroutine smf_set_extra_source (this, q)

    type(source_mf), intent(inout) :: this
    real(r8), intent(in) :: q(:)

    !! Verify that THIS is in the correct state.
    ASSERT(associated(this%index) .and. .not.associated(this%tag))

    ASSERT(size(q) == this%mesh%ncell)

    if (.not.associated(this%xvalue)) allocate(this%xvalue(this%mesh%ncell))
    this%xvalue = q

  end subroutine smf_set_extra_source

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SMF_EVAL
 !!

  subroutine smf_eval (this, t, q)
  
    type(source_mf), intent(inout) :: this
    real(r8), intent(in)  :: t
    real(r8), intent(out) :: q(:)

    integer :: j, n
    logical :: unevaluated
    integer, pointer :: cells(:)
    real(r8) :: args(0:size(this%mesh%x,dim=1))

    !! Verify that THIS is in the correct state.
    ASSERT(associated(this%index) .and. .not.associated(this%tag))

    ASSERT(size(q) == this%mesh%ncell)

    unevaluated = .not.associated(this%fvalue)
    if (unevaluated) then
      allocate(this%fvalue(this%mesh%ncell))
      this%fvalue = this%default
    end if

    if (unevaluated .or. t /= this%tlast) then
      args(0) = t
      do n = 1, this%ngroup
        cells => this%index(this%xgroup(n):this%xgroup(n+1)-1)
        select case (this%hint(n))
        case (SMF_HINT_CONST)
          if (unevaluated) this%fvalue(cells) = this%farray(n)%f%eval(args)
        case (SMF_HINT_X_INDEP)
          this%fvalue(cells) = this%farray(n)%f%eval(args)
        case (SMF_HINT_T_INDEP)
          if (unevaluated) then
            do j = 1, size(cells)
              associate (cnode => this%mesh%cell_node_list_view(cells(j)))
                args(1:) = sum(this%mesh%x(:,cnode),dim=2) / size(cnode)
                this%fvalue(cells(j)) = this%farray(n)%f%eval(args)
              end associate
            end do
          end if
        case default
          do j = 1, size(cells)
            associate (cnode => this%mesh%cell_node_list_view(cells(j)))
              args(1:) = sum(this%mesh%x(:,cnode),dim=2) / size(cnode)
              this%fvalue(cells(j)) = this%farray(n)%f%eval(args)
            end associate
          end do
        end select
      end do
      this%tlast = t
    end if

    !! Return the sum of the function source and extra source.
    if (associated(this%xvalue)) then
      q = this%fvalue + this%xvalue
    else
      q = this%fvalue
    end if

  end subroutine smf_eval

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SMF_DESTROY
 !!

  subroutine smf_destroy (this)

    type(source_mf), intent(inout) :: this

    type(source_mf) :: default

    if (associated(this%fvalue)) deallocate(this%fvalue)
    if (associated(this%xvalue)) deallocate(this%xvalue)
    if (associated(this%xgroup)) deallocate(this%xgroup)
    if (associated(this%index)) deallocate(this%index)
    if (associated(this%hint)) deallocate(this%hint)
    if (associated(this%tag)) deallocate(this%tag)

    this = default  ! assign default initialization values

  end subroutine smf_destroy

end module source_mesh_function
