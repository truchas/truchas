#include "f90_assert.fpp"

module conc_advector_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use bndry_face_func_type
  use parameter_list_type
  implicit none
  private

  type, public :: conc_advector
    private
    type(unstr_mesh), pointer :: mesh => null()  ! unowned reference
    real(r8), pointer :: flux_vol(:,:) => null() ! unowned reference
    type(bndry_face_func) :: c_inflow
    integer :: nfluid
  contains
    procedure :: init
    procedure :: get_advected_scalar
  end type

contains

  subroutine init(this, mesh, n, stat, errmsg)

    use vtrack_driver, only: vtrack_flux_vol_view, vtrack_liq_matid_view
    use flow_driver, only: inflow_plist
    use string_utilities, only: i_to_c

    class(conc_advector), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(in) :: n  ! concentration component
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    this%mesh => mesh
    this%flux_vol => vtrack_flux_vol_view()
    this%nfluid = size(vtrack_liq_matid_view())
    call set_inflow_bc(this, inflow_plist(), 'inflow-conc'//i_to_c(n), stat, errmsg)

  end subroutine

  subroutine get_advected_scalar(this, t, c, dc)

    class(conc_advector), intent(inout) :: this
    real(r8), intent(in)  :: t, c(:)
    real(r8), intent(out) :: dc(:)

    integer :: j, i, i1, i2, n
    real(r8) :: flux_vol, donor_c
    logical :: found

    ASSERT(size(c) == this%mesh%ncell)
    ASSERT(size(dc) >= this%mesh%ncell_onP)

    call this%c_inflow%compute(t) ! update time-dependent inflow concentrations

    do j = 1, this%mesh%ncell_onP
      i1 = this%mesh%xcface(j)
      i2 = this%mesh%xcface(j+1) - 1
      dc(j) = 0.0_r8
      do i = i1, i2 ! sides of cell j
        flux_vol = sum(this%flux_vol(:this%nfluid,i)) ! total outflow through face
        if (flux_vol == 0.0_r8) cycle
        if (flux_vol > 0.0_r8) then
          donor_c = c(j)
        else
          n = this%mesh%cnhbr(i)
          if (n > 0) then
            donor_c = c(n)
          else
            call get_inflow_conc(this, this%mesh%cface(i), found, donor_c)
            if (.not.found) donor_c = c(j)
          end if
        end if
        dc(j) = dc(j) - flux_vol * donor_c
      end do
    end do

  end subroutine get_advected_scalar

  !! This auxiliary subroutine returns the user-specified inflow concentration
  !! for the given boundary face where such data exists (found = true).
  !! NB: I bailed on implementing a binary search in favor of a dumb linear
  !! search. This is only called for boundary faces where inflow is actually
  !! present, and in real applications there are only a small handful of
  !! faces where inflow may occur, so we are only searching a small list a
  !! small number of times.

  subroutine get_inflow_conc(this, face, found, conc)
    class(conc_advector), intent(in) :: this
    integer, intent(in) :: face
    logical, intent(out) :: found
    real(r8), intent(out) :: conc
    integer ::i
    do i = 1, size(this%c_inflow%index)
      if (this%c_inflow%index(i) == face) exit
    end do
    found = (i <= size(this%c_inflow%index))
    if (found) conc = this%c_inflow%value(i)
  end subroutine

  subroutine set_inflow_bc(this, params, varname, stat, errmsg)

    use scalar_func_factories

    class(conc_advector), intent(inout) :: this
    type(parameter_list), intent(inout) :: params
    character(*), intent(in) :: varname
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    integer, allocatable :: setids(:)
    class(scalar_func), allocatable :: f

    call this%c_inflow%init(this%mesh)

    stat = 0
    piter = parameter_list_iterator(params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      if (plist%is_parameter(varname)) then
        call plist%get('face-set-ids', setids, stat, errmsg)
        if (stat /= 0) exit
        call alloc_scalar_func(plist, varname, f, stat, errmsg)
        if (stat /= 0) exit
        call this%c_inflow%add(f, setids, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next
    end do

    if (stat == 0) then
      call this%c_inflow%add_complete
    else
      errmsg = 'FLOW_BC[' // piter%name() // ']: ' // errmsg
      return
    end if

  end subroutine

end module conc_advector_type
