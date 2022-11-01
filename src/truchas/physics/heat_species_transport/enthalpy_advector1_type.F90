#include "f90_assert.fpp"

module enthalpy_advector1_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use enthalpy_advector_class
  use unstr_mesh_type
  use avg_phase_prop_type
  use parameter_list_type
  implicit none
  private

  type, extends(enthalpy_advector), public :: enthalpy_advector1
    type(unstr_mesh), pointer :: mesh => null()
    real(r8), pointer :: flux_vol(:,:) => null()
    type(avg_phase_prop) :: enthalpy
    integer,  allocatable :: inflow_face(:)
    real(r8), allocatable :: inflow_temp(:)
  contains
    procedure :: init
    procedure :: get_advected_enthalpy1
    procedure :: get_advected_enthalpy2
  end type

contains

  subroutine init(this, mesh, stat, errmsg)
    use vtrack_driver, only: vtrack_flux_vol_view, vtrack_liq_matid_view
    use flow_driver, only: inflow_plist
    use material_model_driver, only: matl_model
    class(enthalpy_advector1), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    integer, pointer :: matid(:)
    type(parameter_list), pointer :: plist
    this%mesh => mesh
    this%flux_vol => vtrack_flux_vol_view()
    matid => vtrack_liq_matid_view()
    call this%enthalpy%init('enthalpy', matid, matl_model, stat, errmsg)
    if (stat /= 0) return
    plist => inflow_plist()
    call set_inflow_bc(this, plist, stat, errmsg)
    if (stat /= 0) return
    ASSERT(size(this%flux_vol,dim=1) >= size(matid))
    ASSERT(size(this%flux_vol,dim=2) == size(mesh%cface))
  end subroutine init


  !! Input/output arrays are on-process cells only
  subroutine get_advected_enthalpy1(this, tcell, dq)

    class(enthalpy_advector1), intent(in) :: this
    real(r8), intent(in) :: tcell(:)
    real(r8), intent(out) :: dq(:)

    integer :: i, i1, i2, j, n
    real(r8) :: dq_tot, dq_i, state(1)
    real(r8), allocatable :: tcellx(:)
    logical :: found

    ASSERT(size(dq) == this%mesh%ncell_onP)
    ASSERT(size(tcell) == this%mesh%ncell_onP)

    allocate(tcellx(this%mesh%ncell))
    tcellx(:this%mesh%ncell_onP) = tcell
    call this%mesh%cell_imap%gather_offp(tcellx)

    do j = 1, this%mesh%ncell_onP
      i1 = this%mesh%xcface(j)
      i2 = this%mesh%xcface(j+1) - 1
      dq_tot = 0.0_r8  ! net outflux of heat
      do i = i1, i2 ! sides of cell j
        !! Temperature of the fluxed volume
        if (any(this%flux_vol(:,i) > 0.0)) then ! outflow through face
          state(1) = tcellx(j)
        else if (any(this%flux_vol(:,i) < 0.0)) then ! inflow through face
          n = this%mesh%cnhbr(i)
          if (n > 0) then
            state(1) = tcellx(n)
          else
            call get_inflow_temp(this, this%mesh%cface(i), found, state(1))
            if (.not.found) state(1) = tcellx(j)
          end if
        end if
        call this%enthalpy%compute_value(this%flux_vol(:,i), state, dq_i)
        dq_tot = dq_tot + dq_i
      end do
      dq(j) = -dq_tot
    end do

  end subroutine get_advected_enthalpy1


  !! Input/output arrays are on-process cells only
  subroutine get_advected_enthalpy2(this, tcell, dq, tmin, tmax)

    class(enthalpy_advector1), intent(in) :: this
    real(r8), intent(in) :: tcell(:)
    real(r8), intent(out) :: dq(:), tmin(:), tmax(:)

    integer :: i, i1, i2, j, n
    real(r8) :: dq_tot, dq_i, state(1)
    real(r8), allocatable :: tcellx(:)
    logical :: found

    ASSERT(size(dq) == this%mesh%ncell_onP)
    ASSERT(size(tcell) == this%mesh%ncell_onP)
    ASSERT(size(tmin) >= this%mesh%ncell_onP)
    ASSERT(size(tmax) >= this%mesh%ncell_onP)

    allocate(tcellx(this%mesh%ncell))
    tcellx(:this%mesh%ncell_onP) = tcell
    call this%mesh%cell_imap%gather_offp(tcellx)

    do j = 1, this%mesh%ncell_onP
      tmin(j) = tcellx(j)
      tmax(j) = tcellx(j)
      i1 = this%mesh%xcface(j)
      i2 = this%mesh%xcface(j+1) - 1
      dq_tot = 0.0_r8  ! net outflux of heat
      do i = i1, i2 ! sides of cell j
        !! Temperature of the fluxed volume
        if (any(this%flux_vol(:,i) > 0.0)) then ! outflow through face
          state(1) = tcellx(j)
        else if (any(this%flux_vol(:,i) < 0.0)) then  ! inflow through face
          n = this%mesh%cnhbr(i)
          if (n > 0) then
            state(1) = tcellx(n)
          else
            call get_inflow_temp(this, this%mesh%cface(i), found, state(1))
            if (.not.found) state(1) = tcellx(j)
          end if
          tmin(j) = min(tmin(j),state(1))
          tmax(j) = max(tmax(j),state(1))
        end if
        call this%enthalpy%compute_value(this%flux_vol(:,i), state, dq_i)
        dq_tot = dq_tot + dq_i
      end do
      dq(j) = -dq_tot
    end do

  end subroutine get_advected_enthalpy2

  !! This auxiliary subroutine sets the inflow temperature boundary conditions
  !! specified by the given parameter list. The structure of the list is a list
  !! of sublists. Only those that define an "inflow-temperature" parameter are
  !! significant. Its value and the corresponding value of the required
  !! "face-set-ids" parameter are used to set the inflow material temperature
  !! on a portion of the boundary. Other sublists are ignored. The subroutine
  !! completely redefines the %INFLOW_FACE and %INFLOW_TEMP component arrays.

  subroutine set_inflow_bc(this, params, stat, errmsg)

    use bndry_face_group_builder_type

    class(enthalpy_advector1), intent(inout) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    type(bndry_face_group_builder) :: builder
    integer, allocatable :: setids(:), tlist(:), xgroup(:)
    integer :: j, n
    real(r8) :: temp

    call builder%init(this%mesh)

    piter = parameter_list_iterator(params, sublists_only=.true.)
    n = piter%count()
    allocate(tlist(n))
    n = 0
    do while (.not.piter%at_end())
      plist => piter%sublist()
      if (plist%is_parameter('inflow-temperature')) then
        call plist%get('inflow-temperature', temp, stat=stat, errmsg=errmsg)
        if (stat /= 0) return
        call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
        if (stat /= 0) return
        call builder%add_face_group(setids, stat, errmsg)
        if (stat /= 0) return
        n = n + 1
        tlist(n) = temp
      end if
      call piter%next
    end do

    call builder%get_face_groups(n, xgroup, this%inflow_face)
    if (allocated(this%inflow_temp)) deallocate(this%inflow_temp)
    allocate(this%inflow_temp(size(this%inflow_face)))

    do j = 1, n
      this%inflow_temp(xgroup(j):xgroup(j+1)-1) = tlist(j)
    end do

    !NB: If %INFLOW_FACE and %INFLOW_TEMP are sorted by the former array
    !    we can use a binary search in the following subroutine.

  end subroutine set_inflow_bc

  !! This auxiliary subroutine returns the user-specified inflow temperature
  !! for the given boundary face where such data exists (found = true).
  !! NB: I bailed on implementing a binary search in favor of a dumb linear
  !! search. This is only called for boundary faces where inflow is actually
  !! present, and in real applications there are only a small handful of
  !! faces where inflow may occur, so we are only searching a small list a
  !! small number of times.

  subroutine get_inflow_temp(this, face, found, temp)
    class(enthalpy_advector1), intent(in) :: this
    integer, intent(in) :: face
    logical, intent(out) :: found
    real(r8), intent(out) :: temp
    integer k
    do k = 1, size(this%inflow_face)
     if (this%inflow_face(k) == face) exit
    end do
    found = (k <= size(this%inflow_face))
    if (found) temp = this%inflow_temp(k)
  end subroutine get_inflow_temp

end module enthalpy_advector1_type
