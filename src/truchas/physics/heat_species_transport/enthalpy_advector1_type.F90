#include "f90_assert.fpp"

module enthalpy_advector1_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use enthalpy_advector_class
  use unstr_mesh_type
  implicit none
  private

  type, extends(enthalpy_advector), public :: enthalpy_advector1
    type(unstr_mesh), pointer :: mesh => null()
    real(r8), pointer :: flux_vol(:,:) => null()
    integer, pointer :: matid(:) => null()
  contains
    procedure :: init
    procedure :: get_advected_enthalpy1
    procedure :: get_advected_enthalpy2
  end type

contains

  subroutine init(this, mesh)
    use vtrack_driver, only: vtrack_flux_vol_view, vtrack_liq_matid_view
    class(enthalpy_advector1), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    this%mesh => mesh
    this%flux_vol => vtrack_flux_vol_view()
    this%matid => vtrack_liq_matid_view()
    ASSERT(size(this%flux_vol,dim=1) == size(this%matid))
    ASSERT(size(this%flux_vol,dim=2) == size(mesh%cface))
  end subroutine init


  !! Input/output arrays are on-process cells only
  subroutine get_advected_enthalpy1(this, tcell, dq)

    use material_interop, only: ds_enthalpy_density
    use index_partitioning, only: gather_boundary

    class(enthalpy_advector1), intent(in) :: this
    real(r8), intent(in) :: tcell(:)
    real(r8), intent(out) :: dq(:)

    integer :: i, i1, i2, j, n, m
    real(r8) :: sum, state(1)
    real(r8), allocatable :: tcellx(:)

    ASSERT(size(dq) == this%mesh%ncell_onP)
    ASSERT(size(tcell) == this%mesh%ncell_onP)

    allocate(tcellx(this%mesh%ncell))
    tcellx(:this%mesh%ncell_onP) = tcell
    call gather_boundary(this%mesh%cell_ip, tcellx)

    do j = 1, this%mesh%ncell_onP
      i1 = this%mesh%xcface(j)
      i2 = this%mesh%xcface(j+1) - 1
      sum = 0.0_r8  ! net outflux of heat
      do i = i1, i2 ! sides of cell j
        !! Temperature of the fluxed volume
        if (any(this%flux_vol(:,i) >= 0.0)) then
          state(1) = tcellx(j)
        else
          n = this%mesh%cnhbr(i)
          if (n > 0) then
            state(1) = tcellx(n)
          else
            state(1) = tcellx(j) !FIXME: Want an inflow temperature here
          end if
        end if
        do m = 1, size(this%flux_vol,dim=1)
          if (this%flux_vol(m,i) /= 0) &
              sum = sum + this%flux_vol(m,i)*ds_enthalpy_density(this%matid(m),state)
        end do
      end do
      dq(j) = -sum
    end do

  ! How to handle inflow temps:
  ! 1) side-based array of neighbor temps, overwritten with inflow temps on
  !    on boundary sides. Natural for boundary condition implemented as (side,value)
  !    pairs. Basically what old truchas does.  Could subsequently overwrite side-based
  !    temp array with cell temp. Downside: large temporary array, don't know that an
  !    inflow temp was provided on side with inflow
  ! 2) skip bndry inflow sides on first pass. On second clean-up pass, go through list
  !    of inflow temp sides (BC) and augment the result on sides with inflow. Could
  !    assemble the list of bndry inflow sides on first pass and then cross them off
  !    during second pass to identify bndry inflow sides without an inflow temp.
  ! 3) single pass without pre-computed neighbor temp array. Means we need arbitrary
  !    access to inflow temp BC data by side. Not simple; requires searching as well
  !    as reporting when data doesn't exist.

  end subroutine get_advected_enthalpy1


  !! Input/output arrays are on-process cells only
  subroutine get_advected_enthalpy2(this, tcell, dq, tmin, tmax)

    use material_interop, only: ds_enthalpy_density
    use index_partitioning, only: gather_boundary

    class(enthalpy_advector1), intent(in) :: this
    real(r8), intent(in) :: tcell(:)
    real(r8), intent(out) :: dq(:), tmin(:), tmax(:)

    integer :: i, i1, i2, j, n, m
    real(r8) :: sum, state(1)
    real(r8), allocatable :: tcellx(:)

    ASSERT(size(dq) == this%mesh%ncell_onP)
    ASSERT(size(tcell) == this%mesh%ncell_onP)
    ASSERT(size(tmin) >= this%mesh%ncell_onP)
    ASSERT(size(tmax) >= this%mesh%ncell_onP)

    allocate(tcellx(this%mesh%ncell))
    tcellx(:this%mesh%ncell_onP) = tcell
    call gather_boundary(this%mesh%cell_ip, tcellx)

    do j = 1, this%mesh%ncell_onP
      tmin(j) = tcellx(j)
      tmax(j) = tcellx(j)
      i1 = this%mesh%xcface(j)
      i2 = this%mesh%xcface(j+1) - 1
      sum = 0.0_r8  ! net outflux of heat
      do i = i1, i2 ! sides of cell j
        !! Temperature of the fluxed volume
        if (any(this%flux_vol(:,i) >= 0.0)) then
          state(1) = tcellx(j)
        else
          n = this%mesh%cnhbr(i)
          if (n > 0) then
            state(1) = tcellx(n)
            tmin(j) = min(tmin(j),tcellx(n))
            tmax(j) = max(tmax(j),tcellx(n))
          else
            state(1) = tcellx(j) !FIXME: Want an inflow temperature here
          end if
        end if
        do m = 1, size(this%flux_vol,dim=1)
          if (this%flux_vol(m,i) /= 0) &
              sum = sum + this%flux_vol(m,i)*ds_enthalpy_density(this%matid(m),state)
        end do
      end do
      dq(j) = -sum
    end do

  ! How to handle inflow temps:
  ! 1) side-based array of neighbor temps, overwritten with inflow temps on
  !    on boundary sides. Natural for boundary condition implemented as (side,value)
  !    pairs. Basically what old truchas does.  Could subsequently overwrite side-based
  !    temp array with cell temp. Downside: large temporary array, don't know that an
  !    inflow temp was provided on side with inflow
  ! 2) skip bndry inflow sides on first pass. On second clean-up pass, go through list
  !    of inflow temp sides (BC) and augment the result on sides with inflow. Could
  !    assemble the list of bndry inflow sides on first pass and then cross them off
  !    during second pass to identify bndry inflow sides without an inflow temp.
  ! 3) single pass without pre-computed neighbor temp array. Means we need arbitrary
  !    access to inflow temp BC data by side. Not simple; requires searching as well
  !    as reporting when data doesn't exist.

  end subroutine get_advected_enthalpy2

end module enthalpy_advector1_type
