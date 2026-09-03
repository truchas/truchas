!!
!! FLOW_MATERIAL_MAPPING_TYPE
!!
!! Defines the flow-facing reduction of simulation material distribution. Real
!! fluid phases occupy distinct slots, followed by optional VOID and a lumped
!! SOLID residual slot. The flow phase IDs and their parent material IDs are
!! retained separately because material distribution is indexed by material
!! while flow properties are indexed by phase. Geometric reconstruction uses
!! the separate full PRIORITY permutation.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_material_mapping_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use material_model_type
  use material_distribution_type
  use unstr_2d_mesh_type
  implicit none
  private

  type, public :: flow_material_mapping
    private
    integer, allocatable :: fluid_pid(:), fluid_mid(:)
    character(:), allocatable :: fluid_name(:)
    integer :: void_pid = 0, void_mid = 0
    logical :: has_solid = .false.
    integer, allocatable :: priority(:)
  contains
    procedure :: init
    procedure :: set_priority
    procedure :: num_real_fluid
    procedure :: num_fluid
    procedure :: num_material
    procedure :: get_real_fluid_phase_ids
    procedure :: get_real_fluid_material_ids
    procedure :: get_priority
    procedure :: get_reduced_volume_fractions
    procedure :: put_reduced_volume_fractions
    procedure :: get_phase_volume_fractions
    procedure :: apply_phase_fluxes
  end type

contains

  subroutine init(this, matl_model, stat, errmsg)

    class(flow_material_mapping), intent(out) :: this
    type(material_model), intent(in) :: matl_model
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, j, mid, nchar, first, last, p, nfluid

    allocate(this%fluid_pid(matl_model%nphase_real), this%fluid_mid(matl_model%nphase_real))
    this%fluid_pid = 0
    this%fluid_mid = 0
    j = 0
    do mid = 1, matl_model%nmatl_real
      call matl_model%get_matl_phase_index_range(mid, first, last)
      nfluid = count(matl_model%is_fluid(first:last))
      this%has_solid = this%has_solid .or. any(.not.matl_model%is_fluid(first:last))
      if (nfluid > 1) then
        stat = 1
        errmsg = 'flow currently supports at most one fluid phase per material'
        return
      end if
      if (nfluid == 1) then
        p = first - 1 + findloc(matl_model%is_fluid(first:last), .true., dim=1)
        j = j + 1
        this%fluid_pid(j) = p
        this%fluid_mid(j) = mid
      end if
    end do
    this%fluid_pid = this%fluid_pid(:j)
    this%fluid_mid = this%fluid_mid(:j)
    nchar = 1
    do i = 1, size(this%fluid_pid)
      nchar = max(nchar, len(matl_model%phase_name(this%fluid_pid(i))))
    end do
    allocate(character(nchar) :: this%fluid_name(size(this%fluid_pid)))
    do i = 1, size(this%fluid_pid)
      this%fluid_name(i) = matl_model%phase_name(this%fluid_pid(i))
    end do
    if (matl_model%have_void) then
      this%void_pid = matl_model%void_index
      this%void_mid = matl_model%nmatl
    end if

    allocate(this%priority(this%num_material()))
    this%priority = [(i, i=1,size(this%priority))]

    stat = 0
    errmsg = ''

  end subroutine


  subroutine set_priority(this, params, stat, errmsg)

    class(flow_material_mapping), intent(inout) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: names(:)
    integer :: i

    ASSERT(allocated(this%priority))
    if (.not.params%is_vector('material-priority')) then
      stat = 0
      errmsg = ''
      return
    end if

    call params%get('material-priority', names, stat, errmsg)
    if (stat /= 0) return
    if (size(names) /= size(this%priority)) then
      stat = 1
      errmsg = 'material-priority must name every fluid, VOID, and SOLID slot exactly once'
      return
    end if
    do i = 1, size(names)
      this%priority(i) = slot_index(this, names(i))
      if (this%priority(i) == 0) then
        stat = 1
        errmsg = 'invalid material-priority entry: ' // names(i)
        return
      end if
    end do
    if (.not.all([(count(this%priority == i), i=1,size(this%priority))] == 1)) then
      stat = 1
      errmsg = 'material-priority entries must be unique'
      return
    end if
    stat = 0
    errmsg = ''
  end subroutine


  integer function slot_index(this, name) result(slot)

    class(flow_material_mapping), intent(in) :: this
    character(*), intent(in) :: name
    integer :: k

    slot = 0
    do k = 1, size(this%fluid_pid)
      if (name == this%fluid_name(k)) then
        slot = k
        return
      end if
    end do
    if (this%void_pid /= 0 .and. name == 'VOID') then
      slot = this%num_fluid()
      return
    end if
    if (this%has_solid .and. name == 'SOLID') slot = this%num_material()

  end function


  integer function num_real_fluid(this)
    class(flow_material_mapping), intent(in) :: this

    num_real_fluid = size(this%fluid_pid)
  end function


  integer function num_fluid(this)
    class(flow_material_mapping), intent(in) :: this

    num_fluid = this%num_real_fluid() + merge(1, 0, this%void_pid /= 0)
  end function


  integer function num_material(this)
    class(flow_material_mapping), intent(in) :: this

    num_material = this%num_fluid() + merge(1, 0, this%has_solid)
  end function


  subroutine get_real_fluid_phase_ids(this, phase_ids)
    class(flow_material_mapping), intent(in) :: this
    integer, intent(out) :: phase_ids(:)

    ASSERT(size(phase_ids) == size(this%fluid_pid))
    phase_ids = this%fluid_pid
  end subroutine


  subroutine get_real_fluid_material_ids(this, matl_ids)
    class(flow_material_mapping), intent(in) :: this
    integer, intent(out) :: matl_ids(:)

    ASSERT(size(matl_ids) == size(this%fluid_mid))
    matl_ids = this%fluid_mid
  end subroutine


  subroutine get_priority(this, priority)
    class(flow_material_mapping), intent(in) :: this
    integer, intent(out) :: priority(:)

    ASSERT(size(priority) == size(this%priority))
    priority = this%priority
  end subroutine


  !! Form the reduced flow distribution for the material-level single-phase
  !! contract. SOLID is the residual after the fluid and VOID slots. The
  !! phase-aware variant is GET_PHASE_VOLUME_FRACTIONS.
  subroutine get_reduced_volume_fractions(this, matl_dist, vfrac)
    class(flow_material_mapping), intent(in) :: this
    type(material_distribution), intent(in) :: matl_dist
    real(r8), intent(out) :: vfrac(:,:)

    integer :: m

    ASSERT(size(vfrac,1) == this%num_material())
    ASSERT(size(vfrac,2) >= size(matl_dist%vfrac,2))
    vfrac = 0.0_r8
    do m = 1, this%num_real_fluid()
      vfrac(m,:size(matl_dist%vfrac,2)) = matl_dist%vfrac(this%fluid_mid(m),:)
    end do
    if (this%void_mid /= 0) then
      vfrac(this%num_fluid(),:size(matl_dist%vfrac,2)) = matl_dist%vfrac(this%void_mid,:)
    end if
    if (this%has_solid) then
      vfrac(this%num_material(),:size(matl_dist%vfrac,2)) = &
          1.0_r8 - sum(vfrac(:this%num_fluid(),:size(matl_dist%vfrac,2)), dim=1)
    end if
  end subroutine


  !! Copy transported fluid and VOID slots back to the material distribution.
  !! The lumped SOLID slot is not unpacked; that is unambiguous only for the
  !! current single-phase material contract.
  subroutine put_reduced_volume_fractions(this, vfrac, matl_dist)
    class(flow_material_mapping), intent(in) :: this
    real(r8), intent(in) :: vfrac(:,:)
    type(material_distribution), intent(inout) :: matl_dist

    integer :: m

    ASSERT(size(vfrac,1) == this%num_material())
    ASSERT(size(vfrac,2) >= size(matl_dist%vfrac,2))
    ASSERT(all(abs(sum(vfrac(:,:size(matl_dist%vfrac,2)), dim=1) - 1.0_r8) <= 16.0_r8 * epsilon(1.0_r8)))
    do m = 1, this%num_real_fluid()
      matl_dist%vfrac(this%fluid_mid(m),:) = vfrac(m,:size(matl_dist%vfrac,2))
    end do
    if (this%void_mid /= 0) then
      matl_dist%vfrac(this%void_mid,:) = vfrac(this%num_fluid(),:size(matl_dist%vfrac,2))
    end if
    ASSERT(all(abs(sum(matl_dist%vfrac, dim=1) - 1.0_r8) <= 16.0_r8 * epsilon(1.0_r8)))
  end subroutine


  !! Form the reduced flow distribution from a material distribution and the
  !! current temperature.  A multiphase material contributes the fraction of
  !! its fluid phase to that phase's flow slot; all of its non-fluid phases
  !! contribute to the lumped SOLID slot.  The material distribution itself
  !! remains material-level state and is not modified.
  subroutine get_phase_volume_fractions(this, matl_model, matl_dist, temperature, vfrac)
    class(flow_material_mapping), intent(in) :: this
    type(material_model), intent(in) :: matl_model
    type(material_distribution), intent(in) :: matl_dist
    real(r8), intent(in) :: temperature(:)
    real(r8), intent(out) :: vfrac(:,:)

    integer :: c, m, mid, first, last, nphase, ncell
    real(r8), allocatable :: beta(:)

    ncell = size(matl_dist%vfrac, 2)
    ASSERT(size(matl_dist%vfrac,1) == matl_model%nmatl)
    ASSERT(size(temperature) >= ncell)
    ASSERT(size(vfrac,1) == this%num_material())
    ASSERT(size(vfrac,2) >= ncell)
    allocate(beta(matl_model%nphase_real))
    vfrac = 0.0_r8

    do m = 1, this%num_real_fluid()
      mid = this%fluid_mid(m)
      call matl_model%get_matl_phase_index_range(mid, first, last)
      nphase = last - first + 1
      if (nphase == 1) then
        vfrac(m,:ncell) = matl_dist%vfrac(mid,:)
      else
        do c = 1, ncell
          beta(:nphase) = 0.0_r8
          call matl_model%get_matl_phase_frac(mid, temperature(c), beta(:nphase))
          vfrac(m,c) = matl_dist%vfrac(mid,c) * beta(this%fluid_pid(m)-first+1)
        end do
      end if
    end do
    if (this%void_mid /= 0) vfrac(this%num_fluid(),:ncell) = matl_dist%vfrac(this%void_mid,:)
    if (this%has_solid) vfrac(this%num_material(),:ncell) = &
        1.0_r8 - sum(vfrac(:this%num_fluid(),:ncell), dim=1)
  end subroutine


  !! Apply the divergence of material-transport phase fluxes to the
  !! material-level distribution.  The fluxes, rather than a tracker trial
  !! VOF, are the authoritative result of material transport.  Stationary
  !! solid phases are absent from FLUX_VOLUMES and remain unchanged.
  subroutine apply_phase_fluxes(this, mesh, flux_volumes, matl_dist)
    class(flow_material_mapping), intent(in) :: this
    type(unstr_2d_mesh), intent(in) :: mesh
    real(r8), intent(in) :: flux_volumes(:,:)
    type(material_distribution), intent(inout) :: matl_dist

    integer :: c, m, first_face, last_face

    ASSERT(size(matl_dist%vfrac,2) == mesh%ncell_onP)
    ASSERT(size(flux_volumes,1) == this%num_fluid())
    ASSERT(size(flux_volumes,2) == size(mesh%cface))

    !! A one-material, single-phase fluid has invariant composition.  Avoid
    !! feeding even a small discrete divergence back into its state.
    if (size(matl_dist%vfrac,1) == 1 .and. this%num_real_fluid() == 1 .and. &
        this%void_mid == 0 .and. .not.this%has_solid) return

    do c = 1, mesh%ncell_onP
      first_face = mesh%cstart(c)
      last_face = mesh%cstart(c+1) - 1
      do m = 1, this%num_real_fluid()
        matl_dist%vfrac(this%fluid_mid(m),c) = matl_dist%vfrac(this%fluid_mid(m),c) - &
            sum(flux_volumes(m,first_face:last_face))/mesh%volume(c)
      end do
      if (this%void_mid /= 0) then
        matl_dist%vfrac(this%void_mid,c) = matl_dist%vfrac(this%void_mid,c) - &
            sum(flux_volumes(this%num_fluid(),first_face:last_face))/mesh%volume(c)
      end if
    end do

    !! TODO: Reconcile the moving phase fractions after transport so they
    !! sum to one with the stationary-solid fraction.  This must preserve
    !! small fragments in MATL_DIST rather than adopting a clipped tracker
    !! trial VOF.
  end subroutine

end module flow_material_mapping_type
