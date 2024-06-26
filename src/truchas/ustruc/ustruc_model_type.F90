!!
!! USTRUC_MODEL_TYPE
!!
!! This module defines a derived type that encapsulates the microstructure
!! modeling kernel.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!!  1. This MFD_DISC object is not shared and holds its own large MINV array
!!  component.  That array is only needed on active cells, which for typical
!!  casting applications will be a small fraction of the cells.   This
!!  inefficiency needs to be addressed: pass in a shared MFD_DISC object,
!!  or add specialized capabilities to MFD_DISC itself (FIXME).
!!
!!  2. It is probably more appropriate to instantiate the microstructure
!!  analysis components outside of this object, and pass it in as an argument
!!  to the INIT procedure.  The problem is that the number of analysis points
!!  is determined here, and that is also bound up in the analysis components,
!!  though unnecessarily so.  A redesign may be appropriate here (FIXME).
!!
!!  3. The GET procedures are limited to exactly what the driver uses, and
!!  do not expose the full capabilities of the underlying USTRUC_ANALYSIS%GET
!!  interface.  There are issues here on how best to deal with values on
!!  inactive cells and values on active cells that have invalid data.  The
!!  current implementation is adequate for the present Truchas output needs,
!!  but could be greatly improved.
!!

#include "f90_assert.fpp"

module ustruc_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use mfd_disc_type
  use ustruc_analysis_class
  use parameter_list_type
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  type, public :: ustruc_model
    private
    class(ustruc_analysis), pointer :: analysis => null() ! owned
    type(unstr_mesh), pointer :: mesh => null()  ! reference only -- do not own
    type(mfd_disc) :: disc  ! see Note 1
    logical, allocatable :: mask(:)  ! identifies active mesh cells
    integer :: ncell  ! number of active cells
    integer, allocatable :: map(:)
    !! Model parameters
    real(r8) :: mfrac_min
  contains
    procedure :: init
    procedure :: set_state
    procedure :: update_state
    procedure :: get_comp_list
    procedure :: has
    generic :: get => getr1, getr2
    procedure, private :: getr1
    procedure, private :: getr2
    procedure :: get_map
    procedure :: serialize
    procedure :: deserialize
    final :: delete_ustruc_model
  end type ustruc_model

contains

  subroutine delete_ustruc_model(this)
    type(ustruc_model) :: this
    if (associated(this%analysis)) deallocate(this%analysis)
  end subroutine

  subroutine init(this, mesh, params)

    use ustruc_analysis_factory

    class(ustruc_model), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list) :: params

    integer :: j, stat
    integer, allocatable :: setids(:)
    integer(kind(mesh%cell_set_mask)) :: bitmask
    character(:), allocatable :: errmsg

    this%mesh => mesh

    call this%disc%init(this%mesh, use_new_mfd=.true.) ! See Note 1

    if (params%is_parameter('cell-set-ids')) then
      call params%get('cell-set-ids', setids)
      call this%mesh%get_cell_set_bitmask(setids, bitmask, stat, errmsg)
      if (stat /= 0) call TLS_fatal('USTRUC%INIT: ' // errmsg)
      this%mask = (iand(bitmask, this%mesh%cell_set_mask(:this%mesh%ncell_onP)) /= 0)
    else  ! all cells
      allocate(this%mask(this%mesh%ncell_onP), source=.true.)
    end if

    this%map = pack([(j, j=1,this%mesh%ncell_onP)], this%mask)

    this%ncell = size(this%map)

    call params%get('material-fraction-threshold', this%mfrac_min, default=1.0d-2)
    INSIST(this%mfrac_min <= 1.0_r8)
    INSIST(this%mfrac_min >= 0.0_r8)

    if (params%is_parameter('begin-frac') .and. .not.params%is_parameter('low-temp-phase')) &
        call TLS_fatal('USTRUC%INIT: low-temp-phase must be defined when begin-frac is defined')

    !! For now we instantiate the analysis components here, but this needs to
    !! be done outside and the result passed in.
    this%analysis => new_ustruc_analysis(this%ncell, params)

  end subroutine init

  subroutine set_state(this, t, tcell, tface, liq_vf, sol_vf)
    class(ustruc_model), intent(inout) :: this
    real(r8), intent(in) :: t, tcell(:), tface(:), liq_vf(:), sol_vf(:)
    logical,  allocatable :: invalid(:)
    real(r8), allocatable :: sol_frac(:), temp(:), temp_grad(:,:)
    call get_temp_state(this, tcell, tface, temp, temp_grad)
    call get_sol_frac_state(this, liq_vf, sol_vf, sol_frac, invalid)
    call this%analysis%set_state(t, temp, temp_grad, sol_frac, invalid)
  end subroutine

  subroutine update_state(this, t, tcell, tface, liq_vf, sol_vf)
    class(ustruc_model), intent(inout) :: this
    real(r8), intent(in) :: t, tcell(:), tface(:), liq_vf(:), sol_vf(:)
    logical,  allocatable :: invalid(:)
    real(r8), allocatable :: temp(:), temp_grad(:,:), sol_frac(:)
    call get_temp_state(this, tcell, tface, temp, temp_grad)
    call get_sol_frac_state(this, liq_vf, sol_vf, sol_frac, invalid)
    call this%analysis%update_state(t, temp, temp_grad, sol_frac, invalid)
  end subroutine

  subroutine get_comp_list(this, list)
    class(ustruc_model), intent(in) :: this
    integer, allocatable, intent(out) :: list(:)
    call this%analysis%get_comp_list(list)
  end subroutine

  !! Returns true if one of the analysis components has the named data.
  logical function has(this, name)
    class(ustruc_model), intent(in) :: this
    character(*), intent(in) :: name
    has = this%analysis%has(name)
  end function

  !! I do not know how to handle the logical data with respect to assigning
  !! data on inactive cells -- this really requires knowing what the logical
  !! data means.

  subroutine getr1(this, name, array)
    class(ustruc_model), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(inout) :: array(:)
    real(r8) :: tmp(this%ncell)
    ASSERT(size(array) >= size(this%mask))
    call this%analysis%get(name, tmp)
    call unpack_cc_scalar(this, 0.0_r8, tmp, array)
  end subroutine

  subroutine getr2(this, name, array)
    class(ustruc_model), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(inout) :: array(:,:)
    real(r8) :: tmp(size(array,1),this%ncell)
    ASSERT(size(array,2) >= size(this%mask))
    call this%analysis%get(name, tmp)
    call unpack_cc_vector(this, 0.0_r8, tmp, array)
  end subroutine

  !! These auxillary routines unpack cell-centered data on active cells into
  !! the full cell-centered array over the mesh, assigning specified values
  !! to inactive cells and invalid active cells.  The size of the destination
  !! array must be at least the number of on-process cells.  If it is larger,
  !! say equal to the number of cells, the values on the additional elements
  !! are returned unchanged.

  subroutine unpack_cc_scalar(this, inactive_value, src, dest)
    class(ustruc_model), intent(in) :: this
    real(r8), intent(in) :: inactive_value  ! value to assign to inactive cells
    real(r8), intent(in) :: src(:)          ! values on active cells
    real(r8), intent(inout) :: dest(:)
    integer :: j
    ASSERT(size(src) == this%ncell)
    ASSERT(size(dest) >= size(this%mask))
    do j = 1, size(this%mask)
      if (.not.this%mask(j)) dest(j) = inactive_value
    end do
    dest(this%map) = src
  end subroutine

  subroutine unpack_cc_vector(this, inactive_value, src, dest)
    class(ustruc_model), intent(in) :: this
    real(r8), intent(in) :: inactive_value  ! value to assign to inactive cells
    real(r8), intent(in) :: src(:,:)        ! values on active cells
    real(r8), intent(inout) :: dest(:,:)
    integer :: j
    ASSERT(size(src,1) == size(dest,1))
    ASSERT(size(src,2) == this%ncell)
    ASSERT(size(dest,2) >= size(this%mask))
    do j = 1, size(this%mask)
      if (.not.this%mask(j)) dest(:,j) = inactive_value
    end do
    dest(:,this%map) = src
  end subroutine

  !! This auxiliary subroutine generates the temperature state arrays expected
  !! by the microstructure analysis kernel.  Input temperatures are provided
  !! at all on-process cells and faces, but the output temperature state only
  !! on active cells.
  !!
  !! NB: The gradient is computed using Gauss-Green from the face temperatures.
  !! This is a purely local computation and needs only be done and stored for
  !! the relevant cells, however an existing routine that does it mesh-wide
  !! was used (under control of a mask).  The temporary array GRAD could be
  !! entirely eliminated if this were to be done more smartly (FIXME).

  subroutine get_temp_state(this, tcell, tface, temp, temp_grad)

    class(ustruc_model), intent(in) :: this
    real(r8), intent(in)  :: tcell(:), tface(:)
    real(r8), allocatable, intent(out) :: temp(:), temp_grad(:,:)

    real(r8) :: xtface(this%mesh%nface), grad(3,this%mesh%ncell_onP)

    ASSERT(size(tcell) == this%mesh%ncell_onP)
    ASSERT(size(tface) == this%mesh%nface_onP)

    !! Extract the relevant cell temperatures.
    temp = tcell(this%map)

    !! Compute gradients mesh-wide (under control of mask); see NB above.
    xtface(:this%mesh%nface_onP) = tface  ! off-process extended face temperatures
    call this%mesh%face_imap%gather_offp(xtface)
    call this%disc%compute_cell_grad(xtface, this%mask, grad)

    !! Now extract the relevant gradients.
    temp_grad = grad(:,this%map)

  end subroutine

  !! This auxiliary subroutine generates the solid fraction state arrays
  !! expected by the microstructure analysis kernel.  Input volume fraction
  !! data is provided at all on-process cells, but the output solid fraction
  !! state only on active cells.

  subroutine get_sol_frac_state(this, liq_vf, sol_vf, sol_frac, invalid)

    class(ustruc_model), intent(inout), target :: this
    real(r8), intent(in)  :: liq_vf(:), sol_vf(:)
    real(r8), allocatable, intent(out) :: sol_frac(:)
    logical,  allocatable, intent(out) :: invalid(:)

    integer  :: j
    logical  :: mask(this%mesh%ncell)
    real(r8) :: sfrac(this%mesh%ncell)

    ASSERT(size(liq_vf) == this%mesh%ncell_onP)
    ASSERT(size(sol_vf) == this%mesh%ncell_onP)

    !! Calculate the solid fraction and associated invalid mask arrays.
    do j = 1, this%mesh%ncell_onP
      sfrac(j) = 0.0_r8
      mask(j) = .false.
      if (this%mask(j)) then
        mask(j) = (liq_vf(j) + sol_vf(j) >= this%mfrac_min)
        if (mask(j)) sfrac(j) = sol_vf(j) / (liq_vf(j) + sol_vf(j))
      end if
    end do
    call this%mesh%cell_imap%gather_offp(sfrac)
    call this%mesh%cell_imap%gather_offp(mask)

    !! Extract the relevant solid fraction and invalid mask data.
    sol_frac = sfrac(this%map)
    invalid  = .not.mask(this%map)

  end subroutine get_sol_frac_state

  subroutine get_map(this, map)
    class(ustruc_model), intent(in) :: this
    integer, allocatable, intent(out) :: map(:)
    map = this%map
  end subroutine

  subroutine serialize(this, cid, array)
    use,intrinsic :: iso_fortran_env, only: int8
    class(ustruc_model), intent(in) :: this
    integer, intent(in) :: cid
    integer(int8), allocatable, intent(out) :: array(:,:)
    call this%analysis%serialize(cid, array)
  end subroutine

  subroutine deserialize(this, cid, array)
    use,intrinsic :: iso_fortran_env, only: int8
    class(ustruc_model), intent(inout) :: this
    integer, intent(in) :: cid
    integer(int8), intent(in) :: array(:,:)
    call this%analysis%deserialize(cid, array)
  end subroutine

end module ustruc_model_type
