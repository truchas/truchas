!!
!! USTRUC_MODEL_TYPE
!!
!! This module defines a derived type that encapsulates the microstructure
!! modeling kernel.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2014
!!
!! PROGRAMMING INTERFACE
!!
!!  The derived type USTRUC_MODEL encapsulates the microstructure modeling
!!  kernel.  It has the following type bound procedures.
!!
!!  INIT (MESH, PARAMS) initializes the object.  MESH is a pointer to
!!    the TYPE(UNSTR_MESH) computational mesh.  An internal reference to
!!    the mesh is held by the object.  All input/output state fields are
!!    relative to this mesh.
!!
!!  SET_STATE (THIS, T, TCELL, TFACE, LIQ_VF, SOL_VF) sets the initial
!!    state at time T.  TCELL and TFACE are the cell and face temperatures
!!    associated with the mesh, and LIQ_VF and SOL_VF are the liquid and
!!    solid volume fractions of the subject material on the cells.
!!
!!  UPDATE_STATE (THIS, T, TCELL, TFACE, LIQ_VF, SOL_VF) updates the state
!!    to the specified values.  The interface is identical to SET_STATE.
!!    The difference here is in behavior: the microstructure model is
!!    advanced from the previous state to the new specified state.
!!
!!  GET (NAME, ARRAY) returns the value of the named data in the provided
!!    array.  This is a generic procedure.  ARRAY is a rank-1 or 2 array
!!    of type logical, integer, or real.  NAME is a character string that
!!    identifies the requested data.  The known names are those encoded in
!!    the various concrete implementations of the USTRUC_COMP class which
!!    ultimately handle the request.  The TKR of ARRAY and its shape must
!!    be consistent with the requested data.
!!
!!  HAS (NAME) returns true if NAME identifies known data that the GET
!!    subroutine is able to retrieve.
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
!!  do not expose the full capabilities of the underlying USTRUC_COMP%GET
!!  interface.  There are issues here on how best to deal with values on
!!  inactive cells and values on active cells that have invalid data.  The
!!  current implementation is adequate for the present Truchas output needs,
!!  but could be greatly improved.
!!

#include "f90_assert.fpp"

module ustruc_model_type

  use kinds, only: r8
  use unstr_mesh_type
  use mfd_disc_type
  use ustruc_comp_class
  use parameter_list_type
  use truchas_logging_services
  implicit none
  private

  type, public :: ustruc_model
    private
    class(ustruc_comp), pointer :: comp => null() ! analysis component
    type(unstr_mesh), pointer :: mesh => null()  ! reference only -- do not own
    type(mfd_disc) :: disc  ! see Note 1
    logical, allocatable :: mask(:)  ! identifies active mesh cells
    integer :: ncell  ! number of active cells
    integer, allocatable :: map(:)
    !! Model parameters
    real(r8) :: mfrac_min
    integer, allocatable :: sym_setids(:)
    type(parameter_list) :: grad_params ! gradient solver parameters
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

  subroutine delete_ustruc_model (this)
    type(ustruc_model) :: this
    if (associated(this%comp)) deallocate(this%comp)
  end subroutine delete_ustruc_model

  subroutine init (this, mesh, params)

    use ustruc_comp_factory

    class(ustruc_model), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list) :: params

    integer :: j, stat
    integer, allocatable :: setids(:)
    integer(kind(mesh%cell_set_mask)) :: bitmask
    character(:), allocatable :: errmsg
    real(r8) :: rpar

    this%mesh => mesh

    call this%disc%init (this%mesh, use_new_mfd=.true.) ! See Note 1

    call params%get ('cell-set-ids', setids)
    call this%mesh%get_cell_set_bitmask (setids, bitmask, stat, errmsg)
    if (stat /= 0) call TLS_fatal ('USTRUC%INIT: ' // errmsg)

    allocate(this%mask(mesh%ncell_onP))
    this%mask = (iand(bitmask, this%mesh%cell_set_mask(:this%mesh%ncell_onP)) /= 0)

    this%map = pack([(j, j=1,this%mesh%ncell_onP)], this%mask)

    this%ncell = size(this%map)

    call params%get ('material-fraction-threshold', this%mfrac_min, default=1.0d-2)
    INSIST(this%mfrac_min <= 1.0_r8)
    INSIST(this%mfrac_min >= 0.0_r8)

    call params%get ('symmetry-face-sets', this%sym_setids, default=[integer::])
    
    if (params%is_parameter('grad-abs-tol')) then
      call params%get ('grad-abs-tol', rpar)
      call this%grad_params%set ('grad-abs-tol', rpar)
    end if

    if (params%is_parameter('grad-rel-tol')) then
      call params%get ('grad-rel-tol', rpar)
      call this%grad_params%set ('grad-rel-tol', rpar)
    end if

    !! For now we instantiate the analysis components here, but this needs to
    !! be done outside and the result passed in.
    this%comp => new_ustruc_comp(this%ncell, params)

  end subroutine init

  subroutine set_state (this, t, tcell, tface, liq_vf, sol_vf)
    class(ustruc_model), intent(inout) :: this
    real(r8), intent(in) :: t, tcell(:), tface(:), liq_vf(:), sol_vf(:)
    logical,  allocatable :: invalid(:)
    real(r8), allocatable :: sol_frac(:), sol_frac_grad(:,:), temp(:), temp_grad(:,:)
    call get_temp_state (this, tcell, tface, temp, temp_grad)
    call get_sol_frac_state (this, liq_vf, sol_vf, sol_frac, sol_frac_grad, invalid)
    call this%comp%set_state (t, temp, temp_grad, sol_frac, sol_frac_grad, invalid)
  end subroutine set_state

  subroutine update_state (this, t, tcell, tface, liq_vf, sol_vf)
    use timing_tree
use process_info_module
    class(ustruc_model), intent(inout) :: this
    real(r8), intent(in) :: t, tcell(:), tface(:), liq_vf(:), sol_vf(:)
    logical,  allocatable :: invalid(:)
    real(r8), allocatable :: temp(:), temp_grad(:,:), sol_frac(:), sol_frac_grad(:,:)
    call start_timer ('temperature gradient')
    call get_temp_state (this, tcell, tface, temp, temp_grad)
    call stop_timer ('temperature gradient')
    call start_timer ('solid fraction gradient')
call mem_diag_write ('***A calling get_sol_frac_state')
    call get_sol_frac_state (this, liq_vf, sol_vf, sol_frac, sol_frac_grad, invalid)
    call stop_timer ('solid fraction gradient')
    call start_timer ('analysis')
call mem_diag_write ('***B calling comp%update_state')
    call this%comp%update_state (t, temp, temp_grad, sol_frac, sol_frac_grad, invalid)
    call stop_timer ('analysis')
call mem_diag_write ('***C returning from ustruc_model%update_state')
  end subroutine update_state

  subroutine get_comp_list (this, list)
    class(ustruc_model), intent(in) :: this
    integer, allocatable, intent(out) :: list(:)
    call this%comp%get_comp_list (list)
  end subroutine get_comp_list
  
  !! Returns true if one of the analysis components has the named data.
  logical function has (this, name)
    class(ustruc_model), intent(in) :: this
    character(*), intent(in) :: name
    has = this%comp%has(name)
  end function has

  !! I do not know how to handle the logical data with respect to assigning
  !! data on inactive cells -- this really requires knowing what the logical
  !! data means.

  subroutine getr1 (this, name, array)
    class(ustruc_model), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(inout) :: array(:)
    real(r8) :: tmp(this%ncell)
    ASSERT(size(array) >= size(this%mask))
    call this%comp%get (name, tmp)
    call unpack_cc_scalar (this, 0.0_r8, tmp, array)
  end subroutine getr1

  subroutine getr2 (this, name, array)
    class(ustruc_model), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(inout) :: array(:,:)
    real(r8) :: tmp(size(array,1),this%ncell)
    ASSERT(size(array,2) >= size(this%mask))
    call this%comp%get (name, tmp)
    call unpack_cc_vector (this, 0.0_r8, tmp, array)
  end subroutine getr2

  !! These auxillary routines unpack cell-centered data on active cells into
  !! the full cell-centered array over the mesh, assigning specified values
  !! to inactive cells and invalid active cells.  The size of the destination
  !! array must be at least the number of on-process cells.  If it is larger,
  !! say equal to the number of cells, the values on the additional elements
  !! are returned unchanged.

  subroutine unpack_cc_scalar (this, inactive_value, src, dest)
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
  end subroutine unpack_cc_scalar

  subroutine unpack_cc_vector (this, inactive_value, src, dest)
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
  end subroutine unpack_cc_vector

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

  subroutine get_temp_state (this, tcell, tface, temp, temp_grad)

    use index_partitioning, only: gather_boundary

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
    call gather_boundary (this%mesh%face_ip, xtface)
    call this%disc%compute_cell_grad (xtface, this%mask, grad)

    !! Now extract the relevant gradients.
    temp_grad = grad(:,this%map)

  end subroutine get_temp_state

  !! This auxiliary subroutine generates the solid fraction state arrays
  !! expected by the microstructure analysis kernel.  Input volume fraction
  !! data is provided at all on-process cells, but the output solid fraction
  !! state only on active cells.
  !!
  !! NB1: The gradient is computed using a mimetic finite difference scheme
  !! that computes face fluxes.  This requires a global linear mass matrix
  !! type solve.  This is generally quite an easy solve and underlying solver
  !! parameters have been hardwired within the gradient solver object
  !! (CELL_GRAD), however some of the solver parameters ought to be exposed
  !! and made setable via the input (FIXME).
  !!
  !! NB2: The gradient linear system matrix depends on the calculated mask
  !! (cells where valid solid fraction data exists), which may differ from
  !! one call to the next.  Consequently the solver is created afresh locally.
  !! However, it will be most often the case that the mask is changing very
  !! infrequently, if at all, and there could be significant savings in
  !! reusing the setup solver when possible (FIXME).

  subroutine get_sol_frac_state (this, liq_vf, sol_vf, sol_frac, sol_frac_grad, invalid)

    use cell_grad_type
    use index_partitioning, only: gather_boundary

    class(ustruc_model), intent(inout), target :: this
    real(r8), intent(in)  :: liq_vf(:), sol_vf(:)
    real(r8), allocatable, intent(out) :: sol_frac(:), sol_frac_grad(:,:)
    logical,  allocatable, intent(out) :: invalid(:)

    integer  :: j, stat
    logical  :: mask(this%mesh%ncell)
    real(r8) :: sfrac(this%mesh%ncell), sfrac_grad(3,this%mesh%ncell)
    type(cell_grad), target :: grad
    character(:), allocatable :: errmsg

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
    call gather_boundary (this%mesh%cell_ip, sfrac)
    call gather_boundary (this%mesh%cell_ip, mask)

    !! Extract the relevant solid fraction and invalid mask data.
    sol_frac = sfrac(this%map)
    invalid  = .not.mask(this%map)

    !! Setup the gradient solver; see NB1 and NB2 above.
    call grad%init (this%disc, mask, this%sym_setids, this%grad_params, stat, errmsg)
    if (stat /= 0) call TLS_fatal ('USTRUC_MODEL%INIT: gradient solver setup error: ' // errmsg)

    !! Solve for the solid fraction gradient and extract the relevant data.
    call grad%compute (sfrac, sfrac_grad, stat, errmsg)
    if (stat /= 0) call TLS_fatal ('USTRUC_MODEL%INIT: gradient solver error: ' // errmsg)
    sol_frac_grad = sfrac_grad(:,this%map)

  end subroutine get_sol_frac_state

  subroutine get_map (this, map)
    class(ustruc_model), intent(in) :: this
    integer, allocatable, intent(out) :: map(:)
    map = this%map
  end subroutine

  subroutine serialize (this, cid, array)
    use,intrinsic :: iso_fortran_env, only: int8
    class(ustruc_model), intent(in) :: this
    integer, intent(in) :: cid
    integer(int8), allocatable, intent(out) :: array(:,:)
    call this%comp%serialize (cid, array)
  end subroutine serialize

  subroutine deserialize (this, cid, array)
    use,intrinsic :: iso_fortran_env, only: int8
    class(ustruc_model), intent(inout) :: this
    integer, intent(in) :: cid
    integer(int8), intent(in) :: array(:,:)
    call this%comp%deserialize (cid, array)
  end subroutine deserialize

end module ustruc_model_type
