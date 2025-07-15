!!
!! FLOW_PROPS_TYPE
!!
!! Defines type(flow_props), a container for flow related material properties.
!! All data members are public since they are needed throughout flow but should
!! not be changed, except through `update_XX` subroutines
!!
!!
!! Peter Brady <ptb@lanl.gov>
!! 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  The supported usage pattern of a flow_props instance is:
!!
!!  call fp%read_params
!!  call fp%init
!!  call fp%update_cc(...,initial=.true.)
!!  call fp%update_fc(...,initial=.true.)
!!
!!  do while (flow_simulation)
!!    ...
!!    call fp%update_cc(...)
!!    ...
!!    call fp%update_fc(...)
!!    ...
!!    call fp%accept()
!!    ...
!!  end do
!!
!!
!!
!! PROGRAMMING INTERFACE - MODULE PROCEDURES
!!
!!  READ_FLOW_PROPS_NAMELIST (LUN, PARAMETER_LIST) reads the flow_cutoffs namelist.
!!    This is currently a collection of magic numbers describing when flow should
!!    treat certain things as solid faces and cells.  This is currently called by the
!!    flow driver.
!!
!! PROGRAMMING INTERFACE - TYPE BOUND PROCEDURES
!!
!!  READ_PARAMS(this, PARAMETER_LIST)
!!    extracts `cutvof` and `min face fraction` from parameter list and stores
!!    them in `this`.
!!
!!  INIT(this, FLOW_MESH, REAL DENSITY(:), SCALAR_FUNC_BOX DENSITY_DELTA(:),
!!             SCALAR_FUNC_BOX VISCOSITY (:))
!!    Initializes `this` and allocates all data.  Steals the allocations
!!    of DENSITY_DELTA and VISCOSITY.  Both of these are scalar functions
!!    of temperature.  The i'th index of DENSITY, DENSITY_DELTA, and VISCOSITY
!!    must all refer to the _same_ material
!!
!!  UPDATE_CC(this, REAL VOF(:,:), REAL TEMPERATURE_CC(:) [, LOGICAL INITIAL])
!!    Constructs cell-centered material propererties in this%xxx_cc.
!!    Linear averaging using VOF is used to build the properties (i.e:
!!    rho_cc, mu_cc, ...).  VOF is of the form VOF(nmaterials,
!!    ncells).  VOF(i,:) must be the same material as DENSITY(i),
!!    VISCOSITY(i).  If INITIAL is specified, copy all data to this%xxx_cc_n
!!
!!  UPDATE_FC (this [, LOGICAL INITIAL])
!!    Computes face centered density via volume averaging and face
!!    centered viscosity via harmoic averaging.  This is separate from
!!    this%update_cc since the turbulence model is allowed to modify
!!    the cell centered viscosity.
!!
!!  ACCEPT (this)
!!    Called if the timestep taken by the flow driver is acceptable.
!!    Copies all data to their `_n` counterparts
!!

#include "f90_assert.fpp"

module flow_props_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use flow_domain_types
  use unstr_mesh_type
  use parameter_list_type
  use truchas_logging_services
  use truchas_timers
  use parallel_communication
  use scalar_func_containers
  implicit none
  private

  type, public :: flow_props
    type(unstr_mesh), pointer, private :: mesh => null() ! reference only -- do not own
    real(r8), allocatable :: rho_cc(:), rho_cc_n(:) ! cell centered fluid density
    real(r8), allocatable :: rho_fc(:), rho_fc_n(:) ! face centered fluid density
    real(r8), allocatable :: mu_cc(:), mu_cc_n(:) ! cell centered fluid viscosity
    real(r8), allocatable :: mu_fc(:), mu_fc_n(:) ! face centered fluid viscosity
    real(r8), allocatable :: vof(:), vof_n(:) ! fluid volume fraction (includes void)
    real(r8), allocatable :: vof_novoid(:), vof_novoid_n(:) ! non-void fluid volume fraction
    real(r8), allocatable :: rho_delta_cc(:) ! temperature dependent density delta
    real(r8), allocatable :: solidified_rho(:) ! change in fluid density due to solidification
    real(r8), allocatable :: rho_pre_sol(:) ! volume fraction prior to solidification
    real(r8), allocatable :: void_delta_cc(:) ! change in void volume fraction
    integer, allocatable :: cell_t(:), face_t(:) ! cell and face types (fluid, void or solid)

    integer :: nfluid ! number of mobile materials
    logical :: any_void ! true if there is void in the current timestep
    logical :: any_real_fluid ! true if there is any non-void fluid in the current timestep
    logical :: any_real_fluid_onP ! same as above, but local to this PE
    real(r8) :: cutoff ! fluid volume fractions below this are considered solid
    real(r8) :: min_face_fraction
    ! real(r8) :: cutrho not entirely sure what to do with this
    real(r8) :: minrho
    !
    real(r8), allocatable :: density(:)
    class(scalar_func_box), allocatable :: density_delta(:), viscosity(:)
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: update_cc
    procedure :: update_fc
    procedure :: accept
    procedure :: set_pre_solidification_density
  end type flow_props

contains

  subroutine init(this, mesh, nfluid, density, density_delta, viscosity, params)

    use parameter_list_type

    class(flow_props), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(in) :: nfluid
    real(r8), intent(in) :: density(:)
    class(scalar_func_box), intent(inout) :: density_delta(:), viscosity(:)
    type(parameter_list), intent(inout) :: params

    integer :: nc, fc, s
    type(parameter_list), pointer :: plist

    this%mesh => mesh
    this%nfluid = nfluid
    this%density = density

    plist => params%sublist("cutoffs")
    call plist%get("cutvof", this%cutoff, default=0.01_r8)
    call plist%get("min face fraction", this%min_face_fraction, default=0.001_r8)

    ASSERT(size(viscosity) == size(density_delta))
    allocate(this%viscosity(size(viscosity)))
    allocate(this%density_delta(size(density_delta)))

    do s = 1, size(viscosity)
      call move_alloc(viscosity(s)%f, this%viscosity(s)%f)
      call move_alloc(density_delta(s)%f, this%density_delta(s)%f)
    end do

    nc = this%mesh%ncell
    fc = this%mesh%nface

    allocate(this%rho_cc(nc), this%rho_cc_n(nc), &
        this%rho_fc(fc), this%rho_fc_n(fc), &
        this%mu_cc(nc), this%mu_cc_n(nc), &
        this%mu_fc(fc), this%mu_fc_n(fc), &
        this%vof(nc), this%vof_n(nc), &
        this%vof_novoid(nc), this%vof_novoid_n(nc), &
        this%rho_delta_cc(nc), &
        this%void_delta_cc(nc), &
        this%solidified_rho(nc), &
        this%rho_pre_sol(nc), &
        this%cell_t(nc), this%face_t(fc), &
        stat=s)
    if (s /= 0) call TLS_Fatal("allocation of flow_props failed")

    this%rho_pre_sol = -huge(1.0_r8) ! this forces the initial solidified rho to be 0

  end subroutine init

  subroutine set_initial_state(this, vof, state_cc)
    class(flow_props), intent(inout) :: this
    real(r8), intent(in) :: vof(:,:), state_cc(:,:)
    this%vof_n = 0.0_r8
    this%vof_novoid_n = 0.0_r8
    call update_cc(this, vof, state_cc)
    call update_fc(this)
    call accept(this)
    this%void_delta_cc = 0.0_r8
  end subroutine set_initial_state

  subroutine update_cc(this, vof, state_cc)

    class(flow_props), intent(inout) :: this
    real(r8), intent(in) :: vof(:,:), state_cc(:,:)

    integer :: m, i, j
    real(r8) :: minrho

    call start_timer("update properties")

    minrho = huge(1.0_r8)
    this%any_real_fluid = .false.

    do i = 1, this%mesh%ncell
      this%rho_cc(i) = 0.0_r8
      this%mu_cc(i) = 0.0_r8
      this%rho_delta_cc(i) = 0.0_r8
      do m = 1, size(this%density)
        this%rho_cc(i) = this%rho_cc(i) + vof(m,i)*this%density(m)
        this%mu_cc(i) = this%mu_cc(i) + &
            vof(m,i)*this%viscosity(m)%f%eval(state_cc(:,i))
        this%rho_delta_cc(i) = this%rho_delta_cc(i) + &
            vof(m,i)*this%density_delta(m)%f%eval(state_cc(:,i))
      end do

      this%vof(i) = sum(vof(:this%nfluid,i))
      ! last element of input vof array is void
      this%vof_novoid(i) = sum(vof(:size(this%density),i))
      this%solidified_rho(i) = max(this%rho_pre_sol(i) - this%rho_cc(i), 0.0_r8)

      if (this%vof(i) > 0.0_r8) then
        this%rho_cc(i) = this%rho_cc(i) / this%vof(i)
        this%rho_delta_cc(i) = this%rho_delta_cc(i) / this%vof(i)
        this%mu_cc(i) = this%mu_cc(i) / this%vof(i)
      end if

      if (this%vof(i) < this%cutoff) then ! criteria for solid
        this%cell_t(i) = solid_t
      elseif (this%vof_novoid(i) == 0) then ! criteria for void
        this%cell_t(i) = void_t
      elseif (this%vof(i) > this%vof_novoid(i)) then ! this a mixed fluid/void cell...
        this%cell_t(i) = regular_void_t ! EXPERIMENT
      else ! regular
        this%cell_t(i) = regular_t
      end if

      if (this%vof_novoid(i) > 0.0_r8) then
        minrho = min(minrho, this%rho_cc(i)*this%vof(i) / this%vof_novoid(i))
      end if

      block
        real(r8) :: void, void_n
        void = this%vof(i) - this%vof_novoid(i)
        void_n = this%vof_n(i) - this%vof_novoid_n(i)
        this%void_delta_cc(i) = (void - void_n)!*this%mesh%volume(i)
      end block
    end do

    ! only mark a cell as regular_void if it does not have any neighbors where are void_t
    do i = 1, this%mesh%ncell_onP

      if (this%cell_t(i) /= regular_void_t) cycle

      block
        logical :: void

        associate (cn => this%mesh%cnhbr(this%mesh%xcnhbr(i):this%mesh%xcnhbr(i+1)-1))

          void = .false.
          do j = 1, size(cn)
            if (cn(j) > 0) then
              void = (void .or. this%cell_t(cn(j)) == void_t)
            end if
          end do

          if (void) then
            this%cell_t(i) = regular_t
          end if

        end associate
      end block
    end do

    call this%mesh%cell_imap%gather_offp(this%cell_t)

    this%minrho = global_minval(minrho)
    this%any_void = global_any(this%cell_t == void_t) ! needed for dirichlet boundary conditions
    this%any_real_fluid_onP = any(this%vof_novoid(:this%mesh%ncell_onP) > this%cutoff)
    this%any_real_fluid = global_any(this%any_real_fluid_onP)

    call stop_timer("update properties")

  end subroutine update_cc

  subroutine set_pre_solidification_density(this, vof)

    class(flow_props), intent(inout) :: this
    real(r8), intent(in) :: vof(:,:)

    integer :: i, nfluid

    nfluid = size(this%density)
    do i = 1, this%mesh%ncell_onP
      this%rho_pre_sol(i) = sum(this%density*vof(:nfluid,i))
    end do

  end subroutine set_pre_solidification_density


  subroutine update_fc(this)

    class(flow_props), intent(inout) :: this

    integer :: j
    real(r8) :: w(2), min_face_rho

    min_face_rho = this%minrho*this%min_face_fraction

    ! compute face types
    do j = 1, this%mesh%nface_onP
      associate(cn => this%mesh%fcell(:,j))
        if (cn(2) > 0) then
          if (any(this%cell_t(cn) == void_t) .and. any(this%cell_t(cn) <= regular_t)) then
            this%face_t(j) = regular_void_t
          else
            ! enusre that regular_void_t cells do get inappropriately labeled faces
            this%face_t(j) = max(maxval(this%cell_t(cn)), regular_t)
          end if
        else
          ! enusre that regular_void_t cells do get inappropriately labeled faces
          this%face_t(j) = max(this%cell_t(cn(1)), regular_t)
        end if
      end associate
    end do
    call this%mesh%face_imap%gather_offp(this%face_t)

    ! linear averaged face-centered density
    ! 1) If the face has only one cell neighbor (i.e. a boundary cell)
    ! - the face density is simply the cell density
    ! 2) If the face has shared by a solid cell (where rho_cc == 0) and a fluid cell
    ! - the face density is the fluid cell density
    ! 3) If the face is shared by two void/solid cells
    ! - the face density is 0
    ! 4) any non-zero face density is forced to be at least min_face_rho
    this%rho_fc = 0.0_r8
    do j = 1, this%mesh%nface_onP
      associate(cn => this%mesh%fcell(:,j))
        if (cn(2) == 0) then
          this%rho_fc(j) = this%rho_cc(cn(1))
        else if (any(this%cell_t(cn) <= regular_t)) then
          w = this%mesh%volume(cn)*this%vof(cn)
          this%rho_fc(j) = max(min_face_rho, sum(this%rho_cc(cn)*w)/sum(w))
        end if
      end associate
    end do
    call this%mesh%face_imap%gather_offp(this%rho_fc)

    ! harmonic averaged face viscosity
    ! special cases:
    ! 1) If the face has only one cell neighbor (i.e. a boundary cell)
    ! - the face viscosity is the cell viscosity
    ! 2) If the face is shared by a solid cell (where mu_cc == 0) and a fluid cell
    ! - the face viscosity is the fluid cell viscosity
    ! 3) If the face is shared by a void cell (where mu_cc == 0) and a fluid cell
    ! - the face viscosity is the void cell viscosity (i.e. 0)
    ! 4) If the face is shared by two void/solid cells
    ! - the face viscosity is 0
    this%mu_fc = 0.0_r8
    do j = 1, this%mesh%nface_onP
      associate(cn => this%mesh%fcell(:,j))
        if (cn(2) == 0) then
          this%mu_fc(j) = this%mu_cc(cn(1)) ! boundary
        else if (this%face_t(j) == solid_t) then
          this%mu_fc(j) = maxval(this%mu_cc(cn))
        else if (product(this%mu_cc(cn))  > epsilon(1.0_r8)) then
          this%mu_fc(j) = 2.0_r8*product(this%mu_cc(cn))/sum(this%mu_cc(cn))
        end if
      end associate
    end do
    call this%mesh%face_imap%gather_offp(this%mu_fc)

  end subroutine update_fc

  subroutine accept(this)
    class(flow_props), intent(inout) :: this
    this%rho_cc_n(:) = this%rho_cc(:)
    this%rho_fc_n(:) = this%rho_fc(:)
    this%mu_cc_n(:) = this%mu_cc(:)
    this%mu_fc_n(:) = this%mu_fc(:)
    this%vof_n(:) = this%vof(:)
    this%vof_novoid_n(:) = this%vof_novoid(:)
  end subroutine accept

end module flow_props_type
