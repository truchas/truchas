#include "f90_assert.fpp"

module flow_props_type

  use kinds, only: r8
  use unstr_mesh_type
  use flow_mesh_type
  use parameter_list_type
  use truchas_logging_services
  use truchas_timers
  use index_partitioning
  use parallel_communication
  use phase_property_table
  use scalar_func_containers
  implicit none
  private

  public :: flow_props

  type :: flow_props
    type(flow_mesh), pointer :: mesh => null() ! reference only -- do not own
    real(r8), allocatable :: rho_cc(:), rho_cc_n(:) ! cell centered fluid density
    real(r8), allocatable :: rho_fc(:), rho_fc_n(:) ! face centered fluid density
    real(r8), allocatable :: mu_cc(:), mu_cc_n(:) ! cell centered fluid viscosity
    real(r8), allocatable :: mu_fc(:), mu_fc_n(:) ! face centered fluid viscosity
    real(r8), allocatable :: vof(:), vof_n(:) ! fluid volume fraction (includes void)
    real(r8), allocatable :: vof_novoid(:), vof_novoid_n(:) ! non-void fluid volume fraction
    real(r8), allocatable :: rho_delta_cc(:) ! temperature dependent density deviation
    real(r8), allocatable :: w_face(:,:) ! workspace
    logical :: contains_void
    real(r8) :: cutvof ! fluid volume fractions below this are considered solid
    real(r8) :: min_face_fraction
    real(r8) :: cutrho
    real(r8) :: minrho
    !
    real(r8), allocatable :: density(:)
    class(scalar_func_box), allocatable :: density_delta(:), viscosity(:)
  contains
    procedure :: read_params
    procedure :: init
    procedure :: update
  end type flow_props

contains

  subroutine read_params(this, params)
    class(flow_props), intent(inout) :: this
    type(parameter_list), intent(inout) :: params
    !-
    call params%get("cutvof", this%cutvof, 0.01_r8)
    call params%get("min_face_fraction", this%min_face_fraction, 0.001_r8)
  end subroutine read_params

  subroutine init(this, mesh, density, density_delta, viscosity, contains_void)
    class(flow_props), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: mesh
    real(r8), intent(in) :: density(:)
    class(scalar_func_box), intent(inout) :: density_delta(:), viscosity(:)
    logical, intent(in) :: contains_void
    !-
    integer :: nc, fc, s

    this%mesh => mesh
    this%contains_void = contains_void
    this%density = density

    ASSERT(size(viscosity) == size(density_delta))
    allocate(this%viscosity(size(viscosity)))
    allocate(this%density_delta(size(density_delta)))

    do s = 1, size(viscosity)
      call move_alloc(viscosity(s)%f, this%viscosity(s)%f)
      call move_alloc(density_delta(s)%f, this%density_delta(s)%f)
    end do

    nc = this%mesh%mesh%ncell
    fc = this%mesh%mesh%nface

    allocate(this%rho_cc(nc), this%rho_cc_n(nc), &
        this%rho_fc(fc), this%rho_fc_n(fc), &
        this%mu_cc(nc), this%mu_cc_n(nc), &
        this%mu_fc(fc), this%mu_fc_n(fc), &
        this%vof(nc), this%vof_n(nc), &
        this%vof_novoid(nc), this%vof_novoid_n(nc), &
        this%rho_delta_cc(nc), &
        this%w_face(fc, 2), stat=s)
    if (s /= 0) call TLS_Fatal("allocation of flow_props failed")

  end subroutine init

  subroutine update(this, vof, temperature_cc, initial)
    class(flow_props), intent(inout) :: this
    real(r8), intent(in) :: vof(:,:), temperature_cc(:)
    logical, optional, intent(in) :: initial
    !-
    logical :: ini
    integer :: m, i, j
    real(r8) :: minrho, w, min_face_rho, state(1)
    type(unstr_mesh), pointer :: mesh

    if (present(initial)) then
      ini = initial
    else
      ini = .false.
    end if
    mesh => this%mesh%mesh

    minrho = huge(1.0_r8)
    ! cell-centered quantities
    do i = 1, mesh%ncell
      this%rho_cc(i) = 0.0_r8
      this%mu_cc(i) = 0.0_r8
      this%rho_delta_cc(i) = 0.0_r8
      state(1) = temperature_cc(i)
      do m = 1, size(this%density)
        this%rho_cc(i) = this%rho_cc(i) + vof(m,i)*this%density(m)
        this%mu_cc(i) = this%mu_cc(i) + &
            vof(m,i)*this%viscosity(m)%f%eval(state)
        this%rho_delta_cc = this%rho_delta_cc(i) + &
            vof(m,i)*this%density_delta(m)%f%eval(state)
      end do
      this%vof(i) = sum(vof(:,i))
      if (this%contains_void) then
        this%vof_novoid(i) = sum(vof(:size(this%density),i))
      else
        this%vof_novoid(i) = this%vof(i)
      end if

      if (this%vof(i) > 0.0_r8) then
        this%rho_cc(i) = this%rho_cc(i) / this%vof(i)
        this%rho_delta_cc(i) = this%rho_delta_cc(i) / this%vof(i)
      end if

      if (this%vof_novoid(i) > 0.0_r8) then
        minrho = min(minrho, this%rho_cc(i)*this%vof(i) / this%vof_novoid(i))
      end if
    end do
    this%minrho = global_minval(minrho) ! is this the only place we need a minrho?
    min_face_rho = this%minrho*this%min_face_fraction

    ! linear averaged face-centered density
    ! special cases:
    ! 1) If the face has only one cell neighbor (i.e. a boundary cell)
    ! - the face density is simply the cell density
    ! 2) If the face has shared by a solid cell (where rho_cc == 0) and a fluid cell
    ! - the face density is the fluid cell density
    ! 3) If the face is shared by two solid cells
    ! - the face density is 0
    ! 4) any non-zero face density is forced to be at least min_face_rho
    this%w_face(:,1) = 0.0_r8
    this%rho_fc(:) = 0.0_r8
    do i = 1, mesh%ncell
      associate(fi => mesh%cface(mesh%xcface(i):mesh%xcface(i+1)-1))
        w = mesh%volume(i)*this%vof(i)
        this%w_face(fi,1) = this%w_face(fi,1) + w
        this%rho_fc(fi) = this%rho_fc(fi) + w*this%rho_cc(i)
      end associate
    end do
    do j = 1, mesh%nface ! methinks this loop should be shortend to nface_onP
      if (this%w_face(j,1) > 0.0_r8) then
        this%rho_fc(j) = max(min_face_rho, this%rho_fc(j)/this%w_face(j,1))
      else
        this%rho_fc(j) = 0.0_r8 ! this should already be true but...
      end if
    end do

    ! harmonic averaged face viscosity
    ! special cases:
    ! 1) If the face has only one cell neighbor (i.e. a boundary cell)
    ! - the face viscosity is the cell viscosity
    ! 2) If the face is shared by a solid cell (where mu_cc == 0) and a fluid cell
    ! - the face viscosity is the fluid cell viscosity
    ! 3) If the face is shared by two solid cells
    ! - the face viscosity is 0
    this%w_face(:,:) = 0.0_r8
    do i = 1, mesh%ncell
      associate(fi => mesh%cface(mesh%xcface(i):mesh%xcface(i+1)-1))
        do j = 1, size(fi)
          if (btest(mesh%cfpar(i), pos=j)) then
            this%w_face(fi(j), 1) = this%mu_cc(i)
          else
            this%w_face(fi(j), 2) = this%mu_cc(i)
          end if
        end do
      end associate
    end do
    do j = 1, mesh%nface ! methinks this loop should be shortend to nface_onP
      associate (w1 => this%w_face(j,1), w2 => this%w_face(j,2))
        if (w1 > 0.0_r8 .and. w2 > 0.0_r8) then
          this%mu_fc(j) = 2.0_r8*w1*w2/(w1+w2)
        else if (w1 > 0.0_r8) then
          this%mu_fc(j) = w1
        else if (w2 > 0.0_r8) then
          this%mu_fc(j) = w2
        else
          this%mu_fc(j) = 0.0_r8
        end if
      end associate
    end do

    if (ini) then
      this%rho_cc_n(:) = this%rho_cc(:)
      this%rho_fc_n(:) = this%rho_fc(:)
      this%mu_cc_n(:) = this%mu_cc(:)
      this%mu_fc_n(:) = this%mu_fc(:)
      this%vof_n(:) = this%vof(:)
      this%vof_novoid_n(:) = this%vof_novoid(:)
    end if

  end subroutine update

end module flow_props_type
