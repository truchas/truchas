!!
!! Zach Jibben <zjibben@lanl.gov>
!! August 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_timers
  use truchas_logging_services
  use unstr_mesh_type
  use scalar_func_containers
  use sm_bc_manager_type
  use integration_geometry_type
  use sm_material_model_type
  use viscoplastic_solver_type
  implicit none
  private

  type, public :: sm_model
    private
    !! expose these to precon type
    type(unstr_mesh), pointer, public :: mesh => null() ! unowned reference
    type(integration_geometry), pointer, public :: ig => null()
    type(sm_bc_manager), public :: bc

    type(sm_material_model), pointer :: matl_model => null()
    type(viscoplastic_solver) :: vp_solver
    logical, public :: use_uniform_scaling_factor = .true.

    !! state variables -- expose to driver for checkpointing
    real(r8), allocatable :: strain_plastic(:,:), strain_total(:,:)
    real(r8), allocatable, public :: strain_plastic_old(:,:), strain_thermal_old(:,:), &
        strain_total_old(:,:), strain_pc(:,:)
    real(r8), allocatable :: dstrain_plastic_dt_old(:,:), dstrain_plastic_dt(:,:)

    !! expose to the preconditioner
    real(r8), allocatable, public :: rhs(:,:)
    real(r8), allocatable, public :: lame1_n(:), lame2_n(:), scaling_factor(:)
    real(r8), public :: bnorm3, penalty

    !! expose for computing viz fields
    real(r8), allocatable, public :: lame1(:), lame2(:), strain_thermal(:,:)

    real(r8) :: eff_max_cell_width, tlast, tcurrent
    real(r8), allocatable :: density_c(:), density_n(:)
    real(r8), pointer :: vof(:,:) => null() ! unowned reference
    real(r8), pointer :: temperature(:) => null() ! unowned reference

    !! input parameters
    real(r8), allocatable :: body_force_density(:)

    !! solver parameters -- expose to nlsol
    real(r8), public :: atol, rtol, ftol
  contains
    procedure :: init
    procedure :: init_viscoplastic_rate
    procedure :: apply_state
    procedure :: write_checkpoint
    procedure :: read_checkpoint
    procedure :: update_properties
    procedure :: compute_residual
    procedure :: compute_forces
    procedure :: compute_reference_density
    procedure :: accept_state
    procedure :: viscoplasticity_enabled
    procedure :: compute_viscoplastic_precon
    procedure :: get_plastic_strain
    procedure :: get_plastic_strain_rate
    procedure :: plastic_strain_cell
    procedure :: plastic_strain_rate_cell
    procedure, nopass :: tensor_dot
    procedure, nopass :: strain_tensor
    procedure, private :: compute_total_strain
    final :: sm_model_finalize
  end type sm_model

contains

  elemental subroutine sm_model_finalize(this)
    type(sm_model), intent(inout) :: this
    if (associated(this%matl_model)) deallocate(this%matl_model)
    if (associated(this%ig)) deallocate(this%ig)
  end subroutine sm_model_finalize


  subroutine init(this, mesh, params, nmat, lame1f, lame2f, densityf, reference_density, vp)

    use parameter_list_type
    use viscoplastic_material_model_types

    class(sm_model), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(in) :: nmat
    type(scalar_func_box), allocatable, intent(inout) :: lame1f(:), lame2f(:), densityf(:)
    real(r8), intent(in) :: reference_density(:)
    type(viscoplastic_material_model_box), allocatable, intent(inout) :: vp(:)

    integer :: stat
    real(r8) :: traction, distance
    character(:), allocatable :: errmsg
    type(parameter_list), pointer :: plist => null()

    this%mesh => mesh
    allocate(this%ig)
    call this%ig%init(mesh)
    allocate(this%matl_model)
    call this%matl_model%init(params, nmat, lame1f, lame2f, densityf, reference_density, vp)

    call params%get('body-force-density', this%body_force_density, default=[0d0, 0d0, 0d0])

    ASSERT(size(this%body_force_density) == 3)

    allocate(this%density_c(this%mesh%ncell), this%density_n(this%mesh%nnode_onP), &
        this%strain_total(6,this%ig%npt), &
        this%strain_thermal(6,this%mesh%ncell), this%strain_pc(6,this%mesh%ncell), &
        this%lame1(this%mesh%ncell), this%lame2(this%mesh%ncell), &
        this%lame1_n(this%mesh%nnode), this%lame2_n(this%mesh%nnode), &
        this%rhs(3,this%mesh%nnode_onP), this%scaling_factor(this%mesh%nnode_onP))

    this%eff_max_cell_width = eff_max_cell_width()
    this%scaling_factor = 1

    call params%get('contact-distance', distance, default=1e-7_r8)
    call params%get('contact-normal-traction', traction, default=1e4_r8)
    call params%get('contact-penalty', this%penalty, default=1e3_r8)
    call params%get('abs-displ-tol', this%atol, default=1e100_r8)
    call params%get('rel-displ-tol', this%rtol, default=1e100_r8)
    call params%get('rel-stress-tol', this%ftol, default=1e-10_r8)
    plist => params%sublist('bc')
    call this%bc%init(plist, mesh, this%ig, this%penalty, distance, traction, stat, errmsg)
    if (stat /= 0) call TLS_fatal('Failed to build solid mechanics boundary conditions: '//errmsg)

    !! viscoplasticity
    if (this%matl_model%viscoplasticity_enabled) then
      allocate(this%strain_thermal_old(6,this%mesh%ncell), this%strain_total_old(6,this%ig%npt), &
          this%strain_plastic_old(6,this%ig%npt), this%strain_plastic(6,this%ig%npt), &
          this%dstrain_plastic_dt_old(6,this%ig%npt), this%dstrain_plastic_dt(6,this%ig%npt))
      this%strain_thermal_old = 0
      this%strain_total_old = 0
      this%strain_plastic_old = 0
      this%dstrain_plastic_dt_old = 0
      this%strain_pc = 0 ! phase-change strain is never used
      this%strain_plastic = 0
      this%dstrain_plastic_dt = 0
      this%tlast = 0
      this%tcurrent = 0
      plist => params%sublist('viscoplastic-solver')
      call this%vp_solver%init(plist, this%ig, this%matl_model)
    end if

  contains

    ! Used only for the scaling factor computed in update_properties.
    real(r8) function eff_max_cell_width()
      use parallel_communication, only: global_maxval
      eff_max_cell_width = maxval(this%mesh%volume(:this%mesh%ncell_onP)**(1.0_r8/3.0_r8))
      eff_max_cell_width = global_maxval(eff_max_cell_width)
    end function eff_max_cell_width

  end subroutine init


  !! Compute the reference density, based on input volume fractions. This is
  !! done during initialization.
  subroutine compute_reference_density(this, vof)
    class(sm_model), intent(inout) :: this
    real(r8), intent(in) :: vof(:,:)
    integer :: j
    do j = 1, this%mesh%ncell
      this%density_c(j) = sum(this%matl_model%reference_density * vof(:,j))
    end do
  end subroutine compute_reference_density


  !! Computes the Lame parameters, density, and RHS of the system. This routine
  !! should be called before compute_residual, after which the residual can be
  !! repeatedly computed for any number of different inputs (displacements).
  subroutine update_properties(this, vof, temperature_cc)

    use sm_bc_utilities, only: compute_stress
    use parallel_communication, only: global_maxval

    class(sm_model), intent(inout) :: this
    real(r8), intent(in), target :: vof(:,:), temperature_cc(:)

    integer :: j, n, m, xp, p
    real(r8) :: s, state(1), thermal_stress(6,this%mesh%ncell)
    real(r8) :: density_old, dstrain, max_lame, max_lame1, max_lame2

    ASSERT(size(vof, dim=1) == this%matl_model%nmat)
    ASSERT(size(vof, dim=2) == this%mesh%ncell)
    ASSERT(size(temperature_cc) == this%mesh%ncell)

    call start_timer("properties")
    this%vof => vof
    this%temperature => temperature_cc

    do j = 1, this%mesh%ncell
      density_old = this%density_c(j)
      state(1) = temperature_cc(j)
      this%lame1(j) = 0
      this%lame2(j) = 0
      this%density_c(j) = 0
      do m = 1, this%matl_model%nmat
        if (vof(m,j) == 0) cycle
        this%lame1(j) = this%lame1(j) + vof(m,j) * this%matl_model%lame1f(m)%f%eval(state)
        this%lame2(j) = this%lame2(j) + vof(m,j) * this%matl_model%lame2f(m)%f%eval(state)
        this%density_c(j) = this%density_c(j) + vof(m,j) * this%matl_model%densityf(m)%f%eval(state)
      end do

      if (this%density_c(j) /= 0) then
        density_old = sum(this%matl_model%reference_density * vof(:,j))
        !dstrain = (density_old / this%density_c(j)) ** (1.0_r8 / 3) - 1
        dstrain = log(density_old / this%density_c(j)) / 3
        this%strain_thermal(:3,j) = dstrain
        this%strain_thermal(4:,j) = 0
        call compute_stress(this%lame1(j), this%lame2(j), this%strain_thermal(:,j), &
            thermal_stress(:,j))
      else
        thermal_stress(:,j) = 0
        this%density_c(j) = sum(this%matl_model%reference_density * vof(:,j)) ! for the next iteration
        this%strain_thermal(:,j) = 0
      end if

      this%strain_pc(:,j) = 0 ! TODO: phase-chage strain
    end do

    call cell_to_node(this%mesh, this%density_c, this%density_n)
    call cell_to_node(this%mesh, this%lame1, this%lame1_n)
    call cell_to_node(this%mesh, this%lame2, this%lame2_n)
    call this%mesh%node_imap%gather_offp(this%lame1_n)
    call this%mesh%node_imap%gather_offp(this%lame2_n)

    ! right hand side & scaling factor
    ! The scaling factor is aimed to be an approximation of the ratio stress /
    ! displacement, so that contact & Dirichlet BCs have roughly the same
    ! residual error order of magnitude as stress errors.
    this%bnorm3 = 0
    max_lame1 = global_maxval(abs(this%lame1(:this%mesh%ncell_onP)))
    max_lame2 = global_maxval(abs(this%lame2(:this%mesh%ncell_onP)))
    max_lame = max(max_lame1, max_lame2)

    do n = 1, this%mesh%nnode_onP
      if (this%use_uniform_scaling_factor) then
        this%scaling_factor(n) = max_lame * this%eff_max_cell_width / this%penalty
      else
        this%scaling_factor(n) = max(1.0_r8, this%lame1_n(n), this%lame2_n(n)) &
            * this%eff_max_cell_width / this%penalty
      end if
      this%scaling_factor(n) = max_lame * this%eff_max_cell_width / this%penalty
      !this%scaling_factor(n) = 1

      associate (np => this%ig%npoint(this%ig%xnpoint(n):this%ig%xnpoint(n+1)-1))
        this%rhs(:,n) = -this%body_force_density * this%density_n(n) * this%ig%volume(n)
        do xp = 1, size(np)
          p = np(xp)
          j = this%ig%pcell(p)
          s = merge(-1, 1, this%ig%nppar(xp,n))
          this%rhs(:,n) = this%rhs(:,n) + s * tensor_dot(thermal_stress(:,j), this%ig%n(:,p))
        end do
      end associate

      ! this%bnorm3 = max(this%bnorm3, &
      !     maxval(abs(this%rhs(:,n))) / this%scaling_factor(n), &
      !     this%lame1_n(n) / this%eff_max_cell_width**2 / this%scaling_factor(n), &
      !     this%lame2_n(n) / this%eff_max_cell_width**2 / this%scaling_factor(n))
    end do
    !this%bnorm3 = global_maxval(this%bnorm3)
    this%bnorm3 = 1

    call stop_timer("properties")

  end subroutine update_properties


  subroutine accept_state(this)
    class(sm_model), intent(inout) :: this
    if (this%matl_model%viscoplasticity_enabled) then
      this%strain_total_old = this%strain_total
      this%strain_plastic_old = this%strain_plastic
      this%strain_thermal_old = this%strain_thermal
      this%dstrain_plastic_dt_old = this%dstrain_plastic_dt
      this%tlast = this%tcurrent
    end if
  end subroutine accept_state


  !! Set a viscoplastic solver state, e.g. from a restart file.
  subroutine apply_state(this, tlast, strain_total, strain_thermal, strain_plastic, &
      dstrain_plastic_dt)
    class(sm_model), intent(inout) :: this
    real(r8), intent(in) :: tlast, strain_total(:,:), strain_thermal(:,:), strain_plastic(:,:), &
        dstrain_plastic_dt(:,:)
    this%tlast = tlast
    this%strain_total_old = strain_total
    this%strain_thermal_old = strain_thermal
    this%strain_plastic_old = strain_plastic
    this%dstrain_plastic_dt_old = dstrain_plastic_dt
  end subroutine apply_state


  !! Computes a cell-to-node averaging, weighted by volume It assumes x_cell is
  !! of length ncell, and x_node is of length nnode_onP, and that each of the
  !! cells adjacent to an onP node is available to this rank.
  subroutine cell_to_node(mesh, x_cell, x_node)

    type(unstr_mesh), intent(in) :: mesh
    real(r8), intent(in) :: x_cell(:)
    real(r8), intent(out) :: x_node(:)

    integer :: c, n, xn
    real(r8) :: w_node(mesh%nnode_onP)

    ASSERT(size(x_cell) == mesh%ncell)
    ASSERT(size(x_node) >= mesh%nnode_onP)

    x_node = 0
    w_node = 0
    do c = 1, mesh%ncell
      associate (cn => mesh%cnode(mesh%xcnode(c):mesh%xcnode(c+1)-1))
        do xn = 1, size(cn)
          n = cn(xn)
          if (n > mesh%nnode_onP) cycle
          x_node(n) = x_node(n) + mesh%volume(c)*x_cell(c)
          w_node(n) = w_node(n) + mesh%volume(c)
        end do
      end associate
    end do

    do n = 1, mesh%nnode_onP
      x_node(n) = x_node(n) / w_node(n)
    end do

  end subroutine cell_to_node


  subroutine init_viscoplastic_rate(this)
    class(sm_model), intent(inout) :: this
    if (.not.this%matl_model%viscoplasticity_enabled) return
    this%tcurrent = this%tlast
    call this%vp_solver%compute_plastic_strain(0.0_r8, this%temperature, this%lame1, this%lame2, this%vof, this%strain_pc, &
        this%strain_total_old, this%strain_thermal_old, this%strain_plastic_old, this%dstrain_plastic_dt_old, &
        this%strain_total, this%strain_thermal, this%strain_plastic, this%dstrain_plastic_dt)
    this%strain_plastic_old = 0 !this%strain_plastic
    this%dstrain_plastic_dt_old = this%dstrain_plastic_dt
  end subroutine init_viscoplastic_rate


  !! This evaluates the equations on page 1765 of Bailey & Cross 1995.
  subroutine compute_residual(this, t, displ, r)

    class(sm_model), intent(inout) :: this ! inout to update BCs
    real(r8), intent(in) :: t
    real(r8), intent(inout) :: displ(:,:) ! need to update halo
    real(r8), intent(out) :: r(:,:)

    integer :: n

    ASSERT(size(r,dim=1) == 3 .and. size(r,dim=2) >= this%mesh%nnode_onP)

    call compute_forces(this, t, displ, r)

    call start_timer("residual")
    call this%bc%apply_nontraction(t, this%scaling_factor, displ, r)
    do n = 1, this%mesh%nnode_onP
      r(:,n) = r(:,n) / this%scaling_factor(n)
    end do
    call stop_timer("residual")

  end subroutine compute_residual


  !! This evaluates the equations on page 1765 of Bailey & Cross 1995.
  subroutine compute_forces(this, t, displ, r)

    use sm_bc_utilities, only: compute_stress

    class(sm_model), intent(inout) :: this ! inout to update BCs
    real(r8), intent(in) :: t
    real(r8), intent(inout) :: displ(:,:)
    real(r8), intent(out) :: r(:,:)

    integer :: n, xp, p, j
    real(r8) :: s, stress(6), lhs(3), lhs_strain(6), dt

    ASSERT(size(r,dim=1) == 3 .and. size(r,dim=2) >= this%mesh%nnode_onP)

    call start_timer("residual")

    call this%compute_total_strain(displ, this%strain_total)

    ! viscoplasticity is updated outside the main stress loop so that we
    ! calculate viscoplasticity at integration points associated with
    ! off-process nodes *and* on-process cells. This is purely for visualization.
    if (this%matl_model%viscoplasticity_enabled .and. t /= this%tlast) then
      this%tcurrent = t ! to be read by accept_state
      dt = t - this%tlast
      call this%vp_solver%compute_plastic_strain(dt, this%temperature, this%lame1, this%lame2, &
          this%vof, this%strain_pc, &
          this%strain_total_old, this%strain_thermal_old, &
          this%strain_plastic_old, this%dstrain_plastic_dt_old, &
          this%strain_total, this%strain_thermal, this%strain_plastic, this%dstrain_plastic_dt)
    end if

    call start_timer("stress")
    do n = 1, this%mesh%nnode_onP
      associate (np => this%ig%npoint(this%ig%xnpoint(n):this%ig%xnpoint(n+1)-1))
        lhs = 0
        do xp = 1, size(np)
          p = np(xp)
          j = this%ig%pcell(p)

          lhs_strain = this%strain_total(:,p) ! thermal strain is accounted in the RHS

          if (this%matl_model%viscoplasticity_enabled .and. t /= this%tlast) &
              lhs_strain = lhs_strain - this%strain_plastic(:,p)

          call compute_stress(this%lame1(j), this%lame2(j), lhs_strain, stress)
          s = merge(-1, 1, this%ig%nppar(xp,n))
          lhs = lhs + s * tensor_dot(stress, this%ig%n(:,p))
        end do

        r(:,n) = lhs - this%rhs(:,n)
      end associate
    end do
    call stop_timer("stress")

    call this%bc%apply_traction(t, r)

    call stop_timer("residual")

  end subroutine compute_forces


  !! From the node-centered displacement, compute the integration-point-centered total strain
  subroutine compute_total_strain(this, displ, total_strain)

    class(sm_model), intent(in) :: this
    real(r8), intent(inout) :: displ(:,:) ! inout for off-rank halo update
    real(r8), intent(out) :: total_strain(:,:)

    integer :: p, j, d
    real(r8) :: grad_displ(3,3)

    ASSERT(size(displ, dim=2) == this%mesh%nnode) ! need off-rank halo

    call start_timer("strain")
    call this%mesh%node_imap%gather_offp(displ)

    do j = 1, this%mesh%ncell
      associate (cn => this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1))
        do p = this%ig%xcpoint(j), this%ig%xcpoint(j+1)-1
          !! Compute the derivative of a scalar phi in global coordinates using the
          !! formula on page 1765 of Bailey & Cross 1995--linear interpolation.
          !! The Jacobian inverse multiplication is already performed and stored
          !! in the grad_shape matrix.
          do d = 1, 3
            grad_displ(:,d) = matmul(this%ig%grad_shape(p)%p, displ(d,cn))
          end do
          !grad_displ = matmul(this%ig%grad_shape(p)%p, transpose(displ(:,cn)))
          total_strain(:,p) = this%strain_tensor(grad_displ)
        end do
      end associate
    end do

    call stop_timer("strain")

  end subroutine compute_total_strain


  !! Computes the dot product t_i = a_ij * n_j, given a_ij like the
  !! total_strain: (a_xx, a_yy, a_zz, a_xy, a_xz, a_yz), and assuming
  !! symmetry.
  function tensor_dot(a, n) result(t)
    real(r8), intent(in) :: a(:), n(:)
    real(r8) :: t(3)
    ASSERT(size(a) == 6)
    ASSERT(size(n) == 3)
    t(1) = a(1)*n(1) + a(4)*n(2) + a(5)*n(3)
    t(2) = a(4)*n(1) + a(2)*n(2) + a(6)*n(3)
    t(3) = a(5)*n(1) + a(6)*n(2) + a(3)*n(3)
  end function tensor_dot


  function strain_tensor(grad_displ)
    real(r8), intent(in) :: grad_displ(:,:)
    real(r8) :: strain_tensor(6)
    strain_tensor(1) = grad_displ(1,1) ! exx
    strain_tensor(2) = grad_displ(2,2) ! eyy
    strain_tensor(3) = grad_displ(3,3) ! ezz
    strain_tensor(4) = (grad_displ(1,2) + grad_displ(2,1)) / 2 ! exy
    strain_tensor(5) = (grad_displ(1,3) + grad_displ(3,1)) / 2 ! exz
    strain_tensor(6) = (grad_displ(2,3) + grad_displ(3,2)) / 2 ! eyz
  end function strain_tensor


  pure logical function viscoplasticity_enabled(this)
    class(sm_model), intent(in) :: this
    viscoplasticity_enabled = this%matl_model%viscoplasticity_enabled
  end function

  pure subroutine get_plastic_strain(this, plastic_strain)
    class(sm_model), intent(in) :: this
    real(r8), intent(out), allocatable :: plastic_strain(:,:)
    plastic_strain = this%strain_plastic
  end subroutine

  pure subroutine get_plastic_strain_rate(this, plastic_strain_rate)
    class(sm_model), intent(in) :: this
    real(r8), intent(out), allocatable :: plastic_strain_rate(:,:)
    plastic_strain_rate = this%dstrain_plastic_dt
  end subroutine

  !! Currently just a dumb averaging from integration points to cell center,
  !! purely for visualization.
  !!
  !! NOTE: This routine assumes that plastic strain has been calculated at all
  !! integration points within cells owned by this rank. It is also true that
  !! a rank does not necessarily own all the nodes associated with an owned
  !! cell.
  function plastic_strain_cell(this, j)

    class(sm_model), intent(in) :: this
    integer, intent(in) :: j
    real(r8) :: plastic_strain_cell(6)

    integer :: np, p0, p1

    ASSERT(j <= this%mesh%ncell_onP)
    !ASSERT(all(this%mesh%cnode(this%mesh%xcnode(j):this%mesh%xcnode(j+1)-1) <= this%mesh%nnode_onP))

    if (.not.this%viscoplasticity_enabled()) then
      plastic_strain_cell = 0
      return
    end if

    p0 = this%ig%xcpoint(j)
    p1 = this%ig%xcpoint(j+1)-1
    np = p1 - p0 + 1
    plastic_strain_cell = sum(this%strain_plastic_old(:,p0:p1), dim=2) / np

  end function plastic_strain_cell

  function plastic_strain_rate_cell(this, j, stress) result(g)
    class(sm_model), intent(inout) :: this
    integer, intent(in) :: j
    real(r8), intent(in) :: stress(:)
    real(r8) :: g
    if (this%viscoplasticity_enabled()) then
      g = this%vp_solver%model%strain_rate2(&
          this%temperature(j), this%lame1(j), this%lame2(j), this%vof(:,j), &
          stress)
    else
      g = 0
    end if
  end function plastic_strain_rate_cell


  subroutine compute_viscoplastic_precon(this, precon)
    class(sm_model), intent(inout) :: this
    real(r8), intent(out) :: precon(:,:,:)
    call this%vp_solver%compute_precon(this%temperature, this%lame1, this%lame2, this%vof, &
        this%strain_total_old, this%strain_thermal_old, this%strain_pc, this%strain_plastic_old, &
        precon)
  end subroutine compute_viscoplastic_precon


  !! This is a custom integration-point writer used for restarts.
  subroutine write_checkpoint(this, seq)

    use truchas_h5_outfile, only: th5_seq_group
    use truchas_danu_output_tools

    class(sm_model), intent(in) :: this
    class(th5_seq_group), intent(in) :: seq

    ! NB: The thermal strain is a cell-center field and already output as part of visualization.
    call write_seq_cell_field(seq, reshape(square_array(this%strain_total_old), [6*12, this%mesh%ncell_onP]), 'strain_total', for_viz=.false.)
    call write_seq_cell_field(seq, reshape(square_array(this%strain_plastic_old), [6*12, this%mesh%ncell_onP]), 'strain_plastic', for_viz=.false.)
    call write_seq_cell_field(seq, reshape(square_array(this%dstrain_plastic_dt_old), [6*12, this%mesh%ncell_onP]), 'dstrain_plastic_dt', for_viz=.false.)

  contains

    ! Convert a ragged-array with points for each cell-edge pair into a square array
    function square_array(a) result(b)

      real(r8), intent(in) :: a(:,:)
      !real(r8), intent(out) :: b(:,:,:)
      real(r8) :: b(6, 12, this%mesh%ncell_onP)

      integer :: j, p, e

      ASSERT(size(a,dim=1) == size(b,dim=1))
      ASSERT(size(a,dim=2) >= this%ig%xcpoint(this%mesh%ncell_onP+1)-1)

      b = 0
      do j = 1, this%mesh%ncell_onP
        do p = this%ig%xcpoint(j), this%ig%xcpoint(j+1)-1
          e = p - this%ig%xcpoint(j) + 1
          b(:,e,j) = a(:,p)
        end do
      end do

    end function square_array

  end subroutine write_checkpoint


  subroutine read_checkpoint(this, unit)

    use restart_utilities, only: read_dist_array
    use restart_variables, only: restart_t

    class(sm_model), intent(inout) :: this
    integer, intent(in) :: unit

    real(r8) :: field(6*12,this%mesh%ncell)

    this%tlast = restart_t

    associate (pcell => this%mesh%xcell(:this%mesh%ncell_onp))
      call read_dist_array(unit, field(:,:this%mesh%ncell_onp), pcell, 'READ_STRAIN_TOTAL: error reading strain_total records')
      call this%mesh%cell_imap%gather_offp(field)
      call generate_integrationpt_array(field, this%strain_total_old)

      call read_dist_array(unit, this%strain_thermal_old(:,:this%mesh%ncell_onp), pcell, 'READ_STRAIN_THERMAL: error reading strain_thermal records')
      call this%mesh%cell_imap%gather_offp(this%strain_thermal_old)

      call read_dist_array(unit, field(:,:this%mesh%ncell_onp), pcell, 'READ_STRAIN_PLASTIC: error reading strain_plastic records')
      call this%mesh%cell_imap%gather_offp(field)
      call generate_integrationpt_array(field, this%strain_plastic_old)

      call read_dist_array(unit, field(:,:this%mesh%ncell_onp), pcell, 'READ_DSTRAIN_PLASTIC_DT: error reading dstrain_plastic_dt records')
      call this%mesh%cell_imap%gather_offp(field)
      call generate_integrationpt_array(field, this%dstrain_plastic_dt_old)
    end associate

  contains

    subroutine generate_integrationpt_array(a, b)

      real(r8), intent(in) :: a(:,:)
      real(r8), intent(out) :: b(:,:)

      integer :: j, p, e
      real(r8), allocatable :: a2(:,:,:)

      ASSERT(size(a,dim=2) >= this%mesh%ncell)
      ASSERT(size(b,dim=2) >= this%ig%xcpoint(this%mesh%ncell+1)-1)

      a2 = reshape(a, [6, 12, this%mesh%ncell])

      b = 0
      do j = 1, this%mesh%ncell
        do p = this%ig%xcpoint(j), this%ig%xcpoint(j+1)-1
          e = p - this%ig%xcpoint(j) + 1
          b(:,p) = a2(:,e,j)
        end do
      end do

    end subroutine generate_integrationpt_array

  end subroutine read_checkpoint

end module sm_model_type
