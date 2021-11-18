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
  use unstr_mesh_type
  use scalar_func_containers
  use sm_bc_manager_type
  use integration_geometry_type
  implicit none
  private

  type, public :: sm_model
    private
    ! expose these to precon type
    type(unstr_mesh), pointer, public :: mesh => null() ! unowned reference
    type(integration_geometry), public :: ig
    type(sm_bc_manager), public :: bc

    ! expose to the preconditioner
    real(r8), allocatable, public :: rhs(:,:)
    real(r8), allocatable, public :: lame1_n(:), lame2_n(:), scaling_factor(:)

    ! expose for computing viz fields
    real(r8), allocatable, public :: lame1(:), lame2(:), thermal_strain(:,:)

    integer :: nmat
    real(r8) :: max_cell_halfwidth, penalty
    real(r8), allocatable :: density_c(:), density_n(:), reference_density(:)
    real(r8), allocatable :: delta_temperature(:)
    type(scalar_func_box), allocatable :: lame1f(:), lame2f(:), densityf(:)

    !! input parameters
    real(r8), allocatable :: body_force_density(:)
    real(r8) :: strain_limit

    !! solver parameters -- expose to nlsol
    real(r8), public :: atol, rtol, ftol
  contains
    procedure :: init
    procedure :: size => model_size
    procedure :: update_properties
    procedure :: compute_residual
    procedure :: compute_forces
    procedure :: compute_reference_density
    procedure, nopass :: compute_stress
    procedure, nopass :: tensor_dot
    procedure, private :: compute_total_strain
  end type sm_model

contains

  subroutine init(this, mesh, params, nmat, lame1f, lame2f, densityf, reference_density)

    use parameter_list_type
    use truchas_logging_services

    class(sm_model), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(in) :: nmat
    type(scalar_func_box), allocatable, intent(inout) :: lame1f(:), lame2f(:), densityf(:)
    real(r8), intent(in) :: reference_density(:)

    integer :: stat
    real(r8) :: traction, distance
    character(:), allocatable :: errmsg
    type(parameter_list), pointer :: plist => null()

    ASSERT(nmat > 0)
    ASSERT(size(lame1f) == nmat)
    ASSERT(size(lame2f) == nmat)
    ASSERT(size(densityf) == nmat)
    ASSERT(size(reference_density) == nmat)

    this%mesh => mesh
    call this%ig%init(mesh)

    call params%get('body-force-density', this%body_force_density, default=[0d0, 0d0, 0d0])
    call params%get('strain-limit', this%strain_limit, default=1e-10_r8)

    ASSERT(size(this%body_force_density) == 3)
    if (this%strain_limit < 0) call TLS_fatal('strain-limit must be >= 0')

    allocate(this%density_c(this%mesh%ncell), this%density_n(this%mesh%nnode_onP), &
        this%delta_temperature(this%mesh%ncell), &
        this%thermal_strain(6,this%mesh%ncell), &
        this%lame1(this%mesh%ncell), this%lame2(this%mesh%ncell), &
        this%lame1_n(this%mesh%nnode_onP), this%lame2_n(this%mesh%nnode_onP), &
        this%rhs(3,this%mesh%nnode_onP), &
        this%scaling_factor(this%mesh%nnode))

    this%max_cell_halfwidth = max_cell_halfwidth()
    this%nmat = nmat
    this%reference_density = reference_density
    call move_alloc(lame1f, this%lame1f)
    call move_alloc(lame2f, this%lame2f)
    call move_alloc(densityf, this%densityf)

    call params%get('contact-distance', distance, default=1e-7_r8)
    call params%get('contact-normal-traction', traction, default=1e4_r8)
    call params%get('contact-penalty', this%penalty, default=1e3_r8)
    call params%get('abs-displ-tol', this%atol, default=1e-10_r8)
    call params%get('rel-displ-tol', this%rtol, default=1e-10_r8)
    call params%get('abs-stress-tol', this%ftol, default=1e-10_r8)
    plist => params%sublist('bc')
    call this%bc%init(plist, mesh, this%ig, this%penalty, distance, traction, stat, errmsg)
    if (stat /= 0) call TLS_fatal('Failed to build solid mechanics boundary conditions: '//errmsg)

  contains

    ! compute the maximum cell half-width, used for the
    ! vector scaling factor applied to the system
    real(r8) function max_cell_halfwidth()
      use parallel_communication, only: global_maxval
      real(r8) :: halfwidth(mesh%ncell_onP)
      halfwidth = 0.5_r8 * this%mesh%volume(:this%mesh%ncell_onP)**(1.0_r8/3.0_r8)
      max_cell_halfwidth = global_maxval(halfwidth)
    end function max_cell_halfwidth

  end subroutine init


  !! Compute the reference density, based on input volume fractions. This is
  !! done during initialization.
  subroutine compute_reference_density(this, vof)
    class(sm_model), intent(inout) :: this
    real(r8), intent(in) :: vof(:,:)
    integer :: j
    do j = 1, this%mesh%ncell
      this%density_c(j) = sum(this%reference_density * vof(:,j))
    end do
  end subroutine compute_reference_density


  integer function model_size(this)
    class(sm_model), intent(in) :: this
    model_size = 3*this%mesh%nnode
  end function model_size


  !! Computes the Lame parameters, density, and RHS of the system. This routine
  !! should be called before compute_residual, after which the residual can be
  !! repeatedly computed for any number of different inputs (displacements).
  subroutine update_properties(this, vof, temperature_cc)

    use parallel_communication, only: global_maxval, global_minval
    use index_partitioning, only: gather_boundary

    class(sm_model), intent(inout) :: this
    real(r8), intent(in) :: temperature_cc(:), vof(:,:)

    integer :: j, n, m, xp, p
    real(r8) :: s, state(1), thermal_stress(6,this%mesh%ncell)
    real(r8) :: density_old, dstrain

    ASSERT(size(vof, dim=1) == this%nmat)
    ASSERT(size(vof, dim=2) == this%mesh%ncell)
    ASSERT(size(temperature_cc) == this%mesh%ncell)

    call start_timer("properties")

    do j = 1, this%mesh%ncell
      density_old = this%density_c(j)
      state(1) = temperature_cc(j)
      this%lame1(j) = 0
      this%lame2(j) = 0
      this%density_c(j) = 0
      do m = 1, this%nmat
        if (vof(m,j) == 0) cycle
        this%lame1(j) = this%lame1(j) + vof(m,j) * this%lame1f(m)%f%eval(state)
        this%lame2(j) = this%lame2(j) + vof(m,j) * this%lame2f(m)%f%eval(state)
        this%density_c(j) = this%density_c(j) + vof(m,j) * this%densityf(m)%f%eval(state)
      end do

      if (this%density_c(j) /= 0) then
        density_old = sum(this%reference_density * vof(:,j))
        !dstrain = (density_old / this%density_c(j)) ** (1.0_r8 / 3) - 1
        dstrain = log(density_old / this%density_c(j)) / 3
        this%thermal_strain(:3,j) = dstrain
        this%thermal_strain(4:,j) = 0
        call compute_stress(this%lame1(j), this%lame2(j), this%thermal_strain(:,j), &
            thermal_stress(:,j))
      else
        thermal_stress(:,j) = 0
        this%density_c(j) = sum(this%reference_density * vof(:,j)) ! for the next iteration
        this%thermal_strain(:,j) = 0
      end if
    end do

    call cell_to_node(this%mesh, this%density_c, this%density_n)
    call cell_to_node(this%mesh, this%lame1, this%lame1_n)
    call cell_to_node(this%mesh, this%lame2, this%lame2_n)

    ! right hand side & scaling factor
    do n = 1, this%mesh%nnode_onP
      if (this%lame2_n(n) /= 0) then
        !this%scaling_factor(n) = 1e-3_r8 * this%lame2_n(n) * this%max_cell_halfwidth
        !this%scaling_factor(n) = 1e-3_r8 * this%lame2_n(n) / this%ig%volume(n)
        !this%scaling_factor(n) = this%lame2_n(n) * this%max_cell_halfwidth
        !this%scaling_factor(n) = 2e3_r8 * this%lame2_n(n) * this%max_cell_halfwidth
        this%scaling_factor(n) = this%lame2_n(n) * this%max_cell_halfwidth / this%penalty
      else
        this%scaling_factor(n) = 1
      end if
      !this%scaling_factor(n) = 1 ! DEBUGGING

      associate (np => this%ig%npoint(this%ig%xnpoint(n):this%ig%xnpoint(n+1)-1))
        this%rhs(:,n) = this%body_force_density * this%density_n(n) * this%ig%volume(n)
        do xp = 1, size(np)
          p = np(xp)
          j = this%ig%pcell(p)
          s = merge(-1, 1, this%ig%nppar(xp,n))
          this%rhs(:,n) = this%rhs(:,n) + s * tensor_dot(thermal_stress(:,j), this%ig%n(:,p))
        end do
      end associate
    end do

    call gather_boundary(this%mesh%node_ip, this%scaling_factor)

    !this%contact_penalty = global_maxval(this%lame2_n / this%ig%volume)
    !this%contact_penalty = global_maxval(this%rhs)
    !this%contact_penalty = 1e3_r8
    !this%contact_penalty = global_maxval(this%lame2) * this%max_cell_halfwidth
    !this%contact_penalty = 1

    call stop_timer("properties")

  end subroutine update_properties


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


  !! This evaluates the equations on page 1765 of Bailey & Cross 1995.
  subroutine compute_residual(this, t, displ, r)

    use index_partitioning, only: gather_boundary

    class(sm_model), intent(inout) :: this ! inout to update BCs
    real(r8), intent(in) :: t
    real(r8), intent(inout) :: displ(:,:) ! need to update halo
    real(r8), intent(out) :: r(:,:)

    integer :: n

    ASSERT(size(displ,dim=1) == 3 .and. size(displ,dim=2) >= this%mesh%nnode_onP)
    ASSERT(size(r,dim=1) == 3 .and. size(r,dim=2) >= this%mesh%nnode_onP)

    ! get off-rank halo
    call start_timer("residual")
    call gather_boundary(this%mesh%node_ip, displ)
    call stop_timer("residual")

    call compute_forces(this, t, displ, r)

    call start_timer("residual")
    call this%bc%apply_nontraction(t, this%scaling_factor, displ, r)
    do n = 1, this%mesh%nnode_onP
      r(:,n) = r(:,n) / this%scaling_factor(n)
    end do
    r(:,this%mesh%nnode_onP+1:) = 0 ! clear out halo so it doesn't screw with norms
    call stop_timer("residual")

  end subroutine compute_residual


  !! This evaluates the equations on page 1765 of Bailey & Cross 1995.
  subroutine compute_forces(this, t, displ, r)

    use index_partitioning, only: gather_boundary

    class(sm_model), intent(inout) :: this ! inout to update BCs
    real(r8), intent(in) :: t, displ(:,:)
    real(r8), intent(out) :: r(:,:)

    integer :: n, xp, p, j
    real(r8) :: s, stress(6), total_strain(6,this%ig%npt), lhs(3)

    ASSERT(size(displ,dim=1) == 3 .and. size(displ,dim=2) >= this%mesh%nnode_onP)
    ASSERT(size(r,dim=1) == 3 .and. size(r,dim=2) >= this%mesh%nnode_onP)

    call start_timer("residual")

    call this%compute_total_strain(displ, total_strain)

    call start_timer("stress")
    do n = 1, this%mesh%nnode_onP
      associate (np => this%ig%npoint(this%ig%xnpoint(n):this%ig%xnpoint(n+1)-1))
        lhs = 0
        do xp = 1, size(np)
          p = np(xp)
          j = this%ig%pcell(p)

          call compute_stress(this%lame1(j), this%lame2(j), total_strain(:,p), stress)
          s = merge(-1, 1, this%ig%nppar(xp,n))
          lhs = lhs + s * tensor_dot(stress, this%ig%n(:,p))
        end do

        r(:,n) = lhs - this%rhs(:,n)
      end associate
    end do
    call stop_timer("stress")

    call this%bc%apply_traction(t, r)
    call gather_boundary(this%mesh%node_ip, r)

    call stop_timer("residual")

  end subroutine compute_forces


  !! From the node-centered displacement, compute the integration-point-centered total strain
  subroutine compute_total_strain(this, displ, total_strain)

    class(sm_model), intent(in) :: this
    real(r8), intent(in) :: displ(:,:)
    real(r8), intent(out) :: total_strain(:,:)

    integer :: p, j, d
    real(r8) :: grad_displ(3,3)

    ASSERT(size(displ, dim=2) == this%mesh%nnode) ! need off-rank halo

    call start_timer("strain")

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

          total_strain(1,p) = grad_displ(1,1) ! exx
          total_strain(2,p) = grad_displ(2,2) ! eyy
          total_strain(3,p) = grad_displ(3,3) ! ezz
          total_strain(4,p) = (grad_displ(1,2) + grad_displ(2,1)) / 2 ! exy
          total_strain(5,p) = (grad_displ(1,3) + grad_displ(3,1)) / 2 ! exz
          total_strain(6,p) = (grad_displ(2,3) + grad_displ(3,2)) / 2 ! eyz
        end do
      end associate
    end do

    call stop_timer("strain")

  end subroutine compute_total_strain


  !! Computes the stress from the strain and Lame parameters.
  subroutine compute_stress(lame1, lame2, strain, stress)
    real(r8), intent(in) :: lame1, lame2, strain(:)
    real(r8), intent(out) :: stress(:)
    ASSERT(size(strain) == 6)
    ASSERT(size(stress) == 6)
    stress = 2 * lame2 * strain
    stress(:3) = stress(:3) + lame1 * sum(strain(:3))
  end subroutine compute_stress


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

end module sm_model_type