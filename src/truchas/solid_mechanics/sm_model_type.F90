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
  use sm_bc_type
  use integration_geometry_type
  implicit none
  private

  type, public :: sm_model
    private
    ! expose these to precon type
    type(unstr_mesh), pointer, public :: mesh => null() ! unowned reference
    type(integration_geometry), public :: ig
    type(sm_bc), public :: bc

    integer :: nmat
    real(r8), allocatable :: density_c(:), density_n(:)
    real(r8), allocatable :: delta_temperature(:), lame1(:), lame2(:), rhs(:,:)
    type(scalar_func_box), allocatable :: lame1f(:), lame2f(:), densityf(:)
    logical :: first_step = .true.

    !! input parameters
    real(r8), allocatable :: body_force_density(:)
    real(r8) :: contact_distance, contact_normal_traction, contact_penalty, strain_limit
  contains
    procedure :: init
    procedure :: size => model_size
    procedure :: update_properties
    procedure :: compute_residual
    procedure :: compute_lame_node_parameters
    procedure, nopass :: compute_stress
    procedure, nopass :: tensor_dot
    procedure, private :: compute_total_strain
  end type sm_model

contains

  subroutine init(this, mesh, params, nmat, lame1f, lame2f, densityf)

    use parameter_list_type
    use truchas_logging_services

    class(sm_model), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(in) :: nmat
    type(scalar_func_box), allocatable, intent(inout) :: lame1f(:), lame2f(:), densityf(:)

    integer :: stat
    character(:), allocatable :: errmsg
    type(parameter_list), pointer :: plist => null()

    ASSERT(nmat > 0)
    ASSERT(size(lame1f) == nmat)
    ASSERT(size(lame2f) == nmat)
    ASSERT(size(densityf) == nmat)

    this%mesh => mesh
    call this%ig%init(mesh)

    call params%get('body-force-density', this%body_force_density, default=[0d0, 0d0, 0d0])
    call params%get('contact-distance', this%contact_distance, default=1e-7_r8)
    call params%get('contact-normal-traction', this%contact_normal_traction, default=1e4_r8)
    call params%get('contact-penalty', this%contact_penalty, default=1e3_r8)
    call params%get('strain-limit', this%strain_limit, default=1e-10_r8)

    ASSERT(size(this%body_force_density) == 3)
    if (this%strain_limit < 0) call TLS_fatal('strain-limit must be >= 0')

    allocate(this%density_c(this%mesh%ncell), this%density_n(this%mesh%nnode_onP), &
        this%delta_temperature(this%mesh%ncell), &
        this%lame1(this%mesh%ncell), this%lame2(this%mesh%ncell), this%rhs(3,this%mesh%nnode_onP))
    this%density_c = 0 ! TODO -- set for first time step

    this%nmat = nmat
    call move_alloc(lame1f, this%lame1f)
    call move_alloc(lame2f, this%lame2f)
    call move_alloc(densityf, this%densityf)

    plist => params%sublist('bc')
    call this%bc%init(plist, mesh, stat, errmsg)
    if (stat /= 0) call TLS_fatal('Failed to build solid mechanics boundary conditions: '//errmsg)

  end subroutine init


  integer function model_size(this)
    class(sm_model), intent(in) :: this
    model_size = 3*this%mesh%nnode_onP
  end function model_size


  !! Computes the Lame parameters, density, and RHS of the system. This routine
  !! should be called before compute_residual, after which the residual can be
  !! repeatedly computed for any number of different inputs (displacements).
  subroutine update_properties(this, vof, temperature_cc)

    class(sm_model), intent(inout) :: this
    real(r8), intent(in) :: temperature_cc(:), vof(:,:)

    integer :: j, n, m, xp, p
    real(r8) :: s, state(1), thermal_strain(6), thermal_stress(6, this%mesh%ncell)
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

      if (this%density_c(j) /= 0 .and. .not.this%first_step) then
        dstrain = (density_old / this%density_c(j)) ** (1.0_r8 / 3) - 1
        thermal_strain(:3) = dstrain
        thermal_strain(4:) = 0
        call compute_stress(this%lame1(j), this%lame2(j), thermal_strain, thermal_stress(:,j))
      else
        thermal_stress(:,j) = 0
      end if
    end do

    call cell_to_node(this%mesh, this%density_c, this%density_n)

    ! right hand side
    do n = 1, this%mesh%nnode_onP
      associate (np => this%ig%npoint(this%ig%xnpoint(n):this%ig%xnpoint(n+1)-1))

        this%rhs(:,n) = this%body_force_density * this%density_n(n) * this%ig%volume(n)
        do xp = 1, size(np)
          p = np(xp)
          j = this%ig%pcell(p)

          s = merge(-1, 1, btest(this%ig%nppar(n),xp))
          this%rhs(:,n) = this%rhs(:,n) + s * tensor_dot(thermal_stress(:,j), this%ig%n(:,p))
        end do
      end associate
    end do

    this%first_step = .false.

    call stop_timer("properties")

  end subroutine update_properties


  !! Needed for the preconditioner.
  subroutine compute_lame_node_parameters(this, lame1_n, lame2_n)

    class(sm_model), intent(in) :: this
    real(r8), intent(out) :: lame1_n(:), lame2_n(:)

    call cell_to_node(this%mesh, this%lame1, lame1_n)
    call cell_to_node(this%mesh, this%lame2, lame2_n)

  end subroutine compute_lame_node_parameters


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

    class(sm_model), intent(inout) :: this ! inout to update BCs
    real(r8), intent(in) :: t, displ(:,:)
    real(r8), intent(out) :: r(:)

    integer :: n, xp, p, j, i, f, xn, d
    real(r8) :: s, stress(6), total_strain(6,this%ig%npt), lhs(3)

    ASSERT(size(displ,dim=1) == 3 .and. size(displ,dim=2) == this%mesh%nnode_onP)
    ASSERT(size(r) == 3*this%mesh%nnode_onP)

    call start_timer("residual")

    call this%compute_total_strain(displ, total_strain)

    do n = 1, this%mesh%nnode_onP
      associate (np => this%ig%npoint(this%ig%xnpoint(n):this%ig%xnpoint(n+1)-1))

        lhs = 0
        do xp = 1, size(np)
          p = np(xp)
          j = this%ig%pcell(p)

          call compute_stress(this%lame1(j), this%lame2(j), total_strain(:,p), stress)
          s = merge(-1, 1, btest(this%ig%nppar(n),xp))
          lhs = lhs + s * tensor_dot(stress, this%ig%n(:,p))
        end do

        r(3*(n-1)+1:3*(n-1)+3) = lhs - this%rhs(:,n)
      end associate
    end do

    do d = 1, 3
      ! At traction BCs the stress is defined along the face. Each
      ! face has traction discretized at the integration points, one for each node.
      ! The values array provides the equivalent of traction = tensor_dot(stress, n)
      ! where n is normal of the boundary with magnitude equal to the area of this
      ! integration region (associated with node and face center, not the whole face).
      if (allocated(this%bc%traction(d)%p)) then
        call this%bc%traction(d)%p%compute(t)
        associate (faces => this%bc%traction(d)%p%index, values => this%bc%traction(d)%p%value)
          do i = 1, size(faces)
            f = faces(i)
            do xn = this%mesh%xfnode(f), this%mesh%xfnode(f+1)-1
              n = this%mesh%fnode(xn)
              if (n > this%mesh%nnode_onP) cycle
              r(3*(n-1)+d) = r(3*(n-1)+d) + values(i)
            end do
          end do
        end associate
      end if

      ! Displacement BCs transform the matrix to a diagonal and the right hand
      ! side to the desired solution.
      if (allocated(this%bc%displacement(d)%p)) then
        call this%bc%displacement(d)%p%compute(t)
        associate (faces => this%bc%displacement(d)%p%index, &
            values => this%bc%displacement(d)%p%value)
          do i = 1, size(faces)
            f = faces(i)
            do xn = this%mesh%xfnode(f), this%mesh%xfnode(f+1)-1
              n = this%mesh%fnode(xn)
              if (n > this%mesh%nnode_onP) cycle
              r(3*(n-1)+d) = displ(d,n) - values(i)
            end do
          end do
        end associate
      end if
    end do

    call stop_timer("residual")

  end subroutine compute_residual


  !! From the node-centered displacement, compute the integration-point-centered total strain
  subroutine compute_total_strain(this, displ, total_strain)

    class(sm_model), intent(in) :: this
    real(r8), intent(in) :: displ(:,:)
    real(r8), intent(out) :: total_strain(:,:)

    integer :: p, j, d
    real(r8) :: grad_displ(3,3)

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
          total_strain(4,p) = grad_displ(1,2) + grad_displ(2,1) ! exy
          total_strain(5,p) = grad_displ(1,3) + grad_displ(3,1) ! exz
          total_strain(6,p) = grad_displ(2,3) + grad_displ(3,2) ! eyz
        end do
      end associate
    end do

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
    t(1) = a(1)*n(1) + a(4)*n(2) + a(5)*n(1)
    t(2) = a(2)*n(2) + a(4)*n(1) + a(6)*n(3)
    t(3) = a(3)*n(3) + a(5)*n(2) + a(6)*n(1)
  end function tensor_dot

end module sm_model_type