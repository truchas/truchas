!!
!! FLOW_PROJECTION_TYPE
!!
!! This module provides a type that encapsulates the projection of the face
!! velocities onto a solenoidal field
!!
!! Peter Brady <ptb@lanl.gov>
!! Feb 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines the derived type FLOW_PROJECTION that encapsulates the data
!!  and procedures for setting up and solving the pressure poisson system and
!!  projecting the face velocities onto a solenoidal field.
!!
!!  Interaction with this object occurs through its methods, generally via the flow_driver
!!
!!  READ_PARAMS(PARAMETER_LIST P) - Initializes the parameter list for this object.
!!    This must be call first for the object to behave sensibly.  P is passed on to
!!    a HYPRE_HYBRID instance.
!!
!!  INIT(FLOW_MESH M) - Allocates the data members for this object and initializes the
!!    internal solver
!!
!!  SETUP(FLOW_PROPS FP, REAL(R8) CELL_VELOCITY(:), REAL(R8) DT)
!!    prepares the lhs and rhs of the poisson system. CELL_VELOCITY is an array of cell
!!    centered velocitiesnormal velocities in a "ragged" format indexed directly with xcface.
!!    All face data associated with on process cells must be valid
!!
!!  SOLVE(REAL(R8) SOLUTION(:))
!!    solves the system and returns the solution.  Only the on_process portion of
!!    SOLUTION is valid.
!!
!!  ACCEPT()
!!    Copies new state to the `n` level for use in next time step


module flow_projection_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use flow_mesh_type
  use unstr_mesh_type
  use flow_props_type
  use flow_bc_type
  use fischer_guess_type
  use index_partitioning
  use hypre_hybrid_type
  use pcsr_matrix_type
  use parameter_list_type
  implicit none
  private

  public :: flow_projection

  type :: flow_projection
    private
    type(flow_mesh), pointer :: mesh ! unowned reference
    type(flow_bc), pointer :: bc ! unowned reference
    type(fischer_guess) :: fg
    type(hypre_hybrid) :: solver
    type(parameter_list) :: p
    real(r8), allocatable :: rhs(:)
    real(r8), allocatable :: grad_fc(:,:) ! face centered gradient
    real(r8), allocatable :: grad_p_rho_cc(:) ! dynamic pressure gradient over rho -> cell centered
    real(r8), allocatable :: grad_p_rho_cc_n(:) ! dynamic pressure gradient over rho -> cell centered
    real(r8), allocatable :: delta_p_cc(:)
    real(r8), allocatable :: gravity_head(:,:)
    real(r8), allocatable :: weights_cf(:,:) ! weights for f->c interpolation
  contains
    public
    procedure :: read_params
    procedure :: init
    procedure :: setup
    procedure :: solve
    procedure :: accept
    private
    procedure :: setup_face_velocity
    procedure :: setup_gravity
    procedure :: setup_solver
    procedure :: velocity_fc_correct
    procedure :: pressure_cc_correct
    procedure :: velocity_cc_correct
  end type flow_projection

contains

  subroutine read_params(this, p)
    class(flow_projection), intent(inout) :: this
    type(parameter_list), intent(inout) :: p

    this%p = p
    call this%fg%read_params(p)
  end subroutine read_params

  subroutine init(this, mesh, bc)
    class(flow_projection), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: mesh
    type(flow_bc), target, intent(inout) :: bc
    !-
    integer :: j, i, fi, ni, ri, f0
    type(pcsr_graph), pointer :: g
    type(pcsr_matrix), pointer :: A
    type(ip_desc), pointer :: row_ip
    type(unstr_mesh), pointer :: m
    type(parameter_list), pointer :: p
    real(r8) :: w

    this%mesh => mesh
    this%bc => bc
    m => mesh%mesh
    p = this%p

    allocate(weights_cf(2,m%nface))

    do j = 1, m%nface_onP
      associate (n => mesh%fcell(:,j), fc => mesh%face_centriod(:,j), cc => mesh%cell_centroid(:,j))
        if (n(1) > 0) then
          weights_cf(1,j) = sum((fc(:,j)-cc(:,n(1)))**2)
          weights_cf(2,j) = sum((fc(:,j)-cc(:,n(2)))**2)
          weights_cf(:,j) = weights_cf(:,j) / sum(weights_fc(:,j))
        else
          weights_cf(:,j) = [0.0_r8, 1.0_r8]
        end if
      end associate
    end do
    allocate(this%rhs(m%ncell))

    !! Create a CSR matrix graph for the pressure poisson system.
    allocate(g)
    row_ip => m%cell_ip
    call g%init (row_ip)
    do j = 1, m%ncell_onP
      call g%add_edge(j,j)
      associate (cn => m%cnhbr(m%xcnhbr(j):m%xcnhbr(j+1)-1))
        do i = 1, size(cn)
          if (cn(i) > 0) call g%add_edge(j, i)
        end do
      end associate
    end do
    call g%add_complete

    allocate(A)
    call A%init(g, take_graph=.true.)
    call this%solver%init(A, p)
    call this%fg%init(mesh)

  end subroutine init


  subroutine setup(this, dt, props, body_force, vel_cc, P_cc, vel_fn)
    class(flow_projection), intent(inout) :: this
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: vel_cc(:,:), dt, P_cc(:), body_force(:)
    real(r8), intent(out) :: vel_fn(:)

    call this%setup_gravity(props, body_force)
    call this%setup_face_velocity(dt, props, vel_cc, P_cc, vel_fn)
    call this%setup_solver(dt, props, vel_fn)
  end subroutine setup

  subroutine accept(this)
    class(flow_projection), intent(inout) :: this
    this%grad_p_rho_cc_n = this%grad_p_rho_cc
  end subroutine accept

  subroutine setup_solver(this, dt, props, vel)
    class(flow_projection), intent(inout) :: this
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: vel(:)
    !-
    type(pcsr_matrix), pointer :: A
    real(r8), pointer :: ds(:,:)
    integer :: i, j, fi, ni
    real(r8) :: coeff

    A = this%solver%matrix()
    call A%set_all(0.0_r8)

    ds = flow_gradient_coefficients()

    associate (m => this%mesh%mesh, rho_f => props%rho_fc)
      do j = 1, m%ncell_onP
        associate (nhbr => m%cnhbr(m%xcnhbr(j):m%xcnhbr(j+1)-1), &
            face => m%cface(m%xcface(j):m%xcface(j+1)-1))
          do i = 1, size(nhbr)
            fi = face(i)
            ni = nhbr(i)

            if (rho_f(fi) == 0.0_r8) then
              coeff = 0.0_r8
            elseif (btest(m%cfpar(j), pos=i)) then
              ! note that normal is already weighted by face area
              coeff = dot_product(this%dx_scaled(:,fi), m%normal(:,fi))/rho_f(fi)
            else
              -coeff = dot_product(this%dx_scaled(:,fi), m%normal(:,fi))/rho_f(fi)
            end if
            call A%add_to(j, j, coeff)
            if (ni > 0) call A%add_to(j, ni, -coeff)
          end do
        end associate
      end do
    end associate

    call this%solver%setup()

    !! FIXME: need correction for void/solid... maybe
    associate (m => this%mesh%mesh, rhs => this%rhs)
      do j = 1, m%ncell_onP
        rhs(j) = 0.0_r8
        associate( fn => m%cface(m%xcface(j):m%xcface(j+1)-1))
          do i = 1, size(fn)
            ! face velocity follows the convention of being positive in the inward
            ! face direction as is the normal.  Boundary conditions have already been
            ! applied in computing `vel` so we shouldn't need to use them here again
            if (btest(m%cfpar(j), pos=i)) then
              this%rhs(j) = this%rhs(j) - vel(fn(i))*m%area(fn(i))
            else
              this%rhs(j) = this%rhs(j) - vel(fn(i))*m%area(fn(i))
            end if
          end do
          this%rhs(j) = this%rhs(j) / dt
        end do
      end do
    end associate
  end subroutine setup_solver



  subroutine setup_gravity(this, props, body_force)
    class(flow_projection), intent(inout) :: this
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: body_force(:)
    !-

    integer :: j, i, n
    real(r8) :: rho

    associate( fcell => this%mesh%fcell, m => this%mesh%mesh, g => this%gravity_head, &
        cc => this%mesh%cell_centroid, fc => this%mesh%face_centroid, p => this%props)

      g = 0.0_r8

      do j = 1, m%nface_onP
        do i = 1, 2
          n = fcell(i,j)
          if (n == 0) cycle

          rho = p%rho_cc(n) + p%rho_delta_cc(n)
          g(i,j) = -rho*dot_product(body_force, cc(:,i)-fc(:,j))
        end do
      end do

      call gather_boundary(m%face_ip, g(1,:))
      call gather_boundary(m%face_ip, g(2,:))
    end associate
    ! Truchas has a gravity_limiter... Do we need this?
  end subroutine setup_gravity

  subroutine setup_face_velocity(this, dt, props, vel_cc, p_cc, vel_fn)
    class(flow_projection), intent(inout) :: this
    real(r8), intent(in) :: dt
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: vel_cc(:,:), p_cc(:)
    real(r8), intent(out) :: vel_fn(:)
    !-
    integer :: i, j

    ! cell -> face interpolation
    associate (m => this%mesh%mesh, v => this%vel_cc_star, gpn => this%grad_p_rho_cc_n, &
        w => this%weights_cf)
      do i=1, m%ncell
        v(:,i) = vel_cc(:,i) + dt*gpn(:,i)
      end do
      call interpolate_cf(vel_fn, v, w, this%bc%v_dirichlet, props%inactive_f, 0.0_r8)
    end associate

    ! subtract dynamic pressure grad
    associate (m => this%mesh%mesh, gfc => this%grad_fc, rho_fc => props%rho_fc)

      call gradient_cf(gfc, p_cc, time, this%bc%p_neumann, &
          this%bc%p_dirichlet, props%inactive_f, 0.0_rp, this%gravity_head)

      do j = 1, m%nface_onP
        if (rho_fc(j) > 0.0_r8) then
          vel_fn(j) = vel_fn(j) - dt*dot_product(m%normal(:,j),gfc(:,j))/(m%area(j)*rho_fc(j))
        else
          ! there may be something special with void here....
          vel_fn(j) = 0.0_r8
        end if
      end do
      call gather_boundary(m%face_ip, vel_fn)
    end associate
  end subroutine setup_face_velocity


  subroutine solve(this, dt, props, p_cc, vel_cc, vel_fc)
    class(flow_projection), intent(inout) :: this
    real(r8), intent(in) :: dt
    type(flow_props), intent(in) :: props
    real(r8), intent(inout) :: p_cc(:), vel_cc(:), vel_fc(:)
    !-
    real(r8) :: dp_vof, vof, dec
    integer :: ierr, j, i

    ! solve the pressure poisson system
    this%delta_p_cc = 0.0_r8
    call this%fg%guess(this%rhs, this%delta_p_cc)
    call this%solver%solve(this%rhs, this%delta_p_cc, ierr)
    if (ierr /= 0) call tls_error("projection solve unsuccessful")
    call this%fg%update(this%rhs, this%delta_p_cc, this%solver%matrix())

    ! not sure if this call is necesary
    call gather_boundary(this%mesh%mesh%cell_ip, this%delta_p_cc)

    ! correct face and cell centered quantities
    call this%velocity_fc_correct(dt, props, vel_fc)
    call this%pressure_cc_correct(props, p_cc)
    call this%velocity_cc_correct(dt, props, p_cc, vel_cc)
  end subroutine solve


  subroutine velocity_fc_correct(this, dt, props, vel_fc)
    class(flow_projection), intent(inout) :: this
    real(r8), intent(in) :: dt
    type(flow_props), intent(in) :: props
    real(r8), intent(inout) :: vel_fc(:)
    !-
    integer :: ierr, j, i

    ! face gradient of solution
    call gradient_cf(this%grad_fc, this%delta_p_cc, time, this%bc%p_neumann, &
        this%bc%p_dirichlet, props%inactive_f, 0.0_rp)

    ! divergence free face velocity
    associate (m => this%mesh%mesh)
      do j = 1, m%nface_onP
        if (props%inactive_f(j) > 0) then
          vel_fc(j) = 0.0_rp
        else
          vel_fc(j) = this%vel_fc_star(j) - &
              dt*dot_product(this%grad_fc(:,j),m%normal(:,j)/m%area(j))/props%rho_fc(j)
        end if
      end do

      call this%bc%v_dirichlet%compute(t)
      associate (index => this%bc%v_dirichlet%index, value => this%bc%v_dirichlet%value)
        do i = 1, size(index)
          j = index(i)
          ! is the convention to specify boundary velocity as into or out of domain??????
          ! the following assumes _out_
          if (props%inactive_f(j) == 0) &
              vel_fc(j) = -dot_product(value(:,j), m%normal(:,j))/m%area(j)
        end do
      end associate

      call gather_boundary(m%face_ip, vel_fc)
    end associate

    ! there may be some magic here about setting fluxing velocities on void cells.
    ! not sure if it's applicable when using face based (as opposed to side based) indexing
  end subroutine velocity_fc_correct


  subroutine pressure_cc_correct(this, props, p_cc)
    class(flow_projection), intent(inout) :: this
    real(r8), intent(in) :: dt
    type(flow_props), intent(in) :: props
    real(r8), intent(inout) :: p_cc(:)
    !-
    real(r8) :: dp_vof, vof, dec
    integer :: i

    associate (m => this%mesh%mesh, p => props, dp => this%delta_p_cc)
      ! remove null space from dp if applicable
      if (.not.(this%bc%dirichlet_p .or. p%any_void)) then
        ! maybe this should use inactive_c... needs to be consistenent with how
        ! the operator is constructed
        dp_vof = 0.0_r8
        vof = 0.0_r8
        do i = 1, m%ncell_onP
          vof = vof + p%vof(i)
          dp_vof = dp_vof + p%vof(i)*dp(i)
        end do
        vof = global_sum(vof)
        dp_vof = global_sum(dp_vof)
        dp(1:m%ncell_onP) = dp(1:m%ncell_onP) - dp_vof/vof
      endif
      p_cc = p_cc + dp
      call gather_boundary(m%cell_ip, p_cc)
    end associate
  end subroutine pressure_cc_correct


  subroutine velocity_cc_correct(this, dt, props, p_cc, vel_cc)
    class(flow_projection), intent(inout) :: this
    real(r8), intent(in) :: dt
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: p_cc(:)
    real(r8), intent(inout) :: vel_cc(:,:)
    !-
    integer :: i, j, k, f0, f1

    call gradient_cf(this%grad_fc, p_cc, time, this%bc%p_neumann, &
        this%bc%p_dirichlet, props%inactive_f, 0.0_rp, this%gravity_head)

    associate (m => this%mesh%mesh, gfc => this%grad_fc, rho_fc => props%rho_fc,
      gpn => this%grad_p_rho_cc_n, gp => this%grad_p_rho_cc)

      do j = 1, m%nface_onP
        if (rho_fc(j) > 0.0_r8) then
          gfc(:,j) = gfc(:,j)/rho_fc(j)
        else
          gfc(:,j) = 0.0_r8
        end if
      end do
      call gather_boundary(m%face_ip, gfc)
      call interpolate_fc(gp, gfc)

      do i = 1, m%ncell_onP
        ! is this check actually necessary.. would gpn or gp be non-zero
        ! in either of these cases?
        if (props%rho_cc(i) == 0.0_r8 .or. props%inactive_c(i)) then
          vel_cc(:,i) = 0.0_r8
        else
          vel_cc(:,i) = vel_cc(:,i) + dt*(gpn(:,i)-gp(:,i))
        end if
      end do

    end associate
  end subroutine velocity_cc_correct



end module flow_projection_type
