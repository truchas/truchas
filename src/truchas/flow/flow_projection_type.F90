#include "f90_assert.fpp"
#define ASDF 0
#define Q 42
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
  use flow_domain_types
  use unstr_mesh_type
  use flow_props_type
  use flow_bc_type
  use flow_operators
  use fischer_guess_type
  use index_partitioning
  use parallel_communication
  use hypre_hybrid_type
  use pcsr_matrix_type
  use parameter_list_type
  implicit none
  private

  type, public :: flow_projection
    private
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(flow_bc), pointer :: bc => null() ! unowned reference
    type(fischer_guess) :: fg
    type(hypre_hybrid) :: solver
    real(r8), allocatable :: rhs(:)
    real(r8), allocatable :: grad_fc(:,:) ! face centered gradient
    real(r8), allocatable :: grad_p_rho_cc(:,:)
    real(r8), allocatable :: grad_p_rho_fc(:,:)
    real(r8), allocatable :: delta_p_cc(:)
    real(r8), allocatable :: gravity_head(:,:)
    real(r8), allocatable :: weights_cf(:,:) ! weights for f->c interpolation
    real(r8), allocatable :: vel_cc_star(:,:)
  contains
    procedure :: init
    procedure :: setup
    procedure :: solve
    procedure :: accept
    procedure, private :: setup_face_velocity
    procedure, private :: setup_gravity
    procedure, private :: setup_solver
    procedure :: grad_p_rho
    procedure, private :: velocity_fc_correct
    procedure, private :: pressure_cc_correct
    procedure, private :: velocity_cc_correct
  end type flow_projection

contains

  subroutine init(this, mesh, bc, params)

    use parameter_list_type

    class(flow_projection), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(flow_bc), pointer, intent(in) :: bc
    type(parameter_list), intent(inout) :: params

    integer :: j, i
    type(pcsr_graph), pointer :: g
    type(pcsr_matrix), pointer :: A
    type(ip_desc), pointer :: row_ip
    type(parameter_list), pointer :: plist

    this%mesh => mesh
    this%bc => bc

    allocate(this%grad_fc(3,mesh%nface))
    allocate(this%delta_p_cc(mesh%ncell))
    allocate(this%gravity_head(2,mesh%nface))
    allocate(this%vel_cc_star(3,mesh%ncell))
    allocate(this%weights_cf(2,mesh%nface))
    allocate(this%grad_p_rho_cc(3,mesh%ncell))
    allocate(this%grad_p_rho_fc(3,mesh%nface))


    this%grad_fc = 0.0_r8
    this%delta_p_cc = 0.0_r8
    this%gravity_head = 0.0_r8
    this%vel_cc_star = 0.0_r8
    this%grad_p_rho_cc = 0.0_r8
    this%grad_p_rho_fc = 0.0_r8

    do j = 1, mesh%nface_onP
      associate (n => mesh%fcell(:,j), fc => mesh%face_centroid(:,j), &
          cc => mesh%cell_centroid, w => this%weights_cf(:,j))
        if (n(2) > 0) then
          w(1) = sum((fc(:)-cc(:,n(2)))**2)
          w(2) = sum((fc(:)-cc(:,n(1)))**2)
          w(:) = w(:) / sum(w(:))
        else
          w(:) = [0.0_r8, 1.0_r8]
        end if
      end associate
    end do
    allocate(this%rhs(mesh%ncell))

    !! Create a CSR matrix graph for the pressure poisson system.
    allocate(g)
    row_ip => mesh%cell_ip
    call g%init (row_ip)
    do j = 1, mesh%ncell_onP
      call g%add_edge(j,j)
      associate (cn => mesh%cnhbr(mesh%xcnhbr(j):mesh%xcnhbr(j+1)-1))
        do i = 1, size(cn)
          if (cn(i) > 0) call g%add_edge(j, cn(i))
        end do
      end associate
    end do
    call g%add_complete

    allocate(A)
    call A%init(g, take_graph=.true.)
    plist => params%sublist("corrector")
    ASSERT(plist%count() > 0)
    plist => plist%sublist("solver")
    ASSERT(plist%count() > 0)
    call this%solver%init(A, plist)
    call this%fg%init(mesh, params)
  end subroutine init

  subroutine setup(this, dt, props, grad_p_rho_cc_n, body_force, vel_cc, P_cc, vel_fn, initial)
    class(flow_projection), intent(inout) :: this
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: vel_cc(:,:), dt, P_cc(:), body_force(:), grad_p_rho_cc_n(:,:)
    real(r8), intent(out) :: vel_fn(:)
    logical, optional, intent(in) :: initial
    logical :: init

    init = .false.
    if (present(initial)) init = initial

    call start_timer("setup")
    call this%setup_gravity(props, body_force)
    call this%grad_p_rho(props, P_cc)
    call this%setup_face_velocity(dt, props, grad_p_rho_cc_n, vel_cc, vel_fn)
    ! we pass in dt == 0 to setup initial pass, this causes division by zero
    ! in rhs, so pass in 1.
!!$    if (init) then
!!$      call this%setup_solver(1.0_r8, props, vel_fn)
!!$    else
      call this%setup_solver(dt, props, vel_fn)
    !end if
    call stop_timer("setup")
  end subroutine setup

  subroutine accept(this, grad_p_rho_cc_n)
    class(flow_projection), intent(inout) :: this
    real(r8), intent(inout) :: grad_p_rho_cc_n(:,:)
    grad_p_rho_cc_n = this%grad_p_rho_cc
  end subroutine accept

  subroutine setup_solver(this, dt, props, vel)
    class(flow_projection), intent(inout) :: this
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: vel(:), dt
    !-
    type(pcsr_matrix), pointer :: A
    real(r8), pointer :: ds(:,:), row_val(:)
    integer, pointer :: row_idx(:)
    integer :: i, j, fi, ni
    real(r8) :: coeff
    logical :: pressure_pinned

    A => this%solver%matrix()
    call A%set_all(0.0_r8)

    ds => flow_gradient_coefficients()

    associate (m => this%mesh, rho_f => props%rho_fc, &
        cell_t => props%cell_t, face_t => props%face_t)
      do j = 1, m%ncell_onP
        associate (nhbr => m%cnhbr(m%xcnhbr(j):m%xcnhbr(j+1)-1), &
            face => m%cface(m%xcface(j):m%xcface(j+1)-1))

          if (cell_t(j) /= regular_t) then
            ! solve dummy equations in void/solid cells
            call A%add_to(j, j, 1.0_r8)
            cycle
          end if

          do i = 1, size(nhbr)
            fi = face(i)
            ni = nhbr(i)

            ! ni == 0 implies a domain boundary cell, handle these separately
            if (ni > 0) then
              ! note that normal is already weighted by face area
              coeff = dot_product(ds(:,fi), m%normal(:,fi))/rho_f(fi)/m%volume(j)
              if (face_t(fi) == regular_t) then
                call A%add_to(j, j, coeff)
                call A%add_to(j, ni, -coeff)
              else if (face_t(fi) == regular_void_t) then ! void acts as 0 dirichlet pressure
                call A%add_to(j, j, coeff)
              end if
              ! solid acts as a neumann condition - don't need to explicitly here
            end if
          end do
        end associate
      end do
    end associate


    associate (m => this%mesh, rhs => this%rhs)

      do j = 1, m%ncell_onP
        rhs(j) = 0.0_r8
        if (props%cell_t(j) /= regular_t) cycle

        associate( fn => m%cface(m%xcface(j):m%xcface(j+1)-1))
          do i = 1, size(fn)
            ! face velocity follows the convention of being positive in the face normal
            ! Boundary conditions have already been
            ! applied in computing `vel` so we shouldn't need to use them here again
            if (btest(m%cfpar(j),pos=i)) then ! inward face
              this%rhs(j) = this%rhs(j) + vel(fn(i))*m%area(fn(i))!/m%volume(j)
            else ! outward face
              this%rhs(j) = this%rhs(j) - vel(fn(i))*m%area(fn(i))!/m%volume(j)
            end if
          end do
        end associate
        this%rhs(j) = this%rhs(j) / (dt*m%volume(j))
      end do

      ! handle dirichlet bcs
      ! print *, "DEBUGGING PRESSURE BOUNDARY CONDITIONS"
      associate (faces => this%bc%dp_dirichlet%index, values => this%bc%dp_dirichlet%value, &
          rho_f => props%rho_fc)
        do i = 1, size(faces)
          fi = faces(i)
          if (fi > this%mesh%nface_onP) cycle
          if (props%face_t(fi) /= regular_t) cycle ! no bcs on non-regular faces

          j = this%mesh%fcell(1,fi) ! cell index
          ASSERT(this%mesh%fcell(2,fi) == 0)

          coeff = dot_product(ds(:,fi), m%normal(:,fi))/rho_f(fi)/m%volume(j)
          call A%add_to(j, j, coeff)
          this%rhs(j) = this%rhs(j) + coeff*values(i)

        end do
      end associate

      ! set the pressure at an arbitrary boundary point to 0 if all neumann
      if (.not.(this%bc%pressure_d .or. props%any_void) .and. this%bc%is_p_neumann_fix_PE) then
        pressure_pinned = .false.
        associate (faces => this%bc%p_neumann%index, rho_f => props%rho_fc)
          pressure_fix: do i = 1, size(faces)
            fi = faces(i)
            j = this%mesh%fcell(1,fi) ! cell index
            ASSERT(this%mesh%fcell(2,fi) == 0)

            if (props%face_t(fi) == regular_t) then
              coeff = dot_product(ds(:,fi), m%normal(:,fi))/rho_f(fi)/m%volume(j)
              call A%add_to(j, j, coeff)
              pressure_pinned = .true.
              ! no adjust to rhs requires for 0 dirichlet value
              exit pressure_fix
            end if
          end do pressure_fix
        end associate
        if (.not.pressure_pinned) call TLS_fatal("Cannot pin pressure.  Please fix algorithm.")
      end if
    end associate


    call this%solver%setup()

  end subroutine setup_solver


  subroutine setup_gravity(this, props, body_force)
    class(flow_projection), intent(inout) :: this
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: body_force(:)
    !-

    integer :: j, i, n
    real(r8) :: rho

    associate (g => this%gravity_head, &
        cc => this%mesh%cell_centroid, fc => this%mesh%face_centroid, p => props)

      g = 0.0_r8

      do j = 1, this%mesh%nface_onP
        do i = 1, 2
          n = this%mesh%fcell(i,j)
          if (n == 0) cycle

          rho = p%rho_cc(n) + p%rho_delta_cc(n)
          g(i,j) = -rho*dot_product(body_force, cc(:,n)-fc(:,j))
        end do
      end do

      call gather_boundary(this%mesh%face_ip, g(1,:))
      call gather_boundary(this%mesh%face_ip, g(2,:))
    end associate
    ! Truchas has a gravity_limiter... Do we need this?
  end subroutine setup_gravity

  subroutine setup_face_velocity(this, dt, props, grad_p_rho_cc_n, vel_cc, vel_fn)
    class(flow_projection), intent(inout) :: this
    real(r8), intent(in) :: dt
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: vel_cc(:,:), grad_p_rho_cc_n(:,:)
    real(r8), intent(out) :: vel_fn(:)
    !-
    integer :: i, j

    ! cell -> face interpolation
    vel_fn = 0.0_r8
    associate (v => this%vel_cc_star, gpn => grad_p_rho_cc_n, w => this%weights_cf)
      do i=1, this%mesh%ncell_onP
        v(:,i) = vel_cc(:,i) + dt*gpn(:,i)
      end do
      call gather_boundary(this%mesh%cell_ip, v)
      call interpolate_cf(vel_fn, v, w, this%bc%v_dirichlet, props%face_t, 0.0_r8)
    end associate
    ! subtract dynamic pressure grad
    associate (m => this%mesh, gp_fc => this%grad_p_rho_fc)
      ! enforced in grad_p_rho: this%grad_p_rho_fc == 0 where rho_fc == 0
      do j = 1, m%nface_onP
        vel_fn(j) = vel_fn(j) - dt*dot_product(m%normal(:,j), gp_fc(:,j))/m%area(j)
      end do
      ! should boundary conditions on vel_fn be enforced here again?
      associate (faces => this%bc%v_dirichlet%index, value => this%bc%v_dirichlet%value)
        do j = 1, size(faces)
          i = faces(j)
          vel_fn(i) = dot_product(m%normal(:,i), value(:,j))/m%area(i)
        end do
      end associate
      ! handle zero_normal bcs
      associate (faces => this%bc%v_zero_normal%index)
        do i = 1, size(faces)
          j = faces(i)
          ASSERT(this%mesh%fcell(2,j) == 0)
          vel_fn(j) = 0
        end do
      end associate
      call gather_boundary(m%face_ip, vel_fn)
    end associate
  end subroutine setup_face_velocity

  subroutine solve(this, dt, props, grad_p_rho_cc_n, vel_cc, p_cc, vel_fc, initial)
    class(flow_projection), intent(inout) :: this
    real(r8), intent(in) :: dt, grad_p_rho_cc_n(:,:)
    type(flow_props), intent(in) :: props
    real(r8), intent(inout) :: p_cc(:), vel_cc(:,:), vel_fc(:)
    logical, optional, intent(in) :: initial
    !-
    integer :: ierr, i, in

    ! walk through solve with initial dt = 0 to get divergence free pressure field
!!$    if (present(initial)) then
!!$      if (initial) return
!!$    end if

    call start_timer("solve")

    ! solve the pressure poisson system
    this%delta_p_cc = 0.0_r8
    call this%fg%guess(this%rhs, this%delta_p_cc)

#if ASDF
    write(*,'("PRESSURE RHS[", i3, "]: ",es15.5)') Q, this%rhs(Q)
#endif
    call start_timer("hypre solve")
    call this%solver%solve(this%rhs, this%delta_p_cc, ierr)
    call stop_timer("hypre solve")
    call tls_info('projection solve: ' // this%solver%metrics_string())
    if (ierr /= 0) call tls_error("projection solve unsuccessful")
    call this%fg%update(this%rhs, this%delta_p_cc, this%solver%matrix())
#if ASDF
    write(*,'("DP[",i3,"]: ",es15.5)') Q, this%delta_p_cc(Q)
#endif
    ! not sure if this call is necesary
    call gather_boundary(this%mesh%cell_ip, this%delta_p_cc)

    ! correct face and cell centered quantities
    call this%velocity_fc_correct(dt, props, vel_fc)
    call this%pressure_cc_correct(props, p_cc)
    call this%grad_p_rho(props, p_cc)
    call this%velocity_cc_correct(dt, props, grad_p_rho_cc_n, vel_cc)

    call stop_timer("solve")

  end subroutine solve


  subroutine velocity_fc_correct(this, dt, props, vel_fc)
    class(flow_projection), intent(inout) :: this
    real(r8), intent(in) :: dt
    type(flow_props), intent(in) :: props
    real(r8), intent(inout) :: vel_fc(:)
    !-
    integer :: j, i

    ! face gradient of solution
    call gradient_cf(this%grad_fc, this%delta_p_cc, this%bc%p_neumann, &
        this%bc%dp_dirichlet, props%face_t, 0.0_r8)

    ! divergence free face velocity
    do j = 1, this%mesh%nface_onP
      if (props%face_t(j) > regular_t) then
        vel_fc(j) = 0.0_r8
      else
        vel_fc(j) = vel_fc(j) - &
            dt*dot_product(this%grad_fc(:,j),this%mesh%normal(:,j)/this%mesh%area(j))/props%rho_fc(j)
      end if
    end do

    associate (index => this%bc%v_dirichlet%index, value => this%bc%v_dirichlet%value)
      do i = 1, size(index)
        j = index(i)
        vel_fc(j) = dot_product(value(:,i), this%mesh%normal(:,j))/this%mesh%area(j)
      end do
    end associate

    ! handle zero_normal bcs
    associate (faces => this%bc%v_zero_normal%index)
      do i = 1, size(faces)
        j = faces(i)
        ASSERT(this%mesh%fcell(2,j) == 0)
        vel_fc(j) = 0
      end do
    end associate

    call gather_boundary(this%mesh%face_ip, vel_fc)

    ! there may be some magic here about setting fluxing velocities on void cells.
    ! not sure if it's applicable when using face based (as opposed to side based) indexing
  end subroutine velocity_fc_correct


  subroutine pressure_cc_correct(this, props, p_cc)
    class(flow_projection), intent(inout) :: this
    type(flow_props), intent(in) :: props
    real(r8), intent(inout) :: p_cc(:)
    !-
    real(r8) :: dp_vof, vof, g_vof, g_dp_vof
    integer :: i

    associate (p => props, dp => this%delta_p_cc)
      p_cc = p_cc + dp
      call gather_boundary(this%mesh%cell_ip, p_cc)
    end associate
  end subroutine pressure_cc_correct


  ! compute face and cell centered dynamic_pressure_gradient/rho
  subroutine grad_p_rho(this, props, p_cc)
    class(flow_projection), intent(inout) :: this
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: p_cc(:)
    !-
    integer :: i, j

    associate (rho => props%rho_fc, &
        gp_cc => this%grad_p_rho_cc, gp_fc => this%grad_p_rho_fc)

      call gradient_cf(gp_fc, p_cc, this%bc%p_neumann, &
          this%bc%p_dirichlet, props%face_t, 0.0_r8, this%gravity_head)

      do j = 1, this%mesh%nface_onP
        if (props%face_t(j) <= regular_t) then
          gp_fc(:,j) = gp_fc(:,j)/rho(j)
        else
          gp_fc(:,j) = 0.0_r8
        end if
      end do
      call gather_boundary(this%mesh%face_ip, gp_fc)
      call interpolate_fc(gp_cc, gp_fc, props%face_t, this%bc%p_neumann%index)
    end associate
  end subroutine grad_p_rho


  subroutine velocity_cc_correct(this, dt, props, grad_p_rho_cc_n, vel_cc)
    class(flow_projection), intent(inout) :: this
    real(r8), intent(in) :: dt
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: grad_p_rho_cc_n(:,:)
    real(r8), intent(inout) :: vel_cc(:,:)
    !-
    integer :: i, j

    ! assumes dynamic pressure gradients have already been computed
    do i = 1, this%mesh%ncell_onP
      ! -dt * grad(dP)/rho
      vel_cc(:,i) = vel_cc(:,i) - dt*(this%grad_p_rho_cc(:,i)-grad_p_rho_cc_n(:,i))
    end do
    call gather_boundary(this%mesh%cell_ip, vel_cc)

  end subroutine velocity_cc_correct

end module flow_projection_type
