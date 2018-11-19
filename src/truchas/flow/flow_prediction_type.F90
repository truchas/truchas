#include "f90_assert.fpp"

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

module flow_prediction_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use parameter_list_type
  use index_partitioning
  use flow_domain_types
  use flow_operators
  use unstr_mesh_type
  use flow_props_type
  use flow_bc_type
  use turbulence_models
  use hypre_hybrid_type
  use pcsr_matrix_type
  !use pgslib_module ! barriers for debugging
  implicit none
  private

  type, public :: flow_prediction
    type(unstr_mesh), pointer :: mesh => null()
    type(flow_bc), pointer :: bc => null()
    class(turbulence_model), allocatable :: turb
    !class(porous_drag_model), allocatable :: drag
    type(hypre_hybrid) :: solver(3)
    logical :: inviscid
    logical :: stokes
    real(r8), allocatable :: rhs(:,:), sol(:), rhs1d(:)
    real(r8), allocatable :: grad_fc_vector(:,:,:)
    real(r8) :: viscous_implicitness
    real(r8) :: solidify_implicitness
  contains
    procedure :: init
    procedure :: setup
    procedure :: setup_solver
    procedure :: solve
    procedure :: accumulate_rhs_pressure
    procedure :: accumulate_rhs_viscous_stress
    procedure :: accumulate_rhs_momentum
    procedure :: accumulate_rhs_solidified_rho
    procedure :: accept
  end type flow_prediction

  type graph_container
    type(pcsr_graph), pointer :: g
  end type graph_container

  type matrix_container
    type(pcsr_matrix), pointer :: A
  end type matrix_container


contains

  subroutine init(this, mesh, bc, inviscid, stokes, params)

    use parameter_list_type

    class(flow_prediction), intent(inout) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(flow_bc), intent(in), target :: bc
    logical, intent(in) :: inviscid, stokes
    type(parameter_list), intent(inout) :: params

    type(graph_container) :: g(size(this%solver))
    type(matrix_container) :: A(size(this%solver))
    type(ip_desc), pointer :: row_ip
    type(parameter_list), pointer :: plist
    integer :: i, j, k

    this%mesh => mesh
    this%bc => bc
    this%inviscid = inviscid
    this%stokes = stokes

    plist => params%sublist("options")
    call plist%get('viscous implicitness', this%viscous_implicitness, default=0.5_r8)
    call plist%get('solidfy implicitness', this%solidify_implicitness, default=1.0_r8)
    call alloc_turbulence_model(this%turb, plist, off=this%inviscid)

    allocate(this%rhs(3,mesh%ncell))
    allocate(this%grad_fc_vector(3,3,mesh%nface))
    allocate(this%sol(mesh%ncell))
    allocate(this%rhs1d(mesh%ncell))

    if (this%viscous_implicitness > 0.0_r8 .and. .not.this%inviscid) then
      ! build csr matrix graph for momentum solve, all components can reuse
      ! the same graph/matrix.
      row_ip => mesh%cell_ip
      do k = 1, size(g)
        allocate(g(k)%g)
        call g(k)%g%init(row_ip)
      end do

      do j = 1, mesh%ncell_onP
        do k = 1, size(g)
          call g(k)%g%add_edge(j,j)
        end do

        associate (cn => mesh%cnhbr(mesh%xcnhbr(j):mesh%xcnhbr(j+1)-1))
          do i = 1, size(cn)
            if (cn(i) > 0) then
              do k = 1, size(g)
                call g(k)%g%add_edge(j, cn(i))
              end do
            end if
          end do
        end associate
      end do
      do k = 1, size(g)
        call g(k)%g%add_complete()
      end do

      plist => params%sublist("predictor")
      ASSERT(plist%count() > 0)
      plist => plist%sublist("solver")
      ASSERT(plist%count() > 0)
      do k = 1, size(A)
        allocate(A(k)%A)
        call A(k)%A%init(g(k)%g, take_graph=.true.)
        call this%solver(k)%init(A(k)%A, plist)
      end do
    end if

  end subroutine init

  subroutine setup(this, dt, props, vel_cc)
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: dt, vel_cc(:,:)
    type(flow_props), intent(inout) :: props

    call start_timer("setup")

    this%rhs = 0.0_r8

    ! account for turbulence model
    call this%turb%setup(vel_cc)
    call this%turb%apply(props)
    ! the turbulence model can modify mu_cc so only update
    ! the face-centered material properties after this has occurred
    call props%update_fc()

    call this%setup_solver(dt, props)

    call stop_timer("setup")

  end subroutine setup

  ! solver matrix => lhs of predictor step
  !  (rho - dt*theta_mu*stress_grad)(u*) = rhs
  subroutine setup_solver(this, dt, props)
    class(flow_prediction), intent(inout) :: this
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: dt
    !-
    real(r8), pointer :: ds(:,:)
    integer :: i, j, fi, ni, k
    real(r8) :: coeff, dtv
    type(matrix_container) :: A(size(this%solver))

    ! this should be changed to check all other terms which may be treated implicitly
    if (this%inviscid .or. this%viscous_implicitness == 0.0_r8) return

    do k = 1, size(A)
      A(k)%A => this%solver(k)%matrix()
      call A(k)%A%set_all(0.0_r8)
    end do


    ds => flow_gradient_coefficients()
    dtv = dt*this%viscous_implicitness

    associate (mu_f => props%mu_fc, rho => props%rho_cc, &
        vof => props%vof, cell_t => props%cell_t, face_t => props%face_t)

      do j = 1, this%mesh%ncell_onP
        associate (nhbr => this%mesh%cnhbr(this%mesh%xcnhbr(j):this%mesh%xcnhbr(j+1)-1), &
            face => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))

          if (cell_t(j) /= regular_t) then
            ! solve dummy equations on solid/void
            do k = 1, size(A)
              call A(k)%A%add_to(j,j,1.0_r8)
            end do
            cycle
          end if

          ! inertial terms
          do k = 1, size(A)
            call A(k)%A%add_to(j,j,rho(j))
          end do

          ! viscous terms
          if (.not.this%inviscid .and. this%viscous_implicitness > 0.0_r8) then
            do i = 1, size(nhbr)
              fi = face(i)
              ni = nhbr(i)

              ! ni == 0 implies a domain boundary cell, handle these separately
              if (ni > 0) then
                ! note that normal is already weighted by face area
                coeff = dtv*dot_product(ds(:,fi), this%mesh%normal(:,fi))*mu_f(fi)/this%mesh%volume(j)
                if (face_t(fi) == regular_t) then
                  do k = 1, size(A)
                    call A(k)%A%add_to(j, j, coeff)
                    call A(k)%A%add_to(j, ni, -coeff)
                  end do
                else if (face_t(fi) == solid_t) then
                  do k = 1, size(A)
                    call A(k)%A%add_to(j, j, coeff)
                  end do
                  ! there is no adjustment to the rhs here since it is assumed that the
                  ! velocity in solid cells is zero
                end if
                ! there is no check for a void_t face since that acts as a neumann condition
                ! on velocity which does not contribute to the lhs
              end if
            end do
          end if

          ! porous drag - skip for now

          ! solidification
          if (vof(j) > 0.0_r8) then
            coeff = this%solidify_implicitness*props%solidified_rho(j)/vof(j)
            do k = 1, size(A)
              call A(k)%A%add_to(j, j, coeff)
            end do
          end if
        end associate
      end do

      ! handle dirichlet bcs for viscous terms
      associate (faces => this%bc%v_dirichlet%index, value => this%bc%v_dirichlet%value)
        do i = 1, size(faces)
          fi = faces(i)
          j = this%mesh%fcell(1,fi) ! cell index
          if (j > this%mesh%ncell_onP) cycle ! this cell will be properly handled by a different pe
          if (face_t(fi) /= regular_t) cycle ! no bcs on non-regular faces

          ASSERT(this%mesh%fcell(2,fi) == 0)

          coeff = dtv*dot_product(ds(:,fi), this%mesh%normal(:,fi))*mu_f(fi)/this%mesh%volume(j)
          do k = 1, size(A)
            call A(k)%A%add_to(j, j, coeff)
          end do
          this%rhs(:,j) = this%rhs(:,j) + coeff*value(:,i)
        end do
      end associate

      ! handle zero_normal bcs
      associate (faces => this%bc%v_zero_normal%index)
        do i = 1, size(faces)
          fi = faces(i)
          j = this%mesh%fcell(1,fi) ! cell index
          if (j > this%mesh%ncell_onP) cycle ! this cell will be properly handled by a different pe
          if (face_t(fi) /= regular_t) cycle ! no bcs on non-regular faces

          ASSERT(this%mesh%fcell(2,fi) == 0)

          coeff = dtv*dot_product(ds(:,fi), this%mesh%normal(:,fi))*mu_f(fi)/this%mesh%volume(j)
          do k = 1, size(A)
            call A(k)%A%add_to(j, j, coeff*(this%mesh%normal(k,fi)/this%mesh%area(fi))**2)
          end do
        end do
      end associate
      ! zero gradient of neumann bcs already handled
    end associate

    ! copies data to hypre -> call this last!!!!!
    do k = 1, size(this%solver)
      call this%solver(k)%setup()
    end do

  end subroutine setup_solver

  ! accumulate [-dt * grad(P)^n]
  ! the full contribution to the velocity update is [-dt * grad(P)^n]/rho^n+1
  ! The divisor is handled in the solve routine
  subroutine accumulate_rhs_pressure(this, dt, props, grad_p_rho_cc)
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: dt, grad_p_rho_cc(:,:)
    type(flow_props), intent(in) :: props
    integer :: j

    do j = 1, this%mesh%ncell_onP
      this%rhs(:,j) = this%rhs(:,j) - dt*grad_p_rho_cc(:,j)*props%rho_cc(j)
    end do

  end subroutine accumulate_rhs_pressure

  ! accumulate (1-theta_mu)[ dt * div(mu_f*grad(u_i))^n]
  ! the full contribution to the velocity update is (1-theta_mu)[dt * div(mu_f*grad(u_i))^n]/rho^n+1
  ! The divisor is handled in the solve routine
  subroutine accumulate_rhs_viscous_stress(this, dt, props, vel_cc)
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: dt, vel_cc(:,:)
    type(flow_props), intent(in) :: props
    !
    integer :: i, j
    real(r8) :: dtv

    if (this%inviscid .or. this%viscous_implicitness == 1.0_r8) return

    dtv = dt*(1.0_r8 - this%viscous_implicitness)

    associate (g => this%grad_fc_vector, mu => props%mu_fc)

      ! velocity gradient tensor at faces such that
      ! g(:,1,j) = grad_at_face_j(u_1)
      ! g(:,2,j) = grad_at_face_j(u_2)
      ! g(:,3,j) = grad_at_face_j(u_3)
      call gradient_cf(g, vel_cc, this%bc%v_zero_normal, this%bc%v_dirichlet)
      call gather_boundary(this%mesh%face_ip, g)

      do j = 1, this%mesh%nface
        associate (n1 => this%mesh%fcell(1,j), n2 => this%mesh%fcell(2,j), rhs => this%rhs)
          if (n2 > 0) then ! inward oriented normal
            do i = 1, 3
              rhs(i,n2) = rhs(i,n2) - &
                  dtv*mu(j)*dot_product(this%mesh%normal(:,j),g(:,i,j))/this%mesh%volume(n2)
            end do
          end if
          if (n1 > 0) then ! outward oriented normal
            do i = 1, 3
              rhs(i,n1) = rhs(i,n1) + &
                  dtv*mu(j)*dot_product(this%mesh%normal(:,j),g(:,i,j))/this%mesh%volume(n1)
            end do
          end if
        end associate
      end do
    end associate
  end subroutine accumulate_rhs_viscous_stress

  ! accumulate [-dt*div(rho*u*u)]
  ! flux_volumes is in ragged format where each cell has a potentially unique flux_value
  ! for a given cell face
  ! The flux volume for a given material is dt*A*vel_fn.  A positive value indicates outward flux
  ! Momentum is transported via flux-volumes using a first order upwind scheme
  subroutine accumulate_rhs_momentum(this, props, vel_cc, vel_fn, flux_volumes)

    use f08_intrinsics, only: findloc

    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: flux_volumes(:,:), vel_cc(:,:), vel_fn(:)
    type(flow_props), intent(in) :: props

    integer :: i, j, k, f, f0, f1, nhbr, cn, nfluid
    real(r8) :: v

    ! This is the number of real (non-void) fluids in our system.
    ! Void does not contribute to momentum transport.
    nfluid = size(props%density)

    ! we're going to assume here that the boundary conditions on flux_volumes
    ! and flux_vel have already been set..
    do i = 1, this%mesh%ncell_onP
      f0 = this%mesh%xcface(i)
      f1 = this%mesh%xcface(i+1)-1
      nhbr = this%mesh%xcnhbr(i)
      v = this%mesh%volume(i)

      do j = f0, f1
        k = this%mesh%cface(j)

        ! donor cell upwinding
        if (any(flux_volumes(:,j) > 0.0_r8)) then
          this%rhs(:,i) = this%rhs(:,i) - &
              vel_cc(:,i)*dot_product(flux_volumes(:nfluid,j),props%density(:))/v
        else
          cn = this%mesh%cnhbr(nhbr+(j-f0)) ! neighbor index
          if (cn > 0) &
            this%rhs(:,i) = this%rhs(:,i) - &
                vel_cc(:,cn)*dot_product(flux_volumes(:nfluid,j),props%density(:))/v
          ! inflow boundaries handled in subsequent loops
        end if
      end do
    end do

    ! On dirichlet velocity boundaries, use the boundary condition velocity
    associate (faces => this%bc%v_dirichlet%index, value => this%bc%v_dirichlet%value)
      do k = 1, size(faces)
        f = faces(k)
        i = this%mesh%fcell(1,f) ! cell index
        if (i > this%mesh%ncell_onP) cycle ! this cell will be properly handled by a different pe
        ASSERT(this%mesh%fcell(2,f) == 0)

        ! get the local face id on this cell
        j = this%mesh%xcface(i) - 1 + &
            findloc(this%mesh%cface(this%mesh%xcface(i):this%mesh%xcface(i+1)-1), f)

        ! outflow handled in above loop
        if (.not.any(flux_volumes(:,j) < 0)) cycle

        this%rhs(:,i) = this%rhs(:,i) - &
            value(:,k)*dot_product(flux_volumes(:nfluid,j),props%density) / this%mesh%volume(i)
      end do
    end associate

    ! On pressure dirichlet boundaries we have an implied velocity
    ! neumann condition. The following loop relies on this fact.
    ! On velocity neumann boundaries, assume the upwind velocity is
    ! the same as the local cell velocity. We need tangential components
    ! of the velocity, so the face-based scalar array is insufficient.
    associate (faces => this%bc%p_dirichlet%index)
      do k = 1, size(faces)
        f = faces(k)
        i = this%mesh%fcell(1,f) ! cell index
        if (i > this%mesh%ncell_onP) cycle ! this cell will be properly handled by a different pe
        ASSERT(this%mesh%fcell(2,f) == 0)

        ! get the local face id on this cell
        j = this%mesh%xcface(i) - 1 + &
            findloc(this%mesh%cface(this%mesh%xcface(i):this%mesh%xcface(i+1)-1), f)

        ! outflow handled in above loop
        if (.not.any(flux_volumes(:,j) < 0)) cycle

        this%rhs(:,i) = this%rhs(:,i) - &
            vel_cc(:,i)*dot_product(flux_volumes(:nfluid,j),props%density) / this%mesh%volume(i)
      end do
    end associate

    ! On slip boundaries, the flux volumes are zero, so no momentum flux.

  end subroutine accumulate_rhs_momentum


  subroutine accumulate_rhs_solidified_rho(this, solidified_rho, vel_cc)

    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: solidified_rho(:), vel_cc(:,:)

    integer :: j
    real(r8) :: w

    if (this%solidify_implicitness == 1) return

    w = 1 - this%solidify_implicitness
    do j = 1, this%mesh%ncell_onP
      this%rhs(:,j) = this%rhs(:,j) - w * solidified_rho(j) * vel_cc(:,j)
    end do

  end subroutine accumulate_rhs_solidified_rho


  subroutine solve(this, dt, props, grad_p_rho_cc, flux_volumes, vel_fn, vel_cc)
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: dt, flux_volumes(:,:), grad_p_rho_cc(:,:), vel_fn(:)
    real(r8), intent(inout) :: vel_cc(:,:)
    type(flow_props), intent(in) :: props
    !
    integer :: i, j, ierr
    real(r8) :: coeff
    character(18), parameter :: slabel(3) = 'predictor' // ['1', '2', '3'] // ' solve: '

    call start_timer("solve")
    ! The accumulate routines may add contributions to non-regular cells.  These
    ! contributions are ignored by the solver

    ! Pressure and viscous forces are unaware of solid/fluid boundaries.
    ! These forces are computed first and scaled by vof as a crude approximation.
    ! Once a better solid/fluid model is used, change this numerical hackjob.
    call this%accumulate_rhs_pressure(dt, props, grad_p_rho_cc)
    call this%accumulate_rhs_viscous_stress(dt, props, vel_cc)
    do j = 1, this%mesh%ncell_onP
      this%rhs(:,j) = this%rhs(:,j)*props%vof(j)
    end do

    call this%accumulate_rhs_momentum(props, vel_cc, vel_fn, flux_volumes)
    call this%accumulate_rhs_solidified_rho(props%solidified_rho, vel_cc)

    ! handle surface tension

    ! handle tangential surface tension BC
    associate (faces => this%bc%surface_tension%index, value => this%bc%surface_tension%value)
      do i = 1, size(faces)
        j = this%mesh%fcell(1,faces(i))
        if (j > this%mesh%ncell_onP) cycle
        this%rhs(:,j) = this%rhs(:,j) + dt * value(:,i)
      end do
    end associate

    ! handle porous drag

    associate (rho => props%rho_cc_n, vof => props%vof_n, rhs => this%rhs)
      do j = 1, this%mesh%ncell_onP
        rhs(:,j) = rhs(:,j) + rho(j)*vof(j)*vel_cc(:,j)
      end do
    end associate

    ! skip all the wild limiter hackery for now
    associate (rhs1d => this%rhs1d, rhs => this%rhs, &
        sol => this%sol, cell_t => props%cell_t)
      do i = 1, 3
        ! per-component solve for predictor velocity
        do j = 1, this%mesh%ncell_onP
          sol(j) = vel_cc(i,j)
          if (cell_t(j) /= regular_t) then
            rhs1d(j) = 0.0_r8
          else
            rhs1d(j) = rhs(i,j) / props%vof(j)
          end if
        end do

        if (associated(this%solver(i)%matrix())) then
          call start_timer("hypre solve")
          call this%solver(i)%solve(rhs1d, sol, ierr)
          call stop_timer("hypre solve")
          call tls_info(slabel(i) // this%solver(i)%metrics_string())
          if (ierr /= 0) call tls_error("prediction solve unsuccessful")
          call gather_boundary(this%mesh%cell_ip, sol)
          do j = 1, this%mesh%ncell
            vel_cc(i,j) = sol(j)
          end do
        else
          call gather_boundary(this%mesh%cell_ip, rhs1d)
          do j = 1, this%mesh%ncell
            if (cell_t(j) /= regular_t) then
              vel_cc(i,j) = rhs1d(j)
            else
              coeff = props%rho_cc(j)
              if (props%vof(j) > 0) &
                  coeff = coeff + &
                  this%solidify_implicitness * props%solidified_rho(j) / props%vof(j)
              if (coeff /= 0) then
                vel_cc(i,j) = rhs1d(j) / coeff
              else
                vel_cc(i,j) = 0
              end if
            end if
          end do
        end if

      end do

    end associate

    call stop_timer("solve")

  end subroutine solve

  subroutine accept(this)
    class(flow_prediction), intent(inout) :: this
    call this%turb%accept()
  end subroutine accept
end module flow_prediction_type
