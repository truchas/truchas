#include "f90_assert.fpp"
#define ASDF 0
#define Q 771
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
  use flow_operators
  use flow_mesh_type
  use flow_props_type
  use flow_bc_type
  use turbulence_models
  use hypre_hybrid_type
  use pcsr_matrix_type
  !use pgslib_module ! barriers for debugging
  implicit none
  private

  public :: flow_prediction, read_flow_predictor_namelist

  type :: flow_prediction
    type(flow_mesh), pointer :: mesh => null()
    type(flow_bc), pointer :: bc => null()
    class(turbulence_model), allocatable :: turb
    !class(porous_drag_model), allocatable :: drag
    type(hypre_hybrid) :: solver(3)
    type(parameter_list), pointer :: p => null()
    logical :: inviscid
    logical :: stokes
    real(r8), allocatable :: rhs(:,:), sol(:), rhs1d(:)
    real(r8), allocatable :: grad_fc_vector(:,:,:)
    real(r8) :: viscous_implicitness
    real(r8) :: solidify_implicitness
  contains
    procedure :: read_params
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

  subroutine read_flow_predictor_namelist(lun, p)
    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use flow_input_utils

    integer, intent(in) :: lun
    type(parameter_list), pointer, intent(inout) :: p
    type(parameter_list), pointer :: pp
    integer :: ios
    logical :: found
    character(128) :: iom

    !- hypre related items
    real(r8) :: rel_tol, abs_tol, conv_rate_tol, amg_strong_threshold
    integer :: max_ds_iter, max_amg_iter, gmres_krylov_dim, cg_use_two_norm, logging_level
    integer :: print_level, amg_max_levels, amg_coarsen_type, amg_coarsen_sweeps
    integer :: amg_smoothing_method, amg_smoothing_sweeps, amg_interp_method
    character(128) :: krylov_method, amg_coarsen_method


    namelist /flow_predictor/ rel_tol, &
        abs_tol, conv_rate_tol, max_ds_iter, max_amg_iter, &
        krylov_method, gmres_krylov_dim, cg_use_two_norm, &
        logging_level, print_level, amg_strong_threshold, &
        amg_max_levels, amg_coarsen_method, amg_coarsen_type, &
        amg_smoothing_sweeps, amg_smoothing_method, amg_interp_method

    pp => p%sublist("predictor")
    max_ds_iter = NULL_I
    max_amg_iter = NULL_I
    gmres_krylov_dim = NULL_I
    cg_use_two_norm = NULL_I
    logging_level = NULL_I
    print_level = NULL_I
    amg_max_levels = NULL_I
    amg_coarsen_type = NULL_I
    amg_coarsen_sweeps = NULL_I
    amg_smoothing_method = NULL_I
    amg_smoothing_sweeps = NULL_I
    amg_interp_method = NULL_I
    rel_tol = NULL_R
    abs_tol = NULL_R
    conv_rate_tol = NULL_R
    amg_strong_threshold = NULL_R
    krylov_method = NULL_C
    amg_coarsen_method = NULL_C

    if (is_IOP) then
      rewind lun
      call seek_to_namelist(lun, 'FLOW_PREDICTOR', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('Error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (found) then
      call TLS_info('')
      call TLS_info('Reading FLOW_PREDICTOR namelist ...')
      !! Read the namelist.
      if (is_IOP) then
        read(lun,nml=flow_predictor,iostat=ios,iomsg=iom)
      end if
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading FLOW_PREDICTOR namelist: ' // trim(iom))

      call broadcast(max_ds_iter)
      call broadcast(max_amg_iter)
      call broadcast(gmres_krylov_dim)
      call broadcast(cg_use_two_norm)
      call broadcast(logging_level)
      call broadcast(print_level)
      call broadcast(amg_max_levels)
      call broadcast(amg_coarsen_type)
      call broadcast(amg_coarsen_sweeps)
      call broadcast(amg_smoothing_method)
      call broadcast(amg_interp_method)
      call broadcast(rel_tol)
      call broadcast(abs_tol)
      call broadcast(conv_rate_tol)
      call broadcast(amg_strong_threshold)
      call broadcast(krylov_method)
      call broadcast(amg_coarsen_method)

      pp => pp%sublist("solver")
      call plist_set_if(pp, "rel-tol", rel_tol)
      call plist_set_if(pp, 'abs-tol', abs_tol)
      call plist_set_if(pp, 'conv-rate-tol', conv_rate_tol)
      call plist_set_if(pp, 'max-ds-iter', max_ds_iter)
      call plist_set_if(pp, 'max-amg-iter', max_amg_iter)
      call plist_set_if(pp, 'krylov-method', krylov_method)
      call plist_set_if(pp, 'gmres-krylov-dim', gmres_krylov_dim)
      call plist_set_if(pp, 'cg-use-two-norm', cg_use_two_norm /= 0)
      call plist_set_if(pp, 'logging-level', logging_level)
      call plist_set_if(pp, 'print-level', print_level)
      call plist_set_if(pp, 'amg-strong-threshold', amg_strong_threshold)
      call plist_set_if(pp, 'amg-max-levels', amg_max_levels)
      call plist_set_if(pp, 'amg-coarsen-method', amg_coarsen_method)
      call plist_set_if(pp, 'amg-coarsen-type', amg_coarsen_type)
      call plist_set_if(pp, 'amg-smoothing-sweeps', amg_smoothing_sweeps)
      call plist_set_if(pp, 'amg-smoothing-method', amg_smoothing_method)
      call plist_set_if(pp, 'amg-interp-method', amg_interp_method)
    end if
  end subroutine read_flow_predictor_namelist

  subroutine read_params(this, p)
    class(flow_prediction) :: this
    type(parameter_list), pointer, intent(in) :: p
    !-
    character(:), allocatable :: flow_type
    integer :: stat
    type(parameter_list), pointer :: pp

    this%p => p

    pp => p%sublist("options")

    ! default to crank-nicolson
    call pp%get('viscous implicitness', this%viscous_implicitness, 0.5_r8)
    ! default to fully implicit
    call pp%get('solidfy implicitness', this%solidify_implicitness, 1.0_r8)

    call turbulence_models_read_params(this%turb, p, off=this%inviscid)

    !sub => p%sublist("porous-drag-model")
    !call porous_drag_models_read_params(this%drag, sub, off=this%inviscid)
  end subroutine read_params


  subroutine init(this, mesh, bc, inviscid, stokes)
    class(flow_prediction), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: mesh
    type(flow_bc), pointer, intent(in) :: bc
    logical, intent(in) :: inviscid, stokes
    !-
    type(graph_container) :: g(size(this%solver))
    type(matrix_container) :: A(size(this%solver))
    type(ip_desc), pointer :: row_ip
    type(parameter_list), pointer :: sub
    integer :: i, j, k

    this%mesh => mesh
    this%bc => bc
    this%inviscid = inviscid
    this%stokes = stokes

    associate (m => mesh%mesh)
      allocate(this%rhs(3,m%ncell))
      allocate(this%grad_fc_vector(3,3,m%nface))
      allocate(this%sol(m%ncell))
      allocate(this%rhs1d(m%ncell))

      if (this%viscous_implicitness > 0.0_r8 .and. .not.this%inviscid) then
        ! build csr matrix graph for momentum solve, all components can reuse
        ! the same graph/matrix.
        row_ip => m%cell_ip
        do k = 1, size(g)
          allocate(g(k)%g)
          call g(k)%g%init(row_ip)
        end do

        do j = 1, m%ncell_onP
          do k = 1, size(g)
            call g(k)%g%add_edge(j,j)
          end do

          associate (cn => m%cnhbr(m%xcnhbr(j):m%xcnhbr(j+1)-1))
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

        sub => this%p%sublist("predictor")
        ASSERT(sub%count() > 0)
        sub => sub%sublist("solver")
        ASSERT(sub%count() > 0)
        !sub => this%p%sublist("predictor")%sublist("solver")
        do k = 1, size(A)
          allocate(A(k)%A)
          call A(k)%A%init(g(k)%g, take_graph=.true.)
          call this%solver(k)%init(A(k)%A, sub)
        end do
      end if
    end associate
  end subroutine init

  subroutine setup(this, dt, props, vel_cc, initial)
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: dt, vel_cc(:,:)
    type(flow_props), intent(inout) :: props
    logical, optional, intent(in) :: initial

!!$    if (present(initial)) then
!!$      if (initial) return
!!$    end if

    this%rhs = 0.0_r8

    ! account for turbulence model
    call this%turb%setup(vel_cc)
    call this%turb%apply(props)
    ! the turbulence model can modify mu_cc so only update
    ! the face-centered material properties after this has occurred
    call props%update_fc()

    call this%setup_solver(dt, props)
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

    associate (m => this%mesh%mesh, mu_f => props%mu_fc, rho => props%rho_cc, &
        vof => props%vof)

      do j = 1, m%ncell_onP
        associate (nhbr => m%cnhbr(m%xcnhbr(j):m%xcnhbr(j+1)-1), &
            face => m%cface(m%xcface(j):m%xcface(j+1)-1))

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
                coeff = dtv*dot_product(ds(:,fi), m%normal(:,fi))*mu_f(fi)/m%volume(j)
                do k = 1, size(A)
                  call A(k)%A%add_to(j, j, coeff)
                  call A(k)%A%add_to(j, ni, -coeff)
                end do
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
          j = this%mesh%fcell(2,fi) ! cell index
          if (j > m%ncell_onP) cycle
          ASSERT(this%mesh%fcell(1,fi) == 0)
          ASSERT(j > 0 .and. j <= m%ncell_onP)

          coeff = dtv*dot_product(ds(:,fi), m%normal(:,fi))*mu_f(fi)/m%volume(j)
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
          j = this%mesh%fcell(2,fi) ! cell index
          ASSERT(this%mesh%fcell(1,fi) == 0)

          coeff = dtv*dot_product(ds(:,fi), m%normal(:,fi))*mu_f(fi)/m%volume(j)
          do k = 1, size(A)
            call A(k)%A%add_to(j, j, coeff*(m%normal(k,fi)/m%area(fi))**2)
          end do
        end do
      end associate
      ! zero gradient of neumann bcs already handled


      ! XXX fixup for solid/void cells
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

    associate (m => this%mesh%mesh)
      do j = 1, m%ncell_onP
        this%rhs(:,j) = this%rhs(:,j) - dt*grad_p_rho_cc(:,j)*props%rho_cc_n(j)
      end do
    end associate

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

    associate (m => this%mesh%mesh, g => this%grad_fc_vector, mu => props%mu_fc)
      !write(*,'("Mu_face(:,33):", 6es15.5)') mu(m%cface(m%xcface(1089):m%xcface(1090)-1))

      ! velocity gradient tensor at faces such that
      ! g(:,1,j) = grad_at_face_j(u_1)
      ! g(:,2,j) = grad_at_face_j(u_2)
      ! g(:,3,j) = grad_at_face_j(u_3)
      call gradient_cf(g, vel_cc, this%bc%v_zero_normal, this%bc%v_dirichlet)
      call gather_boundary(m%face_ip, g)

#if ASDF
      associate( faces => m%cface(m%xcface(Q):m%xcface(Q+1)-1))
        write(*,'("grad(u) at [",i4,"] y- ",3es20.12)') faces(1), g(:,1,faces(1))
        write(*,'("grad(u) at [",i4,"] y+ ",3es20.12)') faces(3), g(:,1,faces(3))
        write(*,'("grad(u) at [",i4,"] x- ",3es20.12)') faces(5), g(:,1,faces(5))
        write(*,'("grad(u) at [",i4,"] x+ ",3es20.12)') faces(6), g(:,1,faces(6))
        !write(*,'(a,3es20.12)') 'dx at y+: ', this%mesh%face_centroid(:,faces(3)) - this%mesh%cell_centroid(:,33)
      end associate
#endif
      do j = 1, m%nface
        associate (n1 => this%mesh%fcell(1,j), n2 => this%mesh%fcell(2,j), rhs => this%rhs)
          if (n1 > 0) then ! inward oriented normal
            do i = 1, 3
              rhs(i,n1) = rhs(i,n1) - &
                  dtv*mu(j)*dot_product(m%normal(:,j),g(:,i,j))/m%volume(n1)
            end do
          end if
          if (n2 > 0) then ! outward oriented normal
            do i = 1, 3
              rhs(i,n2) = rhs(i,n2) + &
                  dtv*mu(j)*dot_product(m%normal(:,j),g(:,i,j))/m%volume(n2)
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
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: flux_volumes(:,:), vel_cc(:,:), vel_fn(:)
    type(flow_props), intent(in) :: props
    !
    integer :: i, j, f0, f1, nhbr, cn, k
    real(r8) :: v
#if ASDF
    real(r8) momentum_delta(3)
    momentum_delta = this%rhs(:,Q)
#endif

    ! we're going to assume here that the boundary conditions on flux_volumes
    ! and flux_vel have already been set..
    associate (m => this%mesh%mesh, rhs => this%rhs)
#if ASDF
      write(*, "('vel_fn[',i3,']: ',6es15.5): ") Q, vel_fn(m%cface(m%xcface(Q):m%xcface(Q+1)-1))
      write(*, "('flux_volumes[',i3,']: ',6es15.5): ") Q, flux_volumes(1,m%xcface(Q):m%xcface(Q+1)-1)
#endif
      do i = 1, m%ncell_onP
        f0 = m%xcface(i)
        f1 = m%xcface(i+1)-1
        nhbr = m%xcnhbr(i)

        do j = f0, f1

          k = m%cface(j)
          v = m%volume(i)

          ! donor cell upwinding
          if (any(flux_volumes(:,j) > 0.0_r8)) then
            rhs(:,i) = rhs(:,i) - vel_cc(:,i)*dot_product(flux_volumes(:,j),props%density(:))/v
#if ASDF
            if (i == Q) then
              write(*, "('face offset: ',i4)") j-f0+1
!!$              write(*,"('u_cc[',i4,']: ',es15.5, ' deltaX  : ', 2es15.5)") i, vel_cc(1,i), &
!!$                  -abs(vel_cc(1,i))*dot_product(flux_volumes(:,j),props%density(:))/v

!!$              write(*,"('v_cc[',i4,']: ',es15.5, ' deltaY  : ', 2es15.5)") i, vel_cc(2,i), &
!!$                  -abs(vel_cc(2,i))*dot_product(flux_volumes(:,j),props%density(:))/v
              write(*,"('v_cc[',i4,']: ',es15.5, ' deltaY  : ', 2es15.5)") i, vel_cc(2,i), &
                  -vel_cc(2,i)*dot_product(flux_volumes(:,j),props%density(:))/v
            end if
#endif
          else
            ! neighbor index
            cn = m%cnhbr(nhbr+(j-f0))
            if (cn > 0) then
              rhs(:,i) = rhs(:,i)-vel_cc(:,cn)*dot_product(flux_volumes(:,j),props%density(:))/v
#if ASDF
              if (i == Q) then
!!$                write(*,"('u_cc[',i4,']: ',es15.5, ' deltaX  : ', es15.5)") cn, vel_cc(1,cn), &
!!$                    -abs(vel_cc(1,cn))*dot_product(flux_volumes(:,j),props%density(:))/v
                write(*, "('face offset: ',i4)") j-f0+1
!!$                write(*,"('v_cc[',i4,']: ',es15.5, ' deltaY  : ', es15.5)") cn, vel_cc(2,cn), &
!!$                    -abs(vel_cc(2,cn))*dot_product(flux_volumes(:,j),props%density(:))/v
                write(*,"('v_cc[',i4,']: ',es15.5, ' deltaY  : ', es15.5)") cn, vel_cc(2,cn), &
                    -vel_cc(2,cn)*dot_product(flux_volumes(:,j),props%density(:))/v
              end if
#endif
            else
              ! this must be an inflow boundary so use face velocity.  This seems overly
              ! restrictive and does not account for 'angled' inflow
              rhs(:,i) = rhs(:,i) - vel_fn(k)*m%normal(:,k)/m%area(k) * &
                    dot_product(flux_volumes(:,j),props%density(:))/v
            end if
          end if
        end do
      end do
#if ASDF
      momentum_delta = this%rhs(:,Q)-momentum_delta
      write(*,'("Momentum_Delta(",i4,"):",3es15.5)') Q, Momentum_Delta(:)
#endif
    end associate
  end subroutine accumulate_rhs_momentum


  subroutine accumulate_rhs_solidified_rho(this, props, vel_cc)
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: vel_cc(:,:)
    type(flow_props), intent(in) :: props
    !
    integer :: j
    real(r8) :: w

    w = 1.0_r8-this%solidify_implicitness

    if (this%solidify_implicitness < 1.0_r8) then

      associate (m => this%mesh%mesh, rho => props%solidified_rho, rhs => this%rhs)
        do j = 1, m%ncell_onP
          rhs(:,j) = rhs(:,j) - w*rho(j)*vel_cc(:,j)
        end do
      end associate
    end if
  end subroutine accumulate_rhs_solidified_rho


  subroutine solve(this, dt, props, grad_p_rho_cc, flux_volumes, vel_fn, vel_cc, initial)
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: dt, flux_volumes(:,:), grad_p_rho_cc(:,:), vel_fn(:)
    real(r8), intent(inout) :: vel_cc(:,:)
    type(flow_props), intent(in) :: props
    logical, optional, intent(in) :: initial
    !
    integer :: i, j, ierr

!!$    if (present(initial)) then
!!$      if (initial) return
!!$    end if

    ! Pressure and viscous forces are unaware of solid/fluid boundaries.
    ! These forces are computed first and scaled by vof as a crude approximation.
    ! Once a better solid/fluid model is used, change this numerical hackjob.
    call this%accumulate_rhs_pressure(dt, props, grad_p_rho_cc)
#if ASDF
    write(*,"('post accumulate_rhs_pressure rhs[',i4,']:',3es20.12)") Q, this%rhs(:,Q)
#endif
    call this%accumulate_rhs_viscous_stress(dt, props, vel_cc)
#if ASDF
    write(*,"('post accumulate_rhs_viscous  rhs[',i4,']:',3es20.12)") Q, this%rhs(:,Q)
#endif
    associate (m => this%mesh%mesh, rhs => this%rhs)
      do j = 1, m%ncell_onP
        rhs(:,j) = rhs(:,j)*props%vof(j)
      end do
    end associate

    call this%accumulate_rhs_momentum(props, vel_cc, vel_fn, flux_volumes)
#if ASDF
    write(*,"('post accumulate_rhs_momentum rhs[',i4,']:',3es20.12)") Q, this%rhs(:,Q)
#endif
    call this%accumulate_rhs_solidified_rho(props, vel_cc)

    ! handle surface tension
    ! handle porous drag

    associate (m => this%mesh%mesh, rho => props%rho_cc_n, vof => props%vof_n, rhs => this%rhs)
      do j = 1, m%ncell_onP
        rhs(:,j) = rhs(:,j) + rho(j)*vof(j)*vel_cc(:,j)
      end do
    end associate

    ! skip all the wild limiter hackery for now
    associate (m => this%mesh%mesh, rhs1d => this%rhs1d, rhs => this%rhs, sol => this%sol)
      do i = 1, 3
        ! per-component solve for predictor velocity
        do j = 1, m%ncell_onP
          sol(j) = vel_cc(i,j)
          if (props%inactive_c(j) > 0 .or. props%vof_novoid(j) <= props%cutoff) then
            rhs1d(j) = 0.0_r8
          else
            rhs1d(j) = rhs(i,j) / props%vof(j)
          end if
        end do

        if (associated(this%solver(i)%matrix())) then
          call this%solver(i)%solve(rhs1d, sol, ierr)
          if (ierr /= 0) call tls_error("prediction solve unsuccessful")
          call gather_boundary(m%cell_ip, sol)
          do j = 1, m%ncell
            vel_cc(i,j) = sol(j)
          end do
        else
          ! WARN: Need to replace this with a mass matrix LHS to handle implicit
          !       terms for solidification and porosity.
          call gather_boundary(m%cell_ip, rhs1d)
          do j = 1, m%ncell
            if (props%inactive_c(j) > 0 .or. props%vof_novoid(j) <= props%cutoff) then
              vel_cc(i,j) = rhs1d(j)
            else
              vel_cc(i,j) = rhs1d(j) / props%rho_cc(j)
            end if
          end do
        end if

      end do

    end associate
#if ASDF
    write(*,"('after prediction u[',i4,']:',3es20.12)") Q, vel_cc(:, Q)
#endif
  end subroutine solve

  subroutine accept(this)
    class(flow_prediction), intent(inout) :: this
    call this%turb%accept()
  end subroutine accept
end module flow_prediction_type
