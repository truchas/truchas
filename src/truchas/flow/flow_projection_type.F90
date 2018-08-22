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


module flow_projection_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use flow_mesh_type
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

  public :: flow_projection, read_flow_corrector_namelist

  type :: flow_projection
    private
    type(flow_mesh), pointer :: mesh => null() ! unowned reference
    type(flow_bc), pointer :: bc => null() ! unowned reference
    type(fischer_guess) :: fg
    type(hypre_hybrid) :: solver
    type(parameter_list), pointer :: p => null()
    real(r8), allocatable :: rhs(:)
    real(r8), allocatable :: grad_fc(:,:) ! face centered gradient
    real(r8), allocatable :: grad_p_rho_cc(:,:)
    real(r8), allocatable :: grad_p_rho_fc(:,:)
    real(r8), allocatable :: delta_p_cc(:)
    real(r8), allocatable :: gravity_head(:,:)
    real(r8), allocatable :: weights_cf(:,:) ! weights for f->c interpolation
    real(r8), allocatable :: vel_cc_star(:,:)
  contains
    procedure :: read_params
    procedure :: init
    procedure :: setup
    procedure :: solve
    procedure :: accept
    procedure, private :: setup_face_velocity
    procedure, private :: setup_gravity
    procedure, private :: setup_solver
    procedure, private :: grad_p_rho
    procedure, private :: velocity_fc_correct
    procedure, private :: pressure_cc_correct
    procedure, private :: velocity_cc_correct
  end type flow_projection

contains

  subroutine read_flow_corrector_namelist(lun, p)
    use string_utilities, only: i_to_c
    use parallel_communication, only: is_IOP, broadcast
    use flow_input_utils

    integer, intent(in) :: lun
    type(parameter_list), pointer, intent(inout) :: p
    type(parameter_list), pointer :: pp
    integer :: ios
    logical :: found
    character(128) :: iom

    integer :: fischer_history
    !- hypre related items
    real(r8) :: rel_tol, abs_tol, conv_rate_tol, amg_strong_threshold
    integer :: max_ds_iter, max_amg_iter, gmres_krylov_dim, cg_use_two_norm, logging_level
    integer :: print_level, amg_max_levels, amg_coarsen_type, amg_coarsen_sweeps
    integer :: amg_smoothing_method, amg_smoothing_sweeps, amg_interp_method
    character(128) :: krylov_method, amg_coarsen_method


    namelist /flow_corrector/ fischer_history, rel_tol, &
        abs_tol, conv_rate_tol, max_ds_iter, max_amg_iter, &
        krylov_method, gmres_krylov_dim, cg_use_two_norm, &
        logging_level, print_level, amg_strong_threshold, &
        amg_max_levels, amg_coarsen_method, amg_coarsen_type, &
        amg_smoothing_sweeps, amg_smoothing_method, amg_interp_method

    pp => p%sublist("corrector")

    fischer_history = NULL_I
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
      call seek_to_namelist(lun, 'FLOW_CORRECTOR', found, iostat=ios)
    end if
    call broadcast(ios)
    if (ios /= 0) call TLS_fatal('Error reading input file: iostat=' // i_to_c(ios))

    call broadcast(found)
    if (found) then
      call TLS_info('')
      call TLS_info('Reading FLOW_CORRECTOR namelist ...')
      !! Read the namelist.
      if (is_IOP) then
        read(lun,nml=flow_corrector,iostat=ios,iomsg=iom)
      end if
      call broadcast(ios)
      if (ios /= 0) call TLS_fatal('error reading FLOW_CORRECTOR namelist: ' // trim(iom))

      call broadcast(fischer_history)
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


      call plist_set_if(pp, 'history', fischer_history)
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
  end subroutine read_flow_corrector_namelist


  subroutine read_params(this, p)
    class(flow_projection), intent(inout) :: this
    type(parameter_list), pointer, intent(in) :: p

    this%p => p
    call this%fg%read_params(p)
  end subroutine read_params

  subroutine init(this, mesh, bc)
    class(flow_projection), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: mesh
    type(flow_bc), pointer, intent(in) :: bc
    !-
    integer :: j, i
    type(pcsr_graph), pointer :: g
    type(pcsr_matrix), pointer :: A
    type(ip_desc), pointer :: row_ip
    type(unstr_mesh), pointer :: m
    type(parameter_list), pointer :: sub

    this%mesh => mesh
    this%bc => bc
    m => mesh%mesh

    allocate(this%grad_fc(3,m%nface))
    allocate(this%delta_p_cc(m%ncell))
    allocate(this%gravity_head(2,m%nface))
    allocate(this%vel_cc_star(3,m%ncell))
    allocate(this%weights_cf(2,m%nface))
    allocate(this%grad_p_rho_cc(3,m%ncell))
    allocate(this%grad_p_rho_fc(3,m%nface))


    this%grad_fc = 0.0_r8
    this%delta_p_cc = 0.0_r8
    this%gravity_head = 0.0_r8
    this%vel_cc_star = 0.0_r8
    this%grad_p_rho_cc = 0.0_r8
    this%grad_p_rho_fc = 0.0_r8

    do j = 1, m%nface_onP
      associate (n => mesh%fcell(:,j), fc => mesh%face_centroid(:,j), &
          cc => mesh%cell_centroid, w => this%weights_cf(:,j))
        if (n(1) > 0) then
          w(1) = sum((fc(:)-cc(:,n(1)))**2)
          w(2) = sum((fc(:)-cc(:,n(2)))**2)
          w(:) = w(:) / sum(w(:))
        else
          w(:) = [0.0_r8, 1.0_r8]
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
          if (cn(i) > 0) call g%add_edge(j, cn(i))
        end do
      end associate
    end do
    call g%add_complete

    allocate(A)
    call A%init(g, take_graph=.true.)
    sub => this%p%sublist("corrector")
    ASSERT(sub%count() > 0)
    sub => sub%sublist("solver")
    ASSERT(sub%count() > 0)
    call this%solver%init(A, sub)
    call this%fg%init(mesh)
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

    A => this%solver%matrix()
    call A%set_all(0.0_r8)

    ds => flow_gradient_coefficients()

    associate (m => this%mesh%mesh, rho_f => props%rho_fc)
      do j = 1, m%ncell_onP
        associate (nhbr => m%cnhbr(m%xcnhbr(j):m%xcnhbr(j+1)-1), &
            face => m%cface(m%xcface(j):m%xcface(j+1)-1))
          do i = 1, size(nhbr)
            fi = face(i)
            ni = nhbr(i)

            ! ni == 0 implies a domain boundary cell, handle these separately
            if (rho_f(fi) /= 0.0_r8 .and. ni > 0) then
              ! note that normal is already weighted by face area
              coeff = dot_product(ds(:,fi), m%normal(:,fi))/rho_f(fi)/m%volume(j)
              call A%add_to(j, j, coeff)
              call A%add_to(j, ni, -coeff)
            end if
          end do
        end associate
      end do
    end associate


    !! FIXME: need correction for void/solid... maybe
    associate (m => this%mesh%mesh, rhs => this%rhs)
#if ASDF
      associate( faces => m%cface(m%xcface(Q):m%xcface(Q)-1))
        write(*,'("vel_fn/normal at [",i4,"] y- ",2es20.12)') faces(1), vel(faces(1)), m%normal(2,faces(1))
        write(*,'("vel_fn/normal at [",i4,"] y+ ",2es20.12)') faces(3), vel(faces(3)), m%normal(2,faces(3))
        write(*,'("vel_fn/normal at [",i4,"] x- ",2es20.12)') faces(5), vel(faces(5)), m%normal(1,faces(5))
        write(*,'("vel_fn/normal at [",i4,"] x+ ",2es20.12)') faces(6), vel(faces(6)), m%normal(1,faces(6))
        write(*,'("faces areas: ",4es20.12)') m%area([faces(1),faces(3),faces(5),faces(6)])
      end associate
#endif

      do j = 1, m%ncell_onP
        rhs(j) = 0.0_r8
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
        !if (j == 771) write(*, '("rhs [",i4,"]: ",es20.12)') 771, this%rhs(771)
        this%rhs(j) = this%rhs(j) / (dt*m%volume(j))
      end do

      ! handle dirichlet bcs
      !print *, "DEBUGGING PRESSURE BOUNDARY CONDITIONS"
      associate (faces => this%bc%dp_dirichlet%index, values => this%bc%dp_dirichlet%value, &
          rho_f => props%rho_fc)
        do i = 1, size(faces)
          fi = faces(i)
          j = this%mesh%fcell(2,fi) ! cell index
          ASSERT(this%mesh%fcell(1,fi) == 0)
          if (rho_f(fi) /=  0.0_r8) then
            coeff = dot_product(ds(:,fi), m%normal(:,fi))/rho_f(fi)/this%mesh%mesh%volume(j)
            call A%add_to(j, j, coeff)
            this%rhs(j) = this%rhs(j) + coeff*values(i)
          end if
        end do
      end associate

      ! set the pressure at an arbitrary boundary point to 0 if all neumann
      if (.not.this%bc%pressure_d .and. is_IOP) then
        associate (faces => this%bc%p_neumann%index, rho_f => props%rho_fc)
          i = 1
          fi = faces(i)
          j = this%mesh%fcell(2,fi) ! cell index
          ASSERT(this%mesh%fcell(1,fi) == 0)

          if (rho_f(fi) /=  0.0_r8) then
            coeff = dot_product(ds(:,fi), m%normal(:,fi))/rho_f(fi)/this%mesh%mesh%volume(j)
            call A%add_to(j, j, coeff)
            ! no adjust to rhs requires for 0 dirichlet value
          else
            print *, "FIX THIS: setting ficticous pressure in solid"
          end if
        end associate
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

    associate( fcell => this%mesh%fcell, m => this%mesh%mesh, g => this%gravity_head, &
        cc => this%mesh%cell_centroid, fc => this%mesh%face_centroid, p => props)

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
    associate (m => this%mesh%mesh, v => this%vel_cc_star, gpn => grad_p_rho_cc_n, &
        w => this%weights_cf)
      do i=1, m%ncell_onP
        v(:,i) = vel_cc(:,i) + dt*gpn(:,i)
      end do
      call gather_boundary(m%cell_ip, v)
      call interpolate_cf(vel_fn, v, w, this%bc%v_dirichlet, props%inactive_f, 0.0_r8)
    end associate
    ! subtract dynamic pressure grad
    associate (m => this%mesh%mesh, gp_fc => this%grad_p_rho_fc)
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
      call gather_boundary(m%face_ip, vel_fn)
#if ASDF
      write(*,'("Vel_Fn[",i4,"]: ",6es16.8)') Q, vel_fn(m%cface(m%xcface(Q):m%xcface(Q+1)-1))
#endif
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

    ! solve the pressure poisson system
    this%delta_p_cc = 0.0_r8
    call this%fg%guess(this%rhs, this%delta_p_cc)

#if ASDF
    write(*,'("PRESSURE RHS[", i3, "]: ",es15.5)') Q, this%rhs(Q)
#endif
    call this%solver%solve(this%rhs, this%delta_p_cc, ierr)
    if (ierr /= 0) call tls_error("projection solve unsuccessful")
    call this%fg%update(this%rhs, this%delta_p_cc, this%solver%matrix())
#if ASDF
    write(*,'("DP[",i3,"]: ",es15.5)') Q, this%delta_p_cc(Q)
#endif
    ! not sure if this call is necesary
    call gather_boundary(this%mesh%mesh%cell_ip, this%delta_p_cc)

    ! correct face and cell centered quantities
    call this%velocity_fc_correct(dt, props, vel_fc)
    call this%pressure_cc_correct(props, p_cc)
    call this%grad_p_rho(props, p_cc)
    call this%velocity_cc_correct(dt, props, grad_p_rho_cc_n, vel_cc)
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
        this%bc%dp_dirichlet, props%inactive_f, 0.0_r8)

    ! divergence free face velocity
    associate (m => this%mesh%mesh)
      do j = 1, m%nface_onP
        if (props%inactive_f(j) > 0) then
          vel_fc(j) = 0.0_r8
        else
          vel_fc(j) = vel_fc(j) - &
              dt*dot_product(this%grad_fc(:,j),m%normal(:,j)/m%area(j))/props%rho_fc(j)
        end if
      end do

      associate (index => this%bc%v_dirichlet%index, value => this%bc%v_dirichlet%value)
        do i = 1, size(index)
          j = index(i)
          vel_fc(j) = dot_product(value(:,i), m%normal(:,j))/m%area(j)
        end do
      end associate

      call gather_boundary(m%face_ip, vel_fc)
    end associate

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

    associate (m => this%mesh%mesh, p => props, dp => this%delta_p_cc)
      ! remove null space from dp if applicable
!!$      if (.not.(this%bc%pressure_d .or. p%any_void)) then
!!$        ! maybe this should use inactive_c... needs to be consistenent with how
!!$        ! the operator is constructed
!!$        dp_vof = 0.0_r8
!!$        vof = 0.0_r8
!!$        do i = 1, m%ncell_onP
!!$          vof = vof + p%vof(i)*m%volume(i)
!!$          dp_vof = dp_vof + p%vof(i)*m%volume(i)*dp(i)
!!$        end do
!!$        g_vof = global_sum(vof)
!!$        g_dp_vof = global_sum(dp_vof)
!!$#if ASDF
!!$        print* , "<VOF>: ", g_vof
!!$        print* , "<DP_VOF>: ", g_dp_vof
!!$#endif
!!$        dp(1:m%ncell_onP) = dp(1:m%ncell_onP) - g_dp_vof/g_vof
!!$      endif
      p_cc = p_cc + dp
      call gather_boundary(m%cell_ip, p_cc)
    end associate
#if ASDF
    write(*,'("Corrected P_CC[",i3,"]: ",es15.5)') Q, p_cc(Q)
#endif
  end subroutine pressure_cc_correct


  ! compute face and cell centered dynamic_pressure_gradient/rho
  subroutine grad_p_rho(this, props, p_cc)
    class(flow_projection), intent(inout) :: this
    type(flow_props), intent(in) :: props
    real(r8), intent(in) :: p_cc(:)
    !-
    integer :: i, j

    associate (m => this%mesh%mesh, rho => props%rho_fc, &
        gp_cc => this%grad_p_rho_cc, gp_fc => this%grad_p_rho_fc)

      call gradient_cf(gp_fc, p_cc, this%bc%p_neumann, &
          this%bc%p_dirichlet, props%inactive_f, 0.0_r8, this%gravity_head)

      do j = 1, m%nface_onP
        if (rho(j) > 0.0_r8) then
          gp_fc(:,j) = gp_fc(:,j)/rho(j)
        else
          gp_fc(:,j) = 0.0_r8
        end if
      end do
#if ASDF
      associate (faces => m%cface(m%xcface(771):m%xcface(772)-1))
        write(*,'("grad(p) at [",i4,"] y- ",3es20.12)') faces(1), gp_fc(:,faces(1))
        write(*,'("grad(p) at [",i4,"] y+ ",3es20.12)') faces(3), gp_fc(:,faces(3))
        write(*,'("grad(p) at [",i4,"] x- ",3es20.12)') faces(5), gp_fc(:,faces(5))
        write(*,'("grad(p) at [",i4,"] x+ ",3es20.12)') faces(6), gp_fc(:,faces(6))
      end associate
#endif
      call gather_boundary(m%face_ip, gp_fc)
      call interpolate_fc(gp_cc, gp_fc, props%inactive_f, this%bc%p_neumann%index)
#if ASDF
      write(*,'("grad_p_rho_cc[",i4,"] ",3es20.12)') Q, gp_cc(:,Q)
#endif
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

#if ASDF
    write(*,*) "<<< VELOCITY_CC_CORRECT"
    write(*,*) "dt: ", dt
    write(*,'("vel_cc[",i4,"] ",3es20.12)') Q, vel_cc(:,Q)
    write(*,'("grad_p_rho_cc_n[",i4,"] ",3es20.12)') Q, grad_p_rho_cc_n(:,Q)
    write(*,'("grad_p_rho_cc[",i4,"]   ",3es20.12)') Q, this%grad_p_rho_cc(:,Q)
#endif
    associate (m => this%mesh%mesh)
      ! assumes dynamic pressure gradients have already been computed
      do i = 1, m%ncell_onP
        ! -dt * grad(dP)/rho
        vel_cc(:,i) = vel_cc(:,i) - dt*(this%grad_p_rho_cc(:,i)-grad_p_rho_cc_n(:,i))
      end do
    end associate
#if ASDF
    write(*,'("vel_cc[",i4,"] ",3es20.12)') Q, vel_cc(:,Q)
#endif
  end subroutine velocity_cc_correct



end module flow_projection_type
