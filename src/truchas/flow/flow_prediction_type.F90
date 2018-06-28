module flow_prediction_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use parameter_list_type
  use flow_mesh_type
  use flow_props_type
  use flow_bc_type
  use turbulence_models
  implicit none
  private

  public :: flow_prediction_type

  type :: flow_prediction
    type(flow_mesh), pointer :: mesh
    type(flow_bc), pointer :: bc
    class(turbulence_model), allocatable :: turb
    !class(porous_drag_model), allocatable :: drag
    type(hypre_hybrid) :: solver(3)
    type(parameter_list), pointer :: p
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

  subroutine read_params(this, p)
    class(flow_prediction) :: this
    type(parameter_list), pointer, intent(in) :: p
    !-
    character(:), allocatable :: flow_type
    integer :: stat
    type(parameter_list), pointer :: sub

    this%p => p

    ! default is full navier-stokes equations
    this%inviscid = .false.
    this%stokes =.false.

    call p%get('flow-type', flow_type, stat=stat)
    if(flow_type == 'inviscid') then
      this%inviscid = .true.
    else if(flow_type == 'stokes') then
      this%stokes = .true.
    else
      call TLS_fatal("Unsupported flow-type"//flow_type)
    end if

    ! default to crank-nicolson
    call p%get('viscous-implicitness', this%viscous_implicitness, 0.5_r8)
    ! default to fully implicit
    call p%get('solidfy-implicitness', this%solidify_implicitness, 1.0_r8)

    sub => p%sublist("turbulence-model")
    call turbulence_models_read_params(this%turb, sub, off=this%inviscid)

    !sub => p%sublist("porous-drag-model")
    !call porous_drag_models_read_params(this%drag, sub, off=this%inviscid)
  end subroutine read_params


  subroutine init(this, mesh, bc)
    class(flow_prediction), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: mesh
    type(flow_bc), target, intent(inout) :: bc
    !-
    type(graph_container) :: g(size(this%solver))
    type(matrix_container) :: A(size(this%solver))
    type(ip_desc), pointer :: row_ip
    type(parameter_list), pointer :: sub
    integer :: i, j, k

    this%mesh => mesh
    this%bc => bc

    associate (m => mesh%mesh)
      allocate(this%rhs(3,m%ncell))
      allocate(this%grad_fc_vector(3,3,m%nface))
      allocate(this%sol(m%ncell))
      allocate(this%rhs1d(m%ncell))

      if (this%viscous_implicitness > 0.0_r8) then
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

        do k = 1, size(A)
          allocate(A(k)%A)
          call A(k)%A%init(g(k), take_graph=.true.)
          sub => this%p%sublist("solver")
          call this%solver(k)%init(A(k)%A, sub)
        end do
      end if
    end associate
  end subroutine init

  subroutine setup(this, dt, props, vel_cc)
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: dt, vel_cc(:,:)
    type(flow_props), intent(inout) :: props

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

    !if (this%inviscid .or. this%viscous_implicitness == 0.0_r8) return

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
          ASSERT(this%mesh%fcell(1,fi) == 0)

          coeff = dtv*dot_product(ds(:,fi), m%normal(:,fi))*mu_f(fi)/m%volume(j)
          do k = 1, size(A)
            call A(k)%A%add_to(j, j, coeff)
          end do
          this%rhs(:,j) = this%rhs(:,j) - coeff*value(:,i)
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
    call this%solver%setup()

  end subroutine setup_solver


  subroutine accumulate_rhs_pressure(this, dt, props, grad_p_rho_cc)
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: dt, grad_p_rho_cc(:,:)
    integer :: j

    associate (m => this%mesh%mesh)
      do j = 1, m%ncell_onP
        rhs(:,j) = rhs(:,j) + dt*grad_p_rho_cc(:,j)*props%rho_cc(j)
      end do
    end associate

  end subroutine accumulate_rhs_pressure


  subroutine accumulate_rhs_viscous_stress(this, dt, props, vel_cc)
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: dt, vel_cc(:)
    type(flow_props), intent(in) :: props
    !
    integer :: i, j
    real(r8) :: dtv

    if (this%inviscid .or. this%viscous_implicitness == 1.0_r8) return

    dtv = dt*(1.0_r8 - this%viscous_implicitness)

    associate (m => this%mesh%mesh, g => this%grad_fc_vector, mu => props%mu_fc)
      ! velocity gradient tensor at faces
      call gradient_cf(g, vel_cc, this%bc%v_zero_normal, this%bc%v_dirichlet)
      call gather_boundary(m%face_ip, g)

      do j = 1, m%nface
        associate (n1 => this%mesh%fcell(1,j), n2 => this%mesh%fcell(2,j), rhs => this%rhs)
          if (n1 > 0) then ! inward oriented normal
            do i = 1, 3
              rhs(i,n1) = rhs(i,n1) + &
                  dtv*mu(j)*dot_product(m%normal(:,j),g(:,i,j))/m%volume(n1)
            end do
          end if
          if (n2 > 0) then ! outward oriented normal
            do i = 1, 3
              rhs(i,n2) = rhs(i,n2) - &
                  dtv*mu(j)*dot_product(m%normal(:,j),g(:,i,j))/m%volume(n2)
            end do
          end if
        end associate
      end do
    end associate
  end subroutine accumulate_rhs_viscous_stress


  subroutine accumulate_rhs_momentum(this, props, vel_cc, vel_fc, flux_volumes, flux_vel)
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: flux_volumes(:,:), flux_vel(:), vel_cc(:,:), vel_fc(:,:)
    type(flow_props), intent(in) :: props
    !
    integer :: i, j, f0, f1, nhbr, cn
    real(r8) :: v

    ! we're going to assume here that the boundary conditions on flux_volumes
    ! and flux_vel have already been set..
    associate (m => this%mesh%mesh)
      do i = 1, m%ncell_onP
        f0 = m%xcface(i)
        f1 = m%xcface(i+1)-1
        nhbr = m%xcnhbr(i)

        do j = f0, f1
          k = m%cface(j)
          v = m%volume(i)
          ! donor cell upwinding
          if (any(flux_volumes(:,j) > 0.0_r8)) then
            rhs(:,i) = rhs(:,i) - vel_cc(:,i)* &
                dot_product(flux_volumes(:,j),props%density(:))/v
          else
            ! neighbor index
            cn = m%cnhbr(nhbr+(j-f0+1))
            if (cn > 0) then
              rhs(:,i) = rhs(:,i) - vel_cc(:,cn)* &
                  dot_product(flux_volumes(:,j),props%density(:))/v
            else
              ! this must be an inflow boundary so use face velocity
              ! not so sure if this is the best approach.  I might need a bc for
              ! velocity inflow
              rhs(:,i) = rhs(:,i) - vel_fc(:,k)* &
                  dot_product(flux_volumes(:,j),props%density(:))/v
            end if
          end if
        end do
      end do
    end associate
  end subroutine accumulate_rhs_momentum


  subroutine accumulate_rhs_solidified_rho(this, props, vel_cc)
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: vel_cc(:,:), dt
    type(flow_props), intent(in) :: props
    !
    integer :: j
    real(r8) :: w

    w = 1.0_r8-this%solidify_implicitness

    if (this%solidify_implicitness < 1.0_rp) then

      associate (m => this%mesh%mesh, rho => props%solidified_rho)
        do j = 1, m%ncell_onP
          rhs(:,j) = rhs(:,j) - w*rho(j)*vel_cc(:,j)
        end do
      end associate
    end if
  end subroutine accumulate_rhs_solidified_rho


  subroutine solve(this, dt, props, grad_p_rho_cc, flux_volumes, flux_vel)
    class(flow_prediction), intent(inout) :: this
    real(r8), intent(in) :: dt, flux_volumes(:,:), flux_vel(:)
    type(flow_props), intent(in) :: props
    !
    integer :: i, j, ierr

    ! Pressure and viscous forces are unaware of solid/fluid boundaries.
    ! These forces are computed first and scaled by vof as a crude approximation.
    ! Once a better solid/fluid model is used, change this numerical hackjob.
    call this%accumulate_rhs_pressure(dt, props, grad_p_rho_cc)
    call this%accumulate_rhs_viscous_stress(dt, props, vel_cc)
    associate (m => this%mesh%mesh)
      do j = 1, m%ncell_onP
        rhs(:,j) = rhs(:,j)*props%vof(j)
      end do
    end associate

    call this%accumulate_rhs_momentum()
    call this%accumulate_rhs_solidified_rho()

    ! handle surface tension
    ! handle porous drag

    associate (m => this%mesh%mesh, rho => props%rho_cc_n, vof => props%vof_n)
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

        call this%solver(i)%solve(rhs1d, sol, ierr)
        if (ierr /= 0) call tls_error("prediction solve unsuccessful")
        call gather_boundary(m%cell_ip, sol)
        do j = 1, m%ncell
          vel_cc(i,j) = sol(j)
        end do
      end do
    end associate
  end subroutine solve

  subroutine accept(this)
    class(flow_prediction), intent(inout) :: this
    call this%turb%accept()
  end subroutine accept
end module flow_prediction_type
