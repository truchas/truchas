module flow_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use flow_mesh_type
  use flow_projection_type
  use flow_prediction_type
  use flow_props_type
  use flow_bc_type
  use unstr_mesh_type
  use index_partitioning
  use parallel_communication
  implicit none
  private

  public :: flow

  type :: flow
    type(flow_mesh), pointer :: mesh => null() ! unowned reference
    real(r8), allocatable :: vel_cc(:,:) ! cell-centered velocity (dims, ncell)
    real(r8), allocatable :: vel_fn(:) ! outward oriented face-normal velocity
    real(r8), allocatable :: P_cc(:) ! cell-centered pressure
    real(r8), allocatable :: grad_p_rho_cc_n(:,:) ! dynamic pressure gradient over rho
    real(r8) :: body_force(3)
    type(flow_projection) :: proj
    type(flow_prediction) :: pred
    type(flow_bc), pointer :: bc
    logical :: inviscid, stokes
    real(r8) :: viscous_number
    real(r8) :: courant_number
  contains
    procedure :: read_params
    procedure :: init
    procedure :: step
    procedure :: accept
    procedure :: zero_out_solid_velocities
    procedure :: timestep
    ! views into data
    procedure :: vel_cc_view
    procedure :: vel_fn_view
    procedure :: P_cc_view
  end type flow

contains

  function vel_cc_view(this) result(p)
    class(flow), target, intent(in) :: this
    real(r8), pointer :: p(:,:)
    p => this%vel_cc
  end function vel_cc_view

  function P_cc_view(this) result(p)
    class(flow), target, intent(in) :: this
    real(r8), pointer :: p(:)
    p => this%P_cc
  end function P_cc_view

  function vel_fn_view(this) result(p)
    class(flow), target, intent(in) :: this
    real(r8), pointer :: p(:)
    p => this%vel_fn
  end function vel_fn_view

  subroutine timestep(this, props, dtc, dtv)
    use parallel_communication

    class(flow), intent(in) :: this
    type(flow_props), intent(in) :: props
    real(r8), intent(out) :: dtc, dtv
    real(r8) :: v
    integer :: i, j

    dtc = huge(1.0_r8)
    dtv = huge(1.0_r8)

    ! try something different (simpler) compared to truchas
    associate( m => this%mesh%mesh, vof => props%vof, &
        rho => props%rho_cc, mu => props%mu_cc)

      do j = 1, m%ncell_onP
#ifndef NDEBUG
        print *, "vof(",j,"):",vof(j)
#endif
        if (vof(j) > 0.0_r8) then
          v = (m%volume(j)*vof(j))**(1.0_r8/3.0_r8)
          dtc = min(dtc, v/(epsilon(1.0_r8)+norm2(this%vel_cc(:,j))))
          if (.not.this%inviscid .and. mu(j) > 0.0_r8) &
              dtv = min(dtv, v**2*rho(j)/mu(j))
        end if
      end do
    end associate

    dtc = this%courant_number * global_minval(dtc)
    dtv = this%viscous_number * global_minval(dtv)

  end subroutine timestep

  subroutine read_params(this, p)
    use parameter_list_type
    class(flow), intent(inout) :: this
    type(parameter_list), pointer, intent(in) :: p
    type(parameter_list), pointer :: pp
    real(r8), allocatable :: body_force(:)
    integer :: stat
    !

    pp => p%sublist("options")

    call pp%get('inviscid', this%inviscid, .false.)
    call pp%get('stokes', this%stokes, .false.)
    call pp%get('viscous number', this%viscous_number, 0.1_r8)
    call pp%get('courant number', this%courant_number, 0.5_r8)
    call pp%get('body force', body_force, stat=stat)
    if (stat == 0) then
      this%body_force = body_force
    else
      this%body_force = 0.0_r8
    end if

    allocate(this%bc)
    call this%bc%read_params(p)
    call this%proj%read_params(p)
    call this%pred%read_params(p)

  end subroutine read_params


  subroutine init(this, m, vel_cc, P_cc)
    class(flow), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: m
    real(r8), optional, intent(in) :: vel_cc, P_cc
    !-

    this%mesh => m

    call this%bc%init(m)

    associate (nc => m%mesh%ncell, nf => m%mesh%nface)
      allocate(this%vel_cc(3, nc))
      allocate(this%P_cc(nc))
      allocate(this%vel_fn(nf))
      allocate(this%grad_p_rho_cc_n(3,nc))
    end associate


    if (present(vel_cc)) then
      this%vel_cc = vel_cc
    else
      this%vel_cc = 0.0_r8
    end if

    if (present(P_cc)) then
      this%P_cc = P_cc
    else
      this%P_cc = 0.0_r8
    end if

    this%grad_p_rho_cc_n = 0.0_r8
    this%vel_fn = 0.0_r8

    call this%pred%init(m, this%bc, this%inviscid, this%stokes)
    call this%proj%init(m, this%bc)

  end subroutine init

  subroutine zero_out_solid_velocities(this, props)
    class(flow), intent(inout) :: this
    type(flow_props), intent(in) :: props

    print *, 'zero out all velocities on solid cells and faces'
  end subroutine zero_out_solid_velocities

  subroutine step(this, t, dt, props, flux_volumes, initial)
    class(flow), intent(inout) :: this
    real(r8), intent(in) :: t, dt, flux_volumes(:,:)
    type(flow_props), intent(inout) :: props
    logical, optional, intent(in) ::  initial
    !-
    real(r8) :: p_max, p_min
    integer :: j

    call this%bc%compute(t, dt, initial=initial)

#ifndef NDEBUG
    write(*, *) " >> Pre Predictor"
    do j = 1, this%mesh%mesh%ncell_onP
      write(*,'("v[",i3,"]: ",2es15.5, " | P:",es15.5)') &
          j,this%vel_cc(1:2,j), this%P_cc(j)
    end do
#endif
!!$    do j = 1, this%mesh%mesh%nface_onP
!!$      write(*,'(a,i3,a,3es15.5)') "vel_fn(",j,"):",this%vel_fn(j)
!!$    end do

    call this%pred%setup(dt, props, this%vel_cc, initial=initial)
    call this%pred%solve(dt, props, this%grad_p_rho_cc_n, flux_volumes, this%vel_fn, this%vel_cc, initial=initial)

#ifndef NDEBUG
    write(*, *) " >> Post Predictor"
    do j = 1, this%mesh%mesh%ncell_onP
      write(*,'("v[",i3,"]: ",2es15.5, " | P:",es15.5)') &
          j,this%vel_cc(1:2,j), this%P_cc(j)
    end do
#endif
!!$    do j = 1, this%mesh%mesh%nface_onP
!!$      write(*,'(a,i3,a,3es15.5)') "vel_fn(",j,"):",this%vel_fn(j)
!!$    end do


    call this%proj%setup(dt, props, this%grad_p_rho_cc_n, this%body_force, this%vel_cc, this%P_cc, this%vel_fn, initial=initial)
    call this%proj%solve(dt, props, this%grad_p_rho_cc_n, this%vel_cc, this%P_cc, this%vel_fn, initial=initial)

    p_min = global_minval(this%p_cc)
    p_max = global_maxval(this%p_cc)

#ifndef NDEBUG
    write(*,'(a,es15.5)') "<< time : ", t
    write(*,'(a,2es15.5)') "   [P_min, P_max]: ", p_min, p_max
#endif
  end subroutine step

  subroutine accept(this)
    class(flow), intent(inout) :: this

    call this%pred%accept()
    call this%proj%accept(this%grad_p_rho_cc_n)

  end subroutine accept

end module flow_type
