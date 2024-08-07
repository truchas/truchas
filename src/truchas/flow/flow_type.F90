!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module flow_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use flow_domain_types
  use truchas_logging_services
  use truchas_timers
  use flow_projection_type
  use flow_prediction_type
  use flow_props_type
  use flow_bc_type
  use unstr_mesh_type
  use parallel_communication
  implicit none
  private

  type, public :: flow
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(flow_props), pointer :: props => null() ! unowned reference
    real(r8), allocatable :: vel_cc(:,:) ! cell-centered velocity (dims, ncell)
    real(r8), allocatable :: vel_cc_n(:,:) ! cell-centered velocity (dims, ncell)
    real(r8), allocatable :: vel_fn(:) ! outward oriented face-normal velocity
    real(r8), allocatable :: vel_fn_n(:) ! outward oriented face-normal velocity
    real(r8), allocatable :: P_cc(:) ! cell-centered pressure
    real(r8), allocatable :: grad_p_rho_cc_n(:,:) ! dynamic pressure gradient over rho
    real(r8) :: body_force(3)
    type(flow_projection) :: proj
    type(flow_prediction) :: pred
    type(flow_bc), pointer :: bc
    logical :: inviscid, prescribed
    real(r8) :: viscous_number
    real(r8) :: courant_number
  contains
    procedure :: init
    procedure, private :: set_initial_state_start, set_initial_state_restart
    generic   :: set_initial_state => set_initial_state_start, set_initial_state_restart
    procedure :: step
    procedure :: accept
    procedure :: correct_non_regular_cells
    procedure :: timestep
    ! views into data
    procedure :: vel_cc_view
    procedure :: vel_fn_view
    procedure :: P_cc_view
    procedure :: dump_state
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

  subroutine timestep(this, dt, dt_tag)

    use parallel_communication

    class(flow), intent(in) :: this
    real(r8), intent(out) :: dt
    character(:), allocatable, intent(inout) :: dt_tag

    real(r8) :: v, dtc, dtv
    integer :: j

    dtc = huge(1.0_r8)
    dtv = huge(1.0_r8)

    ! try something different (simpler) compared to truchas
    associate (vof => this%props%vof, rho => this%props%rho_cc, mu => this%props%mu_cc, m => this%mesh)
      do j = 1, m%nface_onP
        dtc = min(dtc, m%face_normal_dist(j)/(epsilon(1.0_r8)+abs(this%vel_fn(j))))
      end do
      dtc = this%courant_number * global_minval(dtc)

      if (.not.this%inviscid .and. this%viscous_number > 0.0_r8) then
        do j = 1, m%ncell_onP
          if (this%props%cell_t(j) <= regular_t .and. mu(j) > 0.0_r8) then
            ! Original simple version
            !v = (m%volume(j)*vof(j))**(1.0_r8/3.0_r8)
            ! Something closer to what Truchas originally had
            associate (cface => m%cface(m%xcface(j):m%xcface(j+1)-1))
              v = minval(m%face_normal_dist(cface))
            end associate
            dtv = min(dtv, v**2*rho(j)/mu(j))
          end if
        end do
        dtv = this%viscous_number * global_minval(dtv)
      end if
    end associate

    dt = min(dtc, dtv)
    if (dt == dtc) then
      dt_tag = 'courant'
    else
      dt_tag = 'viscous'
    end if

  end subroutine timestep

  subroutine init(this, mesh, props, flowbc, prescribed, params)

    use parameter_list_type

    class(flow), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(flow_props), intent(in), target :: props
    type(flow_bc), intent(inout), target :: flowbc
    logical, optional, intent(in) :: prescribed
    type(parameter_list), intent(inout) :: params

    type(parameter_list), pointer :: plist
    real(r8), allocatable :: array(:)

    this%mesh => mesh
    this%props => props
    this%bc => flowbc

    plist => params%sublist('options')
    call plist%get('inviscid', this%inviscid, default=.false.)
    if (.not.this%inviscid) &
        call plist%get('viscous-number', this%viscous_number, default=0.0_r8)
    call plist%get('courant number', this%courant_number, default=0.5_r8)
    call plist%get('body force', array, default=[0.0_r8, 0.0_r8, 0.0_r8])
    this%body_force = array

    associate (nc => mesh%ncell, nf => mesh%nface)
      allocate(this%vel_cc(3, nc))
      allocate(this%vel_cc_n(3, nc))
      allocate(this%P_cc(nc))
      allocate(this%vel_fn(nf))
      allocate(this%vel_fn_n(nf))
      allocate(this%grad_p_rho_cc_n(3,nc))
    end associate

    ! Disable BC initialization and the prediction and
    ! projection solvers if we have prescribed flow.
    if (present(prescribed)) then
      this%prescribed = prescribed
    else
      this%prescribed = .false.
    end if
    if (this%prescribed) return

    call this%pred%init(mesh, this%bc, this%inviscid, params)
    call this%proj%init(mesh, this%bc, params)

    !call this%proj%grad_p_rho

  end subroutine init

  !! This assigns initial values to all the state variables in the case only
  !! external state variables are supplied. The remaining state variables are
  !! computed/derived from them. This is indicated for starting a simulation.

  subroutine set_initial_state_start(this, t, dt, vcell, vof, tcell)

    class(flow), intent(inout) :: this

    real(r8), intent(in) :: t, dt
    real(r8), intent(in) :: vcell(:,:), vof(:,:), tcell(:)

    integer :: i, j
    real(r8) :: vel(3)

    ASSERT(size(vcell,dim=1)==3)
    ASSERT(size(vcell,dim=2)>=this%mesh%ncell_onP)

    call TLS_info('')
    call TLS_info('Computing initial flow pressures and fluxing velocities ...')

    call this%props%set_initial_state(vof, tcell)

    this%vel_cc(:,:this%mesh%ncell_onP) = vcell(:,:this%mesh%ncell_onP)
    call this%mesh%cell_imap%gather_offp(this%vel_cc)

    !! NB: The face velocities should be discretely solenoidal. So the following
    !! naive averaging is only valid for a uniform velocity field. For anything
    !! else we ought to be doing something different.

    do j = 1, this%mesh%nface_onP
      if (this%props%face_t(j) == regular_void_t) then
        ! Cells neighboring void and non-void should take on the value of the non-void cell.
        i = findloc(this%props%cell_t(this%mesh%fcell(:,j)) /= void_t, .true., dim=1)
        vel = this%vel_cc(:,this%mesh%fcell(i,j))
      else
        vel = this%vel_cc(:,this%mesh%fcell(1,j))
        if (this%mesh%fcell(2,j) /= 0) vel = (vel + this%vel_cc(:,this%mesh%fcell(2,j))) / 2
      end if

      this%vel_fn(j) = dot_product(this%mesh%normal(:,j), vel) / this%mesh%area(j)
    end do
    call this%mesh%face_imap%gather_offp(this%vel_fn)

    call this%bc%compute_initial(t)
    associate (faces => this%bc%v_dirichlet%index, vel => this%bc%v_dirichlet%value)
      do i = 1, size(faces)
        j = faces(i)
        this%vel_fn(j) = dot_product(this%mesh%normal(:,j), vel(:,i)) / this%mesh%area(j)
      end do
    end associate

    this%vel_fn_n = this%vel_fn

    this%vel_cc_n = this%vel_cc !TODO: put directly into vel_cc_n and by-pass vel_cc?

    if (this%props%any_real_fluid) then
      call compute_initial_pressure(this, t, dt)
    else
      this%p_cc = 0
      this%grad_p_rho_cc_n = 0
    end if

  end subroutine set_initial_state_start

  !! This assigns initial values to all the state variables in the case they
  !! are all supplied. This is indicated for restarting a simulation from saved
  !! state data. NB: The current restart file format does not contain data for
  !! the cell-centered dynamic pressure gradient, so it is computed.

  subroutine set_initial_state_restart(this, t, pcell, vcell, vface, vof, tcell)

    class(flow), intent(inout) :: this
    real(r8), intent(in) :: t, pcell(:), vcell(:,:), vface(:), tcell(:), vof(:,:)

    ASSERT(size(vcell,dim=1)==3)
    ASSERT(size(vcell,dim=2)>=this%mesh%ncell_onP)
    ASSERT(size(vface)>=this%mesh%nface_onP)

    this%p_cc(:this%mesh%ncell_onP) = pcell(:this%mesh%ncell_onP)
    call this%mesh%cell_imap%gather_offp(this%p_cc)

    this%vel_cc(:,:this%mesh%ncell_onP) = vcell(:,:this%mesh%ncell_onP)
    call this%mesh%cell_imap%gather_offp(this%vel_cc)
    this%vel_cc_n = this%vel_cc

    this%vel_fn(:this%mesh%nface_onP) = vface(:this%mesh%nface_onP)
    call this%mesh%face_imap%gather_offp(this%vel_fn)
    !FIXME? Impose Dirichlet velocity conditions on vel_fn?
    this%vel_fn_n = this%vel_fn

    call this%bc%compute_initial(t)
    call this%props%set_initial_state(vof, tcell)
    call this%proj%get_dyn_press_grad(this%props, this%p_cc, this%body_force, this%grad_p_rho_cc_n)

    !call this%dump_state

  end subroutine set_initial_state_restart

  !! This computes the initial cell-centered pressure field and gradient given
  !! pressure BC and the initial velocity field (both cell and face-centered).
  !! It does this by hijacking flow time step procedures to take an artificial
  !! time step and being selective with the inputs and outputs. The momentum
  !! advection term is currently omitted (in effect) and so this is only
  !! correct for cases where that term is zero at the initial state, which is
  !! the case for a uniform velocity field, for example.

  subroutine compute_initial_pressure(this, t, dt)

    class(flow), intent(inout) :: this
    real(r8), intent(in) :: t, dt

    real(r8), allocatable :: flux_volumes(:,:)

    call this%bc%compute_initial(t)

    !TODO: flow needs a reference to volume tracker or advector object
    allocate(flux_volumes(size(this%props%density),size(this%mesh%cface)))
    flux_volumes = 0.0_r8

    this%vel_cc = this%vel_cc_n
    this%grad_p_rho_cc_n = 0.0_r8
    call this%pred%setup(dt, this%props, this%vel_cc)
    call this%pred%solve(dt, this%props, this%grad_p_rho_cc_n, flux_volumes, this%vel_fn_n, this%vel_cc)

    this%p_cc = 0.0_r8
    call this%proj%setup(dt, this%props, this%grad_p_rho_cc_n, this%body_force, this%vel_cc, this%P_cc, this%vel_fn)
    call this%proj%solve(dt, this%props, this%grad_p_rho_cc_n, this%vel_cc, this%P_cc, this%vel_fn)
    call this%correct_non_regular_cells()

    !! Reject the updated velocities
    this%vel_cc = this%vel_cc_n
    this%vel_fn = this%vel_fn_n

    !! Accept the updated pressures
    !this%p_cc_n = this%p_cc  !FIXME: need to maintain pressure as rejectable state variable
    call this%proj%accept(this%grad_p_rho_cc_n) ! just this%grad_p_rho_cc_n = this%proj%grad_p_rho_cc

  end subroutine compute_initial_pressure


  subroutine correct_non_regular_cells(this)

    class(flow), intent(inout) :: this

    integer :: i

    associate (cell_t => this%props%cell_t, face_t => this%props%face_t, vof_novoid => this%props%vof_novoid)
      do i = 1, this%mesh%ncell
        if (vof_novoid(i) == 0.0_r8) then
          this%vel_cc(:,i) = 0.0_r8
          this%P_cc(i) = 0.0_r8
        end if
      end do

      do i = 1, this%mesh%nface_onP
        if (face_t(i) > regular_t) then
          this%vel_fn(i) = 0.0_r8
        end if
      end do
    end associate

  end subroutine correct_non_regular_cells


  subroutine step(this, t, dt, vof, flux_volumes, tcell)

    class(flow), intent(inout) :: this
    real(r8), intent(in) :: t, dt, vof(:,:), flux_volumes(:,:), tcell(:)

    call this%props%update_cc(vof, tcell)

    if (.not.this%props%any_real_fluid) then
      this%vel_cc = 0
      this%vel_fn = 0
      return
    end if

    call this%bc%compute(t, dt)

    call start_timer("prediction")
    this%vel_cc = this%vel_cc_n
    call this%pred%setup(dt, this%props, this%vel_cc)
    call this%pred%solve(dt, this%props, this%grad_p_rho_cc_n, flux_volumes, this%vel_fn_n, this%vel_cc)
    call stop_timer("prediction")

    call start_timer("projection")
    call this%proj%setup(dt, this%props, this%grad_p_rho_cc_n, this%body_force, this%vel_cc, this%P_cc, this%vel_fn)
    call this%proj%solve(dt, this%props, this%grad_p_rho_cc_n, this%vel_cc, this%P_cc, this%vel_fn)
    call stop_timer("projection")

    call this%correct_non_regular_cells()

  end subroutine step


  subroutine accept(this)

    class(flow), intent(inout) :: this

    this%vel_cc_n = this%vel_cc
    this%vel_fn_n = this%vel_fn

    ! disable the prediction and projection solvers if we have prescribed flow
    if (this%prescribed) return

    call this%pred%accept() ! currently a no-op
    call this%proj%accept(this%grad_p_rho_cc_n) ! just this%grad_p_rho_cc_n = this%proj%grad_p_rho_cc

  end subroutine accept

  subroutine dump_state(this)
    use truchas_logging_services
    class(flow), intent(in) :: this
    integer :: lun
    lun = TLS_debug_unit()
    write(lun,'(/,a)') 'VEL_CC_N'
    write(lun,'(3es20.12)') this%vel_cc_n(:,:this%mesh%ncell_onP)
    write(lun,'(/,a)') 'VEL_FN_N'
    write(lun,'(4es20.12)') this%vel_fn_n(:this%mesh%nface_onP)
    write(lun,'(/,a)') 'P_CC'
    write(lun,'(4es20.12)') this%p_cc(:this%mesh%ncell_onP)
    write(lun,'(/,a)') 'GRAD_P_RHO_CC_N'
    write(lun,'(3es20.12)') this%grad_p_rho_cc_n(:,:this%mesh%ncell_onP)
  end subroutine dump_state

end module flow_type
