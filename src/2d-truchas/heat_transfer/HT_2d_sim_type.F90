!TODO: finish documentation
!!
!! HT_2D_SIM_TYPE
!!
!! This module defines a derived type that encapsulates a heat transfer
!! simulation.  This drives the time integration and generates the output.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! August 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module HT_2d_sim_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use matl_mesh_func_type
  use material_database_type
  use scalar_func_factories
  use mfd_2d_disc_type
  use HT_2d_model_type
  use HT_2d_solver_type
  use xdmf_file_type
  use time_step_sync_type
  use parallel_communication
  use truchas_logging_services
  use timer_tree_type
  implicit none
  private

  type, public:: HT_2d_sim
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(mfd_2d_disc), pointer :: disc => null()
    type(matl_mesh_func), pointer :: mmf => null()
    type(HT_2d_model), pointer :: model => null()
    type(HT_2d_solver), pointer :: solver => null()
    type(xdmf_file) :: xdmf
    !! Integration control
    real(r8) :: t_init
    real(r8) :: tlast, hlast
    real(r8) :: dt_init, dt_min, dt_max
    integer  :: max_try
    real(r8), allocatable :: tout(:)
    type(time_step_sync) :: ts_sync
  contains
    final :: HT_2d_sim_delete
    procedure :: init
    procedure :: run
    procedure :: write_solution
  end type HT_2d_sim

contains

  subroutine HT_2d_sim_delete(this)
    type(HT_2d_sim), intent(inout) :: this
    if (associated(this%mesh)) deallocate(this%mesh)
    if (associated(this%disc)) deallocate(this%disc)
    if (associated(this%model)) deallocate(this%model)
    if (associated(this%solver)) deallocate(this%solver)
  end subroutine HT_2d_sim_delete


  subroutine init(this, params)

    use parameter_list_type
    use unstr_2d_mesh_factory
    use signal_handler, only: init_signal_handler, SIGURG

    class(HT_2d_sim), intent(out) :: this
    type(parameter_list) :: params

    type(parameter_list), pointer :: plist
    type(material_database) :: matl_db
    class(scalar_func), allocatable :: f
    character(:), allocatable :: errmsg, context
    real(r8), allocatable :: temp(:), xmin(:), xmax(:)
    integer, allocatable :: nx(:)
    real(r8) :: eps, rel_tol
    integer :: stat, max_itr

    call start_timer('initialization')
    call TLS_info('Initializing the simulation', TLS_VERB_NOISY)

    !! Catch SIGURG signals.
    call init_signal_handler(SIGURG)

    !! Create the mesh.
    call start_timer('mesh')
    if (params%is_sublist('mesh')) then
      plist => params%sublist('mesh')
      context = 'processing ' // plist%name() // ': '
      call plist%get('xmin', xmin, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call plist%get('xmax', xmax, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call plist%get('nx', nx, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call plist%get('eps', eps, default=0.0_r8, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      this%mesh => new_unstr_2d_mesh(xmin, xmax, nx, eps)
    else
      call TLS_fatal('missing "mesh" sublist parameter')
    end if
    call stop_timer('mesh')

    !! Create the discretization object.
    call start_timer('mfd-discretization')
    allocate(this%disc)
    call this%disc%init(this%mesh)
    call stop_timer('mfd-discretization')

    !TODO: finalize material init
    !! Create the material object.
    ! if (params%is_sublist('material')) then
    !   plist => params%sublist('material')
    !   context = 'processing ' // plist%name() // ': '
    !   allocate(this%mmf)
    !   call load_material_database(matl_db, plist, stat, errmsg)
    !   if (stat /= 0) call TLS_fatal(context//errmsg)
    !   call this%mat%init(plist)
    ! else
    !   call LS_fatal('missing "material" sublist parameter')
    ! end if
    allocate(this%mmf)
    call init_materials(this%mesh, matl_db, this%mmf)

    !! Create the heat conduction model.
    call start_timer('ht-model')
    if (params%is_sublist('ht-model')) then
      plist => params%sublist('ht-model')
      context = 'processing ' // plist%name() // ': '
      allocate(this%model)
      call this%model%init(this%disc, this%mmf, plist, stat, errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
    else
      call TLS_fatal('missing "ht-model" sublist parameter')
    end if
    call stop_timer('ht-model')

    !! Create the heat conduction solver.
    call start_timer('ht-solver')
    if (params%is_sublist('ht-solver')) then
      plist => params%sublist('ht-solver')
      allocate(this%solver)
      call this%solver%init(this%model, plist)
    else
      call TLS_fatal('missing "ht-solver" sublist parameter')
    end if
    call stop_timer('ht-solver')

    !! Create output file.
    call this%xdmf%open('out')
    call this%xdmf%write_mesh(this%mesh)

    !! Simulation control parameters
    if (params%is_sublist('sim-control')) then
      plist => params%sublist('sim-control')
      context = 'processing ' // plist%name() // ': '
      call plist%get('initial-time', this%t_init, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call plist%get('initial-time-step', this%dt_init, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      if (this%dt_init <= 0.0_r8) call TLS_fatal(context//'"initial-time-step" must be > 0.0')
      call plist%get('min-time-step', this%dt_min, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call plist%get('max-time-step', this%dt_max, default=huge(1.0_r8), stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      call plist%get('max-try-at-step', this%max_try, default=10, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      if (this%dt_min > this%dt_init) call TLS_fatal(context//'require "min-time-step" <= "initial-time-step"')
      call plist%get('output-times', this%tout, stat=stat, errmsg=errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
      !TODO: check for strictly increasing values in TOUT, TOUT > t_init, or sort
      !and cull those < t_init.
    else
      call TLS_fatal('missing "sim-control" sublist parameter')
    end if

    this%ts_sync = time_step_sync(4)

    !! Generate the initial temperature field
    call start_timer('initial-state')
    if (params%is_parameter('initial-temperature')) then
      context = 'processing initial-temperature: '
      call alloc_scalar_func(params, 'initial-temperature', f, stat, errmsg)
      if (stat /= 0) call TLS_fatal(context//errmsg)
    else
      call TLS_fatal('missing "initial-temperature" parameter')
    end if
    allocate(temp(this%mesh%ncell_onP))
    ! call this%mesh%compute_mesh_func(f, this%temp)
    !TODO:replace with general temperature computation
    block
      real(r8), allocatable :: Tface(:)
      allocate(Tface(this%mesh%nface_onP))
      call this%mesh%init_cell_centroid
      call this%mesh%init_face_centroid
      call average_integral(this%disc, f, temp, Tface)
    end block

    !! Define the initial heat conduction state
    !TODO: set max_itr, rel_tol as part of param list input
    max_itr = 100
    rel_tol = 1E-6_r8
    call this%solver%set_initial_state(this%t_init, this%dt_init, temp, rel_tol, max_itr)
    call stop_timer('initial-state')

    call stop_timer('initialization')

  end subroutine init


  subroutine run(this, stat, errmsg)

    class(HT_2d_sim), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    real(r8) :: t, hnext, t_write
    character(80) :: string(2)

    call start_timer('integration')

    !! Write the initial solution
    t = this%solver%time()
    call this%write_solution(t)
    t_write = t ! keep track of the last write time

    call TLS_info('')
    write(string(1),'(a,es12.5)') 'Beginning integration at T = ', t
    call TLS_info(string(1))

    hnext = this%dt_init; this%tlast = t; this%hlast = hnext
    do n = 1, size(this%tout)
      call integrate(this, this%tout(n), hnext, t, stat, errmsg)
      if (stat < 0 .and. t == t_write) exit
      call this%write_solution(t)
      t_write = t ! keep track of the last write time
      call this%solver%write_metrics(string)
      call TLS_info('')
      call TLS_info(string(1))
      call TLS_info(string(2))
      if (stat /= 0) exit
    end do

    if (stat > 0) then  ! caught a signal
      call TLS_info('')
      call TLS_info(errmsg // ': current solution written, and now terminating ...')
      stat = 0  ! this is a successful return
      deallocate(errmsg)
    else if (stat < 0) then
      call TLS_info('')
      errmsg = 'unrecoverable integration failure: ' // errmsg
      call TLS_info(errmsg)
    end if

    call TLS_info('')
    write(string(1),'(a,es12.5,a)') 'Completed integration to T = ', t
    call TLS_info(string(1))

    call this%xdmf%close

    call stop_timer('integration')

  end subroutine run

  !! This integrates the system to the target time TOUT.  The final solution
  !! achieved is returned in T and THIS%U.  Nominally this will be at time
  !! TOUT, however it will be at some earlier time when there is a failure in
  !! the time stepping (or a user interrupt).  The input value of HNEXT is the
  !! initial step size to take, and its return value is the suggested next
  !! step size (nominally the value passed to the next call to INTEGRATE).
  !! STAT returns a negative value if a time stepping failure occurs, with an
  !! explanatory message in ERRMSG.  If the process recieves the SIGURG
  !! signal, the procedure returns at the end of the next time step with a
  !! positive value of STAT, with T and THIS%U set accordingly.

  subroutine integrate(this, tout, hnext, t, stat, errmsg)

    use signal_handler, only: read_signal, SIGURG

    class(HT_2d_sim), intent(inout) :: this
    real(r8), intent(in) :: tout
    real(r8), intent(inout) :: hnext
    real(r8), intent(out) :: t
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    logical :: sig_rcvd

    do
      !! Time for next step; nominally TLAST+HNEXT but possibly adjusted
      t = this%ts_sync%next_time(tout, this%tlast, this%hlast, hnext)
      call step(this, t, hnext, stat, errmsg)
      if (stat /= 0) then
        !! Return the last good solution
        t = this%tlast
        return
      end if
      this%hlast = t - this%tlast
      this%tlast = t

      hnext = min(hnext, this%dt_max)

      call read_signal(SIGURG, sig_rcvd)
      if (sig_rcvd) then
        stat = 1
        errmsg = 'received SIGURG signal'
        return
      end if

      if (t == tout) return
    end do

  end subroutine integrate

  !! Take a resilient step.  Nominally this takes a single step from the
  !! current time to time T and returns the solution in THIS%U.  However if
  !! the step attempt fails, the procedure will re-attempt with successively
  !! smaller step sizes, until the step is successful or the number of
  !! attempts exceeds a maximum or the step size gets too small.  Thus the
  !! return value of T may differ from its input value.  HNEXT returns the
  !! suggested next time step.  STAT returns a negative value if the step
  !! was ultimately unsuccessful, with an explanatory message in ERRMSG.

  subroutine step(this, t,  hnext, stat, errmsg)

    class(HT_2d_sim), intent(inout) :: this
    real(r8), intent(inout) :: t
    real(r8), intent(out) :: hnext
    integer,  intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    real(r8) :: tlast

    tlast = this%solver%time()

    do n = 1, this%max_try
      call this%solver%step(t, hnext, stat)
      if (stat == 0) then ! success
        call this%solver%commit_pending_state
        return
      end if
      t = tlast + hnext
      if (t - tlast < this%dt_min) then
        stat = -1
        errmsg = 'next time step is too small'
        return
      end if
    end do

    stat = -2
    errmsg = 'unable to take a time step'

  end subroutine step

  !! This auxiliary subroutine writes the solution to a GMV format viz file.

  subroutine write_solution(this, t)

    class(HT_2d_sim), intent(inout) :: this
    real(r8), intent(in) :: t

    real(r8), allocatable :: Hcell(:), Tcell(:)

    call start_timer('output')

    allocate(Hcell(this%mesh%ncell_onP), Tcell(this%mesh%ncell_onP))

    call this%solver%get_cell_heat_soln(Hcell)
    call this%solver%get_cell_temp_soln(Tcell)

    call this%xdmf%begin_variables(time=t)
    call this%xdmf%write_cell_var(Hcell, 'H')
    call this%xdmf%write_cell_var(Tcell, 'T')
    call this%xdmf%end_variables

    call stop_timer('output')

  end subroutine write_solution


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TODO: remove the following when possible. It's a copy of test_HT_2d_common.F90

  !! Initializes material database and related objects needed by the HT types
  subroutine init_materials(mesh, matl_db, mmf)

    use parameter_list_type
    use parameter_list_json
    use material_model_driver, only: matl_model
    use material_factory, only: load_material_database
    use material_utilities, only: add_enthalpy_prop

    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(material_database), intent(out) :: matl_db
    type(matl_mesh_func), intent(out) :: mmf

    type(parameter_list), pointer :: plist
    integer, allocatable :: matids(:)
    integer :: stat
    character(:), allocatable :: errmsg, string

    string = '{"unobtanium": &
                {"properties":{"conductivity":1.0, &
                               "density":1.0,&
                               "specific-heat":1.0}}}'

    !! Initialize material database/model with single material
    call parameter_list_from_json_string(string, plist, errmsg)
    call load_material_database(matl_db, plist, stat, errmsg)
    if (stat /= 0) call TLS_FATAL(errmsg)
    call matl_model%init(['unobtanium'], matl_db, stat, errmsg)
    if (stat /= 0) call TLS_FATAL(errmsg)

    !! Initialize enthalpy
    call add_enthalpy_prop(matl_model, stat, errmsg)
    if (stat /= 0) call TLS_FATAL(errmsg)

    !! Layout material across the mesh
    matids = matl_model%matl_index(['unobtanium'])
    call mmf%init(mesh)
    call mmf%define_region(mesh%cell_set_id, matids, stat, errmsg)
    if (stat /= 0) call TLS_FATAL(errmsg)
    call mmf%define_complete(stat, errmsg)
    if (stat /= 0) call TLS_FATAL(errmsg)

  end subroutine init_materials


  !! Computes the average integral of a function over on-process faces and cells
  subroutine average_integral(disc, f, ucell, uface)

    type(mfd_2d_disc), intent(in) :: disc
    class(scalar_func), intent(in) :: f
    real(r8), intent(out) :: ucell(:), uface(:)

    ! real(r8) :: x0(2) = [sqrt(3.0_r8)/3.0_r8, -sqrt(3.0_r8)/3.0_r8]  ! Quadrature points
    ! real(r8) :: w(2) = [1.0_r8, 1.0_r8]   ! Quadrature weights
    real(r8) :: x0(3) = [sqrt(3.0_r8/5.0_r8), 0.0_r8, -sqrt(3.0_r8/5.0_r8)]  ! Quadrature points
    real(r8) :: w(3) = [5.0_r8/9.0_r8, 8.0_r8/9.0_r8, 5.0_r8/9.0_r8]   ! Quadrature weights
    real(r8) :: phi(4,size(x0),size(x0))  ! values of DG shape functions
    real(r8) :: N(4)  ! corner normals
    real(r8) :: val, jac, args(3)
    integer :: j, u, v

    ASSERT(size(ucell) == disc%mesh%ncell_onP)
    ASSERT(size(uface) == disc%mesh%nface_onP)

    args(1) = 0.0_r8  ! f is constant in time

    !! Pre-compute shape functions
    do u = 1, size(x0)
      do v = 1, size(x0)
        phi(1,u,v) = 0.25_r8*(1.0_r8-x0(u))*(1.0_r8-x0(v))
        phi(2,u,v) = 0.25_r8*(1.0_r8+x0(u))*(1.0_r8-x0(v))
        phi(3,u,v) = 0.25_r8*(1.0_r8+x0(u))*(1.0_r8+x0(v))
        phi(4,u,v) = 0.25_r8*(1.0_r8-x0(u))*(1.0_r8+x0(v))
      end do
    end do

    !! Compute cell average
    do j = 1, disc%mesh%ncell_onP
      associate (cface => disc%mesh%cface(disc%mesh%cstart(j):disc%mesh%cstart(j+1)-1), &
                 nodes => disc%mesh%x(:,disc%mesh%cnode(disc%mesh%cstart(j):disc%mesh%cstart(j+1)-1)))
        !! Corner normal lengths
        N(1) = cross_product_2D(nodes(:,2)-nodes(:,1), nodes(:,4)-nodes(:,1))
        N(2) = cross_product_2D(nodes(:,3)-nodes(:,2), nodes(:,1)-nodes(:,2))
        N(3) = cross_product_2D(nodes(:,4)-nodes(:,3), nodes(:,2)-nodes(:,3))
        N(4) = cross_product_2D(nodes(:,1)-nodes(:,4), nodes(:,3)-nodes(:,4))

        !! Gaussian quadrature
        val = 0.0_r8
        do u = 1, size(x0)
          do v = 1, size(x0)
            !! Jacobian
            jac = 0.25_r8*dot_product(N, phi(:,u,v))

            !! Quadrature point
            args(2) = dot_product(nodes(1,:), phi(:,u,v))
            args(3) = dot_product(nodes(2,:), phi(:,u,v))

            val = val + (f%eval(args) * w(u) * w(v) * jac)
          end do
        end do
        ucell(j) = val / disc%mesh%volume(j)
      end associate
    end do

    !! Compute face average
    do j = 1, disc%mesh%nface_onP
      associate (nodes => disc%mesh%x(:,disc%mesh%fnode(:,j)))
        !! Jacobian
        jac = 0.5_r8*disc%mesh%area(j)

        !! Gaussian quadrature
        val = 0.0_r8
        do u = 1, size(x0)
          !! Quadrature point
          args(2:) = 0.5_r8*(sum(nodes, dim=2) + (nodes(:,2)-nodes(:,1))*x0(u))
          val = val + (f%eval(args) * w(u) * jac)
        end do
        uface(j) = val / disc%mesh%area(j)
      end associate
    end do

  end subroutine average_integral


  !TODO: add to cell_geometry?
  !! Computes A x B
  pure function cross_product_2D (a, b) result (axb)
    real(r8), intent(in) :: a(2), b(2)
    real(r8) :: axb
    axb = a(1)*b(2) - a(2)*b(1)
  end function cross_product_2D

end module HT_2d_sim_type
