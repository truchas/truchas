!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module EM

  use,intrinsic :: iso_fortran_env, only: rk => real64
  use parallel_communication
  use EM_data_proxy
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  public :: initialize_EM, induction_heating

  logical :: const_prop = .false.

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! INITIALIZE_EM
 !!
 !! This subroutine performs the initializations needed to calculate the Joule
 !! heat.  This consists of initializing the EM data proxy and defining the
 !! initial Joule heat data.  This data may come from the restart file; if not,
 !! it is computed at time T.  Note that the Joule heat data segment of the
 !! restart file must be dealt with, even if electromagnetics is not enabled.
 !!

  subroutine initialize_EM (t)

    use restart_variables, only: restart
    use restart_driver, only: restart_joule_heat
    use EM_properties
    use electromagnetics_namelist

    real(kind=rk), intent(in) :: t

    logical :: jh_defined
    real(kind=rk) :: s
    character(len=10) :: ss

    if (EM_is_on()) then
      call TLS_info ('')
      call TLS_info ('Initializing electromagnetics ...')
      call start_timer('initialization')
      call init_EM_data_proxy ()
      call stop_timer('initialization')
    end if

    !! Process the EM data segment of the restart file if necessary.
    jh_defined = .false.
    if (restart) call restart_joule_heat (jh_defined)

    if (EM_is_off()) return

    call start_timer('joule heat')

    !! Push the initial source and EM material parameters into the EM data proxy.
    call set_source_properties (t)
    call init_material_properties
    block
      use base_mesh_class
      use mesh_manager, only: named_mesh_ptr
      class(base_mesh), pointer :: mesh
      real(rk), allocatable :: value(:)
      mesh => named_mesh_ptr('MAIN')
      allocate(value(mesh%ncell_onP))
      call EM_permittivity(value)
      call set_permittivity(value)
      call EM_permeability(value)
      call set_permeability(value)
      call EM_conductivity(value)
      call set_conductivity(value)
    end block
    const_prop = have_constant_EM_properties()

    if (jh_defined) then

      !! Determine if the restart Joule heat is usable; if not, compute it.
      if (no_source_field()) then
        if (source_has_changed()) then
          call TLS_info ('   Magnetic source field has changed; restart data not usable.')
          call TLS_info ('  No magnetic source field; setting the Joule heat to zero.')
          call zero_joule_power_density ()
        else
          call TLS_info ('   Using the Joule heat data from the restart file.')
        end if
      else if (material_has_changed()) then
        call TLS_info ('   EM material parameters have changed; restart data not usable.')
        call TLS_info ('  Computing the Joule heat ...')
        call compute_joule_heat(params)
      else if (source_has_changed()) then
        if (source_is_scaled(s)) then
          write(ss,fmt='(es9.3)') s
          call TLS_info ('   Magnetic source field was scaled by ' // trim(ss) // '; Joule heat scaled accordingly.')
          call scale_joule_power_density (s)
        else
          call TLS_info ('   Magnetic source field has changed; restart data not usable.')
          call TLS_info ('  Computing the Joule heat ...')
          call compute_joule_heat(params)
        end if
      else
        call TLS_info ('   Using the Joule heat data from the restart file.')
      end if

    else

      !! Compute the initial Joule heat.
      if (no_source_field()) then
        call TLS_info ('  No magnetic source field; setting the Joule heat to zero.')
        call zero_joule_power_density ()
      else
        call TLS_info ('  Computing the Joule heat ...')
        call compute_joule_heat(params)
      end if

    end if

    !! Write the initial Joule heat to the xml output file.
    call danu_write_joule (t)

    call stop_timer('joule heat')
    call TLS_info ('  electromagnetics initialized')

  end subroutine initialize_EM

  !! This auxillary routine ensures that the EM material properties are
  !! defined for every material phase.  Where necessary it assigns constant
  !! default values in keeping with legacy behavior: 0 for conductivity and
  !! the susceptibilities.  This should be called before attempting to
  !! evaluate the properties on the mesh.

  subroutine init_material_properties

    use material_model_driver, only: matl_model
    use material_utilities

    call define_property_default(matl_model, 'electrical-conductivity', 0.0_rk)
    call define_property_default(matl_model, 'electric-susceptibility', 0.0_rk)
    call define_property_default(matl_model, 'magnetic-susceptibility', 0.0_rk)

  end subroutine init_material_properties

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! INDUCTION_HEATING
 !!
 !! Orchestrates the calculation of the time-averaged Joule heat that will be
 !! used over the time interval [t1, t2].  The result is stored in the EM data
 !! proxy for later retrieval with the access function JOULE_POWER_DENSITY.
 !!

  subroutine induction_heating (t1, t2)

    use EM_properties
    use electromagnetics_namelist, only: params

    real(kind=rk), intent(in) :: t1, t2

    real(kind=rk) :: s
    character(len=10) :: ss

    if (EM_is_off()) return
    call start_timer('joule heat')

    !! Push the source and EM material parameters at time T1 into the EM data proxy.
    call set_source_properties (t1)
    if (.not.const_prop) then
      block
        use base_mesh_class
        use mesh_manager, only: named_mesh_ptr
        class(base_mesh), pointer :: mesh
        real(rk), allocatable :: value(:)
        mesh => named_mesh_ptr('MAIN')
        allocate(value(mesh%ncell_onP))
        call EM_permittivity(value)
        call set_permittivity(value)
        call EM_permeability(value)
        call set_permeability(value)
        call EM_conductivity(value)
        call set_conductivity(value)
      end block
    end if

    if (no_source_field()) then ! there is no joule heat ...
      if (source_has_changed()) then
        call TLS_info (' No magnetic source field; setting the Joule heat to zero.')
        call zero_joule_power_density ()
        call danu_write_joule (t1)
      end if
    else if (matl_has_changed()) then
      call TLS_info (' EM material parameters have changed; computing the Joule heat ...')
      call compute_joule_heat(params)
      call danu_write_joule (t1)
    else if (source_has_changed()) then
      if (source_is_scaled(s)) then
        write(ss,fmt='(es9.3)') s
        call TLS_info (' Magnetic source field was scaled by ' // trim(ss) // '; Joule heat scaled accordingly.')
        call scale_joule_power_density (s)
      else
        call TLS_info (' Magnetic source field has changed; computing the Joule heat ...')
        call compute_joule_heat(params)
      end if
      call danu_write_joule (t1)
    end if

    call stop_timer('joule heat')

  contains

    !! A short-circuit version of .not.const_prop .and. material_has_changed()
    !! The latter is not trivial and we do not want the expense if not needed.

    logical function matl_has_changed()
      if (const_prop) then
        matl_has_changed = .false.
      else
        matl_has_changed = material_has_changed()
      end if
    end function

  end subroutine induction_heating

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  DRIVER FOR THE JOULE HEAT SIMULATION -- SERIAL CODE
 !!

  subroutine compute_joule_heat(params)

    use simpl_mesh_type
    use mimetic_discretization
    use MaxwellEddy
    use EM_boundary_data
    use EM_graphics_output
    use parameter_list_type

    type(parameter_list), intent(inout) :: params

    integer :: j, status, n, cg_max_itr, steps_per_cycle, max_cycles, output_level
    real(kind=rk), pointer :: eps(:), mu(:), sigma(:)
    type(simpl_mesh), pointer :: mesh

    type(system), target :: sys
    real(kind=rk), pointer :: efield(:), bfield(:), q(:), q_avg(:), q_avg_last(:)
    real(kind=rk) :: t, dt, error, eps0, mu0, sigma0, escf, bscf, qscf, freq, curr, etasq, delta, cg_red, ss_tol, num_etasq
    logical :: converged, graphics_output
    character(len=256) :: string
    real(kind=rk) :: eps_min, eps_max, mu_min, mu_max, sigma_min, sigma_max

    !call TLS_info (' Beginning Joule Heat Simulation...')

    call start_timer('simulation')

    !! Get the mesh from the EM data proxy and create the discretization.
    mesh => EM_mesh()

    !! Initialize the boundary conditions.
    call cylinder_bv_init(mesh, params)

    eps   => permittivity()
    mu    => permeability()
    sigma => conductivity()

    if (global_any(eps <= 0.0_rk)) call TLS_fatal ('COMPUTE_JOULE_HEAT: Epsilon is not positive')
    if (global_any(mu  <= 0.0_rk)) call TLS_fatal ('COMPUTE_JOULE_HEAT: Mu is not positive')
    if (global_any(sigma <0.0_rk)) call TLS_fatal ('COMPUTE_JOULE_HEAT: Sigma is not nonnegative')

    eps_min = global_minval(eps)
    eps_max = global_maxval(eps)
    mu_min = global_minval(mu)
    mu_max = global_maxval(mu)
    sigma_min = global_minval(sigma, mask=(sigma > 0.0_rk))
    sigma_max = global_maxval(sigma)

    if (is_IOP) then
      write(string,fmt='(3x,2(a,es11.4))') 'Min epsilon=', eps_min,   ', Max epsilon=', eps_max
      call TLS_info (trim(string))
      write(string,fmt='(3x,2(a,es11.4))') 'Min mu=     ', mu_min,    ', Max mu=     ', mu_max
      call TLS_info (trim(string))
      write(string,fmt='(3x,2(a,es11.4))') 'Min sigma=  ', sigma_min, ', Max sigma=  ', sigma_max
      call TLS_info (trim(string))
    end if

    eps0 = get_epsilon_0()
    mu0 = get_mu_0()
    sigma0 = sigma_max

    if (sigma0 <= 0.0_rk) call TLS_fatal ('COMPUTE_JOULE_HEAT: Sigma is uniformly zero!')

    !coil => get_coil()
    freq = source_frequency()
    curr = 1.0_rk ! multiple coils, so no scaling of the current for now.
    !curr = coil%current

    etasq = eps0 * freq / sigma0
    delta = 1.0_rk / sqrt(mu0 * sigma0 * freq)

    if (is_IOP) then
      write(string,fmt='(3x,a,es11.4)') 'DELTA=', delta
      call TLS_info (trim(string))
      write(string,fmt='(3x,a,es11.4)') 'ETASQ=', etasq
      call TLS_info (trim(string))
    end if

    call params%get('num-etasq', num_etasq)
    if (num_etasq > etasq) then
      etasq = num_etasq
      if (is_IOP) then
        write(string,fmt='(3x,a,es11.4)') 'Using input numerical ETASQ value instead:', etasq
        call TLS_info (trim(string))
      end if
    end if

    !! Scale factors: physical var = scale factor * computational var.
    bscf = mu0 * curr
    escf = curr / sigma0
    qscf = curr**2 / sigma0

    call params%get('steps-per-cycle', steps_per_cycle)
    call params%get('ss-stopping-tolerance', ss_tol)

    !! Create and initialize the time-discretized system.
    dt = 1.0_rk / real(steps_per_cycle,kind=rk)
    call sys%init(mesh, eps, mu, sigma/sigma0, etasq, delta, dt, efield_bc%ebedge, params)

    allocate(efield(mesh%nedge), bfield(mesh%nface))
    allocate(q(mesh%ncell), q_avg(mesh%ncell), q_avg_last(mesh%ncell))

    t = 0.0_rk
    efield = 0.0_rk
    bfield = 0.0_rk
    call sys%set_initial_state(t, efield, bfield)

    q = qscf * sys%joule_heat(efield)

    call params%get('graphics-output', graphics_output)
    if (graphics_output) then
      call export_mesh (mesh, eps, mu, sigma)
!NNC!      block
!NNC!        real(r8), allocatable :: probe_point(:,:)
!NNC!        call params%get('probe-points', probe_point)
!NNC!        call initialize_probes (mesh)
!NNC!        call update_probes (t, escf*efield, bscf*bfield, q*mesh%volume)
!NNC!      end block
    end if

    converged = .false.

    call params%get('maximum-source-cycles', max_cycles)
    STEADY_STATE: do n = 1, max_cycles

      if (converged) exit

      if (graphics_output) then
        call initialize_field_output ()
        if (n == 1) call export_fields (mesh, t, escf*efield, bscf*bfield, q)
      end if

      q_avg = 0.0_rk

      SOURCE_CYCLE: do j = 1, steps_per_cycle

        q_avg = q_avg + 0.5_rk * q

        call sys%step(t, efield, bfield, status, set_bv, bndry_src)
        if (global_any(status /= 0)) exit STEADY_STATE

        q = qscf * sys%joule_heat(efield)
        if (graphics_output) then
          call export_fields (mesh, t, escf*efield, bscf*bfield, q)
!NNC!          call update_probes (t, escf*efield, bscf*bfield, q*mesh%volume)
        end if

        q_avg = q_avg + 0.5_rk * q

      end do SOURCE_CYCLE

      !! Time-averaged Joule power density over the last source cycle.
      q_avg = q_avg / real(steps_per_cycle,kind=rk)

      write(string,fmt='(t4,a,i4,2(a,es11.4))') 'Source cycle', n, &
        ': |Q|_max=', global_maxval(q_avg), ', Q_total=', &
        global_dot_product(q_avg(:mesh%ncell_onP), abs(mesh%volume(:mesh%ncell_onP)))
      call TLS_info (trim(string))

      if (graphics_output) call finalize_field_output (q_avg)

      !! Check for 'convergence' to the steady-state periodic solution.
      if (n > 1) then
        error = global_maxval(abs(q_avg-q_avg_last)) / global_maxval(abs(q_avg))
        !converged = (error < get_ss_stopping_tolerance())
        converged = (error < ss_tol)
      end if

      q_avg_last = q_avg

    end do STEADY_STATE

    !! Check for time step failure.
    if (global_any(status /= 0)) then ! finalize the output and bail.
      if (graphics_output) then
!NNC!        call write_probes (em_output)
      end if
      call TLS_fatal ('COMPUTE_JOULE_HEAT: EM time step failure')
    end if

    if (graphics_output) then
!NNC!      call write_probes (em_output)
    end if

    !! Check for convergence to steady state.  We continue in any case.
    if (.not.converged) then
      call TLS_warn ('COMPUTE_JOULE_HEAT: Not converged to steady-state; proceding anyway.')
    end if

    call stop_timer('simulation')

    !! Store the the computed Joule heat in the EM data proxy for later retrieval.
    call set_joule_power_density (q_avg)

    call TLS_info ('  Joule heat computation completed.')

  end subroutine compute_joule_heat

end module EM
