!!
!! TRUCHAS_DANU_OUTPUT
!!
!! Neil N. Carlson <nnc@lanl.gov> Apr 2012
!!
!! Prototype module for developing new HDF5 output using the Danu package.
!! This is organized very much like Tbrook_Utilities so that we can more
!! easily drop this in as a replacement -- just the first step in reworking the
!! output.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module truchas_danu_output

  use truchas_danu_output_data
  use truchas_danu_output_tools
  use,intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR, c_associated
  use parallel_communication
  use truchas_logging_services
  use truchas_h5_outfile, only: th5_mesh_group, th5_seq_group
  use kinds, only: r8
  implicit none
  private

  public :: TDO_open, TDO_close
  public :: TDO_write_default_mesh
  public :: TDO_start_simulation
  public :: TDO_write_timestep

  type(th5_mesh_group), save :: out_mesh
  type(th5_seq_group), save :: seq

  public :: outfile ! others may want to write here

contains

  subroutine TDO_open ()
    use truchas_env, only: output_file_name
    call outfile%open (output_file_name('h5'), io_group_size, is_IOP)
  end subroutine TDO_open

  subroutine TDO_close ()
    call outfile%close()
  end subroutine TDO_close

  subroutine TDO_write_default_mesh

    use kinds, only: r8
    use legacy_mesh_api, only: ndim, nvc, ncells, ncells_tot, nnodes, nnodes_tot
    use legacy_mesh_api, only: vertex, mesh, unpermute_mesh_vector, unpermute_vertex_vector, mesh_has_cblockid_data
    use output_control,  only: part
    use truchas_logging_services

    integer :: k
    real(r8), allocatable :: x(:,:)
    integer, allocatable :: cnode(:,:)

    !! Create the mesh entry.
    call outfile%add_unstr_mesh_group('DEFAULT', nvc, ndim, out_mesh)

    !! Write the node coordinates.
    allocate(x(3,nnodes))
    x(1,:) = vertex%coord(1); x(2,:) = vertex%coord(2); x(3,:) = vertex%coord(3)
    call out_mesh%write_coordinates(nnodes_tot, x)
    deallocate(x)

    !! Write the cell connectivity.
    allocate(cnode(nvc,ncells))
    do k = 1, nvc
      cnode(k,:) = mesh%ngbr_vrtx_orig(k)
    end do
    call out_mesh%write_connectivity(ncells_tot, cnode)
    deallocate(cnode)

    ! I don't know where this should go
    call TDO_start_simulation

    !! Right now the following is being written to the non-series section
    !! of the simulation.  Soon we will write these in the mesh section.

    !! Mapping from internal serial cell numbering to external numbering.
    call sim%write_dist_array('CELLMAP', unpermute_mesh_vector, ncells_tot)

    !! Cell block IDs.
    call sim%write_dist_array('BLOCKID', mesh%cblockid, ncells_tot)

    !! Cell partition assignment.
    call sim%write_dist_array('CELLPART', spread(this_PE,dim=1,ncopies=ncells), ncells_tot)

    !! Mapping from internal serial node numbering to external numbering.
    call sim%write_dist_array('NODEMAP', unpermute_vertex_vector, nnodes_tot)

    !! Node partition assignment.
    call sim%write_dist_array('NODEPART', spread(this_PE,dim=1,ncopies=nnodes), nnodes_tot)

    call sim%write_repl_data('NUMPROCS', nPE)

    !! Parts for movement
    if (size(part) > 0) call sim%write_repl_data('part1', part)

  end subroutine TDO_write_default_mesh

  subroutine TDO_start_simulation

    use physics_module, only: species_transport, number_of_species

    call outfile%add_sim_group('MAIN', sim)

    !! Hard-wired for the moment. Need to accomodate EM simulation
    call sim%add_mesh_link('DEFAULT')

    if (species_transport) then
      call sim%write_attr('NUM_SPECIES', number_of_species)
    end if

  end subroutine TDO_start_simulation

  subroutine TDO_write_timestep

    use time_step_module, only: t, dt, cycle_number
    use physics_module, only: heat_transport, species_transport
    use EM_data_proxy, only: EM_is_on
    use solid_mechanics_input, only: solid_mechanics
    use gap_output, only: set_gap_element_output
    use ustruc_driver, only: ustruc_output
    use flow_driver, only: flow_enabled
    use output_control, only: part_path, write_mesh_partition

    integer :: stat
    real(r8) :: r(3)

    call sim%next_seq_group(cycle_number, t, seq)
    call seq%write_attr('time step', dt)

    !! Part movement
    if (associated(part_path)) then
      call part_path%set_segment(t)
      call part_path%get_position(t, r)
      call seq%write_attr('translate_part1', r)
    end if

    !! Prior to writing, overwrite (bogus) field data on gap elements with
    !! something reasonable.  This doesn't belong here and should be moved.
    !! Currently commented out, because the TBrook output does this.  When
    !! the TBrook output is disabled, this should be uncommented.
    call set_gap_element_output

    !! Cell density, temperature, enthalpy, phase volume fractions.
    call write_common_data

    !! Flow-related fields.
    if (flow_enabled()) call write_new_flow_data

    !! Heat transfer fields (other than temperature and enthalpy).
    if (heat_transport) call write_heat_transfer_data

    !! Induction heating fields.
    if (EM_is_on()) call write_EM_data

    !! Solid mechanics fields.
    if (solid_mechanics) call write_solid_mech_data

    !! Species fields.
    if (species_transport) call write_species_data

    !! Microstructure analysis data (if enabled)
    call ustruc_output (seq)

  contains

    subroutine write_common_data

      use parameter_module, only: nmat
      use legacy_mesh_api, only: ndim, ncells, cell
      use zone_module, only: zone
      use matl_module, only: gather_vof

      integer :: j, m, stat
      real(r8), allocatable :: rho(:), vof(:,:), xc(:,:)
      character(8), allocatable :: name(:)

      !! Average cell density
      !! In TBU_WriteTimeStepData, DAK switched (2/4/09) to using the value
      !! computed like that done by the call to density below, instead of just
      !! using the value in zone%rho; not sure why, as the computation is the
      !! same, but perhaps the zone%rho value is stale?  NNC, 8/9/2012.
      allocate(rho(ncells))
      call write_seq_cell_field (seq, zone%rho, 'Z_RHO', for_viz=.true., viz_name='Density')
      deallocate(rho)

      !! Cell temperature
      call write_seq_cell_field (seq, zone%temp, 'Z_TEMP', for_viz=.true., viz_name='T')

      !! Average cell enthalpy density
      call write_seq_cell_field (seq, zone%enthalpy, 'Z_ENTHALPY', for_viz=.true., viz_name='Enthalpy')

      !! Phase volume fractions
      if (nmat > 1) then
        allocate(vof(nmat,ncells), name(nmat))
        do m = 1, nmat
          call gather_vof (m, vof(m,:))
          write(name(m),'(a,i4.4)') 'VOF', m  !TODO: incorporate material name
        end do
        call write_seq_cell_field (seq, vof, 'VOF', for_viz=.true., viz_name=name)
        deallocate(vof, name)
      end if

      !! Cell centroids
      !allocate(xc(ndim,ncells))
      !do j = 1, ncells
      !  xc(:,j) = cell(j)%centroid
      !end do
      !call write_seq_cell_field (seq, xc, 'CENTROID', for_viz=.true., viz_name=['XC', 'YC', 'ZC'])
      !deallocate(xc)

      !! This is a stop-gap because the Truchas paraview reader ignores this
      !! data written with the mesh. It is only written to the first snapshot.
      if (write_mesh_partition) then
        call write_seq_cell_field(seq, spread(real(this_PE,r8),dim=1,ncopies=ncells), 'CELLPART', for_viz=.true., viz_name='rank')
        write_mesh_partition = .false.
      end if

    end subroutine write_common_data

    subroutine write_new_flow_data
      use legacy_mesh_api, only: ndim, ncells, nfc
      use flow_driver

      real(r8), pointer :: vec_cc(:,:), scalar_cc(:)
      real(r8) :: fluxing_velocity(nfc,ncells)

      vec_cc => flow_vel_cc_view()
      call write_seq_cell_field(seq, vec_cc(:,1:ncells), 'Z_VC', for_viz=.true., viz_name=['U','V','W'])

      scalar_cc => flow_P_cc_view()
      call write_seq_cell_field (seq, scalar_cc(1:ncells), 'Z_P', for_viz=.true., viz_name='P')

      call get_legacy_flux_vel(fluxing_velocity)
      call write_seq_cell_field(seq, fluxing_velocity, 'Face_Vel', for_viz=.false.)

      !call flow_driver_dump_state

!!$
!!$      !! Cell-centered fluid Courant number.
!!$      call write_seq_cell_field (seq, courant, 'COURANT', for_viz=.true.)
!!$
!!$      !! Cell-centered divergence (the volume error).
!!$      allocate(div(ncells))
!!$      call divergence (div)
!!$      call write_seq_cell_field (seq, div, 'Volume_Error', for_viz=.true., viz_name='vol_err')
!!$      deallocate(div)
!!$
!!$      !! Cell-centered fluid density delta.
!!$      if (boussinesq_approximation) then
!!$        allocate(drho(ncells))
!!$        call get_density_delta (zone%temp, drho)
!!$        call write_seq_cell_field (seq, drho, 'del-rho', for_viz=.true., viz_name='delrho')
!!$        deallocate(drho)
!!$      end if

    end subroutine write_new_flow_data

    subroutine write_heat_transfer_data

      use legacy_mesh_api, only: ndim, ncells
      use zone_module, only: zone
      use time_step_module, only: dt
      use diffusion_solver, only: ds_get_temp_grad

      real(r8) :: dTdt(ncells), gradT(ndim,ncells)

      dTdt = (zone%temp - zone%temp_old) / dt
      call write_seq_cell_field (seq, dTdt, 'dTdt', for_viz=.true., viz_name='dT/dt')

      call ds_get_temp_grad (gradT)
      call write_seq_cell_field (seq, gradT, 'Grad_T', for_viz=.true., viz_name=['dT/dx','dT/dy','dT/dz'])

    end subroutine write_heat_transfer_data

    subroutine write_EM_data

      use EM_data_proxy, only: joule_power_density

      real(r8), pointer :: q(:)

      q => joule_power_density()
      call write_seq_cell_field (seq, q, 'Joule_P', for_viz=.true.)

    end subroutine write_EM_data

    subroutine write_solid_mech_data

      use parameter_module, only: ncomps
      use legacy_mesh_api, only: ndim, nnodes, ncells
      use solid_mechanics_output, only: get_sm_displacement, get_sm_thermal_strain,get_sm_rhs, &
          get_sm_rotation_magnitude, get_sm_pc_strain, get_smech_cell_total_strain, &
          get_smech_cell_elastic_stress, get_smech_cell_plastic_strain, &
          get_smech_cell_plastic_strain_rate, smech_num_int_pts, get_smech_ip_total_strain, &
          get_smech_ip_elastic_stress, get_smech_ip_plastic_strain, &
          get_smech_ip_plastic_strain_rate, sm_node_gap_isize, get_sm_node_gap, &
          get_sm_node_norm_trac
      use mech_bc_data_module, only: interface_list

      integer :: n, nipc, isize
      character(32) :: name
      real(r8), allocatable :: scratch1(:), scratch2(:,:),  node_gap(:), node_norm_trac(:)

      !! This is restart data only.
      allocate(scratch1(ncells))
      allocate(scratch2(ncomps,ncells))
      nipc = smech_num_int_pts()
      do n = 1, nipc
        call get_smech_ip_total_strain(n,scratch2)
        write(name,'(a,i2.2)') 'TOTAL_STRAIN_', n
        call write_seq_cell_field (seq, scratch2, trim(name), for_viz=.false.)
        call get_smech_ip_elastic_stress(n,scratch2)
        write(name,'(a,i2.2)') 'ELASTIC_STRESS_', n
        call write_seq_cell_field (seq, scratch2, trim(name), for_viz=.false.)
        call get_smech_ip_plastic_strain(n,scratch2)
        write(name,'(a,i2.2)') 'PLASTIC_STRAIN_', n
        call write_seq_cell_field (seq, scratch2, trim(name), for_viz=.false.)
        call get_smech_ip_plastic_strain_rate(n,scratch1)
        write(name,'(a,i2.2)') 'PLASTIC_STRAIN_RATE_', n
        call write_seq_cell_field (seq, scratch1, trim(name), for_viz=.false.)
      end do
      deallocate(scratch1)
      deallocate(scratch2)

      !! More restart-only data
      allocate(scratch2(ndim,nnodes))
      call get_sm_rhs(scratch2)
      call write_seq_node_field (seq, scratch2, 'RHS', for_viz=.false.)
      deallocate(scratch2)

      !! Restart and viz data (cell based).
      allocate(scratch2(ncomps,ncells))
      call get_smech_cell_elastic_stress(scratch2)
      call write_seq_cell_field (seq, scratch2, 'sigma', for_viz=.true., &
          viz_name=['sigxx', 'sigyy', 'sigzz', 'sigxy', 'sigxz', 'sigyz'])
      call get_smech_cell_total_strain(scratch2)
      call write_seq_cell_field (seq, scratch2, 'epsilon', for_viz=.true., &
          viz_name=['epsxx', 'epsyy', 'epszz', 'epsxy', 'epsxz', 'epsyz'])
      call get_smech_cell_plastic_strain(scratch2)
      call write_seq_cell_field (seq, scratch2, 'e_plastic', for_viz=.true., &
          viz_name=['eplxx', 'eplyy', 'eplzz', 'eplxy', 'eplxz', 'eplyz'])
      call get_sm_thermal_strain(scratch2)
      call write_seq_cell_field (seq, scratch2, 'epstherm', for_viz=.true., &
          viz_name=['epsthxx', 'epsthyy', 'epsthzz', 'epsthxy', 'epsthxz', 'epsthyz'])
      call get_sm_pc_strain(scratch2)
      call write_seq_cell_field (seq, scratch2, 'epspc', for_viz=.true., &
          viz_name=['epspcxx', 'epspcyy', 'epspczz', 'epspcxy', 'epspcxz', 'epspcyz'])
      deallocate(scratch2)
      allocate(scratch1(ncells))
      call get_smech_cell_plastic_strain_rate(scratch1)
      call write_seq_cell_field (seq, scratch1, 'epsdot', for_viz=.true.)

      !! Restart and viz data (node based)
      allocate(scratch2(ndim,nnodes))
      call get_sm_displacement(scratch2)
      call write_seq_node_field (seq, scratch2, 'Displacement', for_viz=.true., &
          viz_name=['Dx', 'Dy', 'Dz'])
      deallocate(scratch2)


      !! Viz data only.
      call get_sm_rotation_magnitude(scratch1)
      call write_seq_cell_field (seq, scratch1, 'Rotation', for_viz=.true.)

      deallocate(scratch1)

      !! Gap displacements and scaled forces.  For each 'interface' (likely a
      !! gap element block?) there is an entire node-based field for the gap
      !! displacement and normal traction, with non-zero values only along the
      !! actual interface.  This REALLY needs to be re-implemented using a more
      !! efficient design.  This is viz data only.

      ! Both node_gap and node_norm_trac are nnodes x isize
      isize = sm_node_gap_isize()
      allocate(node_gap(nnodes))
      allocate(node_norm_trac(nnodes))
      do n = 1, isize
        call get_sm_node_gap(n,node_gap)
        call get_sm_node_norm_trac(n,node_norm_trac)
        if (global_any((node_gap /= 0.0_r8) .or. (node_norm_trac /= 0.0_r8))) then
          write(name,'(a,i2.2)') 'GAP_', interface_list(n)
          call write_seq_node_field (seq, node_gap, trim(name), for_viz=.true.)
          write(name,'(a,i2.2)') 'NTRAC_', interface_list(n)
          call write_seq_node_field (seq, node_norm_trac, trim(name), for_viz=.true.)
        end if
      end do
      deallocate(node_gap)
      deallocate(node_norm_trac)

    end subroutine write_solid_mech_data

    subroutine write_species_data

      use legacy_mesh_api, only: ncells
      use diffusion_solver_data, only: num_species
      use diffusion_solver, only: ds_get_phi
      use string_utilities, only: i_to_c

      integer :: n
      real(r8) :: array(ncells)

      do n = 1, num_species
        call ds_get_phi (n, array)
        call write_seq_cell_field (seq, array, 'phi'//i_to_c(n), for_viz=.true.)
      end do

    end subroutine write_species_data

  end subroutine TDO_write_timestep

end module truchas_danu_output
