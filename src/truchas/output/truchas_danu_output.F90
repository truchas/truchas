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
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module truchas_danu_output

  use truchas_danu_output_data
  use truchas_danu_output_tools
  use,intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR, c_associated
  use parallel_communication
  use truchas_logging_services
  use danu_module
  use kinds, only: r8
  implicit none
  private
  
  public :: TDO_open, TDO_close
  public :: TDO_write_default_mesh
  public :: TDO_start_simulation
  public :: TDO_write_timestep
  
  type(c_ptr), save :: mid = C_NULL_PTR ! Danu mesh id
  type(c_ptr), save :: seq_id = C_NULL_PTR ! Danu sequence id
  
  public :: fid ! others may want to write here
  
contains

  subroutine TDO_open ()
    use truchas_env, only: output_file_name
    integer :: stat
    if (is_IOP) then
      call TLS_info ('DANU: Opening h5 output file')
      call output_file_create (output_file_name('h5'), fid, stat)
      INSIST(c_associated(fid))
    end if
    call broadcast (stat)
    INSIST(stat == DANU_SUCCESS)
  end subroutine TDO_open
  
  subroutine TDO_close ()
    if (c_associated(fid)) then
      call TLS_info ('DANU: Closing h5 output file')
      call output_file_close (fid)
      if (c_associated(fid)) call TLS_warn ('DANU: Danu fid is still associated')
    end if
  end subroutine TDO_close
  
  subroutine TDO_write_default_mesh
  
    use kinds, only: r8
    use parameter_module, only: ndim, nvc, ncells, ncells_tot, nnodes, nnodes_tot
    use mesh_module, only: vertex, mesh, unpermute_mesh_vector, unpermute_vertex_vector, mesh_has_cblockid_data
    use truchas_logging_services
    
    integer :: stat, k
    real(r8), pointer :: x(:), y(:), z(:)
    integer, pointer :: cnode(:,:), iarray(:)
    
    !! Create the mesh entry.
    call TLS_info ('DANU: adding default mesh entry')
    if (is_IOP) then
      INSIST(c_associated(fid))
      call mesh_add_unstructured (fid, 'DEFAULT', nvc, ndim, mid, stat)
    end if
    call broadcast (stat)
    INSIST(stat == DANU_SUCCESS)
    
    !! Write the node coordinates.
    INSIST(ndim == 3)
    call allocate_collated_array (x, nnodes_tot)
    call collate (x, vertex(:)%coord(1))
    call allocate_collated_array (y, nnodes_tot)
    call collate (y, vertex(:)%coord(2))
    call allocate_collated_array (z, nnodes_tot)
    call collate (z, vertex(:)%coord(3))
    call TLS_info ('DANU: writing mesh node coordinates')
    if (is_IOP) then
      call mesh_write_coordinates (mid, nnodes_tot, x, y, z, stat)
    end if
    deallocate(x, y, z)
    call broadcast (stat)
    INSIST(stat == DANU_SUCCESS)
    
    !! Write the cell connectivity.
    call allocate_collated_array (cnode, nvc, ncells_tot)
    do k = 1, nvc
      call collate (cnode(k,:), mesh%ngbr_vrtx_orig(k)) ! then internal serial numbering
    end do
    if (is_IOP) then
      call mesh_write_connectivity (mid, ncells_tot, cnode, stat)
    end if
    deallocate(cnode)
    call broadcast (stat)
    INSIST(stat == DANU_SUCCESS)
    
    ! I don't know where this should go
    call TDO_start_simulation
    
    !! Right now the following is being written to the non-series section
    !! of the simulation.  Soon we will write these in the mesh section.
    
    !! Mapping from internal serial cell numbering to external numbering.
    call allocate_collated_array (iarray, ncells_tot)
    call collate (iarray, unpermute_mesh_vector)
    if (is_iOP) call data_write (sid, 'CELLMAP', iarray, stat)
    call broadcast (stat)
    INSIST(stat == DANU_SUCCESS)
    
    !! Cell block IDs.
    if (mesh_has_cblockid_data) then
      call collate (iarray, mesh%cblockid)
      if (is_iOP) call data_write (sid, 'BLOCKID', iarray, stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
    end if
    
    !! Cell partition assignment.
    if (nPE > 1) then
      call collate (iarray, spread(this_PE, dim=1, ncopies=ncells))
      if (is_IOP) call data_write (sid, 'CELLPART', iarray, stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
    end if
    deallocate(iarray)
    
    !! Mapping from internal serial node numbering to external numbering.
    call allocate_collated_array (iarray, nnodes_tot)
    call collate (iarray, unpermute_vertex_vector)
    if (is_iOP) call data_write (sid, 'NODEMAP', iarray, stat)
    
    !! Node partition assignment.
    if (nPE > 1) then
      call collate (iarray, spread(this_PE, dim=1, ncopies=nnodes))
      if (is_IOP) call data_write (sid, 'NODEPART', iarray, stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
    end if
    deallocate(iarray)
    
    if (is_IOP) call data_write (sid, 'NUMPROCS', nPE, stat)
    call broadcast (stat)
    INSIST(stat == DANU_SUCCESS)
    
  end subroutine TDO_write_default_mesh
  
  subroutine TDO_start_simulation
  
    use physics_module, only: species_transport, heat_species_transport, number_of_species

    integer :: stat
  
    call TLS_info ('DANU: adding main simulation entry')
    if (is_IOP) then
      INSIST(c_associated(fid))
      INSIST(.not.c_associated(sid))
      call simulation_add (fid, 'MAIN', sid, stat)
    end if
    call broadcast (stat)
    INSIST(stat == DANU_SUCCESS)

    !! Hard-wired for the moment. Need to accomodate EM simulation 
    if (is_IOP) then
      INSIST(c_associated(fid))
      INSIST(c_associated(sid))
      call simulation_link_mesh(fid,sid,'DEFAULT',stat)
    end if
    call broadcast (stat)
    INSIST(stat == DANU_SUCCESS)
    
    if (species_transport .or. heat_species_transport) then
      if (is_IOP) call attribute_write (sid, 'NUM_SPECIES', number_of_species, stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
    end if
    
  end subroutine TDO_start_simulation
  
  subroutine TDO_write_timestep
  
    use time_step_module, only: t, dt, cycle_number
    use fluid_data_module, only: fluid_flow
    use physics_module, only: heat_transport, species_transport, heat_species_transport
    use EM_data_proxy, only: EM_is_on
    use solid_mechanics_input, only: solid_mechanics
    use gap_output, only: set_gap_element_output
    use ustruc_driver, only: ustruc_output
    
    integer :: stat
  
    if (is_IOP) then
      INSIST(c_associated(sid))
      call sequence_next_id (sid, cycle_number, t, seq_id, stat)
    end if
    call broadcast (stat)
    INSIST(stat == DANU_SUCCESS)
    if (is_IOP) call attribute_write (seq_id, 'time step', dt, stat)
    call broadcast (stat)
    INSIST(stat == DANU_SUCCESS)
    
    !! Prior to writing, overwrite (bogus) field data on gap elements with
    !! something reasonable.  This doesn't belong here and should be moved.
    !! Currently commented out, because the TBrook output does this.  When
    !! the TBrook output is disabled, this should be uncommented.
    call set_gap_element_output

    !! Cell density, temperature, enthalpy, phase volume fractions.
    call write_common_data
    
    !! Flow-related fields.
    if (fluid_flow) call write_fluid_flow_data
    
    !! Heat transfer fields (other than temperature).
    if (heat_transport .or. heat_species_transport) call write_heat_transfer_data
    
    !! Induction heating fields.
    if (EM_is_on()) call write_EM_data
    
    !! Solid mechanics fields.
    if (solid_mechanics) call write_solid_mech_data
    
    !! Species fields.
    if (species_transport .or. heat_species_transport) call write_species_data

    !! Microstructure analysis data (if enabled)
    call ustruc_output (seq_id)
    
  contains
  
    subroutine write_common_data

      use parameter_module, only: ndim, ncells, nmat
      use mesh_module, only: cell
      use zone_module, only: zone
      use property_module, only: get_density, get_user_material_id
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
      call get_density (zone%temp, rho)
      call write_seq_cell_field (seq_id, rho, 'Z_RHO', for_viz=.true., viz_name='Density')
      deallocate(rho)
    
      !! Cell temperature
      call write_seq_cell_field (seq_id, zone%temp, 'Z_TEMP', for_viz=.true., viz_name='T')

      !! Average cell enthalpy density
      call write_seq_cell_field (seq_id, zone%enthalpy, 'Z_ENTHALPY', for_viz=.true., viz_name='Enthalpy')

      !! Phase volume fractions
      if (nmat > 1) then
        allocate(vof(nmat,ncells), name(nmat))
        do m = 1, nmat
          call gather_vof (m, vof(m,:))
          write(name(m),'(a,i4.4)') 'VOF', get_user_material_id(m)
        end do
        call write_seq_cell_field (seq_id, vof, 'VOF', for_viz=.true., viz_name=name)
        deallocate(vof, name)
      end if
      
      !! Cell centroids
      allocate(xc(ndim,ncells))
      do j = 1, ncells
        xc(:,j) = cell(j)%centroid
      end do
      call write_seq_cell_field (seq_id, xc, 'CENTROID', for_viz=.true., viz_name=['XC', 'YC', 'ZC'])
      deallocate(xc)

    end subroutine write_common_data
  
    subroutine write_fluid_flow_data
    
      use parameter_module, only: ndim, ncells
      use zone_module, only: zone
      use fluid_data_module, only: fluxing_velocity, courant, boussinesq_approximation
      use property_module, only: get_density_delta
      use diagnostics_module, only: divergence
      
      integer :: n
      real(r8), allocatable :: vcell(:,:), div(:), drho(:)
    
      !! Cell-centered fluid velocity.
      allocate(vcell(ndim,ncells))
      do n = 1, ndim
        vcell(n,:) = zone%vc(n) ! work around flawed data structure
      end do
      call write_seq_cell_field (seq_id, vcell, 'Z_VC', for_viz=.true., viz_name=['U','V','W'])
      deallocate(vcell)
      
      !! Cell-centered fluid pressure.
      call write_seq_cell_field (seq_id, zone%p, 'Z_P', for_viz=.true., viz_name='P')
      
      !! Face fluxing velocities.
      call write_seq_cell_field (seq_id, fluxing_velocity, 'Face_Vel', for_viz=.false.)
      
      !! Cell-centered fluid Courant number.
      call write_seq_cell_field (seq_id, courant, 'COURANT', for_viz=.true.)
      
      !! Cell-centered divergence (the volume error).
      allocate(div(ncells))
      call divergence (div)
      call write_seq_cell_field (seq_id, div, 'Volume_Error', for_viz=.true., viz_name='vol_err')
      deallocate(div)
      
      !! Cell-centered fluid density delta.
      if (boussinesq_approximation) then
        allocate(drho(ncells))
        call get_density_delta (zone%temp, drho)
        call write_seq_cell_field (seq_id, drho, 'del-rho', for_viz=.true., viz_name='delrho')
        deallocate(drho)
      end if

    end subroutine write_fluid_flow_data
    
    subroutine write_heat_transfer_data
    
      use parameter_module, only: ndim, ncells
      use zone_module, only: zone
      use time_step_module, only: dt
      use diffusion_solver, only: ds_get_temp_grad
    
      real(r8) :: dTdt(ncells), gradT(ndim,ncells)
      
      dTdt = (zone%temp - zone%temp_old) / dt
      call write_seq_cell_field (seq_id, dTdt, 'dTdt', for_viz=.true., viz_name='dT/dt')
      
      call ds_get_temp_grad (gradT)
      call write_seq_cell_field (seq_id, gradT, 'Grad_T', for_viz=.true., viz_name=['dT/dx','dT/dy','dT/dz'])

    end subroutine write_heat_transfer_data
    
    subroutine write_EM_data
    
      use EM_data_proxy, only: joule_power_density
      
      real(r8), pointer :: q(:)
      
      q => joule_power_density()
      call write_seq_cell_field (seq_id, q, 'Joule_P', for_viz=.true.)

    end subroutine write_EM_data
  
    subroutine write_solid_mech_data
    
      use parameter_module, only: ndim, nnodes, ncomps, ncells
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
        call write_seq_cell_field (seq_id, scratch2, name, for_viz=.false.)
        call get_smech_ip_elastic_stress(n,scratch2)
        write(name,'(a,i2.2)') 'ELASTIC_STRESS_', n
        call write_seq_cell_field (seq_id, scratch2, name, for_viz=.false.)
        call get_smech_ip_plastic_strain(n,scratch2)
        write(name,'(a,i2.2)') 'PLASTIC_STRAIN_', n
        call write_seq_cell_field (seq_id, scratch2, name, for_viz=.false.)
        call get_smech_ip_plastic_strain_rate(n,scratch1)
        write(name,'(a,i2.2)') 'PLASTIC_STRAIN_RATE_', n
        call write_seq_cell_field (seq_id, scratch1, name, for_viz=.false.)
      end do
      deallocate(scratch1)
      deallocate(scratch2)
      
      !! More restart-only data
      allocate(scratch2(ndim,nnodes))
      call get_sm_rhs(scratch2)
      call write_seq_node_field (seq_id, scratch2, 'RHS', for_viz=.false.)
      deallocate(scratch2)
      
      !! Restart and viz data (cell based).
      allocate(scratch2(ncomps,ncells))
      call get_smech_cell_elastic_stress(scratch2)
      call write_seq_cell_field (seq_id, scratch2, 'sigma', for_viz=.true., &
          viz_name=['sigxx', 'sigyy', 'sigzz', 'sigxy', 'sigxz', 'sigyz'])
      call get_smech_cell_total_strain(scratch2)
      call write_seq_cell_field (seq_id, scratch2, 'epsilon', for_viz=.true., &
          viz_name=['epsxx', 'epsyy', 'epszz', 'epsxy', 'epsxz', 'epsyz'])
      call get_smech_cell_plastic_strain(scratch2)
      call write_seq_cell_field (seq_id, scratch2, 'e_plastic', for_viz=.true., &
          viz_name=['eplxx', 'eplyy', 'eplzz', 'eplxy', 'eplxz', 'eplyz'])
      call get_sm_thermal_strain(scratch2)
      call write_seq_cell_field (seq_id, scratch2, 'epstherm', for_viz=.true., &
          viz_name=['epsthxx', 'epsthyy', 'epsthzz', 'epsthxy', 'epsthxz', 'epsthyz'])
      call get_sm_pc_strain(scratch2)  
      call write_seq_cell_field (seq_id, scratch2, 'epspc', for_viz=.true., &
          viz_name=['epspcxx', 'epspcyy', 'epspczz', 'epspcxy', 'epspcxz', 'epspcyz'])
      deallocate(scratch2)
      allocate(scratch1(ncells))
      call get_smech_cell_plastic_strain_rate(scratch1)  
      call write_seq_cell_field (seq_id, scratch1, 'epsdot', for_viz=.true.)
      
      !! Restart and viz data (node based)
      allocate(scratch2(ndim,nnodes))
      call get_sm_displacement(scratch2)
      call write_seq_node_field (seq_id, scratch2, 'Displacement', for_viz=.true., &
          viz_name=['Dx', 'Dy', 'Dz'])
      deallocate(scratch2)  

 
      !! Viz data only.
      call get_sm_rotation_magnitude(scratch1)
      call write_seq_cell_field (seq_id, scratch1, 'Rotation', for_viz=.true.)
     
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
          call write_seq_node_field (seq_id, node_gap, name, for_viz=.true.)
          write(name,'(a,i2.2)') 'NTRAC_', interface_list(n)
          call write_seq_node_field (seq_id, node_norm_trac, name, for_viz=.true.)
        end if
      end do
      deallocate(node_gap)
      deallocate(node_norm_trac)

    end subroutine write_solid_mech_data
    
    subroutine write_species_data
    
      use parameter_module, only: ncells
      use diffusion_solver_data, only: num_species
      use diffusion_solver, only: ds_get_phi
      use string_utilities, only: i_to_c
      
      integer :: n
      real(r8) :: array(ncells)
      
      do n = 1, num_species
        call ds_get_phi (n, array)
        call write_seq_cell_field (seq_id, array, 'phi'//i_to_c(n), for_viz=.true.)
      end do

    end subroutine write_species_data
    
  end subroutine TDO_write_timestep

end module truchas_danu_output
