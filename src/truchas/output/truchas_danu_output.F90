!!
!! TRUCHAS_DANU_OUTPUT
!!
!! Neil N. Carlson <nnc@lanl.gov> Apr 2012
!!
!! Prototype module for developing new HDF5 output using the Danu package.
!! This is organized very much like Tbrook_Utilities so that we can more
!! easily drop this in as a replacement -- just the first step in reworking the
!! output.  Compiles to an empty module unless the macro USE_DANU is defined.
!!

#include "f90_assert.fpp"

module truchas_danu_output

#ifdef USE_DANU

  use truchas_danu_output_data
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
  
#ifdef PATHSCALE_COMPILER_WORKAROUND
  type(c_ptr), save :: mid ! Danu mesh id
  type(c_ptr), save :: seq_id ! Danu sequence id
#else
  !moved to truchas_danu_output_data! type(c_ptr), save :: fid = C_NULL_PTR ! h5 file id
  type(c_ptr), save :: mid = C_NULL_PTR ! Danu mesh id
  !moved to truchas_danu_output_data! type(c_ptr), save :: sid = C_NULL_PTR ! Danu simulation id
  type(c_ptr), save :: seq_id = C_NULL_PTR ! Danu sequence id
#endif
  
  public :: fid ! others may want to write here
  
  interface write_cell_field
    module procedure write_cell_field_r0, write_cell_field_r1
  end interface
  
  interface write_node_field
    module procedure write_node_field_r0, write_node_field_r1
  end interface
  
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
    use solid_mechanics_data, only: solid_mechanics
    !use gap_output, only: set_gap_element_output
    
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
    !call set_gap_element_output

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
    
  contains
  
    subroutine write_common_data

      use parameter_module, only: ncells, nmat
      use zone_module, only: zone
      use property_module, only: get_density, get_user_material_id
      use matl_module, only: gather_vof

      integer :: m, stat
      real(r8), allocatable :: rho(:), vof(:,:)
      character(8), allocatable :: name(:)

      !! Average cell density
      !! In TBU_WriteTimeStepData, DAK switched (2/4/09) to using the value
      !! computed like that done by the call to density below, instead of just
      !! using the value in zone%rho; not sure why, as the computation is the
      !! same, but perhaps the zone%rho value is stale?  NNC, 8/9/2012. 
      allocate(rho(ncells))
      call get_density (zone%temp, rho)
      call write_cell_field (rho, 'Z_RHO', for_viz=.true., viz_name='Density')
      deallocate(rho)
    
      !! Cell temperature
      call write_cell_field (zone%temp, 'Z_TEMP', for_viz=.true., viz_name='T')

      !! Average cell enthalpy density
      call write_cell_field (zone%enthalpy, 'Z_ENTHALPY', for_viz=.true., viz_name='Enthalpy')

      !! Phase volume fractions
      if (nmat > 1) then
        allocate(vof(nmat,ncells), name(nmat))
        do m = 1, nmat
          call gather_vof (m, vof(m,:))
          write(name(m),'(a,i4.4)') 'VOF', get_user_material_id(m)
        end do
        call write_cell_field (vof, 'VOF', for_viz=.true., viz_name=name)
        deallocate(vof, name)
      end if

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
      call write_cell_field (vcell, 'Z_VC', for_viz=.true., viz_name=['U','V','W'])
      deallocate(vcell)
      
      !! Cell-centered fluid pressure.
      call write_cell_field (zone%p, 'Z_P', for_viz=.true., viz_name='P')
      
      !! Face fluxing velocities.
      call write_cell_field (fluxing_velocity, 'Face_Vel', for_viz=.false.)
      
      !! Cell-centered fluid Courant number.
      call write_cell_field (courant, 'COURANT', for_viz=.true.)
      
      !! Cell-centered divergence (the volume error).
      allocate(div(ncells))
      call divergence (div)
      call write_cell_field (div, 'Volume_Error', for_viz=.true., viz_name='vol_err')
      deallocate(div)
      
      !! Cell-centered fluid density delta.
      if (boussinesq_approximation) then
        allocate(drho(ncells))
        call get_density_delta (zone%temp, drho)
        call write_cell_field (drho, 'del-rho', for_viz=.true., viz_name='delrho')
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
      call write_cell_field (dTdt, 'dTdt', for_viz=.true., viz_name='dT/dt')
      
      call ds_get_temp_grad (gradT)
      call write_cell_field (gradT, 'Grad_T', for_viz=.true., viz_name=['dT/dx','dT/dy','dT/dz'])

    end subroutine write_heat_transfer_data
    
    subroutine write_EM_data
    
      use EM_data_proxy, only: joule_power_density
      
      real(r8), pointer :: q(:)
      
      q => joule_power_density()
      call write_cell_field (q, 'Joule_P', for_viz=.true.)

    end subroutine write_EM_data
  
    subroutine write_solid_mech_data
    
      use parameter_module, only: ndim, nnodes
      use node_operator_module, only: nipc
      use solid_mechanics_data, only: smech_cell, smech_ip, rhs, displacement, thermal_strain, &
                                      pc_strain, rotation_magnitude, node_gap, node_norm_trac
      use mech_bc_data_module, only: interface_list
      
      integer :: n
      character(32) :: name
      real(r8), allocatable :: tmp(:,:)
      
      !! This is restart data only.
      do n = 1, nipc
        write(name,'(a,i2.2)') 'TOTAL_STRAIN_', n
        call write_cell_field (smech_ip(n)%total_strain, name, for_viz=.false.)
        write(name,'(a,i2.2)') 'ELASTIC_STRESS_', n
        call write_cell_field (smech_ip(n)%elastic_stress, name, for_viz=.false.)
        write(name,'(a,i2.2)') 'PLASTIC_STRAIN_', n
        call write_cell_field (smech_ip(n)%plastic_strain, name, for_viz=.false.)
        write(name,'(a,i2.2)') 'PLASTIC_STRAIN_RATE_', n
        call write_cell_field (smech_ip(n)%plastic_strain_rate, name, for_viz=.false.)
      end do
      
      !! More restart data.  Rhs is ndim-vector node-based data but
      !! stored in a flat array.  We've got to unflatten it, argh!
      allocate(tmp(ndim,nnodes))
      do n = 1, ndim
        tmp(n,:) = rhs(n::ndim)
      end do
      call write_node_field (tmp, 'RHS', for_viz=.false.)
      deallocate(tmp)
      
      !! Restart and viz data.
      call write_cell_field (smech_cell%elastic_stress, 'sigma', for_viz=.true., &
          viz_name=['sigxx', 'sigyy', 'sigzz', 'sigxy', 'sigxz', 'sigyz'])
      call write_cell_field (smech_cell%total_strain, 'epsilon', for_viz=.true., &
          viz_name=['epsxx', 'epsyy', 'epszz', 'epsxy', 'epsxz', 'epsyz'])
      call write_cell_field (smech_cell%plastic_strain, 'e_plastic', for_viz=.true., &
          viz_name=['eplxx', 'eplyy', 'eplzz', 'eplxy', 'eplxz', 'eplyz'])
      call write_cell_field (smech_cell%plastic_strain_rate, 'epsdot', for_viz=.true.)
      call write_cell_field (thermal_strain, 'epstherm', for_viz=.true., &
          viz_name=['epsthxx', 'epsthyy', 'epsthzz', 'epsthxy', 'epsthxz', 'epsthyz'])
      call write_cell_field (pc_strain, 'epspc', for_viz=.true., &
          viz_name=['epspcxx', 'epspcyy', 'epspczz', 'epspcxy', 'epspcxz', 'epspcyz'])
      call write_node_field (displacement, 'Displacement', for_viz=.true., &
          viz_name=['Dx', 'Dy', 'Dz'])
      
      !! Viz data only.
      call write_cell_field (rotation_magnitude, 'Rotation', for_viz=.true.)
      
      !! Gap displacements and scaled forces.  For each 'interface' (likely a
      !! gap element block?) there is an entire node-based field for the gap
      !! displacement and normal traction, with non-zero values only along the
      !! actual interface.  This REALLY needs to be re-implemented using a more
      !! efficient design.  This is viz data only.
      
      do n = 1, size(node_gap,2)
        if (global_any((node_gap(:,n) /= 0.0_r8) .or. (node_norm_trac(:,n) /= 0.0_r8))) then
          write(name,'(a,i2.2)') 'GAP_', interface_list(n)
          call write_node_field (node_gap(:,n), name, for_viz=.true.)
          write(name,'(a,i2.2)') 'NTRAC_', interface_list(n)
          call write_node_field (node_norm_trac(:,n), name, for_viz=.true.)
        end if
      end do
      
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
        call write_cell_field (array, 'phi'//i_to_c(n), for_viz=.true.)
      end do

    end subroutine write_species_data
    
  end subroutine TDO_write_timestep
  
  subroutine write_cell_field_r0 (ldata, name, for_viz, viz_name)
  
    use parameter_module, only: ncells, ncells_tot

    real(r8), intent(in) :: ldata(:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name

    integer :: stat
    real(r8), allocatable :: gdata(:)
    type(c_ptr) :: dataset_id
    
    INSIST(size(ldata) == ncells)
    
    if (nPE == 1) then
      call simulation_data_write (seq_id, name, ldata, stat)
    else
      if (is_IOP) then
        allocate(gdata(ncells_tot))
      else
        allocate(gdata(0))
      end if
      call collate (gdata, ldata)
      if (is_IOP) call simulation_data_write (seq_id, name, gdata, stat)
      call broadcast (stat)
    end if
    
    INSIST(stat == DANU_SUCCESS)
    
    if (for_viz) then
      if (is_IOP) call simulation_open_data (seq_id, name, dataset_id, stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      if (is_IOP) call attribute_write (dataset_id, 'FIELDTYPE', 'CELL', stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      if (is_IOP) then
        if (present(viz_name)) then
          call attribute_write (dataset_id, 'FIELDNAME', viz_name, stat)
        else
          call attribute_write (dataset_id, 'FIELDNAME', name, stat)
        end if
      end if
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
    end if
    
  end subroutine write_cell_field_r0
  
  subroutine write_cell_field_r1 (ldata, name, for_viz, viz_name)
  
    use parameter_module, only: ncells, ncells_tot
    use string_utilities, only: i_to_c

    real(r8), intent(in) :: ldata(:,:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name(:)

    integer :: stat, n
    real(r8), allocatable :: gdata(:,:)
    type(c_ptr) :: dataset_id
    
    INSIST(size(ldata,dim=2) == ncells)
    
    if (nPE == 1) then
      call simulation_data_write (seq_id, name, ldata, stat)
    else
      if (is_IOP) then
        allocate(gdata(size(ldata,1),ncells_tot))
      else
        allocate(gdata(size(ldata,1),0))
      end if
      call collate (gdata, ldata)
      if (is_IOP) call simulation_data_write (seq_id, name, gdata, stat)
      call broadcast (stat)
    end if
    
    INSIST(stat == DANU_SUCCESS)
    
    if (for_viz) then
      if (is_IOP) call simulation_open_data (seq_id, name, dataset_id, stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      if (is_IOP) call attribute_write (dataset_id, 'FIELDTYPE', 'CELL', stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      INSIST(present(viz_name))
      INSIST(size(viz_name) == size(ldata,1))
      ! Argh, no array attributes!  So we'll do it this way for now.
      !if (is_IOP) call attribute_write (dataset_id, 'FIELDNAME', viz_name, stat)
      !call broadcast (stat)
      !INSIST(stat == DANU_SUCCESS)
      do n = 1, size(viz_name)
        if (is_IOP) call attribute_write (dataset_id, 'FIELDNAME'//i_to_c(n), viz_name(n), stat)
        call broadcast (stat)
        INSIST(stat == DANU_SUCCESS)
      end do
    end if
    
  end subroutine write_cell_field_r1
  
  subroutine write_node_field_r0 (ldata, name, for_viz, viz_name)
  
    use parameter_module, only: nnodes, nnodes_tot

    real(r8), intent(in) :: ldata(:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name

    integer :: stat
    real(r8), allocatable :: gdata(:)
    type(c_ptr) :: dataset_id
    
    INSIST(size(ldata) == nnodes)
    
    if (nPE == 1) then
      call simulation_data_write (seq_id, name, ldata, stat)
    else
      if (is_IOP) then
        allocate(gdata(nnodes_tot))
      else
        allocate(gdata(0))
      end if
      call collate (gdata, ldata)
      if (is_IOP) call simulation_data_write (seq_id, name, gdata, stat)
      call broadcast (stat)
    end if
    
    INSIST(stat == DANU_SUCCESS)
    
    if (for_viz) then
      if (is_IOP) call simulation_open_data (seq_id, name, dataset_id, stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      if (is_IOP) call attribute_write (dataset_id, 'FIELDTYPE', 'NODE', stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      if (is_IOP) then
        if (present(viz_name)) then
          call attribute_write (dataset_id, 'FIELDNAME', viz_name, stat)
        else
          call attribute_write (dataset_id, 'FIELDNAME', name, stat)
        end if
      end if
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
    end if
    
  end subroutine write_node_field_r0
  
  subroutine write_node_field_r1 (ldata, name, for_viz, viz_name)
  
    use parameter_module, only: nnodes, nnodes_tot
    use string_utilities, only: i_to_c

    real(r8), intent(in) :: ldata(:,:)
    character(*), intent(in) :: name
    logical, intent(in) :: for_viz
    character(*), intent(in), optional :: viz_name(:)

    integer :: stat, n
    real(r8), allocatable :: gdata(:,:)
    type(c_ptr) :: dataset_id
    
    INSIST(size(ldata,dim=2) == nnodes)
    
    if (nPE == 1) then
      call simulation_data_write (seq_id, name, ldata, stat)
    else
      if (is_IOP) then
        allocate(gdata(size(ldata,1),nnodes_tot))
      else
        allocate(gdata(size(ldata,1),0))
      end if
      call collate (gdata, ldata)
      if (is_IOP) call simulation_data_write (seq_id, name, gdata, stat)
      call broadcast (stat)
    end if
    
    INSIST(stat == DANU_SUCCESS)
    
    if (for_viz) then
      if (is_IOP) call simulation_open_data (seq_id, name, dataset_id, stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      if (is_IOP) call attribute_write (dataset_id, 'FIELDTYPE', 'NODE', stat)
      call broadcast (stat)
      INSIST(stat == DANU_SUCCESS)
      INSIST(present(viz_name))
      INSIST(size(viz_name) == size(ldata,1))
      ! Argh, no array attributes!  So we'll do it this way for now.
      !if (is_IOP) call attribute_write (dataset_id, 'FIELDNAME', viz_name, stat)
      !call broadcast (stat)
      !INSIST(stat == DANU_SUCCESS)
      do n = 1, size(viz_name)
        if (is_IOP) call attribute_write (dataset_id, 'FIELDNAME'//i_to_c(n), viz_name(n), stat)
        call broadcast (stat)
        INSIST(stat == DANU_SUCCESS)
      end do
    end if
    
  end subroutine write_node_field_r1
  
#endif

end module truchas_danu_output
