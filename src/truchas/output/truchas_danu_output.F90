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

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_danu_output_data
  use truchas_danu_output_tools
  use parallel_communication
  use truchas_logging_services
  use truchas_h5_outfile, only: th5_mesh_group, th5_seq_group
  implicit none
  private

  public :: TDO_write_mesh
  public :: TDO_write_timestep

  type(th5_mesh_group), save :: out_mesh
  type(th5_seq_group), save :: seq

  public :: outfile ! others may want to write here

contains

  subroutine TDO_write_mesh(mesh)

    use unstr_mesh_type
    use degen_hex_topology, only: OLD_TET_NODE_MAP, OLD_PYR_NODE_MAP, OLD_PRI_NODE_MAP
    use output_control,  only: part
    use truchas_logging_services
    use truchas_env, only: output_file_name, overwrite_output

    type(unstr_mesh), intent(in) :: mesh

    integer :: j
    integer, allocatable :: hnode(:,:), cblockid(:)
    logical :: exists
    integer :: ndim, nvc, ncells, ncells_tot, nnodes, nnodes_tot

    if (is_IOP .and. .not. overwrite_output) then
      inquire(file=output_file_name('h5'), exist=exists)
        if (exists) then
          call TLS_panic("must specify `-f` flag to overwrite `" // output_file_name('h5')//"`")
        end if
    endif
    call outfile%open (output_file_name('h5'), io_group_size, is_IOP)

    ndim = 3
    nvc  = 8
    ncells = mesh%ncell_onP
    ncells_tot = mesh%cell_imap%global_size
    nnodes = mesh%nnode_onP
    nnodes_tot = mesh%node_imap%global_size

    !! Create the mesh entry.
    call outfile%add_unstr_mesh_group('DEFAULT', nvc, ndim, out_mesh)

    !! Write the node coordinates.
    call out_mesh%write_coordinates(nnodes_tot, mesh%x(:,:nnodes))

    !! Write the cell connectivity as degenerate hexes
    allocate(hnode(nvc,ncells))
    do j = 1, ncells
      associate (cnode => mesh%cnode(mesh%xcnode(j):mesh%xcnode(j+1)-1))
        select case (size(cnode))
        case (4) ! tet
          hnode(:,j) = cnode(OLD_TET_NODE_MAP)
        case (5) ! pyramid
          hnode(:,j) = cnode(OLD_PYR_NODE_MAP)
        case (6) ! prism
          hnode(:,j) = cnode(OLD_PRI_NODE_MAP)
        case (8) ! hex
          hnode(:,j) = cnode
        case default
          INSIST(.false.)
        end select
        hnode(:,j) = mesh%node_imap%global_index(hnode(:,j))  ! map to global node IDs
      end associate
    end do
    call out_mesh%write_connectivity(ncells_tot, hnode)
    deallocate(hnode)

    ! I don't know where this should go
    call TDO_start_simulation

    !! Right now the following is being written to the non-series section
    !! of the simulation.  Soon we will write these in the mesh section.

    !! Mapping from internal serial cell numbering to external numbering.
    call sim%write_dist_array('CELLMAP', mesh%xcell(:mesh%ncell_onP), ncells_tot)

    !! Cell block IDs.
    allocate(cblockid(ncells))
    do j = 1, ncells
     associate (bitmask => mesh%cell_set_mask(j))
       INSIST(popcnt(bitmask) == 1)
       cblockid(j) = mesh%cell_set_id(trailz(bitmask))
     end associate
    end do
    call sim%write_dist_array('BLOCKID', cblockid, ncells_tot)
    deallocate(cblockid)

    !! Cell partition assignment.
    call sim%write_dist_array('CELLPART', spread(this_PE,dim=1,ncopies=ncells), ncells_tot)

    !! Mapping from internal serial node numbering to external numbering.
    call sim%write_dist_array('NODEMAP', mesh%xnode(:mesh%nnode_onP), nnodes_tot)

    !! Node partition assignment.
    call sim%write_dist_array('NODEPART', spread(this_PE,dim=1,ncopies=nnodes), nnodes_tot)

    call sim%write_repl_data('NUMPROCS', nPE)

    !! Parts for movement
    if (size(part) > 0) call sim%write_repl_data('part1', part)

    call outfile%close()

  end subroutine TDO_write_mesh

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

    use unstr_mesh_type
    use mesh_manager, only: unstr_mesh_ptr
    use time_step_module, only: t, dt, cycle_number
    use physics_module, only: heat_transport, species_transport
    use ih_driver, only: ih_enabled
    use ustruc_driver, only: ustruc_output
    use flow_driver, only: flow_enabled
    use solid_mechanics_driver, only: solid_mechanics_enabled
    use output_control, only: part_path, write_mesh_partition
    use truchas_env, only: output_file_name

    integer :: ncells
    real(r8) :: r(3)
    type(unstr_mesh), pointer :: mesh

    mesh => unstr_mesh_ptr('MAIN')
    INSIST(associated(mesh))
    ncells = mesh%ncell_onP

    call outfile%reopen (output_file_name('h5'), io_group_size, is_IOP)

    call sim%next_seq_group(cycle_number, t, seq)
    call seq%write_attr('time step', dt)

    !! Part movement
    if (allocated(part_path)) then
      call part_path%get_position(t, r)
      call seq%write_attr('translate_part1', r)
    end if

    !! Cell density, temperature, enthalpy, phase volume fractions.
    call write_common_data

    !! Flow-related fields.
    if (flow_enabled()) call write_new_flow_data

    !! Heat transfer fields (other than temperature and enthalpy).
    if (heat_transport) call write_heat_transfer_data

    !! Induction heating fields.
    if (ih_enabled()) call write_ih_data

    !! Solid mechanics fields.
    if (solid_mechanics_enabled()) call write_solid_mechanics_data

    !! Species fields.
    if (species_transport) call write_species_data

    !! Microstructure analysis data (if enabled)
    call ustruc_output (seq)

    call outfile%close()

  contains

    subroutine write_common_data

      use zone_module, only: zone
      use legacy_matl_api, only: nmat, gather_vof
      use material_model_driver, only: matl_model

      integer :: m
      real(r8), allocatable :: vof(:,:)
      character(32), allocatable :: name(:)

      !! Average cell density
      call write_seq_cell_field (seq, zone%rho, 'Z_RHO', for_viz=.true., viz_name='Density')

      !! Cell temperature
      call write_seq_cell_field (seq, zone%temp, 'Z_TEMP', for_viz=.true., viz_name='T')

      !! Average cell enthalpy density
      call write_seq_cell_field (seq, zone%enthalpy, 'Z_ENTHALPY', for_viz=.true., viz_name='Enthalpy')

      !! Phase volume fractions
      if (nmat > 1) then
        allocate(vof(nmat,ncells), name(nmat))
        do m = 1, nmat
          call gather_vof (m, vof(m,:))
          write(name(m),'(a)') matl_model%phase_name(m)
        end do
        call write_seq_cell_field (seq, vof, 'VOF', for_viz=.true., viz_name=name)
        deallocate(vof, name)
      end if

      !! This is a stop-gap because the Truchas paraview reader ignores this
      !! data written with the mesh. It is only written to the first snapshot.
      if (write_mesh_partition) then
        call write_seq_cell_field(seq, spread(real(this_PE,r8),dim=1,ncopies=ncells), 'CELLPART', for_viz=.true., viz_name='rank')
        write_mesh_partition = .false.
      end if

    end subroutine write_common_data

    subroutine write_new_flow_data
      use flow_driver

      real(r8), pointer :: vec_cc(:,:), scalar_cc(:)
      real(r8) :: fluxing_velocity(6,ncells)

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

      use zone_module, only: zone
      use time_step_module, only: dt
      use diffusion_solver, only: ds_get_temp_grad

      real(r8) :: dTdt(ncells), gradT(3,ncells)

      dTdt = (zone%temp - zone%temp_old) / dt
      call write_seq_cell_field (seq, dTdt, 'dTdt', for_viz=.true., viz_name='dT/dt')

      call ds_get_temp_grad (gradT)
      call write_seq_cell_field (seq, gradT, 'Grad_T', for_viz=.true., viz_name=['dT/dx','dT/dy','dT/dz'])

    end subroutine write_heat_transfer_data

    subroutine write_ih_data

      use ih_driver, only: joule_power_density

      real(r8), pointer :: q(:)

      q => joule_power_density()
      call write_seq_cell_field (seq, q, 'Joule_P', for_viz=.true.)

    end subroutine

    subroutine write_solid_mechanics_data

      use solid_mechanics_driver

      real(r8), allocatable :: displ(:,:), thermal_strain(:,:), total_strain(:,:), &
          elastic_stress(:,:), rotation(:), gap_displacement(:), gap_normal_traction(:), &
          plastic_strain(:,:), plastic_strain_rate(:)

      call solid_mechanics_compute_viz_fields(displ, thermal_strain, total_strain, &
          elastic_stress, rotation, gap_displacement, gap_normal_traction, &
          plastic_strain, plastic_strain_rate)

      call write_seq_node_field(seq, displ, 'Displacement', for_viz=.true., &
          viz_name=['Dx', 'Dy', 'Dz'])
      call write_seq_cell_field(seq, total_strain, 'epsilon', for_viz=.true., &
          viz_name=['epsxx', 'epsyy', 'epszz', 'epsxy', 'epsxz', 'epsyz'])
      call write_seq_cell_field(seq, thermal_strain, 'epstherm', for_viz=.true., &
          viz_name=['epsthxx', 'epsthyy', 'epsthzz', 'epsthxy', 'epsthxz', 'epsthyz'])
      call write_seq_cell_field (seq, elastic_stress, 'sigma', for_viz=.true., &
          viz_name=['sigxx', 'sigyy', 'sigzz', 'sigxy', 'sigxz', 'sigyz'])
      call write_seq_cell_field(seq, rotation, 'Rotation', for_viz=.true.)
      !! Note: the legacy solver also output phase change strain on cells, but this
      !! field seemed to never be set.

      !! NB: These gap fields are nonzero only where there is a gap BC.
      !! They hold data at nodes, but only for one gap condition. At
      !! nodes where more than one gap intersect, each node only holds
      !! one value; i.e., there is only one displacement vizualised even
      !! though there are multiple displacements across multiple gaps at
      !! that node. This could be addressed by dumping separate
      !! displacement and traction fields for every gap BC (wasteful),
      !! or a new output type which associates a small field with given
      !! sidesets.
      call write_seq_node_field(seq, gap_displacement, 'Gap Displacement', for_viz=.true.)
      call write_seq_node_field(seq, gap_normal_traction, 'Gap Normal Traction', for_viz=.true.)

      if (solid_mechanics_viscoplasticity_enabled()) then
        call write_seq_cell_field(seq, plastic_strain, 'e_plastic', for_viz=.true., &
          viz_name=['eplxx', 'eplyy', 'eplzz', 'eplxy', 'eplxz', 'eplyz'])
        call write_seq_cell_field(seq, plastic_strain_rate, 'epsdot', for_viz=.true.)
        call solid_mechanics_write_checkpoint(seq) ! restart-only data
      end if

    end subroutine write_solid_mechanics_data

    subroutine write_species_data

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
