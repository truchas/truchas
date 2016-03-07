!!
!! SOLID_MECHANICS_OUTPUT  
!!
!! The purpose of this module is to provide interfaces to data required for
!! routines outside of the solid mechanics package. This module is primarily
!! used by output modules such as long_edit_moudle.F90 to write out data to
!! the screen and files. These modules still use the old mesh module and this 
!! module provides the means to translate between old and DIST_MESH.
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
!! PUBLIC Datatypes
!!
!! CELL_MECH_INVARIANT
!!
!! An array of user defined types that holds pointers to stress and strain tensor 
!! quanities. Only used in  edit_module.F90 and long_edit_module.F90. Should 
!! be replaced.
!!
!! PUBLIC Subroutines and functions
!!
!! GET_SM_THERMAL_STRAIN              
!! GET_SM_DISPLACEMENT                
!! GET_SM_ROTATION_MAGNITUDE          
!! GET_SM_PC_STRAIN  
!! GET_SM_RHS
!!
!! Public subroutines that access data defined in the SOLID_MECHANICS_DATA
!! module. Created to allow modules outside of solid mechanics to access data
!! that is internally defined on DIST_MESH.
!!
!! GET_SMECH_CELL_*
!!
!! A number of data accessor subroutines that return data held in SMECH_CELL
!! SMECH_CELL is an array of defined type, MECH_DATA, defined on cells
!! in the solid_mechanics_data.F90 file. 
!! In the original code, the output subroutines outside of solid mechanics
!! would import this variable into the local namespace to dump out data. 
!! With the new DIST_MESH, modules outside of solid mechanics must use these
!! routines to access this information.
!!
!! GET_SMECH_IP_TOTAL_STRAIN
!! GET_SMECH_IP_ELASTIC_STRESS
!! GET_SMECH_IP_PLASTIC_STRAIN
!! GET_SMECH_IP_PLASTIC_STRAIN_RATE  
!!
!! A number of data accessor subroutines that return data held in SMECH_IP
!! SMECH_IP is an array of defined type, MECH_DATA, defined on integration
!! points in the solid_mechanics_data.F90 file. 
!! In the original code, the output subroutines outside of solid mechanics
!! would import this variable into the local namespace to dump out data. 
!! With the new DIST_MESH, modules outside of solid mechanics must use these
!! routines to access this information.
!!
!! SET_SM_ROTATION_MAGNITUDE
!! SET_SMECH_CELL_TOTAL_STRAIN
!! SET_SMECH_CELL_ELASTIC_STRESS
!! SET_SMECH_CELL_PLASTIC_STRAIN
!! SET_SMECH_CELL_PLASTIC_STRAIN_RATE
!!
!! In gap_module.F90, the subroutine set_gap_element_output alters theses
!! data sets. These subroutines were added to oreserve the current output
!! state of the solid mechanics tests in the test suite. Should be removed
!! at some point in the future.
!!
!! GET_SM_NODE_GAP
!! GET_SM_NODE_NORM_TRAC
!!
!! Accessor subroutines for the NODE_GAP and NODE_NORM_TRAC node-based arrays
!! defined in solid_mechanics_data.F90. Only returns array slices. Must pass
!! in an index to request a specific slice. Calling subroutines should call 
!! SM_NODE_GAP_ISIZE to define the size of the second dimension of NODE_GAP
!! and NODE_NORM_TRAC. 
!!
!! SM_NODE_GAP_ISIZE
!!
!! The array NODE_GAP from the SOLID_MECHANICS_DATA module is
!! rank 2 and the second index size is based on the number of
!! non-zero entries in INTERFACE_LIST defined in the 
!! MECH_BC_DATA_MODULE. Call this function before GET_SM_NODE_GAP
!! or GET_SM_NODE_NORM_TRAC to pass the correct sized array required
!! by these routines.
!!
!! SMECH_NUM_INT_PTS
!!
!! Returns the number of integration points parameter defined in the
!! NODE_OPERATOR_MODULE. Need to call this to determine the size of
!! SMech_IP array.
!!

#include "f90_assert.fpp"

module solid_mechanics_output
      use kinds, only: r8
      implicit none
      private

      ! This data type and instance are only for use in long edit.  We should not be
      ! setting up a data structure this way (an array of types).
      public :: CELL_MECH_INVARIANT
      type CELL_MECH_INVARIANT
         real(r8) :: mises_stress
         real(r8) :: eff_plastic_strain
         real(r8) :: mean_stress
         real(r8) :: volumetric_strain
      end type CELL_MECH_INVARIANT

      ! Public data accessors
      public :: get_sm_thermal_strain,              &
                get_sm_displacement,                &
                get_sm_rhs,                         &
                get_sm_rotation_magnitude,          &
                set_sm_rotation_magnitude,          &
                get_sm_pc_strain,                   &
                get_smech_cell_total_strain,        &
                set_smech_cell_total_strain,        &
                get_smech_cell_elastic_stress,      &
                set_smech_cell_elastic_stress,      &
                get_smech_cell_plastic_strain,      &
                set_smech_cell_plastic_strain,      &
                get_smech_cell_plastic_strain_rate, &
                set_smech_cell_plastic_strain_rate, &
                smech_num_int_pts,                  &
                get_smech_ip_total_strain,          &
                get_smech_ip_elastic_stress,        &
                get_smech_ip_plastic_strain,        &
                get_smech_ip_plastic_strain_rate,   &
                sm_node_gap_isize,                  &
                get_sm_node_gap,                    &
                get_sm_node_norm_trac

      ! Public procedures
      public :: read_sm_data, skip_sm_data

      ! Public functions
      public :: stress_strain_invariants

contains

  subroutine get_sm_thermal_strain(output)
    use solid_mechanics_data, only: Thermal_Strain
    use parameter_module, only: ncomps
    use legacy_mesh_api, only: ncells
    real(r8), intent(out) :: output(:,:)
    ASSERT(size(output,1) == ncomps)
    ASSERT(size(output,2) == ncells)
    output = Thermal_Strain
  end subroutine get_sm_thermal_strain              

  subroutine get_sm_displacement(output)
    use solid_mechanics_data, only: Displacement
    use legacy_mesh_api, only: nnodes
    real(r8), intent(out) :: output(:,:)
    ASSERT(size(output,1) == 3)
    ASSERT(size(output,2) == nnodes)
    output = Displacement
  end subroutine get_sm_displacement

  subroutine get_sm_rhs(output)
    use solid_mechanics_data, only: RHS
    use legacy_mesh_api, only: nnodes
    real(r8), intent(out) :: output(:,:)
    real(r8), pointer :: rhs2(:,:)
    ASSERT(size(output,1) == 3)
    ASSERT(size(output,2) == nnodes)
    rhs2(1:3,1:nnodes) => RHS
    output = rhs2
  end subroutine get_sm_rhs
               
  subroutine get_sm_rotation_magnitude(output)
    use solid_mechanics_data, only: ROTATION_MAGNITUDE
    use legacy_mesh_api, only: ncells
    real(r8), intent(out) :: output(:)
    ASSERT(size(output) == ncells)
    output = ROTATION_MAGNITUDE
  end subroutine get_sm_rotation_magnitude
 
  subroutine set_sm_rotation_magnitude(output)
    use solid_mechanics_data, only: ROTATION_MAGNITUDE
    use legacy_mesh_api, only: ncells
    real(r8), intent(in) :: output(:)
    ASSERT(size(output) == ncells)
    ROTATION_MAGNITUDE = output
  end subroutine set_sm_rotation_magnitude
 
  subroutine get_sm_pc_strain(output)
    use solid_mechanics_data, only: PC_STRAIN
    use parameter_module, only: ncomps
    use legacy_mesh_api, only: ncells
    real(r8), intent(out) :: output(:,:)
    ASSERT(size(output,1) == ncomps)
    ASSERT(size(output,2) == ncells)
    output = PC_STRAIN
  end subroutine get_sm_pc_strain

  subroutine get_smech_cell_total_strain(output)
    use solid_mechanics_data, only: SMECH_CELL
    use parameter_module, only: ncomps
    use legacy_mesh_api, only: ncells
    real(r8), intent(out) :: output(:,:)
    ASSERT(size(output,1) == ncomps)
    ASSERT(size(output,2) == ncells)
    output = SMECH_CELL%TOTAL_STRAIN
  end subroutine get_smech_cell_total_strain

  subroutine set_smech_cell_total_strain(output)
    use solid_mechanics_data, only: SMECH_CELL
    use parameter_module, only: ncomps
    use legacy_mesh_api, only: ncells
    real(r8), intent(in) :: output(:,:)
    ASSERT(size(output,1) == ncomps)
    ASSERT(size(output,2) == ncells)
    SMECH_CELL%TOTAL_STRAIN = output
  end subroutine set_smech_cell_total_strain

  subroutine get_smech_cell_elastic_stress(output)
    use solid_mechanics_data, only: SMECH_CELL
    use parameter_module, only: ncomps
    use legacy_mesh_api, only: ncells
    real(r8), intent(out) :: output(:,:)
    ASSERT(size(output,1) == ncomps)
    ASSERT(size(output,2) == ncells)
    output = SMECH_CELL%ELASTIC_STRESS
  end subroutine get_smech_cell_elastic_stress

  subroutine set_smech_cell_elastic_stress(output)
    use solid_mechanics_data, only: SMECH_CELL
    use parameter_module, only: ncomps
    use legacy_mesh_api, only: ncells
    real(r8), intent(in) :: output(:,:)
    ASSERT(size(output,1) == ncomps)
    ASSERT(size(output,2) == ncells)
    SMECH_CELL%ELASTIC_STRESS = output
  end subroutine set_smech_cell_elastic_stress

  subroutine get_smech_cell_plastic_strain(output)
    use solid_mechanics_data, only: SMECH_CELL
    use parameter_module, only: ncomps
    use legacy_mesh_api, only: ncells
    real(r8), intent(out) :: output(:,:)
    ASSERT(size(output,1) == ncomps)
    ASSERT(size(output,2) == ncells)
    output = SMECH_CELL%PLASTIC_STRAIN
  end subroutine get_smech_cell_plastic_strain

  subroutine set_smech_cell_plastic_strain(output)
    use solid_mechanics_data, only: SMECH_CELL
    use parameter_module, only: ncomps
    use legacy_mesh_api, only: ncells
    real(r8), intent(in) :: output(:,:)
    ASSERT(size(output,1) == ncomps)
    ASSERT(size(output,2) == ncells)
    SMECH_CELL%PLASTIC_STRAIN = output
  end subroutine set_smech_cell_plastic_strain

  subroutine get_smech_cell_plastic_strain_rate(output)
    use solid_mechanics_data, only: SMECH_CELL
    use legacy_mesh_api, only: ncells
    real(r8), intent(out) :: output(:)
    ASSERT(size(output) == ncells)
    output = SMECH_CELL%PLASTIC_STRAIN_RATE
  end subroutine get_smech_cell_plastic_strain_rate
 
  subroutine set_smech_cell_plastic_strain_rate(output)
    use solid_mechanics_data, only: SMECH_CELL
    use legacy_mesh_api, only: ncells
    real(r8), intent(in) :: output(:)
    ASSERT(size(output) == ncells)
    SMECH_CELL%PLASTIC_STRAIN_RATE = output
  end subroutine set_smech_cell_plastic_strain_rate

  integer function smech_num_int_pts() result(n)
    use node_operator_module, only: nipc
    n = nipc
  end function smech_num_int_pts

  subroutine get_smech_ip_total_strain(idx,output)
    use solid_mechanics_data, only: SMECH_IP
    use parameter_module, only: ncomps
    use legacy_mesh_api, only: ncells
    integer,  intent(in)  :: idx
    real(r8), intent(out) :: output(:,:)
    ASSERT(size(output,1) == ncomps)
    ASSERT(size(output,2) == ncells)
    output = SMECH_IP(idx)%TOTAL_STRAIN
  end subroutine get_smech_ip_total_strain

  subroutine get_smech_ip_elastic_stress(idx,output)
    use solid_mechanics_data, only: SMECH_IP
    use parameter_module, only: ncomps
    use legacy_mesh_api, only: ncells
    integer,  intent(in)  :: idx
    real(r8), intent(out) :: output(:,:)
    ASSERT(size(output,1) == ncomps)
    ASSERT(size(output,2) == ncells)
    output = SMECH_IP(idx)%ELASTIC_STRESS
  end subroutine get_smech_ip_elastic_stress

  subroutine get_smech_ip_plastic_strain(idx,output)
    use solid_mechanics_data, only: SMECH_IP
    use parameter_module, only: ncomps
    use legacy_mesh_api, only: ncells
    integer,  intent(in)  :: idx
    real(r8), intent(out) :: output(:,:)
    ASSERT(size(output,1) == ncomps)
    ASSERT(size(output,2) == ncells)
    output = SMECH_IP(idx)%PLASTIC_STRAIN
  end subroutine get_smech_ip_plastic_strain

  subroutine get_smech_ip_plastic_strain_rate(idx,output)
    use solid_mechanics_data, only: SMECH_IP
    use legacy_mesh_api, only: ncells
    integer,  intent(in)  :: idx
    real(r8), intent(out) :: output(:)
    ASSERT(size(output) == ncells)
    output = SMECH_IP(idx)%PLASTIC_STRAIN_RATE
  end subroutine get_smech_ip_plastic_strain_rate


  integer function sm_node_gap_isize() result(isize)
    use solid_mechanics_data, only: Node_Gap
    isize = size(Node_Gap,2)
  end function sm_node_gap_isize

  subroutine get_sm_node_gap(idx,output)
    use solid_mechanics_data, only: Node_Gap
    integer, intent(in) :: idx
    real(r8), dimension(:), intent(out) :: output
    output = Node_Gap(:,idx)
  end subroutine get_sm_node_gap

  subroutine get_sm_node_norm_trac(idx,output)
    use solid_mechanics_data, only: Node_Norm_Trac
    integer, intent(in) :: idx
    real(r8), dimension(:), intent(out) :: output
    output = Node_Norm_Trac(:,idx)
  end subroutine get_sm_node_norm_trac

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! READ_SM_DATA
 !!
 !! Neil N. Carlson <nnc@lanl.gov>
 !! 20 Apr 2005
 !!
 !! This subroutine reads the solid mechanics data from a restart file opened
 !! (and pre-positioned) on UNIT, and initializes the module structures
 !! SMECH_IP, SMECH_CELL, RHS, THERMAL_STRAIN, PC_STRAIN, and DISPLACEMENT
 !! with this data (properly distributed and permuted).  VERSION is the version
 !! number of the restart file format.
 !!
 !! NB: It is assumed that storage for the module structures being initialized
 !! has already been suitably allocated.
 !!

  subroutine read_sm_data (unit, version)

    use restart_utilities, only: read_dist_array
    use node_operator_module, only: nipc
    use legacy_mesh_api, only: pcell => unpermute_mesh_vector, pnode => unpermute_vertex_vector
    use solid_mechanics_data, only: displacement, pc_strain, thermal_strain, &
                                     rhs, smech_ip, smech_cell, mech_data
    use legacy_mesh_api, only: ndim                              

    integer, intent(in) :: unit, version

    integer :: n

    ASSERT( associated(smech_ip) )
    ASSERT( size(smech_ip) == nipc )

    do n = 1, nipc
      call read_mech_data (smech_ip(n))
    end do
    call read_mech_data (smech_cell)

    do n = 1, ndim  ! data layout same as a rank-2 array dimensions (ndim,nnodes)
      call read_dist_array (unit, rhs(n::ndim), pnode, 'READ_SM_DATA: error reading RHS records')
    end do
    call read_dist_array (unit, thermal_strain, pcell, 'READ_SM_DATA: error reading THERMAL_STRAIN records')
    call read_dist_array (unit, pc_strain,      pcell, 'READ_SM_DATA: error reading PC_STRAIN records')
    call read_dist_array (unit, displacement,   pnode, 'READ_SM_DATA: error reading DISPLACEMENT records')

  contains

    subroutine read_mech_data (mech)
      type(mech_data), intent(inout) :: mech  ! has pointer components
      call read_dist_array (unit, mech%total_strain,   pcell, 'READ_SM_DATA: error reading TOTAL_STRAIN records')
      call read_dist_array (unit, mech%elastic_stress, pcell, 'READ_SM_DATA: error reading ELASTIC_STRESS records')
      call read_dist_array (unit, mech%plastic_strain, pcell, 'READ_SM_DATA: error reading PLASTIC_STRAIN records')
      call read_dist_array (unit, mech%plastic_strain_rate, pcell, 'READ_SM_DATA: error reading PLASTIC_STRAIN_RATE records')
    end subroutine read_mech_data

  end subroutine read_SM_data

  subroutine skip_sm_data (unit, version)
    use restart_utilities, only: skip_records
    integer, intent(in) :: unit, version
    call skip_records (unit, 265, 'SKIP_SM_DATA: error skipping the solid mechanics data')
  end subroutine skip_SM_data

  !! This function was originally in solid_mechanics_module
  !! It is used by the edit modules EDIT_MODULE and LONG_EDIT_DATA_TYPES to 
  !! define the CELL_MECH_DATA defined types declared in these modules

  function stress_strain_invariants () result(stress_strain)
    !=============================================================================
    !
    ! Calculate invariants of stress and strain tensors at cell centers, which are
    ! returned in Stress_Strain.  Right now this should only be called for long 
    ! edit output.
    ! 
    !=============================================================================
    use parameter_module, only: ncomps
    use legacy_mesh_api, only: ncells
    use truchas_logging_services, only: TLS_panic

    integer :: status

    type(CELL_MECH_INVARIANT), pointer, dimension(:) :: Stress_Strain

    ! local variables
    real(r8), allocatable :: scratch(:,:)

    allocate(stress_strain(ncells),stat = status)
    if (status /= 0) call TLS_panic ( 'STRESS_STRAIN_INVARIANTS: allocation error: Stress_Strain')

    ! Scratch buffer for the accessor functions
    allocate(scratch(ncomps,ncells))


    ! Von Mises stress
    !Stress_Strain(:)%mises_stress = sqrt(((SMech_Cell%Elastic_Stress(1,:) - SMech_Cell%Elastic_Stress(2,:))**2 + &
    !     (SMech_Cell%Elastic_Stress(2,:) - SMech_Cell%Elastic_Stress(3,:))**2 + &
    !     (SMech_Cell%Elastic_Stress(3,:) - SMech_Cell%Elastic_Stress(1,:))**2 + &
    !     6.0 * (SMech_Cell%Elastic_Stress(4,:)**2 + SMech_Cell%Elastic_Stress(5,:)**2 + SMech_Cell%Elastic_Stress(6,:)**2))/2.)
    call get_smech_cell_elastic_stress(scratch)
    Stress_Strain(:)%mises_stress = sqrt(((scratch(1,:) - scratch(2,:))**2 + &
         (scratch(2,:) - scratch(3,:))**2 + &
         (scratch(3,:) - scratch(1,:))**2 + &
         6.0 * (scratch(4,:)**2 + scratch(5,:)**2 + scratch(6,:)**2))/2.)
    Stress_Strain(:)%mean_stress = (scratch(1,:) + scratch(2,:) + &
                                    scratch(3,:))/3.0
    ! Not sure how to define this yet
    !Stress_Strain(:)%eff_plastic_strain = sqrt(((SMech_Cell%Plastic_Strain(1,:) - SMech_Cell%Plastic_Strain(2,:))**2 + &
    !                                       (SMech_Cell%Plastic_Strain(2,:) - SMech_Cell%Plastic_Strain(3,:))**2 + &
    !                                       (SMech_Cell%Plastic_Strain(3,:) - SMech_Cell%Plastic_Strain(1,:))**2) * 2.0/9.0 + &
    !                                      (SMech_Cell%Plastic_Strain(4,:)**2 + SMech_Cell%Plastic_Strain(5,:)**2 + &
    !                                       SMech_Cell%Plastic_Strain(6,:)**2) * 4.0/3.0)
    call get_smech_cell_plastic_strain(scratch)
    Stress_Strain(:)%eff_plastic_strain = sqrt(((scratch(1,:) - scratch(2,:))**2 + &
                                           (scratch(2,:) - scratch(3,:))**2 + &
                                           (scratch(3,:) - scratch(1,:))**2) * 2.0/9.0 + &
                                          (scratch(4,:)**2 + scratch(5,:)**2 + &
                                           scratch(6,:)**2) * 4.0/3.0)

    ! Mean stress and strain
    !Stress_Strain(:)%mean_stress = (SMech_Cell%Elastic_Stress(1,:) + SMech_Cell%Elastic_Stress(2,:) + &
    !                                SMech_Cell%Elastic_Stress(3,:))/3.0
    !Stress_Strain(:)%volumetric_strain =  SMech_Cell%Total_Strain(1,:) + SMech_Cell%Total_Strain(2,:) + &
    !                                      SMech_Cell%Total_Strain(3,:)
    call get_smech_cell_total_strain(scratch)
    Stress_Strain(:)%volumetric_strain =  scratch(1,:) + scratch(2,:) + &
                                          scratch(3,:)

    ! Clean up
    deallocate(scratch)

  end function stress_strain_invariants
  !
  !

end module solid_mechanics_output
