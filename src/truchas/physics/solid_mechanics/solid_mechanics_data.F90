#include "f90_assert.fpp"

Module SOLID_MECHANICS_DATA
  !-----------------------------------------------------------------------------
  ! Purpose:
  !
  !  Data for solid mechanics routines
  !
  !             
  ! Authors:  Dave Korzekwa (dak@lanl.gov)
  !-----------------------------------------------------------------------------
  !
  use kinds, only: r8
  use parameter_module, only: string_len, ncomps, ndim
  implicit none
  private

  ! Public data
  Public :: solid_mechanics,                  &
            solid_mechanics_body_force,     &
            Solid_Mask,                       &
            gap_cell_mask,                    &
            mech_interface,                   &
            !
            Thermal_Strain,                   &
            Thermal_Strain_Inc,               &
            tm_dens_old,                       &
            tm_dens,                          &
            PC_Strain,                        &
            Elastic_Strain,                   &
            Rotation_Magnitude,               &
            !
            SMech_IP,                          &
            SMech_IP_Old,                      &
            SMech_Cell,                        &
            SMech_Cell_Old,                    &
            !
            Eff_Stress_Cell_old,              &
            Eff_Stress_IP_old,                &
            !
            Dev_Stress_Cell,                  &
            Dev_Stress_IP,                    &
            !
            Lame1,                            &
            Lame2,                            &
            CTE,                              &
            Displacement,                     &
            CELL_MECH_INVARIANT,              &
            MECH_DATA,                        &
            displacement_linear_solution,     &
            displacement_nonlinear_solution,  &
            Ubik_DISPLACEMENT,                &
            thermo_elastic_iterations,        &
            viscoplastic_iterations,          &
            thermo_elastic_precond_iter,      &
            NK_DISPLACEMENT,                  &
            RHS,                              &
            Src,                              &
            Viscoplastic_Model_Forms,         &
            maxvpforms,                       &
            plasticity,                       &
            Node_Gap,                         &
            Node_Norm_Trac,                   &
            contact_distance,                 &
            contact_norm_trac,                   &
            contact_penalty,                  &
            strain_limit,                     &
            cscale

  !
  ! public procedures
  public :: read_SM_data, skip_SM_data, GAP_NODE_DISPLACEMENT

  INTERFACE COLLATE
     MODULE PROCEDURE COLLATE_MECH
  END INTERFACE
  
  ! Linear and non-linear solution names
  character(string_len), save :: displacement_linear_solution
  character(string_len), save :: displacement_nonlinear_solution
  !
  ! Ubik_user element number to use for energy linear solve
  integer, save :: Ubik_DISPLACEMENT
  ! Linear iteration count
  integer, save :: thermo_elastic_iterations
  ! Nonlinear iteration count
  integer, save :: viscoplastic_iterations
  ! Precondition iteration count
  integer, save :: thermo_elastic_precond_iter
  !
  ! NKuser element number to use for energy solution
  integer, save :: NK_DISPLACEMENT
  ! This data type and instance are only for use in long edit.  We should not be
  ! setting up a data structure this way (an array of types).
  type CELL_MECH_INVARIANT
     real(r8) :: mises_stress
     real(r8) :: eff_plastic_strain
     real(r8) :: mean_stress
     real(r8) :: volumetric_strain
  end type CELL_MECH_INVARIANT
  !
  ! Needed in processing the namelist input
  integer, parameter                     :: maxvpforms = 4
  character(80), dimension(maxvpforms), save :: Viscoplastic_Model_Forms
  ! This data type is used to store the solid mechanics state at integration points and 
  ! cell centers
  type MECH_DATA
     real(r8), pointer, dimension(:,:) :: Total_Strain => null()
     real(r8), pointer, dimension(:,:) :: Elastic_Stress => null()
     real(r8), pointer, dimension(:,:) :: Plastic_Strain => null()
     real(r8), pointer, dimension(:)   :: Plastic_Strain_Rate => null()
  end type MECH_DATA
  !
  ! Various stresses and strains

  type(MECH_DATA), save, pointer, dimension(:) :: SMech_IP
  type(MECH_DATA), save, pointer, dimension(:) :: SMech_IP_Old
  type(MECH_DATA), save                        :: SMech_Cell
  type(MECH_DATA), save                        :: SMech_Cell_Old

  real(r8), Save, Pointer, Dimension(:,:)   :: Thermal_Strain
  real(r8), Save, Pointer, Dimension(:,:)   :: Thermal_Strain_Inc
  real(r8), Save, Pointer, Dimension(:)     :: tm_dens_old
  real(r8), Save, Pointer, Dimension(:)     :: tm_dens
  real(r8), Save, Pointer, Dimension(:,:)   :: PC_Strain
  real(r8), Save, Pointer, Dimension(:,:)   :: Elastic_Strain
  real(r8), Save, Pointer, Dimension(:)     :: Rotation_Magnitude
  !
  real(r8), Save, Pointer, Dimension(:)     :: Eff_Stress_Cell_old
  real(r8), Save, Pointer, Dimension(:,:)   :: Eff_Stress_IP_old
  !
  real(r8), Save, Pointer, Dimension(:,:)   :: Dev_Stress_Cell
  real(r8), Save, Pointer, Dimension(:,:,:) :: Dev_Stress_IP
  !
  ! Thermo-elastic properties
  real(r8), Save, Pointer, Dimension(:) :: Lame1
  real(r8), Save, Pointer, Dimension(:) :: Lame2
  real(r8), Save, Pointer, Dimension(:) :: CTE
  ! Nodal displacements
  real(r8), Save, Pointer, Dimension(:,:) :: Displacement
  ! Gap displacements and scaled forces
  real(r8), Save, Pointer, Dimension(:,:) :: Node_Gap => NULL()
  real(r8), Save, Pointer, Dimension(:,:) :: Node_Norm_Trac => NULL()
  ! Right hand side of system of equations
  real(r8), Save, Pointer, Dimension(:) :: RHS
  ! Source terms needed for contact equations
  real(r8), Save, Pointer, Dimension(:) :: Src
  ! Physics flags
  logical, Save                          :: solid_mechanics
  logical, Save                          :: plasticity
  logical, Save                          :: solid_mechanics_body_force
  logical, Save, Pointer, dimension(:)   :: Solid_Mask, gap_cell_mask
  ! Interface BC flag
  logical, Save                          :: mech_interface
  ! Contact parameters
  real(r8), Save                          :: contact_distance, contact_norm_trac, contact_penalty
  ! Maximum plsatic strain increment
  real(r8), Save                          :: strain_limit
  ! System scale factor
  real(r8), save, pointer, dimension(:)   :: cscale

CONTAINS

  FUNCTION MECH_COLLATE (Mech)
    !==================================================================
    ! Purpose(s):
    !   Collate a distributed mech data type into a single large mech data on IO PE
    !==================================================================
    use parallel_info_module, only: p_info
    use parameter_module,     only: ncells_tot

    ! Arguments
    type(MECH_DATA), intent(IN) :: Mech
    type(MECH_DATA), pointer    :: Mech_Collate
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (p_info%IOP) then
       ALLOCATE(Mech_Collate)
       ALLOCATE(Mech_Collate%Total_Strain(ncomps,ncells_tot))
       ALLOCATE(Mech_Collate%Elastic_Stress(ncomps,ncells_tot))
       ALLOCATE(Mech_Collate%Plastic_Strain(ncomps,ncells_tot))
       ALLOCATE(Mech_Collate%Plastic_Strain_Rate(ncells_tot))
    else
       ALLOCATE(Mech_Collate)
       ALLOCATE(Mech_Collate%Total_Strain(ncomps,0))
       ALLOCATE(Mech_Collate%Elastic_Stress(ncomps,0))
       ALLOCATE(Mech_Collate%Plastic_Strain(ncomps,0))
       ALLOCATE(Mech_Collate%Plastic_Strain_Rate(0))
    end if
    
    call COLLATE(Mech_Collate, Mech)

  END FUNCTION MECH_COLLATE

  SUBROUTINE COLLATE_MECH (Collated_Mech, Local_Mech)
    !==================================================================
    ! Purpose(s):
    !   Collate a distributed mech data type into a single large mech on IO PE
    !==================================================================
    use pgslib_module, only: PGSLib_COLLATE

    ! Arguments
    type(MECH_DATA), intent(IN)  :: Local_Mech
    type(MECH_DATA), intent(OUT) :: Collated_Mech
    
    ! Local variables
    integer :: n

    do n = 1,ncomps
       call PGSLib_COLLATE (Collated_Mech%Total_Strain(n,:),     Local_Mech%Total_Strain(n,:))
       call PGSLib_COLLATE (Collated_Mech%Elastic_Stress(n,:),   Local_Mech%Elastic_Stress(n,:))
       call PGSLib_COLLATE (Collated_Mech%Plastic_Strain(n,:),   Local_Mech%Plastic_Strain(n,:))
    end do
    call PGSLib_COLLATE (Collated_Mech%Plastic_Strain_Rate,      Local_Mech%Plastic_Strain_Rate)

  END SUBROUTINE COLLATE_MECH

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

  subroutine read_SM_data (unit, version)

    use restart_utilities, only: read_dist_array
    use node_operator_module, only: nipc
    use mesh_module, only: pcell => unpermute_mesh_vector, pnode => unpermute_vertex_vector

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

  subroutine skip_SM_data (unit, version)
    use restart_utilities, only: skip_records
    integer, intent(in) :: unit, version
    call skip_records (unit, 265, 'SKIP_SM_DATA: error skipping the solid mechanics data')
  end subroutine skip_SM_data

  SUBROUTINE PERMUTE_MECH (Permuted_Mech, Orig_Mech, Permuter, SCOPE)
    !==================================================================
    ! Purpose(s):
    !   Permute mech according the Permuter vector
    !==================================================================
    use parallel_scope

    ! Arguments
    type(MECH_DATA), intent(IN   ) :: Orig_Mech
    type(MECH_DATA), intent(INOUT) :: Permuted_Mech
    integer, dimension(:), intent(IN) :: Permuter
    type (PL_SCOPE), OPTIONAL, intent(IN   ) :: SCOPE

    ! Local variables
    type (PL_SCOPE) :: Desired_Scope

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Default scope is global
    if (PRESENT(SCOPE)) then
       Desired_Scope = SCOPE
    else
       Desired_Scope = GLOBAL_SCOPE
    end if

    if (DESIRED_SCOPE == GLOBAL_SCOPE) then
       call PERMUTE_MECH_GLOBAL(Permuted_Mech, Orig_Mech, Permuter)
    end if

    if (DESIRED_SCOPE == LOCAL_SCOPE) then
       call PERMUTE_MECH_LOCAL (Permuted_Mech, Orig_Mech, Permuter)
    end if

  end SUBROUTINE PERMUTE_MECH

  SUBROUTINE PERMUTE_MECH_GLOBAL (Permuted_Mech, Orig_Mech, Permuter)
    !==================================================================
    ! Purpose(s):
    !   Permute mech according the Permuter vector, global version
    !==================================================================
    use pgslib_module,    only: PGSLib_Permute,    &
                                PGSLIB_Deallocate_Trace, &
                                PGSLib_GS_Trace

    ! Arguments
    type(MECH_DATA), intent(IN   ) :: Orig_Mech
    type(MECH_DATA), intent(INOUT) :: Permuted_Mech
    integer, dimension(:), intent(IN) :: Permuter
    
    ! Local variables
    integer :: n
    type (PGSLib_GS_Trace), POINTER :: Mech_Trace
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    NULLIFY(Mech_Trace)

    do n = 1,ncomps
       call PGSLib_PERMUTE (DEST   = Permuted_Mech%Total_Strain(n,:),          &
                         SOURCE = Orig_Mech%Total_Strain(n,:),              &
                         INDEX  = Permuter,       &
                         TRACE  = Mech_Trace)
       call PGSLib_PERMUTE (DEST   = Permuted_Mech%Elastic_Stress(n,:),      &
                         SOURCE = Orig_Mech%Elastic_Stress(n,:),              &
                         INDEX  = Permuter,       &
                         TRACE  = Mech_Trace)
       call PGSLib_PERMUTE (DEST   = Permuted_Mech%Plastic_Strain(n,:),         &
                         SOURCE = Orig_Mech%Plastic_Strain(n,:),                 &
                         INDEX  = Permuter,       &
                         TRACE  = Mech_Trace)
    end do
    call PGSLib_PERMUTE (DEST   = Permuted_Mech%Plastic_Strain_Rate,     &
                         SOURCE = Orig_Mech%Plastic_Strain_Rate,             & 
                         INDEX  = Permuter,       &
                         TRACE  = Mech_Trace)

    ! Done with the trace
    call PGSLib_DEALLOCATE_TRACE (Mech_Trace)

  END SUBROUTINE PERMUTE_MECH_GLOBAL

  SUBROUTINE PERMUTE_MECH_LOCAL (Permuted_Mech, Orig_Mech, Permuter)
    !==================================================================
    ! Purpose(s):
    !   Permute mech according the Permuter vector, local version
    !   The Permuter vector refers to local indices.
    !   The input and output vectors must have the same size.
    !==================================================================

    ! Arguments
    type(MECH_DATA), intent(IN   ) :: Orig_Mech
    type(MECH_DATA), intent(INOUT) :: Permuted_Mech
    integer, dimension(:), intent(IN) :: Permuter
    
    ! Local variables
    integer :: cell, n
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do cell = 1, SIZE(Permuter)
       do n = 1,ncomps
          Permuted_Mech%Total_Strain(n,Permuter(cell)) = Orig_Mech%Total_Strain(n,cell)
          Permuted_Mech%Elastic_Stress(n,Permuter(cell)) = Orig_Mech%Elastic_Stress(n,cell)
          Permuted_Mech%Plastic_Strain(n,Permuter(cell)) = Orig_Mech%Plastic_Strain(n,cell)
       end do
       Permuted_Mech%Plastic_Strain_Rate(Permuter(cell)) = Orig_Mech%Plastic_Strain_Rate(cell)
    end do

  END SUBROUTINE PERMUTE_MECH_LOCAL


  SUBROUTINE GAP_NODE_DISPLACEMENT ()
    !
    ! Purpose: Transfer gap normal displacement data from Node_Displacement_BC
    ! to an nnodes array for plotting
    !
    use mech_bc_data_module
    use parameter_module, only: nnodes

    ! Local variables
    integer :: inode, nnum, n, isize, iint
    
    if (.not. (ASSOCIATED(Node_Gap))) then
       isize = 0
       do n = 1, SIZE(Interface_List)
          if (Interface_List(n) /= 0) isize = isize + 1
       end do
       allocate (Node_Gap(nnodes,isize))
       allocate (Node_Norm_Trac(nnodes,isize))
    end if

    Node_Gap = 0.0
    Node_Norm_Trac = 0.0

    do inode = 1, SIZE(Node_Displacement_BC%Node)
       nnum = Node_Displacement_BC%Node(inode)
       ! For now, just calculate one value for nodes with more than one gap
       select case (Node_Displacement_BC%Combination(inode))
       case(ONE_NORM_CONST)
          iint = Node_Displacement_BC%Interface(1,inode)
          Node_Gap(nnum, iint) = Node_Displacement_BC%Gap_Disp(inode,1)
          Node_Norm_Trac(nnum, iint) = Node_Displacement_BC%Normal_Traction(inode,1)
       case(TWO_NORM_CONST)
          do n = 1,2
             iint = Node_Displacement_BC%Interface(n,inode)
             Node_Gap(nnum, iint) = Node_Displacement_BC%Gap_Disp(inode,n)
             Node_Norm_Trac(nnum, iint) = Node_Displacement_BC%Normal_Traction(inode,n)
          end do

       case(THREE_NORM_CONST)
          do n = 1,3
             iint = Node_Displacement_BC%Interface(n,inode)
             Node_Gap(nnum, iint) = Node_Displacement_BC%Gap_Disp(inode,n)
             Node_Norm_Trac(nnum, iint) = Node_Displacement_BC%Normal_Traction(inode,n)
          end do
       case(ONE_D_ONE_NC)
          iint = Node_Displacement_BC%Interface(2,inode)
          Node_Gap(nnum, iint) = Node_Displacement_BC%Gap_Disp(inode,2)
          Node_Norm_Trac(nnum, iint) = Node_Displacement_BC%Normal_Traction(inode,2)

       case(TWO_D_ONE_NC)
          iint = Node_Displacement_BC%Interface(3,inode)
          Node_Gap(nnum, iint) = Node_Displacement_BC%Gap_Disp(inode,3)
          Node_Norm_Trac(nnum, iint) = Node_Displacement_BC%Normal_Traction(inode,3)

       case(ONE_D_TWO_NC)
          do n = 2,3
             iint = Node_Displacement_BC%Interface(n,inode)
             Node_Gap(nnum, iint) = Node_Displacement_BC%Gap_Disp(inode,n)
             Node_Norm_Trac(nnum, iint) = Node_Displacement_BC%Normal_Traction(inode,n)
          end do

       end select      
    end do

  end SUBROUTINE GAP_NODE_DISPLACEMENT

end Module SOLID_MECHANICS_DATA
