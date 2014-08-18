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
  !use parameter_module, only: string_len, ncomps, ndim
  use parameter_module, only: string_len
  use solid_mechanics_mesh, only: ncomps, ndim
  implicit none
  private

  ! Public data
  Public :: Solid_Mask,                       &
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
            MECH_DATA,                        &
            thermo_elastic_precond_iter,      &
            RHS,                              &
            Src,                              &
            plasticity,                       &
            Node_Gap,                         &
            Node_Norm_Trac,                   &
            cscale

  !
  ! public procedures
  public :: GAP_NODE_DISPLACEMENT

  INTERFACE COLLATE
     MODULE PROCEDURE COLLATE_MECH
  END INTERFACE
  
    !
     ! Precondition iteration count
  integer, save :: thermo_elastic_precond_iter
  !
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
  logical, Save                          :: plasticity
  logical, Save, Pointer, dimension(:)   :: Solid_Mask, gap_cell_mask
  ! Interface BC flag
  logical, Save                          :: mech_interface
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
