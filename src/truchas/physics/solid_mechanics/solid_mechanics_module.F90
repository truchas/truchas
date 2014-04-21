#include "f90_assert.fpp"

Module SOLID_MECHANICS_MODULE
  !-----------------------------------------------------------------------------
  ! Purpose:
  !
  !   Calculate the thermomechanical response of all relevant solid materials 
  !   present in a given problem
  !
  ! Public Interfaces:
  !
  !   Call SOLID_MECHANICS_ALLOCATE ()
  !   Call THERMO_MECHANICS ()2
  !   Call STRESS_STRAIN_INVARIANTS ()
  !   Call SOLID_MECH_INIT ()
  !
  ! Contains:  SOLID_MECHANICS_ALLOCATE ()
  !            SOLID_MECH_INIT
  !            THERMO_MECHANICS
  !               MATERIAL_PROPERTIES ()
  !	          MATERIAL_DISPLACEMENTS ()
  !               RHS_THERMO_MECH()
  !               MECH_PRECONDITION)
  !                  PRECON_DISPLACEMENT_GRADIENTS
  !                     TET_GRADIENT_FACTORS()
  !               ELAS_VP_RESIDUAL ()
  !               VP_MATVEC 
  !            STRESS_STRAIN_INVARIANTS
  !             
  ! Authors:  Dave Korzekwa (dak@lanl.gov), Mark Schraad (schraad@lanl.gov)
  !-----------------------------------------------------------------------------
  use kinds, only: r8
  Use var_vector_module
  use solid_mechanics_data
  use nonlinear_solution, only: NK_SOLUTION_FIELD
  use truchas_logging_services
  Implicit None
  Private

  ! Public procedures
  Public :: SOLID_MECHANICS_ALLOCATE,   &
            THERMO_MECHANICS,           &
            STRESS_STRAIN_INVARIANTS,   &
            SOLID_MECH_INIT

  ! Reference temperature for thermal stress/strain
  real(r8), Save, Pointer, Dimension(:)   :: Ref_Temp
  ! Gradient data for preconditioners
  type(real_var_vector), save, pointer, Dimension(:)                :: M1
  type(real_var_vector), save, pointer, Dimension(:)                :: M2
  !
  ! Instantiate an NK space-time data type.
  type(NK_SOLUTION_FIELD),  save    :: VP_Solution_Field

  !-----------------------------------------------------------------------------

Contains

  Subroutine SOLID_MECHANICS_ALLOCATE ()
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    Allocate space for and initialize the material property, displacement,
    !    strain, and stress arrays
    !
    !---------------------------------------------------------------------------
    Use parameter_module, Only: ncells, ndim, nnodes, ncomps
    Use node_operator_module, Only: nipc
    use timing_tree

    Integer :: ip

    !---------------------------------------------------------------------------

    ! If not calculating solid mechanics, then skip the rest
    If (.not. solid_mechanics) Return
    call start_timer("Solid Mechanics")

    Allocate(Thermal_Strain(ncomps,ncells))
    Allocate(Thermal_Strain_Inc(ncomps,ncells))
    Allocate(tm_dens_old(ncells))
    Allocate(tm_dens(ncells))
    Allocate(PC_Strain(ncomps,ncells))
    Allocate(Elastic_Strain(ncomps,ncells))
    Allocate(Rotation_Magnitude(ncells))

    Allocate(SMech_Cell%Total_Strain(ncomps,ncells))
    Allocate(SMech_Cell%Elastic_Stress(ncomps,ncells))
    Allocate(SMech_Cell%Plastic_Strain(ncomps,ncells))
    Allocate(SMech_Cell%Plastic_Strain_Rate(ncells))

    Allocate(SMech_Cell_Old%Total_Strain(ncomps,ncells))
    Allocate(SMech_Cell_Old%Elastic_Stress(ncomps,ncells))
    Allocate(SMech_Cell_Old%Plastic_Strain(ncomps,ncells))
    Allocate(SMech_Cell_Old%Plastic_Strain_Rate(ncells))

    Allocate(SMech_IP(nipc))
    Allocate(SMech_IP_Old(nipc))

    do ip = 1,nipc
       Allocate(SMech_IP(ip)%Total_Strain(ncomps,ncells))
       Allocate(SMech_IP(ip)%Elastic_Stress(ncomps,ncells))
       Allocate(SMech_IP(ip)%Plastic_Strain(ncomps,ncells))
       Allocate(SMech_IP(ip)%Plastic_Strain_Rate(ncells))

       Allocate(SMech_IP_Old(ip)%Total_Strain(ncomps,ncells))
       Allocate(SMech_IP_Old(ip)%Elastic_Stress(ncomps,ncells))
       Allocate(SMech_IP_Old(ip)%Plastic_Strain(ncomps,ncells))
       Allocate(SMech_IP_Old(ip)%Plastic_Strain_Rate(ncells))
    end do

    Allocate(Eff_Stress_Cell_old(ncells))
    Allocate(Eff_Stress_IP_old(ncells,nipc))

    Allocate(Dev_Stress_Cell(ncomps,ncells))
    Allocate(Dev_Stress_IP(ncomps,ncells,nipc))

    Allocate(Lame1(ncells))
    Allocate(Lame2(ncells))
    Allocate(CTE(ncells))

    Allocate(Displacement(ndim,nnodes))

    Allocate(RHS(ndim*nnodes))
    Allocate(Src(ndim*nnodes))
    Allocate(Ref_Temp(ncells))
    Allocate(Solid_Mask(ncells))
    Allocate(gap_cell_mask(ncells))

    Thermal_Strain(:,:)                    = 0.0_r8
    Thermal_Strain_Inc(:,:)                = 0.0_r8
    tm_dens_old(:)                         = 0.0_r8
    tm_dens(:)                             = 0.0_r8
    PC_Strain(:,:)                         = 0.0_r8
    Elastic_Strain(:,:)                    = 0.0_r8
    Rotation_Magnitude(:)                  = 0.0_r8

    SMech_Cell%Total_Strain(:,:)           = 0.0_r8
    SMech_Cell%Elastic_Stress(:,:)         = 0.0_r8
    SMech_Cell%Plastic_Strain(:,:)         = 0.0_r8
    SMech_Cell%Plastic_Strain_Rate(:)      = 0.0_r8

    SMech_Cell_Old%Total_Strain(:,:)       = 0.0_r8
    SMech_Cell_Old%Elastic_Stress(:,:)     = 0.0_r8
    SMech_Cell_Old%Plastic_Strain(:,:)     = 0.0_r8
    SMech_Cell_Old%Plastic_Strain_Rate(:)  = 0.0_r8

    do ip = 1,nipc
       SMech_IP(ip)%Total_Strain(:,:)          = 0.0_r8
       SMech_IP(ip)%Elastic_Stress(:,:)        = 0.0_r8
       SMech_IP(ip)%Plastic_Strain(:,:)        = 0.0_r8
       SMech_IP(ip)%Plastic_Strain_Rate(:)     = 0.0_r8

       SMech_IP_Old(ip)%Total_Strain(:,:)      = 0.0_r8
       SMech_IP_Old(ip)%Elastic_Stress(:,:)    = 0.0_r8
       SMech_IP_Old(ip)%Plastic_Strain(:,:)    = 0.0_r8
       SMech_IP_Old(ip)%Plastic_Strain_Rate(:) = 0.0_r8
    end do

    Eff_Stress_Cell_old(:)                 = 0.0_r8
    Eff_Stress_IP_old(:,:)                 = 0.0_r8

    Dev_Stress_Cell(:,:)                   = 0.0_r8
    Dev_Stress_IP(:,:,:)                   = 0.0_r8

    Lame1(:)                               = 0.0_r8
    Lame2(:)                               = 0.0_r8
    CTE(:)                                 = 0.0_r8

    Displacement(:,:)                      = 0.0_r8

    RHS(:)                                 = 0.0_r8
    Src(:)                                 = 0.0_r8
    Ref_Temp(:)                            = 0.0_r8

    call stop_timer("Solid Mechanics")

  End Subroutine SOLID_MECHANICS_ALLOCATE
  !
  !-----------------------------------------------------------------------------
  !
  Subroutine SOLID_MECH_INIT
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !  Calculate initial displacements, stresses and strains and strain rates 
    ! resulting from the initial or restart temperature field and BCs
    !
    !---------------------------------------------------------------------------

    Use parameter_module,     Only: ndim, nnodes, ncells, nvc, nfc, nmat, mat_slot
    use node_operator_module, only: cv_init, CV_Internal, nipc
    use node_op_setup_module, only: ALLOCATE_CONTROL_VOLUME, CELL_CV_FACE, BOUNDARY_CV_FACE
    use timing_tree
    use UbikSolve_module
    use viscoplasticity,      only: VISCOPLASTICITY_INIT, MATERIAL_STRESSES, MATERIAL_STRAINS, VISCOPLASTIC_STRAIN_RATE
    use linear_module,        only: LINEAR_GRAD
    Use gs_module,            Only: EN_GATHER, EN_SUM_SCATTER
    use restart_variables,    only: restart, have_solid_mechanics_data, ignore_solid_mechanics
    use truchas_logging_services
    use fluid_data_module,    only: isImmobile
    use matl_module,          only: Matl
    use mesh_module,            only: Cell, Mesh, GAP_ELEMENT_1
    use solid_mech_constraints, only: FACE_GAP_INITIALIZE, FACE_GAP_UPDATE
    use pgslib_module,        only: PGSLib_Global_MAXVAL
    Use zone_module,          Only: Zone

    logical :: have_initial_state
    real(r8), Dimension(ndim*nnodes)      :: Solution
    real(r8), Dimension(nvc,ncells,ndim)  :: U
    real(r8), Dimension(ndim,ncells)      :: Arow
    real(r8), Dimension(ndim,ndim,ncells) :: Ugrad
    integer :: idim, inodes, ip, i, j, k, status, icell, imat, islot
    real(r8), pointer, dimension(:) :: Lame2_Node, Nvol, L2tmp
    real(r8), dimension(nfc) :: htemp

    character(len=200) :: errmsg
    character(128) :: message

    !---------------------------------------------------------------------------

    ! If not calculating solid mechanics, then skip the rest
    If (.not. solid_mechanics) Return
    
    ! Allocate and calculate control volume structures
    
    if (.not. cv_init) then
       call start_timer("Initialization")
       call start_timer("Solid Mech Init")

       call ALLOCATE_CONTROL_VOLUME()
       call CELL_CV_FACE()
       call BOUNDARY_CV_FACE()
       cv_init = .true.
       
       call stop_timer("Solid Mech Init")
       call stop_timer("Initialization")
    end if

    ! Determine if plastic flow is to be included
    call viscoplasticity_init (plasticity)

    call define_tm_density_property (status, errmsg)
    if (status /= 0) call TLS_fatal ('SOLID_MECH_INIT: ' // trim(errmsg))
    
    gap_cell_mask = (mesh(:)%Cell_Shape >= GAP_ELEMENT_1)
    
    ! Calculate the volume-averaged, cell-centered solid material properties
    Call MATERIAL_PROPERTIES ()
    
    ! Scale factor calculated for each node based on the surface area of the control volume and
    ! the second Lame' constant

    allocate (cscale(ndim*nnodes))
    
    ! Get second elastic constant at the nodes
    allocate(Lame2_Node(nnodes), stat=status)
    if (status /= 0) call TLS_panic ( 'MECH_PRECONDITION: allocation error: Lame2_Node')
    allocate(Nvol(nnodes), stat=status)
    if (status /= 0) call TLS_panic ( 'MECH_PRECONDITION: allocation error: Nvol')
    allocate(L2tmp(ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'MECH_PRECONDITION: allocation error: L2tmp')
    
    ! Get nodal values for Lame'2
    call EN_SUM_SCATTER (Nvol, Cell%Volume)
    L2tmp = Lame2*Cell%Volume
    where (gap_cell_mask) L2tmp = 0.0
    call EN_SUM_SCATTER (Lame2_Node, L2tmp)
    Lame2_Node = Lame2_Node/Nvol

    do i = 1, nfc
       htemp(i) = MAXVAL(Cell%Halfwidth(i))
    end do
    
    ! SCale factor for each node
    do idim = 1, ndim
       cscale(idim:(ndim*(nnodes-1)+idim):ndim) = Lame2_node * 0.001 * PGSLib_Global_MAXVAL(htemp)
    end do

    where (cscale(:) == 0.0) cscale = 1.0

    deallocate(L2tmp)
    deallocate(Nvol)
    deallocate(Lame2_Node)

!    cscale = PGSLib_Global_MAXVAL(Lame2) * 0.0001

    ! Initialize mask to find all cells that contain any solid (immobile) material
    Solid_Mask = .false.
    do icell = 1, ncells
       do imat = 1,nmat
          if (isImmobile(imat)) then
             do islot = 1, mat_slot
                if (Matl(islot)%Cell(icell)%Id == imat) then
                   Solid_Mask(icell) = .true.
                end if
             end do
          end if
       end do
    end do

    !! Calculate the initial thermo-elastic state if necessary.
    have_initial_state = restart .and. have_solid_mechanics_data .and. .not.ignore_solid_mechanics
    INITIAL_STATE: if (.not.have_initial_state) then
       call start_timer("Solid Mechanics")

       ! Solve for the initial elastic stress, strain and strain rate fields
       call TLS_info ('')
       call TLS_info (' Calculating initial thermo-elastic state.', advance=.false.)

       call compute_cell_property ('TM reference density', zone%temp, tm_dens_old)
       where (gap_cell_mask) tm_dens_old = 1.0

       ! Get initial displacements
       Call MATERIAL_DISPLACEMENTS()

       Do inodes = 1, nnodes
          Do idim = 1, ndim
             Solution((inodes - 1)*ndim + idim) = Displacement(idim,inodes)
          End Do
       End Do

       ! Calculate elastic stress tensor for all integration points in all cells

       do idim = 1, ndim
          call EN_GATHER (U(:,:,idim),Solution(idim:(ndim*(nnodes-1)+idim):ndim))
       end do

       IP_STRESS_LOOP: do ip = 1,nipc
          ! Calculate displacement gradient and strain components
          Ugrad(:,:,:) = 0.0
          do i=1,ndim
             call LINEAR_GRAD(ncells,CV_Internal%Face_Coord(:,ip,:),U(:,:,i),Arow(:,:))
             do j = 1,ndim
                do k=1,ndim
                   Ugrad(i,j,:)= Ugrad(i,j,:) + Arow(k,:) * CV_Internal%Face_Ijac(k,j,ip,:)
                end do
             end do
          end do
          ! Total strains
          SMech_IP(ip)%Total_Strain(1,:) = Ugrad(1,1,:)
          SMech_IP(ip)%Total_Strain(2,:) = Ugrad(2,2,:)
          SMech_IP(ip)%Total_Strain(3,:) = Ugrad(3,3,:)
          SMech_IP(ip)%Total_Strain(4,:) = 0.5 * (Ugrad(1,2,:) + Ugrad(2,1,:))
          SMech_IP(ip)%Total_Strain(5,:) = 0.5 * (Ugrad(1,3,:) + Ugrad(3,1,:))
          SMech_IP(ip)%Total_Strain(6,:) = 0.5 * (Ugrad(3,2,:) + Ugrad(2,3,:))
          Elastic_Strain = SMech_IP(ip)%Total_Strain(:,:) - Thermal_Strain - PC_Strain
          call MATERIAL_STRESSES(Elastic_Strain, SMech_IP(ip)%Elastic_Stress(:,:))
       end do IP_STRESS_LOOP
       
       ! Initial effective stress.
       do ip = 1,nipc
          Eff_Stress_IP_old(:,ip) = sqrt(((SMech_IP(ip)%Elastic_Stress(1,:) - SMech_IP(ip)%Elastic_Stress(2,:))**2 + &
               (SMech_IP(ip)%Elastic_Stress(2,:) - SMech_IP(ip)%Elastic_Stress(3,:))**2 + &
               (SMech_IP(ip)%Elastic_Stress(3,:) - SMech_IP(ip)%Elastic_Stress(1,:))**2 + &
               6.0 * (SMech_IP(ip)%Elastic_Stress(4,:)**2 + SMech_IP(ip)%Elastic_Stress(5,:)**2 + &
               SMech_IP(ip)%Elastic_Stress(6,:)**2))/2.)
       end do

       ! Initial strain rate
       do ip = 1,nipc
          call viscoplastic_strain_rate (Eff_Stress_IP_old(:,ip), zone%temp, SMech_IP(ip)%Plastic_Strain_Rate)
       end do

       ! Initialize cell centered stresses, strains and strain rates for output
       ! Calculate the cell-centered solid material total strain field
       Call MATERIAL_STRAINS (Solution, SMech_Cell%Total_Strain,Rotation=Rotation_Magnitude)

       ! Calculate the cell-centered solid material elastic strain field, noting 
       ! that the total strains are decomposed into elastic and thermal strains
       Elastic_Strain(:,:) = SMech_Cell%Total_Strain(:,:) - Thermal_Strain(:,:) - PC_Strain(:,:)

       ! Calculate the cell-centered elastic stress field
       Call MATERIAL_STRESSES (Elastic_Strain, SMech_Cell%Elastic_Stress)

       ! Cell centered effective stress
       Eff_Stress_Cell_old(:) = sqrt(((SMech_Cell%Elastic_Stress(1,:) - SMech_Cell%Elastic_Stress(2,:))**2 + &
            (SMech_Cell%Elastic_Stress(2,:) - SMech_Cell%Elastic_Stress(3,:))**2 + &
            (SMech_Cell%Elastic_Stress(3,:) - SMech_Cell%Elastic_Stress(1,:))**2 + &
            6.0 * (SMech_Cell%Elastic_Stress(4,:)**2 + SMech_Cell%Elastic_Stress(5,:)**2 + &
            SMech_Cell%Elastic_Stress(6,:)**2))/2.)

       ! Cell centered plastic strain rate
       call viscoplastic_strain_rate (Eff_Stress_Cell_old, zone%temp, SMech_Cell%Plastic_Strain_Rate)

       call TLS_info (' Done.')
       
       call stop_timer("Solid Mechanics")
 
!       thermo_elastic_iterations   = NKuser(NK_DISPLACEMENT)%linear_tot
!       viscoplastic_iterations = NKuser(NK_DISPLACEMENT)%Newton_tot
       call TLS_info ('')
       write(message,'(5x,i0,a)') thermo_elastic_iterations, ' Thermo-elastic iterations (linear)'
       call TLS_info (message)
       write(message,'(5x,i0,a)') viscoplastic_iterations, ' Thermo-elastic iterations (nonlinear)'
       call TLS_info (message)

    else

       call compute_cell_property ('TM density', zone%temp, tm_dens)
       where (gap_cell_mask) tm_dens = 1.0
    end if INITIAL_STATE

    ! Put gap displacements into a nodal array
    call GAP_NODE_DISPLACEMENT ()
    ! Set up face and node gap arrays and connectivity
    call FACE_GAP_INITIALIZE()
    ! Calculate gap displacements for htc_gap bcs
    call FACE_GAP_UPDATE()

    status=0

  end Subroutine SOLID_MECH_INIT

  !-----------------------------------------------------------------------------

  Subroutine THERMO_MECHANICS ()
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    Control the solid mechanics calculation cycle
    !
    !---------------------------------------------------------------------------
    !
    Use parameter_module,     Only: ndim, nnodes, ncells, ncomps, nmat, mat_slot
    Use node_operator_module, Only: nipc
    use timing_tree
    use viscoplasticity,      only: MATERIAL_STRESSES, MATERIAL_STRAINS, &
                                    VISCOPLASTIC_STRAIN_RATE, PLASTIC_STRAIN_INCREMENT,&
                                    DEVIATORIC_STRESS
    use fluid_data_module,    only: isImmobile
    use matl_module,          only: Matl
    use solid_mech_constraints, only: FACE_GAP_UPDATE
    use zone_module, only: zone

    integer :: idim, inodes, ip, icell, imat, islot
    real(r8), Dimension(ndim*nnodes) :: Temp
    real(r8), Dimension(ncomps,ncells) :: Pl_Inc_Cell, Dummy

    !---------------------------------------------------------------------------

    ! If not calculating solid mechanics, then skip the rest
    If (.not. solid_mechanics) Return
    
    call start_timer("Solid Mechanics")

    ! Update solid material mask
    Solid_Mask = .false.
    do icell = 1, ncells
       do imat = 1,nmat
          if (isImmobile(imat)) then
             MAT_SLOT_LOOP: do islot = 1, mat_slot
                if (Matl(islot)%Cell(icell)%Id == imat) then
                   Solid_Mask(icell) = .true.
                end if
             end do MAT_SLOT_LOOP
          end if
       end do
    end do

    ! Set stresses, strains and plastic strain rate at the beginning of the time step

    do ip = 1,nipc
       SMech_IP_Old(ip)%Plastic_Strain_Rate = SMech_IP(ip)%Plastic_Strain_Rate
       SMech_IP_Old(ip)%Elastic_Stress      = SMech_IP(ip)%Elastic_Stress
       SMech_IP_Old(ip)%Plastic_Strain      = SMech_IP(ip)%Plastic_Strain
       SMech_IP_Old(ip)%Total_Strain        = SMech_IP(ip)%Total_Strain
    end do
    !
    SMech_Cell_Old%Plastic_Strain_Rate = SMech_Cell%Plastic_Strain_Rate
    SMech_Cell_Old%Elastic_Stress = SMech_Cell%Elastic_Stress
    SMech_Cell_Old%Plastic_Strain = SMech_Cell%Plastic_Strain
    SMech_Cell_Old%Total_Strain = SMech_Cell%Total_Strain
    !
    do ip = 1,nipc
       Eff_Stress_IP_old(:,ip) = sqrt(((SMech_IP_Old(ip)%Elastic_Stress(1,:) - SMech_IP_Old(ip)%Elastic_Stress(2,:))**2 + &
            (SMech_IP_Old(ip)%Elastic_Stress(2,:) - SMech_IP_Old(ip)%Elastic_Stress(3,:))**2 + &
            (SMech_IP_Old(ip)%Elastic_Stress(3,:) - SMech_IP_Old(ip)%Elastic_Stress(1,:))**2 + &
            6.0 * (SMech_IP_Old(ip)%Elastic_Stress(4,:)**2 + SMech_IP_Old(ip)%Elastic_Stress(5,:)**2 + &
            SMech_IP_Old(ip)%Elastic_Stress(6,:)**2))/2.)
    end do

    ! Update tm_dens_old
    tm_dens_old(:) = tm_dens(:)
    tm_dens(:) = 0.0
    
    do ip = 1,nipc
       call viscoplastic_strain_rate (Eff_Stress_IP_old(:,ip), zone%temp, SMech_IP_Old(ip)%Plastic_Strain_Rate)
    end do

    ! Deviatoric stress components for all integration points - used to set the direction of the 
    ! strain increment for the time step

    do ip = 1,nipc
       do icell = 1,ncells
          call DEVIATORIC_STRESS(SMech_IP_Old(ip)%Elastic_Stress(:,icell),Dev_Stress_IP(:,icell,ip))
       end do
    end do

    ! Cell centered effective stress
    Eff_Stress_Cell_old(:) = sqrt(((SMech_Cell_Old%Elastic_Stress(1,:) - SMech_Cell_Old%Elastic_Stress(2,:))**2 + &
         (SMech_Cell_Old%Elastic_Stress(2,:) - SMech_Cell_Old%Elastic_Stress(3,:))**2 + &
         (SMech_Cell_Old%Elastic_Stress(3,:) - SMech_Cell_Old%Elastic_Stress(1,:))**2 + &
         6.0 * (SMech_Cell_Old%Elastic_Stress(4,:)**2 + SMech_Cell_Old%Elastic_Stress(5,:)**2 + &
         SMech_Cell_Old%Elastic_Stress(6,:)**2))/2.)

    ! Deviatoric stress components for cell centered quantities

    do icell = 1,ncells
       call DEVIATORIC_STRESS(SMech_Cell_Old%Elastic_Stress(:,icell),Dev_Stress_Cell(:,icell))
    end do

    ! Calculate the volume-averaged, cell-centered solid material properties
    Call MATERIAL_PROPERTIES ()

    ! End initialization

    ! Calculate the displacements elastic stresses and plastic strains and rate at the end of the time step.
    Call MATERIAL_DISPLACEMENTS ()

    ! Store the displacement vector in a temporary one-dimensional array
    Do inodes = 1, nnodes
       Do idim = 1, ndim
          Temp((inodes - 1) * ndim + idim) = Displacement(idim,inodes)
       End Do
    End Do

    ! Put gap displacements into a nodal array
    call GAP_NODE_DISPLACEMENT ()

    ! Calculate the cell-centered total strain field
    Call MATERIAL_STRAINS (Temp, SMech_Cell%Total_Strain,Rotation=Rotation_Magnitude)

    ! Get plastic strain increment at cell centroid
    Call PLASTIC_STRAIN_INCREMENT(Pl_Inc_Cell,                         &
                                  SMech_Cell%Plastic_Strain_Rate,      &
                                  SMech_Cell%Elastic_Stress,           &
                                  Dummy,                               &
                                  Dev_Stress_Cell,                     &
                                  SMech_Cell_Old%Elastic_Stress,       &
                                  SMech_Cell_Old%Plastic_Strain_Rate,  &
                                  Eff_Stress_Cell_old,                 &
                                  SMech_Cell%Total_Strain,             &
                                  SMech_Cell_Old%Total_Strain,         &
                                  SMech_Cell_Old%Plastic_Strain)

    ! Update total plastic strain
    SMech_Cell%Plastic_Strain(:,:) = SMech_Cell%Plastic_Strain(:,:) + Pl_Inc_Cell


    ! Calculate the cell-centered elastic strain field, noting  that the total strains
    ! are decomposed into elastic, plastic and thermal strains
    Elastic_Strain(:,:) = SMech_Cell%Total_Strain(:,:) - Thermal_Strain(:,:) - PC_Strain(:,:) - &
                          SMech_Cell%Plastic_Strain(:,:)

    ! Calculate the cell-centered elastic stress field
    Call MATERIAL_STRESSES (Elastic_Strain, SMech_Cell%Elastic_Stress)

    ! Calculate gap displacements for htc_gap bcs
    call FACE_GAP_UPDATE()

    call stop_timer("Solid Mechanics")

  End Subroutine THERMO_MECHANICS
  !
  !-----------------------------------------------------------------------------
  !
  Subroutine MATERIAL_PROPERTIES ()
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    Calculate the volume-averaged, cell-centered solid material properties
    !    that use the conic temperature dependence
    !
    !---------------------------------------------------------------------------
    use zone_module, only: zone
    
    call compute_cell_property ('Lame1', zone%temp, Lame1)
    call compute_cell_property ('Lame2', zone%temp, Lame2)

  End Subroutine MATERIAL_PROPERTIES

  !
  !-----------------------------------------------------------------------------
  !
  Subroutine MATERIAL_DISPLACEMENTS ()
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    Solve the thermomechanical equations of static equilibrium for the 
    !    vertex-centered solid material displacement field
    !
    !---------------------------------------------------------------------------
    !
    Use linear_solution, Only: Ubik_user, PRECOND_NONE, PRECOND_TM_SSOR, PRECOND_TM_DIAG
    Use parameter_module, Only: ndim, nnodes
    Use preconditioners, Only: PRECONDITION
    use UbikSolve_module
    use nonlinear_solution, only: Nonlinear_Solve, NKuser,                     &
                                  NK_GET_SOLUTION_FIELD, NK_INITIALIZE, NK_FINALIZE
    use string_utilities, only: i_to_c

    real(r8), Dimension(ndim*nnodes) :: Solution

    ! Local variables
    integer :: idim, inodes, status, precon_type
    ! Solution time and space information
    integer, parameter :: FUTURE            = 2
    integer, parameter :: PRESENT           = 1
    integer, parameter :: PAST              = 1

    !---------------------------------------------------------------------------

    ! Initialize the solution of the linear system of equations
    Do inodes = 1, nnodes
       Do idim = 1, ndim
          Solution((inodes - 1)*ndim + idim) = Displacement(idim,inodes)
       End Do
    End Do

    ! Get right hand side from thermal stresses
    Call RHS_THERMO_MECH()

    ! Solve the nonlinear system of equations for the vertex-centered displacement field

    ! Call preconditioner if required
    precon_type = Ubik_user(NKuser(NK_DISPLACEMENT)%linear_solver_index)%precond

    ! Check for allowable preconditioners for solid mechanics
    select case (precon_type)
    case (PRECOND_TM_SSOR, PRECOND_TM_DIAG)
       call MECH_PRECONDITION ()
    case (PRECOND_NONE)
       ! Do nothing
    case DEFAULT
       call TLS_fatal ('MATERIAL_DISPLACEMENTS: invalid displacement preconditioner: ' // i_to_c(precon_type))
    end select


    ! Initialize Newton Krylov data
    Call NK_INITIALIZE(Solution_field=VP_Solution_Field, VECTORSIZE = (ndim*nnodes), UNKNOWNS_PER_ELEMENT = 1, &
                          TIME_LEVELS = 2, SOLUTION_OLD = Solution, SOLUTION_CURRENT = Solution)

    ! Call the nonlinear solver.
    call Nonlinear_Solve (NK_DATA=VP_Solution_Field, RESIDUAL = ELAS_VP_RESIDUAL, MATVEC = VP_MATVEC,            &
                PRECONDITIONER = PRECONDITION, NLS = NKuser(NK_DISPLACEMENT), &
                LS = Ubik_user(NKuser(NK_DISPLACEMENT)%linear_solver_index), STATUS = status)

    ! Get the NK solution field and copy it into the data structure.
    call NK_GET_SOLUTION_FIELD (VP_Solution_Field, Solution, FUTURE)

    ! Destroy/deallocate the HT NK solution data structure.
    call NK_FINALIZE (VP_Solution_Field)

    ! Store the solution of the linear system of equations in the displacement vector
    Do inodes = 1, nnodes
       Do idim = 1, ndim
          Displacement(idim,inodes) = Solution((inodes - 1)*ndim + idim)
       End Do
    End Do
    ! Store the iteration counts
    thermo_elastic_iterations   = NKuser(NK_DISPLACEMENT)%linear_tot
    viscoplastic_iterations = NKuser(NK_DISPLACEMENT)%Newton_tot

  End Subroutine MATERIAL_DISPLACEMENTS
  !
  !-----------------------------------------------------------------------------
  !
  Subroutine RHS_THERMO_MECH ()
    !---------------------------------------------------------------------------
    !  Purpose:
    !   Compute the thermal stress terms (and volume change from phase change 
    !   terms) for the right hand side of the linear system or the constant part
    !   of the nonlinear residual.
    !---------------------------------------------------------------------------
    Use mesh_module,                 Only: Cell_Edge, Cell
    Use parameter_module,            Only: ndim, nnodes, ncells, nvc, ncomps, nvf
    Use zone_module,                 Only: Zone
    use node_operator_module,        only: CV_Internal, CV_Boundary, nipc, nbface, Nodal_Volume
    Use gs_module,                   Only: EN_SUM_Scatter
    use bc_data_types
    use mech_bc_data_module
    use viscoplasticity,             only: MATERIAL_STRESSES
    use solid_mech_constraints,      only: RHS_DISPLACEMENT_CONSTRAINTS
    use body_data_module,            only: Body_Force

    ! Local variables
    type (BC_Operator), POINTER :: Operator
    type (BC_Atlas),    POINTER :: Atlas
    real(r8), pointer, dimension(:,:) :: BC_Value_List

    integer :: idim, inode, icomps, icell, ibcop, status, ibface, atlas_size, ip
    real(r8), Dimension(nvc,ncells)    :: RHS_Temp
    real(r8), Dimension(ncells)        :: Dstrain
    real(r8), Dimension(nnodes)        :: SM_Body_Force
    real(r8), Dimension(ncomps,ncells) :: Thermal_Stress
    integer, parameter :: NON_MECH_BCS = 8
    real(r8), pointer, dimension(:) :: Rho_Node, Nvol, Rho_tmp

    ! Calculate the change in the cell-centered thermal strain field

    Thermal_Strain_Inc = 0.0_r8
    Thermal_Stress     = 0.0_r8
    Dstrain            = 0.0_r8

    call compute_cell_property ('TM density', zone%temp, tm_dens)
    where (gap_cell_mask) tm_dens = 1.0
    
    ! NNC, May 2012.  What's wanted here is (rho_old/rho)**(1/3) - 1.
    ! The Dstrain computed below is asymptotically correct for rho_old/rho
    ! nearly 1.  Does it give better accuracy (less cancellation error)
    ! than the straightforward calculation?
    
    ! Isotropic dilatation from both thermal expansion and phase change
    ! log strain - we may want this for finite strain accuracy...
    Dstrain(:) = log(tm_dens_old(:)/tm_dens(:))/3.0
    ! but for now it may be more consistent to use (delta l)/l... No, accumulates errors
    ! Dstrain(:) = exp(Dstrain(:)) - 1.0
    where (.not. Solid_Mask)
       Dstrain(:) = 0.0
    end where

    do icomps = 1,ndim
       Thermal_Strain_Inc(icomps,:) = Dstrain(:)
    end do

    ! Calculate the change in the cell-centered thermal stress field
    Thermal_Strain(:,:) = Thermal_Strain(:,:) + Thermal_Strain_Inc(:,:)
    Call MATERIAL_STRESSES (Thermal_Strain, Thermal_Stress)

    ! Calculate right hand side thermal stress terms for control volume faces internal to cells.
    ! For the first vertex of each edge: RHS = RHS + Stress_ii * normal_i * area
    ! For the second vertex of each edge: RHS = RHS + Stress_ii * (-normal_i) * area

    NDIM_LOOP: do idim = 1,ndim
       RHS_Temp(:,:) = 0.0_r8
       NIPC_LOOP: do ip=1,nipc
          do icell=1,ncells
             RHS_Temp(Cell_Edge(1,ip),icell) = RHS_Temp(Cell_Edge(1,ip),icell) + &
                   Thermal_Stress(idim,icell) * &
                   CV_Internal%Face_Normal(idim,ip,icell) * CV_Internal%Face_Area(ip,icell)
             RHS_Temp(Cell_Edge(2,ip),icell) = RHS_Temp(Cell_Edge(2,ip),icell) - &
                   Thermal_Stress(idim,icell) * &
                   CV_Internal%Face_Normal(idim,ip,icell) * CV_Internal%Face_Area(ip,icell)
          end do
       end do NIPC_LOOP

       ! For traction BCs substitute the specified boundary traction for all other
       ! contributions for the CV face on the boundary.  This can go on the RHS 
       ! since it does not depend on the displacements.
       !
       ! Do tractions first and sum-scatter along with the internal CV faces to get 
       ! contributions of all faces - take care of displacement BCs later.
       !
       ! This must be done every time step
       select case(idim)
       case(1)
          ibcop = BC_X_TRACTION_OP
       case(2)
          ibcop = BC_Y_TRACTION_OP
       case(3)
          ibcop = BC_Z_TRACTION_OP
       end select

       Operator => BC_Spec_Get_Operator(Displacement_BC,ibcop)
       Atlas => BC_Op_Get_Atlas(Operator)
       atlas_size = DATA_SIZE(Atlas)
       if (atlas_size > 0) then
          BC_Value_List => BC_Get_Values(Atlas)
          !
          ! idim = 1,2,3 corresponds to x,y,z tractions
          do ibface = 1, nbface(ibcop)
             do inode = 1, nvf
                RHS_Temp(CV_Boundary(ibcop)%Face_Node(inode,ibface),CV_Boundary(ibcop)%Cell(ibface)) = &
                  RHS_Temp(CV_Boundary(ibcop)%Face_Node(inode,ibface),CV_Boundary(ibcop)%Cell(ibface)) - &
                  BC_Value_List(1,ibface) * &
                  CV_Boundary(ibcop)%Face_Area(inode,ibface)
             end do
          end do
       end if
       ! Tractions specified normal to the boundary
       Operator => BC_Spec_Get_Operator(Displacement_BC,BC_NORM_TRACTION_OP)
       Atlas => BC_Op_Get_Atlas(Operator)
       atlas_size = DATA_SIZE(Atlas)
       if (atlas_size > 0) then
          BC_Value_List => BC_Get_Values(Atlas)
          ibcop =  BC_NORM_TRACTION_OP
          do ibface = 1, nbface(ibcop)
             do inode = 1, nvf
                RHS_Temp(CV_Boundary(ibcop)%Face_Node(inode,ibface),CV_Boundary(ibcop)%Cell(ibface)) = &
                  RHS_Temp(CV_Boundary(ibcop)%Face_Node(inode,ibface),CV_Boundary(ibcop)%Cell(ibface)) - &
                  BC_Value_List(1,ibface) * CV_Boundary(ibcop)%Face_Normal(idim,inode,ibface) * &
                  CV_Boundary(ibcop)%Face_Area(inode,ibface)
             end do
          end do
       end if

       !
       ! For now, enforce zero displacements for materials with zero stiffness(?)
       ! This will probably cause problems when we have materials in the cell solidifying, etc.
       Where(Lame2(:) == 0)
          RHS_Temp(1,:) = 0.0_r8
          RHS_Temp(2,:) = 0.0_r8
          RHS_Temp(3,:) = 0.0_r8
          RHS_Temp(4,:) = 0.0_r8
          RHS_Temp(5,:) = 0.0_r8
          RHS_Temp(6,:) = 0.0_r8
          RHS_Temp(7,:) = 0.0_r8
          RHS_Temp(8,:) = 0.0_r8
       End Where

       call EN_SUM_Scatter(RHS(idim:(ndim*(nnodes-1)+idim):ndim),RHS_Temp)

       NULLIFY (Operator, Atlas, BC_Value_List)
    end do NDIM_LOOP

    ! Add body forces - uses Body_Force in the physics namelist
    if (solid_mechanics_body_force) then
    ! Get nodal averaged cell density for body force
    ! Use straight volume averaging to avoid gap element problems 
       allocate(Rho_Node(nnodes), stat=status)
       if (status /= 0) call TLS_panic ( 'RHS_THERMO_MECH: allocation error: Rho_Node')
       allocate(Nvol(nnodes), stat=status)
       if (status /= 0) call TLS_panic ( 'RHS_THERMO_MECH: allocation error: Nvol')
       allocate(Rho_tmp(ncells), stat=status)
       if (status /= 0) call TLS_panic ( 'RHS_THERMO_MECH: allocation error: Rho_tmp')
       
       call EN_SUM_SCATTER (Nvol, Cell%Volume)
       Rho_tmp = Zone(:)%Rho*Cell%Volume
       call EN_SUM_SCATTER (Rho_Node, Rho_tmp)
       Rho_Node = Rho_Node/Nvol

       do idim = 1, ndim
          SM_Body_Force = Rho_Node * Nodal_Volume * Body_Force(idim)
          RHS(idim:ndim*(nnodes - 1) + idim:ndim) = RHS(idim:ndim*(nnodes-1) + idim:ndim) - &
               SM_Body_Force
       end do

       deallocate(Rho_Node, Nvol, Rho_tmp)
    end if

    RHS = RHS/cscale    

    ! Save the current RHS vector which contains only source terms and traction BCs
    ! Src will be needed by the contact BC routines.  Change the sign for consistency
    ! with the equations f_j + r_j = 0, where r_j = source terms.
    Src = -RHS

    CALL RHS_DISPLACEMENT_CONSTRAINTS()

  END Subroutine RHS_THERMO_MECH
!

  !
  SUBROUTINE MECH_PRECONDITION()
    !=============================================================================
    !
    ! Construct preconditioning matrix, P, for thermo-elastic preconditioning 
    ! operator.  The matrix M, of second derivatives of displacements (strain
    ! gradients), which depends on mesh geometry only, is calculated and 
    ! stored.  This is the displacement derivatives part of the preconditioning
    ! operator, consisting of u_i,jj and u_j,ij in the linear elastic constitutive 
    ! equation:
    !   
    ! Lame1 * u_i,jj + (Lame1 + Lame2) * u_j,ij = RHS
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================
    use preconditioners,        only: TM_P, TM_P_Map
    use parameter_module,       only: ndim, nnodes, ncells
    use node_operator_module,   only: mech_precond_init
    use mesh_module,            only: Vertex_Ngbr_All, Vertex_Ngbr_All_Orig, Cell, Mesh, &
                                      GAP_ELEMENT_1
    use mech_bc_data_module
    use timing_tree
    Use gs_module,              Only: EN_SUM_Scatter
    use solid_mech_constraints, only: MECH_PRECOND_DISP_CONSTRAINTS

    ! Preconditioning matrix.  TM_P will be pointed at this.
    type(real_var_vector), pointer, save, dimension(:) :: A_Elas
    ! Connectivity array for A_Elas.  TM_P_Map will be pointed at this.
    type(int_var_vector), pointer, save, dimension(:) :: A_Conn
    ! Node centered elastic constants
    real(r8), pointer, dimension(:) :: Lame1_Node, Lame2_Node , Nvol, L1tmp, L2tmp
    integer :: idim, jdim, lnode,  inode, jnode, nmax
    integer :: status
    integer, dimension(nnodes) :: NN_Sizes
    integer, dimension(nnodes*ndim) :: C_Sizes
    ! Pointers to var_vector data
    real(r8), pointer, dimension(:) :: A_Vec
    real(r8), pointer, dimension(:) :: M1_Vec
    real(r8), pointer, dimension(:) :: M2_Vec
    integer, pointer, dimension(:) :: A_C_Vec
    integer, pointer, dimension(:) :: NN_Vec

    if (.not. mech_precond_init) then
       ! This if block is really part of initialization and is only done once
       
       call stop_timer ("Solid Mechanics")
       call start_timer("Initialization")
       call start_timer("Solid Mech Init")

       call PRECON_DISPLACEMENT_GRADIENTS()

       ! A_Elas and A_Conn are created once and saved
       allocate(A_Elas(ndim*nnodes), stat=status)
       if (status /= 0) call TLS_panic ( 'MECH_PRECONDITION: allocation error: A_Elas')
       allocate(A_Conn(ndim*nnodes), stat=status)
       if (status /= 0) call TLS_panic ( 'MECH_PRECONDITION: allocation error: A_Conn')
       NN_Sizes = SIZES(Vertex_Ngbr_All_Orig)

       ! The size of the connectivity array for a given node is ndim * (number of neighbors + 1)
       do idim = 1,ndim
          C_Sizes(idim:(ndim*(nnodes-1)+idim):ndim) = (NN_Sizes(:) + 1) * ndim
       end do
       call CREATE(A_Elas,C_Sizes)
       call CREATE(A_Conn,C_Sizes)

       ! Connectivity array does not change
       A_C_Vec => FLATTEN(A_Conn)
       A_C_Vec = 0
       ! Construct connectivity array
       ! Elements 4:(nmax*3) are neighbor node displacement components in x,y,z order
       do inode = 1,nnodes
          nmax = SIZES(Vertex_Ngbr_All(inode))
          NN_Vec => FLATTEN(Vertex_Ngbr_All(inode))
          do idim = 1, ndim
             A_C_Vec => FLATTEN(A_Conn(ndim*(inode-1)+idim))
             do lnode = 1, nmax
                jnode = NN_Vec(lnode)
                do jdim = 1,ndim
                   if(jnode > 0) then
                      A_C_Vec(ndim*lnode+jdim) = ndim*(jnode-1)+jdim
                   else
                      ! This connectivity should match that in NN_Gather_BoundaryData in 
                      ! the preconditioner routines
                      A_C_Vec(ndim*lnode+jdim) = -(ndim*(-jnode-1)+jdim)
                   end if
                end do
             end do
          end do
       end do

       ! Fill the first 3 elements of the connectivity array
       do inode = 1,nnodes    
          ! If the DOF is an x displacement component, element 2 is y, 3 is z
          A_C_Vec => FLATTEN(A_Conn(ndim*(inode-1)+1))
          A_C_Vec(1) = ndim * (inode - 1) + 1
          A_C_Vec(2) = ndim * (inode - 1) + 2
          A_C_Vec(3) = ndim * (inode - 1) + 3
          ! If the DOF is a y displacement component, element 2 is x, 3 is z
          A_C_Vec => FLATTEN(A_Conn(ndim*(inode-1)+2))
          A_C_Vec(1) = ndim * (inode - 1) + 2
          A_C_Vec(2) = ndim * (inode - 1) + 1
          A_C_Vec(3) = ndim * (inode - 1) + 3
          ! If the DOF is a z displacement component, element 2 is y, 3 is x
          A_C_Vec => FLATTEN(A_Conn(ndim*(inode-1)+3))
          A_C_Vec(1) = ndim * (inode - 1) + 3
          A_C_Vec(2) = ndim * (inode - 1) + 2
          A_C_Vec(3) = ndim * (inode - 1) + 1
       end do

       mech_precond_init = .true.
       !
       call stop_timer ("Solid Mech Init")
       call stop_timer ("Initialization")
       call start_timer("Solid Mechanics")
       
       !NNC: call TIMER_START(TIMER_CYCLE)
    end if

    call start_timer ("Precondition")
    
    ! Initialize
    A_Vec => FLATTEN(A_Elas)
    A_Vec = 0.0

    ! Get elastic constants at the nodes
    allocate(Lame1_Node(nnodes), stat=status)
    if (status /= 0) call TLS_panic ( 'MECH_PRECONDITION: allocation error: Lame1_Node')
    allocate(Lame2_Node(nnodes), stat=status)
    if (status /= 0) call TLS_panic ( 'MECH_PRECONDITION: allocation error: Lame2_Node')
    allocate(Nvol(nnodes), stat=status)
    if (status /= 0) call TLS_panic ( 'MECH_PRECONDITION: allocation error: Nvol')
    ! 
    allocate(L1tmp(ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'MECH_PRECONDITION: allocation error: L1tmp')
    allocate(L2tmp(ncells), stat=status)
    if (status /= 0) call TLS_panic ( 'MECH_PRECONDITION: allocation error: L2tmp')

    ! To avoid weighting problems with gap elements (~0 volume), use a straight
    ! volume weighted average for the elastic constants ("It's only a preconditioner.")
    call EN_SUM_SCATTER (Nvol, Cell%Volume)
    L1tmp = Lame1*Cell%Volume
    L2tmp = Lame2*Cell%Volume
    where (Mesh(:)%Cell_Shape >= GAP_ELEMENT_1)
       L1tmp = 0.0
       L2tmp = 0.0
    end where

    call EN_SUM_SCATTER (Lame1_Node, L1tmp)
    call EN_SUM_SCATTER (Lame2_Node, L2tmp)

    Lame1_Node = Lame1_Node/Nvol
    Lame2_Node = Lame2_Node/Nvol

    ! Fill preconditioning matrix by multiplying displacement gradient factors in M by 
    ! elastic constants Lame1 and Lame2
    do inode = 1,nnodes    
       nmax = SIZES(Vertex_Ngbr_All(inode))
       do idim = 1, ndim
          A_Vec => FLATTEN(A_Elas(ndim*(inode-1)+idim))
          M1_Vec => FLATTEN(M1(ndim*(inode-1)+idim))
          M2_Vec => FLATTEN(M2(ndim*(inode-1)+idim))
          if (Lame1_Node(inode) > 1.0e-6) then
             do lnode = 1, (nmax + 1)
                do jdim = 1,ndim
                   A_Vec(ndim*(lnode-1)+jdim) = (Lame1_Node(inode) * M1_Vec(ndim*(lnode-1)+jdim) + &
                        Lame2_Node(inode) * M2_Vec(ndim*(lnode-1)+jdim)) / cscale(ndim*(inode-1)+idim)
                end do
             end do
          else
          ! Set displacements to zero for empty or fluid filled cells
             A_Vec(1) = 1.0
             do lnode = 2,ndim
                A_Vec(lnode) = 0.0
             end do
             do lnode = 2, (nmax + 1)
                do jdim = 1,ndim
                   A_Vec(ndim*(lnode-1)+jdim) = 0.0
                end do
             end do
          end if
       end do
    end do

    ! Enforce displacement constraints

    call MECH_PRECOND_DISP_CONSTRAINTS(A_Elas)

    ! Nullify var_vector pointers
    NULLIFY(A_C_Vec, A_Vec, NN_Vec, M1_Vec, M2_Vec)
    ! Deallocate material property arrays
    deallocate(Lame1_Node)
    deallocate(Lame2_Node)
    deallocate(Nvol)
    deallocate(L1tmp)
    deallocate(L2tmp)

    ! Associate preconditioner matrix and connectivity

    TM_P => A_Elas
    TM_P_Map => A_Conn
   
    call stop_timer("Precondition")
    
  END SUBROUTINE MECH_PRECONDITION

  !-----------------------------------------------------------------------------
  !
  SUBROUTINE PRECON_DISPLACEMENT_GRADIENTS()
    !=============================================================================
    !
    !  Construct strain gradient matrix, M, for thermo-elastic preconditioning 
    !  operator
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================
    use parameter_module, only: ndim, nvc, ncells_tot, nnodes, nnodes_tot
    use mesh_module, only: Mesh, Vertex, CELL_TET,&
         CELL_PYRAMID, CELL_PRISM, CELL_HEX, MESH_COLLATE_VERTEX, VERTEX_COLLATE, &
          Vertex_Ngbr_All_Orig, VERTEX_DATA, GAP_ELEMENT_1, GAP_ELEMENT_3
    use node_operator_module, only: CV_Internal, nipc
    use parallel_info_module
    use pgslib_module, only: PGSLib_DIST, PGSLib_COLLATE, PGSLib_Global_Sum

    integer, parameter :: nnt = ndim + 1

    ! Local variables
    real(r8), dimension(ndim) :: Tcen, Fsign 
    integer, dimension(ndim) :: Cvf
    integer, dimension(nnt) :: n
    integer, dimension(nnt-1) :: Tnode
    integer :: inode, idim, lnode, neq, i, j, mindex_i, mindex_j, icell, ivrtx, ip, flat_size
    integer :: status
    ! Interpolation function coefficients for tets
    real(r8), Dimension(ndim,nnt) :: A
    ! Coordinates of tet vertices
    real(r8), Dimension(nnt,ndim) :: Tc
    ! Used for reordering M
    real(r8) :: M_Temp
    ! Local node-node connectivity data 
    integer, dimension(nnodes) :: NN_Sizes
    integer, pointer, dimension(:) :: Node_Ngbr
    ! Collated mesh connectivity data, vertex mapping only (nvc,ncells_tot)
    integer, pointer, dimension(:,:) :: Mesh_Tot_Vertex => NULL()
    ! Collated mesh cell shape data (ncells_tot)
    integer, pointer, dimension(:) :: Cell_Shape_Tot
    ! Collated vertex data (nnodes_tot)
    type(VERTEX_DATA), pointer, dimension(:) :: Vertex_Tot => NULL()
    ! M var vector sizes
    integer, pointer, dimension(:) :: M_Sizes_Tot
    integer, pointer, dimension(:) :: M_Sizes
    ! Collated NN sizes data
    integer, pointer, dimension(:) :: NN_Sizes_Tot
    ! Collated NN vector data
    integer, pointer, dimension(:) :: Node_Ngbr_Tot
    ! Collated NN var vector (nnodes_tot)
    type(int_var_vector), pointer, dimension(:) :: Vertex_Ngbr_Tot
    ! M1 and M2 var vectors for collated mesh
    type(real_var_vector), pointer, dimension(:) :: M1_Tot
    type(real_var_vector), pointer, dimension(:) :: M2_Tot

    integer, pointer, dimension(:) :: NN_Vec
    real(r8), pointer, dimension(:) :: M1_Tot_Vec
    real(r8), pointer, dimension(:) :: M2_Tot_Vec
    real(r8), pointer, dimension(:) :: M_Data_Tot
    real(r8), pointer, dimension(:) :: M_Data
    real(r8), allocatable, dimension(:,:) :: CV_Area_Tot
    real(r8), allocatable, dimension(:,:,:) :: CV_Normal_Tot

    !Allocate temporary arrays

    ! Collate vertex cell connectivity, vertex coordinates and cell shape

    Mesh_Tot_Vertex => MESH_COLLATE_VERTEX (Mesh)
    Vertex_Tot => VERTEX_COLLATE(Vertex)
    if (p_info%IOP) then
       allocate(Cell_Shape_Tot(ncells_tot), stat = status)
    else
       allocate(Cell_Shape_Tot(0), stat = status)
    end if
    if (status /= 0) call TLS_panic ( 'PRECON_DISPLACEMENT_GRADIENT: allocation error: Cell_Shape_Tot')
    call PGSLib_COLLATE (Cell_Shape_Tot,Mesh%Cell_Shape)

    ! Collate control volume face areas and normals
    if (p_info%IOP) then
       allocate(CV_Area_Tot(nipc,ncells_tot), stat = status)
    else
       allocate(CV_Area_Tot(nipc,0), stat = status)
    end if
    if (status /= 0) call TLS_panic ( 'PRECON_DISPLACEMENT_GRADIENT: allocation error: CV_Area_Tot')
    if (p_info%IOP) then
       allocate(CV_Normal_Tot(ndim,nipc,ncells_tot), stat = status)
    else
       allocate(CV_Normal_Tot(ndim,nipc,0), stat = status)
    end if
    if (status /= 0) call TLS_panic ( 'PRECON_DISPLACEMENT_GRADIENT: allocation error: CV_Normal_Tot')

    do ip = 1,nipc
       call PGSLib_COLLATE (CV_Area_Tot(ip,:),CV_Internal%Face_Area(ip,:))
       do idim=1,ndim
          call PGSLib_COLLATE (CV_Normal_Tot(idim,ip,:),CV_Internal%Face_Normal(idim,ip,:))
       end do
    end do
    ! Collate the node neighbors

    ! Get collated sizes 
    NN_Sizes = SIZES(Vertex_Ngbr_All_Orig)
    if (p_info%IOP) then
       allocate(NN_Sizes_Tot(nnodes_tot), stat = status)
    else
       allocate(NN_Sizes_Tot(0), stat = status)
    end if
    if (status /= 0) call TLS_panic ( 'PRECON_DISPLACEMENT_GRADIENT: allocation error: NN_Sizes_Tot')
    call PGSLib_COLLATE (NN_Sizes_Tot, NN_Sizes)

    ! Get collated node neighbors

    Node_Ngbr => FLATTEN(Vertex_Ngbr_All_Orig)
    flat_size = PGSLib_Global_Sum(SIZE(Node_Ngbr))
    if (p_info%IOP) then
       allocate(Node_Ngbr_Tot(flat_size), stat = status)
    else
       allocate(Node_Ngbr_Tot(0), stat = status)
    end if
    if (status /= 0) call TLS_panic ( 'PRECON_DISPLACEMENT_GRADIENT: allocation error: Node_Ngbr_Tot')
    call PGSLib_COLLATE (Node_Ngbr_Tot, Node_Ngbr) 
    if (p_info%IOP) then
       allocate(Vertex_Ngbr_Tot(nnodes_tot), stat = status)
    else
       allocate(Vertex_Ngbr_Tot(0), stat = status)
    end if
    if (status /= 0) call TLS_panic ( 'PRECON_DISPLACEMENT_GRADIENT: allocation error: Vertex_Ngbr_Tot')
    call CREATE(Vertex_Ngbr_Tot,NN_Sizes_Tot)
    call STUFF(Vertex_Ngbr_Tot,Node_Ngbr_Tot)

    ! Create var vector to hold static geometry data for constructing the preconditioning
    ! matrix (also a var_vector)
    if (p_info%IOP) then
       allocate(M_Sizes_Tot(nnodes_tot * ndim), stat=status)
    else
       allocate(M_Sizes_Tot(0), stat=status)
    end if
    if (status /= 0) call TLS_panic ( 'PRECON_DISPLACEMENT_GRADIENT: allocation error: M_Sizes_Tot')

    ! Each node has ndim displacement components

    if (p_info%IOP) then
       do idim = 1,ndim
          M_Sizes_Tot(idim:(ndim*(nnodes_tot-1)+idim):ndim) = (NN_Sizes_Tot(:) + 1) * ndim
       end do
    end if

    if (p_info%IOP) then
       allocate(M1_Tot(nnodes_tot * ndim), stat = status)
    else
       allocate(M1_Tot(0), stat = status)
    end if
    if (status /= 0) call TLS_panic ( 'PRECON_DISPLACEMENT_GRADIENT: allocation error: M1_Tot')
    if (p_info%IOP) then
       allocate(M2_Tot(nnodes_tot * ndim), stat = status)
    else
       allocate(M2_Tot(0), stat = status)
    end if
    if (status /= 0) call TLS_panic ( 'PRECON_DISPLACEMENT_GRADIENT: allocation error: M2_Tot')
    call CREATE(M1_Tot,M_Sizes_Tot)
    call CREATE(M2_Tot,M_Sizes_Tot)

    M1_Tot_Vec => FLATTEN(M1_Tot)
    M1_Tot_Vec = 0.0
    M2_Tot_Vec => FLATTEN(M2_Tot)
    M2_Tot_Vec = 0.0

    ! Loop over cells, getting matrix coefficients for each node

    CELL_LOOP: do icell = 1, SIZE(Mesh_Tot_Vertex,2)
       ! For each vertex, generate gradient factors A(), vector to the centroid Dx,
       ! squared distance to the centroid D2
       if (Cell_Shape_Tot(icell) >= GAP_ELEMENT_1) CYCLE CELL_LOOP
       VERTEX_LOOP: do ivrtx = 1,nvc
          A(:,:) = 0.0
          TNode(:) = 0
          inode = Mesh_Tot_Vertex(ivrtx,icell)
          ! Get vertex numbers, control volume faces and vector polarity 
          ! for this tet.  Accounting for degenerate cells makes this 
          ! kind of messy.  The control volume faces and normals do not
          ! change, but the node numbers do.
          select case (ivrtx)
          case(1)  
             n(2) = 2; n(3) = 4; n(4) = 5
             if (Cell_Shape_Tot(icell) == CELL_TET) n(2) = 3
             Cvf(1) = 1; Cvf(2) = 4; Cvf(3) = 8
             Fsign(1) =  1.0; Fsign(2) = -1.0;Fsign(3) =  1.0
          case(2)
             n(2) = 3; n(3) = 1; n(4) = 6
             if (Cell_Shape_Tot(icell) == CELL_TET) n(3) = 4
             Cvf(1) = 2; Cvf(2) = 1; Cvf(3) = 5
             Fsign(1) =  1.0; Fsign(2) = -1.0;Fsign(3) =  1.0
          case(3)
             n(2) = 4; n(3) = 2; n(4) = 7
             Cvf(1) = 3; Cvf(2) = 2; Cvf(3) = 6
             Fsign(1) =  1.0; Fsign(2) = -1.0;Fsign(3) =  1.0
          case(4)
             n(2) = 1; n(3) = 3; n(4) = 8
             Cvf(1) = 4; Cvf(2) = 3; Cvf(3) = 7
             Fsign(1) =  1.0; Fsign(2) = -1.0;Fsign(3) =  1.0
          case(5)
             n(2) = 6; n(3) = 1; n(4) = 8
             if (Cell_Shape_Tot(icell) /= CELL_HEX) then
                select case (Cell_Shape_Tot(icell))
                case(CELL_PRISM)
                   n(4) = 4
                case(CELL_PYRAMID)
                   n(2) = 2; n(4) = 4
                case(CELL_TET)
                   n(2) = 2; n(3) = 4; n(4) = 3
                ! Gap elements with prism type connectivity (maybe we should not use cell shape
                ! for gap element identification?)
                case(GAP_ELEMENT_3)
                   if (Mesh_Tot_Vertex(5,icell) == Mesh_Tot_Vertex(8,icell)) n(4) = 4
                end select
             end if
             Cvf(1) = 9; Cvf(2) = 8; Cvf(3) = 12
             Fsign(1) =  1.0; Fsign(2) = -1.0;Fsign(3) = -1.0
          case(6)
             n(2) = 7; n(3) = 2; n(4) = 5
             if (Cell_Shape_Tot(icell) /= CELL_HEX) then
                select case (Cell_Shape_Tot(icell))
                case(CELL_PRISM)
                   n(2) = 3
                case(CELL_PYRAMID)
                   n(2) = 3; n(4) = 1
                case(CELL_TET)
                   n(2) = 3; n(4) = 4
                case(GAP_ELEMENT_3)
                   if (Mesh_Tot_Vertex(5,icell) == Mesh_Tot_Vertex(8,icell)) n(2) = 3
                end select
             end if
             Cvf(1) = 10; Cvf(2) = 5; Cvf(3) = 9
             Fsign(1) =  1.0; Fsign(2) = -1.0;Fsign(3) = -1.0
          case(7)
             n(2) = 8; n(3) = 3; n(4) = 6
             if (Cell_Shape_Tot(icell) /= CELL_HEX) then
                select case (Cell_Shape_Tot(icell))
                case(CELL_PRISM)
                   n(4) = 2
                case(CELL_PYRAMID)
                   n(2) = 4; n(4) = 2
                case(CELL_TET)
                   n(2) = 4; n(4) = 2
                case(GAP_ELEMENT_3)
                   if (Mesh_Tot_Vertex(5,icell) == Mesh_Tot_Vertex(8,icell)) n(4) = 2
                end select
             end if
             Cvf(1) = 11; Cvf(2) = 6; Cvf(3) = 10
             Fsign(1) =  1.0; Fsign(2) = -1.0;Fsign(3) = -1.0
          case(8)
             n(2) = 5; n(3) = 4; n(4) = 7
             if (Cell_Shape_Tot(icell) /= CELL_HEX) then
                select case (Cell_Shape_Tot(icell))
                case(CELL_PRISM)
                   n(2) = 1
                case(CELL_PYRAMID)
                   n(2) = 1; n(4) = 3
                case(CELL_TET)
                   n(2) = 1; n(4) = 3
                case(GAP_ELEMENT_3)
                   if (Mesh_Tot_Vertex(5,icell) == Mesh_Tot_Vertex(8,icell)) n(2) = 1
                end select
             end if
             Cvf(1) = 12; Cvf(2) = 7; Cvf(3) = 11
             Fsign(1) =  1.0; Fsign(2) = -1.0;Fsign(3) = -1.0
          end select
          ! Connectivity indices for this tet
          TET_NODE_LOOP: do lnode = 2, nnt
             do i = 1,SIZES(Vertex_Ngbr_Tot(inode))
                NN_Vec => FLATTEN(Vertex_Ngbr_Tot(inode))
                if (NN_Vec(i) == Mesh_Tot_Vertex(n(lnode),icell)) then
                   Tnode(lnode-1) = i + 1
                   exit
                end if
             end do
          end do TET_NODE_LOOP
          if (any(Tnode(:) == 0)) then
             call TLS_fatal ('PRECON_DISPLACEMENT_GRADIENT: error finding connectivity for mech preconditioner')
          end if
          ! Vertex coordinates for this tet only
          Tc(1,:) = Vertex_Tot(Mesh_Tot_Vertex(ivrtx,icell))%Coord(:)
          Tc(2,:) = Vertex_Tot(Mesh_Tot_Vertex(n(2),icell))%Coord(:)
          Tc(3,:) = Vertex_Tot(Mesh_Tot_Vertex(n(3),icell))%Coord(:)
          Tc(4,:) = Vertex_Tot(Mesh_Tot_Vertex(n(4),icell))%Coord(:)
          ! Get gradient factors for this tet
          ! Set A to zero for gap elements
          Call TET_GRADIENT_FACTORS(Tc, A(:,:), Tcen(:))
          ! n(1:4) are indicies for the connectivity matrix relating global nodes to
          ! nodes 1-4 for a single tet.
          n(1) = 1
          n(2) = Tnode(1)
          n(3) = Tnode(2)
          n(4) = Tnode(3)
          IP_LOOP: do ip = 1,ndim
             FORCE_COMP_LOOP: do i = 1,ndim
                ! neq corresponds to dof i in sigma_ij
                neq = (inode - 1) * ndim + i
                M1_Tot_Vec => FLATTEN(M1_Tot(neq))
                M2_Tot_Vec => FLATTEN(M2_Tot(neq))
                ! Loop over all nodes of the tet
                do lnode = 1, nnt
                   ! The array M is in order u1,v1,w1,u2,v2,w2... , starting at 1
                   mindex_i = ndim * (n(lnode) - 1) + i
                   ! Summing over j, where T_i = Lame_1 u_jj n_i + Lame_2 * (u_i,j + u_j,i) n_j
                   do j = 1,ndim
                      mindex_j = ndim * (n(lnode) - 1) + j
                      M1_Tot_Vec(mindex_j) = M1_Tot_Vec(mindex_j) + A(j,lnode) * CV_Normal_Tot(i,Cvf(ip),icell) * &
                           Fsign(ip) * CV_Area_Tot(Cvf(ip),icell)
                      M2_Tot_Vec(mindex_i) = M2_Tot_Vec(mindex_i) + A(j,lnode) * CV_Normal_Tot(j,Cvf(ip),icell) * &
                           Fsign(ip) * CV_Area_Tot(Cvf(ip),icell)
                      M2_Tot_Vec(mindex_j) = M2_Tot_Vec(mindex_j) + A(i,lnode) * CV_Normal_Tot(j,Cvf(ip),icell) * &
                           Fsign(ip) * CV_Area_Tot(Cvf(ip),icell)
                   end do
                end do
             end do FORCE_COMP_LOOP
          end do IP_LOOP
       end do VERTEX_LOOP
    end do CELL_LOOP

    ! Switch terms of M1 and M2 to put the diagonal term in the first vector location for v and w equations
    do i = 1,(SIZE(M1_Tot)/ndim)
       M1_Tot_Vec => FLATTEN(M1_Tot(ndim*(i-1)+2))
       M2_Tot_Vec => FLATTEN(M2_Tot(ndim*(i-1)+2))
       M_Temp = M1_Tot_Vec(1)
       M1_Tot_Vec(1) = M1_Tot_Vec(2)
       M1_Tot_Vec(2) = M_Temp
       M_Temp = M2_Tot_Vec(1)
       M2_Tot_Vec(1) = M2_Tot_Vec(2)
       M2_Tot_Vec(2) = M_Temp
       if (ndim == 3) then
          M1_Tot_Vec => FLATTEN(M1_Tot(ndim*i))
          M2_Tot_Vec => FLATTEN(M2_Tot(ndim*i))
          M_Temp = M1_Tot_Vec(1)
          M1_Tot_Vec(1) = M1_Tot_Vec(3)
          M1_Tot_Vec(3) = M_Temp
          M_Temp = M2_Tot_Vec(1)
          M2_Tot_Vec(1) = M2_Tot_Vec(3)
          M2_Tot_Vec(3) = M_Temp
       end if
    end do

    ! Now create local M1 and M2 and distribute the data to them

    allocate(M_Sizes(nnodes*ndim), stat=status)
    if (status /= 0) call TLS_panic ( 'PRECON_DISPLACEMENT_GRADIENT: allocation error: M_Sizes')
    allocate(M1(nnodes*ndim), stat=status)
    if (status /= 0) call TLS_panic ( 'PRECON_DISPLACEMENT_GRADIENT: allocation error: M1')
    allocate(M2(nnodes*ndim), stat=status)
    if (status /= 0) call TLS_panic ( 'PRECON_DISPLACEMENT_GRADIENT: allocation error: M2')

    do idim = 1,ndim
       M_Sizes(idim:(ndim*(nnodes-1)+idim):ndim) = (NN_Sizes(:) + 1) * ndim
    end do
    allocate(M_Data(SUM(M_Sizes)), stat=status)
    if (status /= 0) call TLS_panic ( 'PRECON_DISPLACEMENT_GRADIENT: allocation error: M_Data')
    call CREATE(M1, M_Sizes)
    call CREATE(M2, M_Sizes)
    M_Data_Tot => Flatten(M1_Tot)
    call PGSLib_Dist(M_Data,M_Data_Tot)
    call STUFF(M1,M_Data)
    M_Data_Tot => Flatten(M2_Tot)
    call PGSLib_Dist(M_Data,M_Data_Tot)
    call STUFF(M2,M_Data)

    ! Clean up...

    ! Deallocate temporary arrays
    deallocate(Cell_Shape_Tot)
    deallocate(CV_Area_Tot)
    deallocate(CV_Normal_Tot)
    deallocate(M_Sizes_Tot)
    deallocate(NN_Sizes_Tot)
    deallocate(Node_Ngbr_Tot)
    ! Destroy temporary var_vectors
    call DESTROY(Vertex_Ngbr_Tot)
    call DESTROY(M1_Tot)
    call DESTROY(M2_Tot)
    ! Deallocate temporary var_vectors
    deallocate(Vertex_Ngbr_Tot)
    deallocate(M1_Tot)
    deallocate(M2_Tot)
    deallocate(Mesh_Tot_Vertex)
    deallocate(Vertex_Tot)

    !
    nullify(Node_Ngbr)
    nullify(NN_Sizes_Tot)
    nullify(NN_Vec, M1_Tot_Vec, M2_Tot_Vec)

  END SUBROUTINE PRECON_DISPLACEMENT_GRADIENTS

  !
  SUBROUTINE TET_GRADIENT_FACTORS(Tc, A, Tcen)
    !=============================================================================
    !
    !  Calculate interpolation function coefficients A_jl for a linear tetrahedron,
    !  such that:
    !
    !           4
    !         =====  
    !         \\
    !  u_i,j = \\  A_jl u_il      where u_il are displacements at the nodes l
    !          //
    !         //
    !         =====
    !          l=1
    !
    !  adapted from The Finite Element Method, O.C. Zienkiewicz, Third Edition 
    !
    ! Author(s): Dave Korzekwa, LANL (dak@lanl.gov)
    !=============================================================================
    ! Tc are the coordinates of the nodes
    ! Tcen are the coordinates of the centroid of the tet

    use parameter_module, only: ndim

    real(r8), Dimension(ndim+1,ndim) :: Tc
    real(r8), Dimension(ndim,ndim+1) :: A
    real(r8), Dimension(ndim) :: Tcen
    real(r8) :: sixv, det, a1, a2, b1, b2, c1, c2
    integer :: idim, inode

    ! This mess is the determinant of
    !   
    !   | 1 x1 y1 z1 |
    !   | 1 x2 y2 z2 |
    !   | 1 x3 y3 z3 |
    !   | 1 x4 y4 z4 |
    !
    ! as expanded by Maxima:
    !
    ! - tc11 (tc43 yc32 + tc23 (tc42 - yc32) - tc22 (tc43 - tc33) - tc33 tc42)
    !
    ! + tc21 (tc43 yc32 - tc33 tc42) - tc13 (- tc41 yc32 - tc21 (tc42 - yc32)
    !
    ! + tc31 tc42 + tc22 (tc41 - tc31)) + tc23 (tc31 tc42 - tc41 yc32)
    !
    ! + tc12 (- tc21 (tc43 - tc33) + tc31 tc43 + tc23 (tc41 - tc31) - tc33 tc41)
    !
    ! - tc22 (tc31 tc43 - tc33 tc41)

    sixv = -Tc(1,1)*(Tc(4,3)*Tc(3,2)+Tc(2,3)*(Tc(4,2)-Tc(3,2))-Tc(2,2)*(Tc(4,3)-Tc(3,3))- &
         Tc(3,3)*Tc(4,2)) + Tc(2,1)*(Tc(4,3)*Tc(3,2)-Tc(3,3)*Tc(4,2)) - Tc(1,3)*(-Tc(4,1)*Tc(3,2)- &
         Tc(2,1)*(Tc(4,2)-Tc(3,2))+Tc(3,1)*Tc(4,2)+Tc(2,2)*(Tc(4,1)-Tc(3,1))) + &
         Tc(2,3)*(Tc(3,1)*Tc(4,2)-Tc(4,1)*Tc(3,2)) + Tc(1,2)*(-Tc(2,1)*(Tc(4,3)-Tc(3,3))+ &
         Tc(3,1)*Tc(4,3)+Tc(2,3)*(Tc(4,1)-Tc(3,1))-Tc(3,3)*Tc(4,1)) - Tc(2,2)*(Tc(3,1)*Tc(4,3)- &
         Tc(3,3)*Tc(4,1))

    DIM_LOOP: do idim = 1,ndim
       NODE_LOOP: do inode = 1,ndim+1
          select case (idim)
          case (1)
             select case (inode)
             case (1)
                a1=Tc(2,2); a2=Tc(2,3); b1=Tc(3,2); b2=Tc(3,3); c1=Tc(4,2); c2=Tc(4,3)
             case (2)
                a1=Tc(3,2); a2=Tc(3,3); b1=Tc(4,2); b2=Tc(4,3); c1=Tc(1,2); c2=Tc(1,3)
             case (3)
                a1=Tc(4,2); a2=Tc(4,3); b1=Tc(1,2); b2=Tc(1,3); c1=Tc(2,2); c2=Tc(2,3)
             case (4)
                a1=Tc(1,2); a2=Tc(1,3); b1=Tc(2,2); b2=Tc(2,3); c1=Tc(3,2); c2=Tc(3,3)
             end select
             ! - aa1 (cc2 - bb2) + bb1 cc2 + aa2 (cc1 - bb1) - bb2 cc1
             det = -a1 * (c2 - b2) + b1 * c2 + a2 * (c1 - b1) - b2 * c1
          case(2)
             select case (inode)
             case (1)
                a1=Tc(2,1); a2=Tc(2,3); b1=Tc(3,1); b2=Tc(3,3); c1=Tc(4,1); c2=Tc(4,3)
             case (2)
                a1=Tc(3,1); a2=Tc(3,3); b1=Tc(4,1); b2=Tc(4,3); c1=Tc(1,1); c2=Tc(1,3)
             case (3)
                a1=Tc(4,1); a2=Tc(4,3); b1=Tc(1,1); b2=Tc(1,3); c1=Tc(2,1); c2=Tc(2,3)
             case (4)
                a1=Tc(1,1); a2=Tc(1,3); b1=Tc(2,1); b2=Tc(2,3); c1=Tc(3,1); c2=Tc(3,3)
             end select
             !aa1 (cc2 - bb2) - bb1 cc2 + bb2 cc1 + aa2 (bb1 - cc1)
             det = a1 * (c2 - b2) - b1 * c2 + b2 * c1 + a2 * (b1 - c1)
          case(3)
             select case (inode)
             case (1)
                a1=Tc(2,1); a2=Tc(2,2); b1=Tc(3,1); b2=Tc(3,2); c1=Tc(4,1); c2=Tc(4,2)
             case (2)
                a1=Tc(3,1); a2=Tc(3,2); b1=Tc(4,1); b2=Tc(4,2); c1=Tc(1,1); c2=Tc(1,2)
             case (3)
                a1=Tc(4,1); a2=Tc(4,2); b1=Tc(1,1); b2=Tc(1,2); c1=Tc(2,1); c2=Tc(2,2)
             case (4)
                a1=Tc(1,1); a2=Tc(1,2); b1=Tc(2,1); b2=Tc(2,2); c1=Tc(3,1); c2=Tc(3,2)
             end select
             !bb1 cc2 + aa1 (bb2 - cc2) - bb2 cc1 - aa2 (bb1 - cc1)
             det = b1 * c2 + a1 * (b2 - c2) - b2 * c1 - a2 * (b1 - c1)
          end select
          A(idim,inode) = (-1)**inode * det / sixv
       end do NODE_LOOP
       Tcen(idim) = SUM(Tc(:,idim))/4.0
    end do DIM_LOOP
  END SUBROUTINE TET_GRADIENT_FACTORS

  !
  FUNCTION STRESS_STRAIN_INVARIANTS () RESULT(STRESS_STRAIN)
    !=============================================================================
    !
    ! Calculate invariants of stress and strain tensors at cell centers, which are
    ! returned in Stress_Strain.  Right now this should only be called for long 
    ! edit output.
    ! 
    !=============================================================================
    use parameter_module, only: ncells

    integer :: status

    type(CELL_MECH_INVARIANT), pointer, dimension(:) :: Stress_Strain

    allocate(stress_strain(ncells),stat = status)
    if (status /= 0) call TLS_panic ( 'STRESS_STRAIN_INVARIANTS: allocation error: Stress_Strain')

    ! Von Mises stress
    Stress_Strain(:)%mises_stress = sqrt(((SMech_Cell%Elastic_Stress(1,:) - SMech_Cell%Elastic_Stress(2,:))**2 + &
         (SMech_Cell%Elastic_Stress(2,:) - SMech_Cell%Elastic_Stress(3,:))**2 + &
         (SMech_Cell%Elastic_Stress(3,:) - SMech_Cell%Elastic_Stress(1,:))**2 + &
         6.0 * (SMech_Cell%Elastic_Stress(4,:)**2 + SMech_Cell%Elastic_Stress(5,:)**2 + SMech_Cell%Elastic_Stress(6,:)**2))/2.)
    ! Not sure how to define this yet
    Stress_Strain(:)%eff_plastic_strain = sqrt(((SMech_Cell%Plastic_Strain(1,:) - SMech_Cell%Plastic_Strain(2,:))**2 + &
                                           (SMech_Cell%Plastic_Strain(2,:) - SMech_Cell%Plastic_Strain(3,:))**2 + &
                                           (SMech_Cell%Plastic_Strain(3,:) - SMech_Cell%Plastic_Strain(1,:))**2) * 2.0/9.0 + &
                                          (SMech_Cell%Plastic_Strain(4,:)**2 + SMech_Cell%Plastic_Strain(5,:)**2 + &
                                           SMech_Cell%Plastic_Strain(6,:)**2) * 4.0/3.0)
    ! Mean stress and strain
    Stress_Strain(:)%mean_stress = (SMech_Cell%Elastic_Stress(1,:) + SMech_Cell%Elastic_Stress(2,:) + &
                                    SMech_Cell%Elastic_Stress(3,:))/3.0
    Stress_Strain(:)%volumetric_strain =  SMech_Cell%Total_Strain(1,:) + SMech_Cell%Total_Strain(2,:) + &
                                          SMech_Cell%Total_Strain(3,:)
  end FUNCTION STRESS_STRAIN_INVARIANTS
  !
  !-----------------------------------------------------------------------------
  !
  Subroutine ELAS_VP_RESIDUAL (X_old, X, Residual)
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !   Calculate Ax + B(x) - RHS = 0, where A represents a discrete operator for 
    !   the thermomechanical equations of static equilibrium, B(x) represents the 
    !   viscoplastic strain correction to the elastic stress, and RHS represents 
    !   source terms that are due to thermal strains and volume changes from phase 
    !   change, etc.  This is in the form of a vector of forces calculated from
    !   a surface integral over the control volumes.
    !
    !---------------------------------------------------------------------------
    Use mesh_module,            Only: Mesh, Cell_Edge, GAP_ELEMENT_1
    Use parameter_module,       Only: ndim, ncells, nnodes, nvc, ncomps
    use node_operator_module,   only: CV_Internal, nipc, stress_reduced_integration
    use linear_module,          only: LINEAR_GRAD
    Use gs_module,              Only: EN_GATHER, EN_SUM_Scatter
    use mech_bc_data_module
    Use UbikSolve_module
    Use viscoplasticity,        only: PLASTIC_STRAIN_INCREMENT
    use solid_mech_constraints, only: DISPLACEMENT_CONSTRAINTS

    ! Argument list
    real(r8), Dimension(:), Intent(IN) :: X, X_old
    real(r8), Dimension(:), Intent(OUT) :: Residual
    integer :: status

    ! Local variables
    integer :: ip, icell, idim, icomps
    integer :: ix, iy, iz, i, j, k
    real(r8), Dimension(nvc,ncells)         :: Y_Temp
    real(r8), Dimension(ndim*nnodes)        :: Y
    real(r8), Dimension(nvc,ncells,ndim)    :: U
    real(r8), Dimension(ndim,ncells)        :: Arow
    real(r8), Dimension(ndim,ndim,ncells)   :: Ugrad
    real(r8), Dimension(ncomps,ncells)      :: Ip_Pl_Inc
    real(r8), Dimension(ncomps,ncells,nipc) :: Ipstress
    real(r8), Dimension(nnodes)             :: Lame2_Sum

    !---------------------------------------------------------------------------

    ! If using reduced integration, then calculate only one stress per cell.  Otherwise
    ! calculate the elastic stress at the centroid of each controlvolume face within the cell.
    if (stress_reduced_integration) then 
       ! This option is currently deprecated and not functional
       call TLS_fatal ('ELAS_VP_RESIDUAL: stress_reduced_integration is not currently a valid option')
    else
       ! Calculate stress tensors for all integration points in all cells
       U = 0.0
       do idim=1,ndim
          call EN_GATHER (U(:,:,idim),X(idim:(ndim*(nnodes-1)+idim):ndim))
       end do
       IP_STRESS_LOOP: do ip = 1,nipc
          ! Calculate displacement gradient and strain components
          Ugrad(:,:,:) = 0.0_r8
          do i=1,ndim
             call LINEAR_GRAD(ncells,CV_Internal%Face_Coord(:,ip,:),U(:,:,i),Arow(:,:))
             do j = 1,ndim
                do k=1,ndim
                   Ugrad(i,j,:)= Ugrad(i,j,:) + Arow(k,:) * CV_Internal%Face_Ijac(k,j,ip,:)
                end do
             end do
          end do
          ! Total strains
          SMech_IP(ip)%Total_Strain(1,:) = Ugrad(1,1,:)
          SMech_IP(ip)%Total_Strain(2,:) = Ugrad(2,2,:)
          SMech_IP(ip)%Total_Strain(3,:) = Ugrad(3,3,:)
          SMech_IP(ip)%Total_Strain(4,:) = (Ugrad(1,2,:) + Ugrad(2,1,:)) / 2
          SMech_IP(ip)%Total_Strain(5,:) = (Ugrad(1,3,:) + Ugrad(3,1,:)) / 2
          SMech_IP(ip)%Total_Strain(6,:) = (Ugrad(3,2,:) + Ugrad(2,3,:)) / 2
          ! Set total strains to zero in gap elements

          do icell = 1,ncells
             do icomps = 1,ncomps
                if (Mesh(icell)%Cell_Shape >= GAP_ELEMENT_1) then
                   SMech_IP(ip)%Total_Strain(icomps,icell) = 0.0
                end if
             end do
          end do
          ! Calculate elastic stresses, plastic strain rates and left hand side stress (Ipstress)

          ! Inputs:
          !         Dev_Stress_IP(:,:,ip)
          !         SMech_IP_Old(ip)%Elastic_Stress(:,:)
          !         SMech_IP_Old(ip)%Plastic_Strain_Rate(:)
          !         Eff_Stress_IP_old(:,ip)
          !         SMech_IP(ip)%Total_Strain(:,:)
          !         SMech_IP_Old(ip)%Total_Strain(:,:)
          !         SMech_IP_Old(ip)%Plastic_Strain(:,:))

          ! Outputs:
          !        Ip_Pl_Inc                            Plastic strain increment for one integration point
          !        SMech_IP(ip)%Plastic_Strain_Rate(:)  Plastic strain rate
          !        SMech_IP(ip)%Elastic_Stress(:,:)     Elastic stress
          !        Ipstress(:,:,ip)                     LHS stress (elastic stress + thermal ans PC stresses)
          
          Call PLASTIC_STRAIN_INCREMENT(Ip_Pl_Inc,                               &
                                        SMech_IP(ip)%Plastic_Strain_Rate(:),     &
                                        SMech_IP(ip)%Elastic_Stress(:,:),        &
                                        Ipstress(:,:,ip),                        &
                                        Dev_Stress_IP(:,:,ip),                   &
                                        SMech_IP_Old(ip)%Elastic_Stress(:,:),    &
                                        SMech_IP_Old(ip)%Plastic_Strain_Rate(:), &
                                        Eff_Stress_IP_old(:,ip),                 &
                                        SMech_IP(ip)%Total_Strain(:,:),          &
                                        SMech_IP_Old(ip)%Total_Strain(:,:),      &
                                        SMech_IP_Old(ip)%Plastic_Strain(:,:))
          ! Update plastic strain
          SMech_IP(ip)%Plastic_Strain(:,:) = SMech_IP_Old(ip)%Plastic_Strain(:,:) + Ip_Pl_Inc

       end do IP_STRESS_LOOP
    end if

    ! Initialize Y
    Y(:) = 0.0_r8
    ! Accumulate Y, looping over cells for each integration point
    NDIM_LOOP: do idim = 1,ndim
       Y_Temp(:,:) = 0.0_r8
       select case (idim)
          case (1)
             ix=1; iy = 4; iz = 5
          case (2)
             ix=4; iy = 2; iz = 6
          case(3)
             ix=5; iy = 6; iz = 3
          end select
       NIPC_LOOP: do ip=1,nipc
          do icell=1,ncells
             ! For first vertex of the edge that uses this integration point:
             !     Y = Y + sigma_ij * n_j * area
             Y_Temp(Cell_Edge(1,ip),icell) = Y_Temp(Cell_Edge(1,ip),icell) + &
                  (Ipstress(ix,icell,ip)* CV_Internal%Face_Normal(1,ip,icell) + &
                  Ipstress(iy,icell,ip)* CV_Internal%Face_Normal(2,ip,icell) + &
                  Ipstress(iz,icell,ip)* CV_Internal%Face_Normal(3,ip,icell)) * &
                  CV_Internal%Face_Area(ip,icell)
             !
             ! For the second vertex:
             !     Y = Y + sigma_ij * (-n_j) * area
             Y_Temp(Cell_Edge(2,ip),icell) = Y_Temp(Cell_Edge(2,ip),icell) - &
                  (Ipstress(ix,icell,ip)* CV_Internal%Face_Normal(1,ip,icell) + &
                  Ipstress(iy,icell,ip)* CV_Internal%Face_Normal(2,ip,icell) + &
                  Ipstress(iz,icell,ip)* CV_Internal%Face_Normal(3,ip,icell)) * &
                  CV_Internal%Face_Area(ip,icell)
          end do
       end do NIPC_LOOP

       call EN_SUM_Scatter(Y(idim:(ndim*(nnodes-1)+idim):ndim),Y_Temp)
    end do NDIM_LOOP

    ! If a node is connected only to cells with zero shear strength,
    ! set displacement to zero
    call EN_SUM_Scatter(Lame2_Sum, Lame2)
    do idim = 1,ndim
       Where(Lame2_Sum(:) == 0.0)
          Y(idim:(ndim*(nnodes-1)+idim):ndim) = X(idim:(ndim*(nnodes-1)+idim):ndim)
       end Where
    end do

    Y = Y/cscale

    ! Add displacement constraints (contact and interface constraints) here.
    
    call DISPLACEMENT_CONSTRAINTS(X,Y)

    Residual = Y - RHS
    status = 0

  End Subroutine ELAS_VP_RESIDUAL

  !
  Subroutine VP_MATVEC (X_vec, Y, status)
    !=======================================================================
    ! Purpose:
    !
    !   Compute the matrix-vector multiply y = J*v = F(h + eps) - F(h)]/eps
    !   where J is the Jacobian needed for the matrix-free Newton-Krylov method
    !   and F(h) is the nonlinear function for which solutions F(h)=0 are
    !   being sought.
    !
    !=======================================================================
    use lnorm_module,       only: L1NORM, L2NORM
    use nonlinear_solution, only: P_Residual, P_Future, p_control, P_Past
    use parameter_module,   only: ndim, nnodes, nnodes_tot
    use UbikSolve_module

    ! arguments
    type (Ubik_vector_type), intent(INOUT) :: X_vec
    real(r8), dimension(:), target, intent(INOUT) :: Y
    integer, intent(OUT) :: status

    ! Local Variables
    real(r8), dimension(ndim * nnodes) :: Perturbed_Residual, Perturbed_X
    real(r8) :: pert, vnorm
    real(r8), dimension(:), pointer :: X
    !
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    !
    X => Ubik_values_ptr(X_vec)

    ! Define scalar pert and perturbed vector.
    pert = L1NORM(P_Future)
    vnorm = L2NORM(X)
    if (vnorm == 0.0_r8) vnorm = 1.0_r8

    pert = p_control%eps_NK * pert / nnodes_tot / vnorm + p_control%eps_NK
    Perturbed_X = X * pert + P_Future

    ! Define perturbed residual: (this is the only line a user should change!)
    call ELAS_VP_RESIDUAL (P_Past, Perturbed_X, Perturbed_Residual)

    ! Perform approximate "matvec"  (note sign!)
    Y = (Perturbed_Residual - P_Residual) / pert  

    status = 0

  END SUBROUTINE VP_MATVEC
  
 !!
 !! DEFINE_TM_DENSITY_PROPERTY
 !!
 !! Defines the temperature-dependent density function for each immobile
 !! material phase using the specified reference density and temperature
 !! values and the linear CTE function.
 !!
  
  subroutine define_tm_density_property (stat, errmsg)
  
    use parameter_module, only: nmat
    use material_interop, only: void_material_index, material_to_phase
    use fluid_data_module, only: isImmobile
    use phase_property_table
    use scalar_func_class
    use scalar_func_tools, only: is_const
    use tm_density, only: alloc_tm_density_func
    
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg
    character(:), allocatable :: errm
    
    integer :: m, ref_dens_id, ref_temp_id, cte_id, dens_id, phase_id
    real(r8) :: dens0, temp0, state(0)
    class(scalar_func), allocatable :: temp0_fun, dens0_fun, cte_fun, rho_fun
    
    ASSERT(ppt_has_property('TM reference density'))
    ref_dens_id = ppt_property_id('TM reference density')
    
    ASSERT(ppt_has_property('TM reference temperature'))
    ref_temp_id = ppt_property_id('TM reference temperature')
    
    ASSERT(ppt_has_property('TM linear CTE'))
    cte_id = ppt_property_id('TM linear CTE')
    
    call ppt_add_property ('TM density', dens_id)
    
    do m = 1, nmat
      if (m == void_material_index .or. .not.isImmobile(m)) cycle
      phase_id = material_to_phase(m)
      ASSERT(phase_id > 0)
      
      !! Get the value of the (constant) reference density.
      call ppt_get_phase_property (phase_id, ref_dens_id, dens0_fun)
      ASSERT(allocated(dens0_fun))
      ASSERT(is_const(dens0_fun))
      dens0 = dens0_fun%eval(state)
      
      !! Get the value of the (constant) reference temperature.
      call ppt_get_phase_property (phase_id, ref_temp_id, temp0_fun)
      ASSERT(allocated(temp0_fun))
      ASSERT(is_const(temp0_fun))
      temp0 = temp0_fun%eval(state)
      
      !! Get the CTE function.
      call ppt_get_phase_property (phase_id, cte_id, cte_fun)
      ASSERT(allocated(cte_fun))
      
      !! Create the temperature-dependent density function.
      call alloc_tm_density_func (rho_fun, dens0, temp0, cte_fun, stat, errm)
      if (stat /= 0) then
        errmsg = 'problem with the "TM linear CTE" property: ' // trim(errm)
        return
      end if
      
      !! Assign the density function as the TM density property.
      if (ppt_has_phase_property(phase_id, dens_id)) then
        stat = -1
        errmsg = 'found conflicting "TM density" property for phase "' // &
                 trim(ppt_phase_name(phase_id)) // '"'
        return
      end if
      call ppt_assign_phase_property (phase_id, dens_id, rho_fun)
    end do
    
    stat = 0
    errmsg = ''

  end subroutine define_tm_density_property
  
 !!
 !! COMPUTE_CELL_PROPERTY
 !!
 !! This computes the specified property for the cells on the original Truchas
 !! mesh. The property is one given in a PHASE namelist PROPERTY_NAME variable,
 !! or one created internally by Truchas, and the associated property value is
 !! either constant or a function of temperature only.  The routine essentially
 !! handles the material averaging of the property over a cell using the volume
 !! fraction data from MATL.  The void material (one with a MATERIAL namelist
 !! density of zero) contributes nothing to the property value; e.g., zero is
 !! returned for an entirely void cell.
 !! 
  
  subroutine compute_cell_property (prop, temp, value)
  
    use parameter_module, only: ncells, nmat
    use fluid_data_module, only: isImmobile
    use phase_property_table
    use material_interop, only: void_material_index, material_to_phase
    use matl_module, only: gather_vof
    use scalar_func_class
    use scalar_func_tools, only: is_const
    
    character(*), intent(in) :: prop
    real(r8), intent(in) :: temp(:)
    real(r8), intent(out) :: value(:)
    
    integer :: m, j, phase_id, prop_id
    real(r8) ::vofm (ncells), state(1), mval
    class(scalar_func), allocatable :: prop_fun
    
    ASSERT(size(temp) == ncells)
    ASSERT(size(value) == ncells)
    
    ASSERT(ppt_has_property(prop))
    prop_id = ppt_property_id(prop)
    
    value = 0.0_r8
    do m = 1, nmat
      if (m == void_material_index .or. .not.isImmobile(m)) cycle
      phase_id = material_to_phase(m)
      ASSERT(phase_id > 0)
      call ppt_get_phase_property (phase_id, prop_id, prop_fun)
      ASSERT(allocated(prop_fun))
      call gather_vof (m, vofm)
      if (is_const(prop_fun)) then
        mval = prop_fun%eval(state)  ! state is ignored, but required
        value = value + mval*vofm
      else
        do j = 1, ncells
          if (vofm(j) > 0.0_r8) then
            state(1) = temp(j)
            mval = prop_fun%eval(state)
            value(j) = value(j) + mval*vofm(j)
          end if
        end do
      end if
    end do

  end subroutine compute_cell_property

End Module SOLID_MECHANICS_MODULE
