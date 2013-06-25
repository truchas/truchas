MODULE DIAGNOSTICS_MODULE
  !=======================================================================
  ! Purpose:
  !
  !    Data and procedures for diagnostic quantities for output
  !
  ! Contains: DIAGNOSTICS
  !           PROBES_FIELDS
  !           PROBES_POSITIONS
  !           DIVERGENCE
  !
  ! Author(s): Kin Lam (klam@lanl.gov)
  !            Sharen Cummins (scummins@lanl.gov)
  !=======================================================================
  implicit none
  private

  ! public procedures
  public :: DIAGNOSTICS, PROBES_FIELDS, PROBES_POSITIONS, DIVERGENCE, &
            get_global_cell, get_global_node

  logical, save, public :: solidification_diagnostics = .false., &
                           alloy_solid_diagnostics    = .false.
  
CONTAINS

 !!
 !! NNC, 2 May 2012 -- As part of removing old HT, DIAGNOSTICS, which dealt with
 !! solidification times, has been gutted leaving just a stub.  We need to add
 !! back the appropriate analysis for the new heat transfer solver.
 !!

  SUBROUTINE DIAGNOSTICS ()
  
    logical, save :: first_time = .true.

    ! Decide whether to do various diagnostics calculations
    ! Create (and initialize) or nullify pointer arrays
    if (first_time) then
       ! If problem does not involve phase change, no solidification diagnostics
       solidification_diagnostics = .false.
       alloy_solid_diagnostics    = .false.
       if (solidification_diagnostics) then
          !TODO: Do something here based on the new phase change.
       end if
       first_time = .false.
       return   ! don't do diagnostics calculations the first time
    end if

    if (solidification_diagnostics) then
      !TODO: do something for the new phase change
    end if

  END SUBROUTINE DIAGNOSTICS

  SUBROUTINE PROBES_POSITIONS
  !=======================================================================
  ! Purpose:
  !
  !    For all probes:
  !    Calculates the nearest CELL index and resulting CELL centroids to the probe location
  !    Calculates the nearest NODE index and resulting NODE coords to the probe location
  !    Initialise the probe's scalar,vector,tensor fields 

  ! Author(s): Sharen Cummins (scummins@lanl.gov)
  !=======================================================================

    use kind_module,       only: real_kind,int_kind
    use mesh_module,       only: Cell, Vertex, UnPermute_Mesh_Vector, UnPermute_Vertex_Vector
    use PGSLIB_module,     only: pgslib_global_sum
    use probe_module,      only: probes
    use parameter_module,  only: nprobes, ndim

    implicit none

    ! Local variables.
    integer(KIND=int_kind) :: i, j, icell, inode, cellindex, nodeindex
    real(KIND=real_kind)   :: multiplicity, cellcoords(ndim), nodecoords(ndim)

    do i=1,nprobes
       call get_global_cell(probes(i)%coords,icell,multiplicity)
       
       if (icell > -1) then
          cellcoords = Cell(icell)%Centroid
          cellindex  = UnPermute_Mesh_Vector(icell)
       else
          cellindex  = 0
          cellcoords = 0.0
       end if

       probes(i)%cell%index        = pgslib_global_sum(cellindex) !i.e in global numbering
       do j=1,ndim
          probes(i)%cell%coords(j) = pgslib_global_sum(cellcoords(j))
       end do

       call get_global_node(probes(i)%coords,inode,multiplicity)

       if (inode > -1) then
          nodecoords = Vertex(inode)%Coord
          nodeindex  = UnPermute_Vertex_Vector(inode)
       else
          nodeindex  = 0
          nodecoords = 0.0
       end if

       probes(i)%node%index  = pgslib_global_sum(nodeindex) !i.e in global numbering
       do j=1,ndim
          probes(i)%node%coords(j) = pgslib_global_sum(nodecoords(j))
       end do

    end do

    !now get the initial field values at the probe locations
    call PROBES_FIELDS

  END SUBROUTINE PROBES_POSITIONS


  SUBROUTINE PROBES_FIELDS
  !=======================================================================
  ! Purpose:
  !
  !    Calculate scalar,vector,tensor field quantities at probe locations, for all probes
  !    For CELL fields, the result is the field value at the nearest CELL
  !    For NODAL fields, the result is the field value at the nearest NODE
  !
  ! Author(s): Sharen Cummins (scummins@lanl.gov)
  !=======================================================================

    use kind_module,          only: real_kind, int_kind, log_kind
    use zone_module,          only: Zone
    use fluid_data_module,    only: fluid_flow, boussinesq_approximation
    use property_module,      only: get_density_delta
    use time_step_module,     only: dt
    use discrete_op_module,   only: GRADIENT_CELL
    use matl_module,          only: GATHER_VOF
    use EM_data_proxy,        only: joule_power_density, EM_is_on
    use solid_mechanics_data, only: SMech_Cell, Thermal_Strain, PC_Strain, &
                                    Displacement, Node_Gap, Node_Norm_Trac, solid_mechanics
    use probe_module,         only: probes
    use parameter_module,     only: string_len, ncells, ndim, nmat, nprobes
    use constants_module,     only: zero

    implicit none

    ! Local variables.
    integer(KIND=int_kind)    :: i, j, k, n, m, count, theindex, vfieldsize, tfieldsize
    character(LEN=string_len) :: themeshspace

    type afield
       real(KIND=real_kind), dimension(:), pointer   :: field
    end type afield
    
    type(afield), dimension(:),   pointer               :: scalarfields
    type(afield), dimension(:,:), pointer               :: vectorfields
    type(afield), dimension(:,:), pointer               :: tensorfields
    real(KIND=real_kind), dimension(:,:),   pointer     :: tmp_2darray 
    real(KIND=real_kind), dimension(ncells),target      :: tmp,tmp2,tmp3,tmp4
    real(KIND=real_kind), dimension(nmat,ncells),target :: vof
    real(KIND=real_kind), dimension(ncells)             :: tmp_1darray

    real(real_kind),  dimension(ncells) :: mVOF
    logical(log_kind),dimension(ncells) :: mMask

    scalarfields => NULL()
    vectorfields => NULL()
    tensorfields => NULL()
    tmp_2darray  => NULL()

    !Scalar probe fields are stored first in the following order:
    !    a) cell fields  : [Rho,Delta_Rho,Temp,dT_dt,|Grad(T)|,Enthalpy,P,PC,Vof,Joule,Epsdot]
    !    b) nodal fields : [GapDisp,GapForce]

    !Vector probe fields are stored next in the following order:
    !    a) cell fields  : [Vc]
    !    b) nodal fields : [Disp]

    !Tensor probe fields are stored next in the following order:
    !    a) cell fields  : [Elastic Stress, Total Strain, Thermal Strain, PC Strain]
    !    b) nodal fields : []

    !***********************************************************************************
    !
    !NOTE:This ordering must remain consistent with the ordering provided in the 
    !initialization subroutine 'PROBE_INIT' in src/setup/initialize/init_module.F90
    !
    !***********************************************************************************

    if (nprobes > 0) then

       if (SIZE(probes(1)%ScalarVarLU) > 0) then

          ALLOCATE(scalarfields(SIZE(probes(1)%ScalarVarLU))) 
          ALLOCATE(tmp_2darray(ndim,ncells))

          scalarfields(1)%field => Zone(:)%Rho

          ! calculate delta-rho 
          tmp = 0 
          if ( fluid_flow .and. boussinesq_approximation) then 
             call get_density_delta (Zone%Temp, tmp)
          endif
          scalarfields(2)%field => tmp
          scalarfields(3)%field => Zone(:)%Temp
          tmp2                  = (Zone(:)%Temp-Zone(:)%Temp_Old)/dt
          scalarfields(4)%field => tmp2
          call GRADIENT_CELL(tmp_2darray, Zone(:)%Temp) 
          tmp_1darray  = zero
          do n = 1,ndim 
             tmp_1darray = tmp_1darray + tmp_2darray(n,:)**2 
          end do
          tmp3                  = SQRT(tmp_1darray(:))
          scalarfields(5)%field => tmp3 
          scalarfields(6)%field => Zone(:)%Enthalpy
          scalarfields(7)%field => Zone(:)%P
          count = 8
          do j=1,nmat
             call GATHER_VOF(j,vof(j,:))
             scalarfields(count)%field      => vof(j,:)
             count = count + 1
          end do
          if (EM_is_on()) then
             tmp4                           = joule_power_density()
             scalarfields(count)%field      => tmp4
             count = count + 1
          end if
          if (solid_mechanics) then
             scalarfields(count)%field      => SMech_Cell%Plastic_Strain_Rate(:)
             count = count + 1
          end if

          !get all node centered scalar fields
          if (solid_mechanics) then
             scalarfields(count)%field      => Node_Gap(:,1)
             count = count + 1
             scalarfields(count)%field      => Node_Norm_Trac(:,1)
             count = count + 1
          end if

          DEALLOCATE(tmp_2darray)

       end if

       if (SIZE(probes(1)%VectorVarLU) > 0) then

          vfieldsize    = SIZE(probes(1)%VectorVarLU(1)%field)
          ALLOCATE(vectorfields(SIZE(probes(1)%VectorVarLU),vfieldsize))

          !get all cell-centered vector fields
          count         = 1
          do j=1,vfieldsize
             vectorfields(count,j)%field    => Zone(:)%Vc(j)
          end do
          count = count + 1     

          !get all node-centered vector fields
          if (solid_mechanics) then
             do j=1,vfieldsize
                vectorfields(count,j)%field => Displacement(j,:) 
             end do
             count = count + 1 
          end if

       end if

       if (SIZE(probes(1)%TensorVarLU) > 0) then
          tfieldsize    = SIZE(probes(1)%TensorVarLU(1)%field)

          ALLOCATE(tensorfields(SIZE(probes(1)%TensorVarLU),tfieldsize))

          !get all cell-centered tensor fields
          count = 1

          if (solid_mechanics) then
             do j=1,tfieldsize
                tensorfields(count,j)%field => SMech_Cell%Elastic_Stress(j,:)
             end do
             count = count + 1
             do j=1,tfieldsize
                tensorfields(count,j)%field => SMech_Cell%Total_Strain(j,:)
             end do
             count = count + 1
             do j=1,tfieldsize
                tensorfields(count,j)%field => Thermal_Strain(j,:)
             end do
             count = count + 1
             do j=1,tfieldsize
                tensorfields(count,j)%field => PC_Strain(j,:)
             end do
             count = count + 1
          end if
          !get all node-centered tensor fields....

       end if

    end if


    do i=1,nprobes
       
       do j=1,SIZE(probes(i)%ScalarVarLU)
          themeshspace  = probes(i)%ScalarVarLU(j)%meshspace
          call get_probe_field(themeshspace,probes(i)%coords,scalarfields(j)%field,theindex,probes(i)%ScalarVarLU(j)%field)
       end do

       do j=1,SIZE(probes(i)%VectorVarLU)
          themeshspace  = probes(i)%VectorVarLU(j)%meshspace
          do k=1,vfieldsize
             call get_probe_field(themeshspace,probes(i)%coords,vectorfields(j,k)%field,theindex,probes(i)%VectorVarLU(j)%field(k))
          end do
       end do

       do j=1,SIZE(probes(i)%TensorVarLU)
          themeshspace  = probes(i)%TensorVarLU(j)%meshspace
          do k=1,tfieldsize
             call get_probe_field(themeshspace,probes(i)%coords,tensorfields(j,k)%field,theindex,probes(i)%TensorVarLU(j)%field(k))
          end do
       end do

    end do

    if (ASSOCIATED(scalarfields)) DEALLOCATE(scalarfields)
    if (ASSOCIATED(vectorfields)) DEALLOCATE(vectorfields)
    if (ASSOCIATED(tensorfields)) DEALLOCATE(tensorfields)


  END SUBROUTINE PROBES_FIELDS


  SUBROUTINE DIVERGENCE (D)

    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the cell-centered divergence (D) for diagnostic purposes
    !
    !=======================================================================
    use constants_module,  only: zero
    use fluid_data_module, only: fluidRho
    use fluid_data_module, only: Fluxing_Velocity
    use kind_module,       only: int_kind, real_kind
    use mesh_module,       only: Cell 
    use parameter_module,  only: nfc, ncells
    use time_step_module,  only: dt

    implicit none

    ! Arguments
    real(KIND = real_kind),   dimension(ncells), intent(OUT) :: D

    ! Local Variables
    integer(KIND = int_kind) :: f 

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Loop over faces and accumulate the divergence
    ! from these face-centered velocities
    D = zero
    do f = 1,nfc
       D = D + Fluxing_Velocity(f,:)*Cell%Face_Area(f)
    end do

    ! Normalize by cell volume
    D = D*dt/Cell%Volume

    ! Zero out divergences in void cells
    where (FluidRho == zero) D = zero

    ! Eliminate noise
!   D = MERGE(zero, D, ABS(D) <= alittle)

    return

  END SUBROUTINE DIVERGENCE

  subroutine get_global_cell (coordinates, icell_o, multiplicity)

    !=======================================================================
    ! Purpose(s):
    !
    !   get the global cell number that is closest to a given coordinate set
    !
    !=======================================================================


    use kind_module,          only: int_kind, real_kind
    use mesh_module,          only: Cell
    use parameter_module,     only: ndim, ncells
    use PGSLIB_module,        only: PGSLib_Global_SUM, pgslib_global_minval
    use constants_module,           only: one, zero

    implicit none

    ! Argument List

    integer(int_kind), intent(out)  :: icell_o
    real(real_kind), intent(out)    :: multiplicity
    real(real_kind), intent(in)     :: coordinates(ndim)

    ! Local Variables

    real(real_kind)                 :: tmp(ncells)
    real(real_kind)                 :: Minimum_distance 
    real(real_kind)                 :: Minimum_distance_global
    real(real_kind)                 :: Minimum_distance_global_tot
    integer(int_kind)               :: icell(1)
    
    ! Find cell at minimum distance from the point

    tmp(:) = (coordinates(1) - Cell(:)%Centroid(1))**2 +        &
         (coordinates(2) - Cell(:)%Centroid(2))**2 +            &
         (coordinates(3) - Cell(:)%Centroid(3))**2 

    Minimum_distance = minval(tmp)
    icell            = minloc(tmp)

    Minimum_distance_global = pgslib_global_minval(tmp)

    if ( Minimum_distance_global == Minimum_distance) then
       icell_o          = icell(1)
    else
       icell_o          = -1
       Minimum_distance = zero
    endif

    if (Minimum_distance_global == 0) then
       multiplicity                = one
    else
       Minimum_distance_global_tot = pgslib_global_sum(Minimum_distance)
       multiplicity                = Minimum_distance_global_tot / Minimum_distance_global 
    endif

  END SUBROUTINE get_global_cell 


  subroutine get_global_node (coordinates, inode_o, multiplicity)

    !=======================================================================
    ! Purpose(s):
    !
    !   get the global node number that is closest to a given coordinate set
    !
    !=======================================================================


    use kind_module,          only: int_kind, real_kind
    use mesh_module,          only: Vertex
    use parameter_module,     only: ndim, nnodes
    use PGSLIB_module,        only: PGSLib_Global_SUM, pgslib_global_minval
    use constants_module,           only: one, zero

    implicit none

    ! Argument List

    integer(int_kind), intent(out)  :: inode_o
    real(real_kind), intent(out)    :: multiplicity
    real(real_kind), intent(in)     :: coordinates(ndim)

    ! Local Variables

    real(real_kind)                 :: tmp(nnodes)
    real(real_kind)                 :: Minimum_distance 
    real(real_kind)                 :: Minimum_distance_global
    real(real_kind)                 :: Minimum_distance_global_tot
    integer(int_kind)               :: inode(1)
    
    ! Find node at minimum distance from the point

    tmp(:) = (coordinates(1) - Vertex(:)%Coord(1))**2 +        &
         (coordinates(2) - Vertex(:)%Coord(2))**2 +            &
         (coordinates(3) - Vertex(:)%Coord(3))**2 

    Minimum_distance = minval(tmp)
    inode            = minloc(tmp)

    Minimum_distance_global = pgslib_global_minval(tmp)

    if ( Minimum_distance_global == Minimum_distance) then
       inode_o          = inode(1)
    else
       inode_o          = -1
       Minimum_distance = zero
    endif

    if (Minimum_distance_global == 0) then
       multiplicity     = one
    else
       Minimum_distance_global_tot = pgslib_global_sum(Minimum_distance)
       multiplicity                = Minimum_distance_global_tot / Minimum_distance_global 
    endif

  END SUBROUTINE get_global_node

  SUBROUTINE get_probe_field (themeshspace, coordinates, field, findex, response)

    !=======================================================================
    ! Purpose(s):
    !
    !   get the field value at the index 'index' - applicable for CELL and NODE fields 
    !
    !=======================================================================


    use kind_module,                only: int_kind, real_kind
    use PGSLIB_module,              only: pgslib_global_sum
    use parameter_module,           only: ndim, string_len 
    use constants_module,           only: zero

    implicit none

    ! Argument List

    character(LEN=string_len), intent(in)             :: themeshspace
    real(real_kind),   intent(in), dimension(ndim)    :: coordinates
    real(real_kind),   intent(in), dimension(:)       :: field
    integer(int_kind), intent(out)                    :: findex
    real(real_kind),   intent(out)                    :: response

    ! Local Variables

    real(real_kind)   :: local_response, multiplicity
    integer(int_kind) :: index, theindex

    if (themeshspace == 'cell') then
       call get_global_cell(coordinates,index,multiplicity)
    else
       call get_global_node(coordinates,index,multiplicity)
    end if

    if (index > -1) then

       local_response = field(index)
       theindex       = index

    else
 
       local_response = zero
       theindex       = 0

    endif

    response =  pgslib_global_sum (local_response) / multiplicity
    findex   =  pgslib_global_sum (theindex)

  END SUBROUTINE get_probe_field
  
END MODULE DIAGNOSTICS_MODULE
