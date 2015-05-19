#include "f90_assert.fpp"

! set to 1 to turn on local debugging code
#define DEBUG_MAT_SLOT 0
MODULE INIT_MODULE
!=======================================================================
! Purpose(s):
!
!   Define procedures that initialize the Matl and Zone derived types.
!   The initialization is based on input variables read in from the
!   input file.
!
! Public Interface(s):
!
!   * call INITIAL ()
!
!     Initialize the Matl and Zone types.
!
! Contains: INITIAL
!           BC_INIT
!           MATL_INIT
!           ZONE_INIT
!           CHECK_VOF
!           VELOCITY_INIT
!
! Author(s): Douglas B. Kothe (dbk@lanl.gov)
!            Larry J. Cox (ljcox@lanl.gov
!
!=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  ! Public Procedures
  public :: INITIAL

  ! For testing:
  public :: check_vof_aux

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE INITIAL (t, dt)
    !=======================================================================
    ! Purpose(s):
    !
    !   Initialize the Matl, Zone, and Alloy (if necessary) derived types.
    !   Initialize all BC data structures and information.
    !   Initialize all probe structures
    !
    !=======================================================================
    use fluid_data_module,      only: Void_Material_Exists,     &
                                      Void_Material_Index,      &
                                      Void_Material_Count, fluid_flow
    use fluid_utilities_module, only: FLUID_INIT
    use interfaces_module,      only: nbody
    use overwrite_module,       only: OVERWRITE_BC, OVERWRITE_MATL,       &
                                      OVERWRITE_VEL, OVERWRITE_ZONE
    use parameter_module,       only: ncells, nmat
    use restart_variables,      only: restart
    use restart_driver,         only: restart_matlzone, restart_solid_mechanics, restart_species
    use property_module,        only: DENSITY_MATERIAL
    use zone_module,            only: Zone
    use var_vector_module
    use gs_module
    use solid_mechanics_input,  only: solid_mechanics
    use solid_mechanics_module, only: SOLID_MECH_INIT
    use vof_init,               only: VOF_INITIALIZE
    use diffusion_solver_data,  only: ds_enabled, num_species, &
                                      ds_sys_type, DS_SPEC_SYS, DS_TEMP_SYS, DS_TEMP_SPEC_SYS
    use diffusion_solver,       only: ds_init, ds_set_initial_state
    use material_interop,       only: generate_material_mappings
    use probe_output_module,    only: probe_init
    use physics_module,         only: heat_transport, heat_species_transport

    real(r8), intent(in) :: t, dt

    ! Local Variables
    integer :: m, stat
    real(r8) :: density
    logical :: found

    real(r8), dimension(nbody,ncells) :: Hits_Vol
    real(r8), dimension(nbody,ncells) :: volume_fractions
    real(r8), allocatable :: phi(:,:)
    character(200) :: errmsg
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    !! Initialize the MATERIAL_INTEROP module which provides tools for moving
    !! between Truchas materials and the new material software model introduced
    !! with the diffusion solver.
    call generate_material_mappings (stat, errmsg)
    if (stat /= 0) call TLS_fatal ('INITIAL: ' // trim(errmsg))

    ! Determine if a void material exists in this calculation and save its index if one does
    ! NOTE:  There is an inherent assumption that there exists only one void material in any problem
    Void_Material_Exists = .false.
    Void_Material_Index  = 0
    Void_Material_Count  = 0
    MATERIALS: do m = 1,nmat
       density = DENSITY_MATERIAL(m)
       if (density == 0.0_r8) then
          Void_Material_Exists = .true.
          Void_Material_Count = Void_Material_Count + 1
          Void_Material_Index(Void_Material_Count) = m
       end if
    end do MATERIALS


    ! Set up Thermodynamic Reference Temperatures and Enthalpies
    !call THERMO_REFERENCES

    ! vof_initialize returns the volume of each body in each cell
    ! hits_vol is historically used to hold this information
    ! I think this is converted to volume fractions in matl_init.
    ! That's hopelessly confusing and will be corrected later.

    call TLS_info ('')
    call TLS_info ('Computing initial volume fractions ... ')

    volume_fractions = vof_initialize()
    hits_vol = volume_fractions   ! temporary compatibility

    ! Either read Zone and Matl from a restart file or initialize them
    if (restart) then

       call restart_matlzone ()
       call restart_solid_mechanics ()

    else

       ! Initialize Matl%Cell and Zone%Vc.
       call MATL_INIT (Hits_Vol)

       call TLS_info ('  Initial volume fractions computed.')

    end if

    call property_init

    ! Initialize Zone%Rho and Zone%Temp
    call ZONE_INIT (Hits_Vol)

    ! Initialize probe%name, probe%description, probe%coords
    call PROBE_INIT ()

    ! [sriram] End of original restart if condition

    ! Initialize BC quantities.
    ! Must be called before FLUID_INIT as Pressure_BC is required
    ! for the flow restart when a non-ortho operator is used.
    call BC_INIT()

    !! NNC, December 2012.  Flow initialization is really messed up.  It does
    !! necessary stuff even when flow is inactive.  The work of ensuring the
    !! relevant properties are defined belongs there but for now it's here so
    !! it doesn't get lost in the mess that is fluid_init.  It needs to be done
    !! always because fluid_init uses its results regardless of whether flow is on.
    call flow_property_init
    call FLUID_INIT (t)

    ! Allow arbitrary overwriting of the Matl, Zone, and BC types.
    call OVERWRITE_MATL ()
    call OVERWRITE_ZONE ()
    call OVERWRITE_VEL  ()
    call OVERWRITE_BC   ()

    ! Calculate initial stress-strain field and displacements if solid mechanics is active
    if (solid_mechanics) Call SOLID_MECH_INIT

    ! Get the initial species concentration fields.
    if (ds_enabled .and. num_species > 0) then
      allocate(phi(ncells,num_species))
      if (restart) then
        call restart_species (phi, found)
        if (.not.found) call init_conc (hits_vol, phi)
      else
        call init_conc (hits_vol, phi)
      end if
    else if (restart) then
      call restart_species () ! skip any species data that may be present
    end if

    ! Initialize the diffusion solver.
    if (ds_enabled) then
      call ds_init
      select case (ds_sys_type)
      case (DS_SPEC_SYS)
        call ds_set_initial_state (t, dt, conc=phi)
        deallocate(phi)
      case (DS_TEMP_SYS)
        call ds_set_initial_state (t, dt, temp=zone%temp)
      case (DS_TEMP_SPEC_SYS)
        call ds_set_initial_state (t, dt, temp=zone%temp, conc=phi)
        deallocate(phi)
      end select
    end if

  END SUBROUTINE INITIAL

  subroutine property_init

    use material_table
    use material_utilities
    use physics_module

    integer :: stat
    character(128) :: errmsg
    character(:), allocatable :: errmsg1

    integer, allocatable :: matids(:)

    allocate(matids(mt_num_material()))
    call mt_get_material_ids (matids)

    call required_property_check (matids, 'density', stat, errmsg)
    if (stat /= 0) call TLS_fatal ('PROPERTY_INIT: ' // errmsg)

    call constant_property_check (matids, 'density', stat, errmsg)
    if (stat /= 0) call TLS_fatal ('PROPERTY_INIT: ' // errmsg)

    if (heat_transport .or. heat_species_transport) then
      call required_property_check (matids, 'specific heat', stat, errmsg)
      if (stat /= 0) call TLS_fatal ('PROPERTY_INIT: ' // errmsg)

      call define_enthalpy_density_property (matids, stat, errmsg1)
      if (stat /= 0) call TLS_fatal ('PROPERTY_INIT: ' // errmsg1)
    else ! unused but expected for initialization and output.
      call request_property (matids, 'specific heat', default = 0.0_r8)
      call define_enthalpy_density_property (matids, stat, errmsg1)
      if (stat /= 0) call TLS_fatal ('PROPERTY_INIT: ' // errmsg1)
    end if

  end subroutine property_init

  subroutine flow_property_init

    use fluid_data_module, only: boussinesq_approximation
    use viscous_data_module, only: inviscid
    use property_module, only: request_fluid_property

    if (boussinesq_approximation) then
      call request_fluid_property ('density deviation', default=0.0_r8)
    end if

    if (.not.inviscid) then
      call request_fluid_property ('viscosity')
    end if

  end subroutine flow_property_init

  subroutine init_conc (hits_vol, phi)

    use parameter_module, only: ncells
    use mesh_module, only: mesh, GAP_ELEMENT_1
    use interfaces_module, only: body_phi

    real(r8), intent(in) :: hits_vol(:,:)
    real(r8), intent(out) :: phi(:,:)

    integer :: j, nbody, nconc

    nbody = size(hits_vol,dim=1)
    nconc = size(phi,dim=2)
    do j = 1, ncells
      if (mesh(j)%cell_shape < GAP_ELEMENT_1) then  ! valid HITS_VOL data
        phi(j,:) = matmul(hits_vol(:,j), body_phi(:nbody,:nconc)) / sum(hits_vol(:,j))
      else  ! a gap cell with bogus HITS_VOL data?
        phi(j,:) = 0.0
      end if
    end do

  end subroutine init_conc

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE BC_INIT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Initialize boundary condition quantities in the BC structure.
    !=======================================================================
    use bc_data_module,         only: BC_Type, BC_Name,                                  &
                                      BC_Value, BC_Table, BC_Variable, Inflow_Material,  &
                                      Inflow_Index, nbc_surfaces, BC_Pressure,           &
                                      BC_Conc, BC_Prs, BC_Mat, bndry_vel, BC_Temp, BC_Zero, & !BC_Vel
                                      Conic_XX, Conic_YY,                                &
                                      Conic_ZZ, Conic_XY, Conic_XZ, Conic_YZ, Conic_X,   &
                                      Conic_Y, Conic_Z, Conic_Constant, Conic_Tolerance, &
                                      Bounding_Box, Inflow_Temperature,                  &
                                      Surface_Name,                                      &
                                      Conic_Relation, Surfaces_In_This_BC,               &
                                      Surface_Materials, Srfmatl_Index,                  &
                                      Mesh_Face_Set, Mesh_Surface
    use bc_module,              only: DIRICHLET, FREE_SLIP, DIRICHLET_VEL,               &
                                      SET_DIRICHLET, NEUMANN_VEL,                        &
                                      SET_NO_VEL_BC, SET_FREE_SLIP, SET_DIRICHLET_VEL,   &
                                      SET_INTERNAL_BC, SET_NEUMANN, SET_NEUMANN_VEL,     &
                                      BC, InitializeBoundaryConditions, Prs, Vel
    use bc_operations
    use bc_displacement_init,   only: Initialize_Displacement_BC, Node_Set_BC_Init,         &
                                      append_to_displacement_bc,                            &
                                      Make_Displacement_BC_Atlases, Interface_Surface_Id
    use bc_pressure_init,       only: Initialize_Pressure_BC
    use input_utilities,        only: NULL_I, NULL_R
    use fluid_data_module,      only: fluid_flow
    use gs_module,              only: EE_GATHER
    use mesh_module,            only: Cell, Mesh, DEGENERATE_FACE
    use parameter_module,       only: ncells, ndim, nfc, nvc
    use pgslib_module,          only: PGSLIB_GLOBAL_COUNT, PGSLIB_GLOBAL_SUM
    use projection_data_module, only: dirichlet_pressure
    use property_module,        only: Get_Truchas_Material_Id
    use solid_mechanics_input,  only: solid_mechanics
    use physics_module,         only: heat_transport, heat_species_transport
    use string_utilities,       only: i_to_c
    use vector_func_factories

    ! Local Variables
    integer                              :: bit_position, f, n, nf, p, &
                                                             m, c, mat1, mat2, set1, nssets
    logical, dimension(:),   allocatable :: Mask1, Mask2
    logical, dimension(:,:), allocatable :: Disp_Mask
    logical, dimension(:,:), allocatable :: Found_Faces
    integer, allocatable, dimension(:,:,:,:) :: Ngbr_Mesh_Face_Set
    real(r8), dimension(:),   allocatable :: Coeff_XX
    real(r8), dimension(:),   allocatable :: Coeff_X
    real(r8), dimension(:),   allocatable :: Area
    integer,  dimension(:),   allocatable :: Faces
    real(r8), dimension(:,:), allocatable :: tmp_r2
    integer,  dimension(:,:), allocatable :: Tmp_BC
    real(r8), dimension(:),   allocatable :: Tmp1, Tmp2
    integer :: istatus
    character(128) :: message
    class(vector_func), pointer :: vel_func

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

    ALLOCATE(Mask1(ncells),                      &
             Mask2(ncells),                      &
             Disp_Mask(nfc,ncells),              &
             Found_Faces(nfc,nbc_surfaces),      &
             Coeff_XX(3*(ndim-1)),               &
             Coeff_X(ndim),                      &
             Area(nbc_surfaces),                 &
             Faces(nbc_surfaces),                &
             tmp_r2(nfc,ncells),                 &
             Tmp_BC(nvc,ncells),                 &
             Tmp1(ncells),                       &
             Tmp2(ncells),                       &
             Stat = istatus)
    if (istatus /= 0) call TLS_panic ('BC_INIT: allocation error')

    ! Set the number of side sets to zero to indicate we haven't setup
    !  the parallel data structures yet
    nssets      = 0

    ! Nullify the pointers to BC arrays, and then initialize
    !NULLIFY (BC_Vel, BC_Mat, BC_Prs, BC_Conc, BC_Temp, BC_Zero, BC_Pressure)
    NULLIFY (BC_Mat, BC_Prs, BC_Conc, BC_Temp, BC_Zero, BC_Pressure)

    if (fluid_flow) then
       !ALLOCATE (BC_Vel(ndim,nfc,ncells))
       !BC_Vel = NULL_R
       ALLOCATE (BC_Pressure(nfc,ncells))
       BC_Pressure = 0.0_r8
       ALLOCATE (BC_Zero(nfc,ncells))
       BC_Zero = 0.0_r8
       ALLOCATE (BC_Mat(nfc,ncells))
       BC_Mat = NULL_I
       if (heat_transport .or. heat_species_transport) then
          ALLOCATE (BC_Temp(nfc,ncells))
          BC_Temp = NULL_R
       end if
    end if

    BC_Prs => BC_Pressure

    ! Thermo-mechanics bcs are handled differently using the "new" bc stuff
    ! We still use the face and bc surface loops to get the masks for each surface
    ! specified in the bc namelists
    ! Always initialize displacement BCs so that the code does not crash if
    ! displacement BCs are defined but solid mechanics is not active.  A warning will
    ! displayed.
!    if (solid_mechanics) call Initialize_Displacement_BC
    call Initialize_Displacement_BC

    ! Initialize BC flag variables - heat transfer, concentration, fluid flow
    call InitializeBoundaryConditions (.true., .false., .false.)

    ! Initialize velocity BC's; default is a no-BC condition (for internal
    ! faces);  we set mesh boundary faces to free-slip
    do f = 1,nfc
       bit_position = Vel%Face_bit(f)
       Mask1 = .false.
       where (Mesh%Ngbr_Face(f) == 0) Mask1 = .true.
       call SET_FREE_SLIP (Mask1, BC%Flag, bit_position)
    end do
    dirichlet_pressure = .false.

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! NNC, Jan 2014.  New code to handle time-dependent velocity Dirichlet BC
    !! Here we create and store functions that return the Dirichlet velocity
    !! data specified in the BC namelist using either the BC_Value or BC_Table
    !! arrays.  If anything was specified for BC_Table, we use it to define a
    !! time-dependent tabular function; otherwise we define a constant function
    !! using BC_value.
    call bndry_vel%init (ncells, nbc_surfaces)
    do p = 1, nbc_surfaces
      if (BC_Variable(p) == 'velocity' .and. BC_Type(p) == 'dirichlet') then
        if (any(BC_Table(:,:,p) /= NULL_R)) then  ! time-dependent
          !! work out the size of the table read from the input file.
          do n = 1, size(BC_Table,dim=2)
            if (any(BC_Table(:,n,p) == NULL_R)) exit
          end do
          m = n - 1 ! table size
          if (any(BC_Table(:,m+1:,p) /= NULL_R)) &
              call TLS_fatal (' Badly specified BC_Table array for BC namelist '// i_to_c(p))
          do n = 2, m
            if (BC_Table(1,n,p) <= BC_Table(1,n-1,p)) &
                call TLS_fatal ('BC_Table time values not in ascending order')
          end do
          vel_func => new_tabular_vector_func(BC_Table(1,:m,p), BC_Table(2:,:m,p))
        else  ! constant using BC_Value (first three elements)
          vel_func => new_const_vector_func(BC_Value(:3,p))
        end if
        call bndry_vel%set_vel_func (p, vel_func)
      end if
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Write notice
    call TLS_info ('')
    call TLS_info (' Locating cell faces for which BCs are to be applied ...')

    ! Initialize relevant quantities
    Found_Faces = .true.
    Area        = 0.0_r8
    Faces       = 0

    ! Loop over faces, finding those that need a BC applied to them
    FACE_LOOP: do f = 1,nfc

       if (nbc_surfaces == 0) exit FACE_LOOP

       ! Loop over each BC surface, looking for those cell faces that lie on
       ! on the surface within the specified bounding box.
       ! For those faces satisfying these conditions, apply the BC
       ! associated with that surface.
       SURFACE_LOOP: Do p = 1, nbc_surfaces

          Mask1 = .true.
          ! Assign conic function coefficients
          do c = 1,Surfaces_In_This_BC(p)

             select case(Surface_Name(c,p))

             ! BC surface provided in the mesh file.
             case('from mesh file')

                ! If Ngbr_Mesh_Face is not built, do it now then set nssets to
                !   indicate it is done
                if (nssets .eq. 0) then
                   nssets      = SIZE(Mesh_Face_Set,1)

                   ALLOCATE(Ngbr_Mesh_Face_Set(nssets,nfc,nfc,ncells), &
                        Stat=istatus)
                   if (istatus /= 0) call TLS_panic ('MeshBCFaces: allocation error')

                   ! Gather the face set arrays of all the neighboring cells so that we check for
                   ! shared faces that are in the face sets.  This is used by MeshBCFace, but it
                   ! is initialized here to avoid gathering it for every side set
                   do n = 1,nssets
                      do nf = 1,nfc
                         call EE_GATHER(Ngbr_Mesh_Face_Set(n,nf,:,:), Mesh_Face_Set(n,nf,:))
                      end do
                   end do
                endif

                set1 = Mesh_Surface(c,p)
                call MeshBCFaces(f, set1, Mask2, nssets, Ngbr_Mesh_Face_Set )

                Mask1 = Mask1 .and. Mask2

             ! Node set specified instead of a surface; this is handled elsewhere,
             ! only for solid mechanics.
             case('node set')

             ! Find those faces which are boundaries between two specified materials.
             case('material boundary')

                mat1 = Surface_Materials(1,Srfmatl_Index(c,p),p)
                mat2 = Surface_Materials(2,Srfmatl_Index(c,p),p)

                call BoundaryBetweenMaterials (f, mat1, mat2, Mask2)
                Mask1 = Mask1 .and. Mask2

             ! Find those boundary faces having the specified material.
             case('external material boundary')

                mat1 = Surface_Materials(1,Srfmatl_Index(c,p),p)

                call ExternalMaterialBoundary (f, mat1, Mask2)
                Mask1 = Mask1 .and. Mask2

                ! Find those faces which are on this conic function.
             case('conic')

                Coeff_X(1) = Conic_X(c,p)
                Coeff_X(2) = Conic_Y(c,p)
                Coeff_X(3) = Conic_Z(c,p)

                Coeff_XX(1) = Conic_XX(c,p)
                Coeff_XX(2) = Conic_XY(c,p)
                Coeff_XX(3) = Conic_XZ(c,p)
                Coeff_XX(4) = Conic_YY(c,p)
                Coeff_XX(5) = Conic_YZ(c,p)
                Coeff_XX(6) = Conic_ZZ(c,p)

                   ! Mask off those faces whose centroids do not lie on this BC surface.
                Mask2 = BC_SURFACE (f, Coeff_X, Coeff_XX, Conic_Constant(c,p), &
                     Conic_Tolerance(c,p), Conic_Relation(c,p))
                Mask1 = Mask1 .and. Mask2

             end select

          end do

          ! Further mask off those faces not contained within this bounding box.
          do n = 1,ndim
             Mask1 = Mask1 .and. &
                    Cell%Face_Centroid(n,f) >= Bounding_Box(1,n,p) .and. &
                    Cell%Face_Centroid(n,f) <= Bounding_Box(2,n,p)
          end do

          ! Mask off any degenerate faces.
          Mask1 = Mask1 .and. (Mesh%Ngbr_Cell(f) /= DEGENERATE_FACE)

          ! Treat temperature differently than all others.
          select case (BC_Variable(p))

             case ('pressure')

                bit_position = Prs%Face_bit(f)

                select case (BC_Type(p))

                   case ('neumann')

                      call SET_NEUMANN (Mask1, BC%Flag, bit_position)

                   case ('dirichlet')

                      call SET_DIRICHLET (Mask1, BC%Flag, bit_position)
                      where (Mask1) BC_Prs(f,:) = BC_Value(1,p)
                      dirichlet_pressure = .true.

                      ! if dirichlet pressure, then no velocity BC
                      bit_position = Vel%Face_bit(f)
                      call SET_NO_VEL_BC (Mask1, BC%Flag, bit_position)

                end select

             case ('velocity')

                bit_position = Vel%Face_bit(f)

                select case (BC_Type(p))

                   case ('free-slip')

                      call SET_FREE_SLIP (Mask1, BC%Flag, bit_position)

                   case ('no-slip')

                      call SET_DIRICHLET_VEL (Mask1, BC%Flag, bit_position)
                      !! NNC, Jan 2014.  Time-dependent dirichlet velocity.
                      !ORIG: do n = 1, ndim
                      !ORIG:    where (Mask1) BC_Vel(n,f,:) = 0.0_r8
                      !ORIG: end do
                      call bndry_vel%set_no_slip (f, Mask1)

                   case ('dirichlet')

                      call SET_DIRICHLET_VEL (Mask1, BC%Flag, bit_position)
                      !! NNC, Jan 2014.  Time-dependent dirichlet velocity.
                      !ORIG: do n = 1, ndim
                      !ORIG:    where (Mask1) BC_Vel(n,f,:) = BC_Value(n,p)
                      !ORIG: end do
                      call bndry_vel%set_dirichlet (p, f, Mask1)

                   case ('neumann')

                      CALL SET_NEUMANN_VEL (Mask1, BC%Flag, bit_position)

                   end select

                   ! Handle all displacement bc setup with "new" bc stuff in another routine
                case ('displacement')
                   if (TRIM(Surface_Name(1,p)) == 'node set') then
                      Mask1 = .false.
                      if (f == 1) call Node_Set_BC_Init(1,p)
                   else
                      Disp_Mask = .false.
                      Disp_Mask(f,:) = Mask1
                      call Append_to_Displacement_BC(Disp_Mask,p,f)
                      ! Get multiple side set info for gap interfaces
                      if ((TRIM(BC_Type(p)) == 'normal-constraint') .or. &
                           (TRIM(BC_Type(p)) == 'contact') .or. (TRIM(BC_Type(p)) == 'normal-displacement') &
                           .or. (TRIM(BC_Type(p)) == 'free-interface')) &
                           call Interface_Surface_Id(Disp_Mask,p,f,p)
                   end if
                end select

          ! if velocity or pressure, dirichlet or neumann, assign inflow material,
          ! temperature & concentration, if provided

          if ((BC_Variable(p) == 'velocity' .or. BC_Variable(p) == 'pressure') .and. &
              (BC_Type(p) == 'dirichlet' .or. BC_Type(p) == 'neumann')) then

             where (Mask1) BC_Mat(f,:)  = Get_Truchas_Material_Id(Inflow_Material(Inflow_Index(p)))

             if (heat_transport .or. heat_species_transport) then
                where (Mask1) BC_Temp(f,:) = Inflow_Temperature(Inflow_Index(p))
             end if

             !FIXME: need inflow species concentrations.

          end if

          ! Print out how many cells were affected by this BC
          n = PGSLib_Global_COUNT(Mask1)
          if (n > 0) then

             Faces(p) = Faces(p) + n
             Tmp2     = MERGE (Cell%Face_Area(f), 0.0_r8, Mask1)
             Area(p)  = Area(p) + PGSLib_Global_SUM(Tmp2)
             n = PGSLib_Global_COUNT(Mask1 .and. Mesh%Ngbr_cell(f) /= 0)

             ! If the BC is internal, warn the user and
             ! set the internal BC flag (BC%Internal)
             if (n > 0) then

                write (message, 20) n, f, p
20              format (i6,' interior #',i1,' faces are affected by BC namelist ',i2,'!')
                call TLS_warn (message)
                Mask1 = Mask1 .and. Mesh%Ngbr_cell(f) /= 0
                ! don't do this for temperature or displacement...
                if (BC_Variable(p) /= 'temperature' .and. &
                    BC_Variable(p) /= 'displacement') then
                   call SET_INTERNAL_BC (Mask1, BC%Internal, bit_position)
                end if

             end if

          else

             ! No faces matched this surface; set the flag.
             Found_Faces(f,p) = .false.

          end if

       end do SURFACE_LOOP

    end do FACE_LOOP

    ! Make sure each BC surface matches at least one face, otherwise terminate.
    ! If things are OK, then write out diagnostics.
    BC_SURFACE_CHECK: do p = 1,nbc_surfaces

       if ((.not. ANY(Found_Faces(:,p))) .and. (TRIM(Surface_Name(1,p)) /= 'node set')) then

          write(message,'(a,i0,2a)') 'BC_INIT: found no faces for BC namelist ', p, ' named ', trim(BC_Name(p))
          call TLS_fatal (message)

       else

          if (p == 1) call TLS_info ('')
          write (message, 11) p, Area(p), Faces(p)
11           format (4x,'Boundary conditions in BC namelist ',i2,' will be applied on an area of ', &
                  1pe13.5,' (',i6,' faces)')
          call TLS_info (message)

       end if

    end do BC_SURFACE_CHECK

    ! Make sure that faces don't have both a pressure dirichlet condition and some
    ! velocity BC assigned
    VEL_PRS_CONSISTENCY_CHECK: do f = 1,nfc

       Mask1 = DIRICHLET (BC%Flag, Prs%Face_bit(f))
       Mask2 = FREE_SLIP (BC%Flag, Vel%Face_bit(f)) .or.     &
               DIRICHLET_VEL (BC%Flag, Vel%Face_bit(f)) .or. &
               NEUMANN_VEL (BC%Flag, Vel%Face_bit(f))
       if (ANY (Mask1 .and. Mask2)) then
          call TLS_fatal ('BC_INIT: faces found with both velocity and pressure BCs assigned')
       end if

    end do VEL_PRS_CONSISTENCY_CHECK

    ! Setup the "New" BCs.  Temperature always is used.  Pressure only with new pressure BCs.
    ! Always setup pressure BCs, since need them for the discrete operators.
    call Initialize_Pressure_BC(Pressure_BC)
    if (solid_mechanics) call Make_Displacement_BC_Atlases

    !  Should we DEALLOCATE BC_T at this point?   JMS

    ! Write notice
    call TLS_info ('')
    call TLS_info (' Finished BC initialization.')
    DEALLOCATE(Mask1,        &
               Mask2,        &
               Disp_Mask,    &
               Found_Faces,  &
               Coeff_XX,     &
               Coeff_X,      &
               Area,         &
               Faces,        &
               tmp_r2,       &
               Tmp_BC,       &
               Tmp1,         &
               Tmp2)

     if (nssets.gt.0) DEALLOCATE(Ngbr_Mesh_Face_Set)

  END SUBROUTINE BC_INIT

  FUNCTION BC_SURFACE (face, Coeff_X, Coeff_XX, constant, tolerance, relation)
    !=======================================================================
    ! Purpose(s):
    !
    !=======================================================================
    use parameter_module, only: ndim, ncells
    use mesh_module,      only: Cell

    ! Arguments
    integer, intent(IN) :: face
    real(r8), dimension(ndim), intent(IN) :: Coeff_X
    real(r8), dimension(3*(ndim-1)), intent(IN) :: Coeff_XX
    real(r8), intent(IN) :: constant
    real(r8), intent(IN) :: tolerance
    character(*), intent(IN) :: relation

    ! Function Return
    logical, dimension(ncells) :: BC_SURFACE

    ! Local variables
    integer :: m, n, k
    real(r8), dimension(ncells) :: F

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Constant term.
    F = Constant

    ! Linear terms.
    do n = 1,ndim
       F = F + Coeff_X(n)*Cell%Face_Centroid(n,face)
    end do

    ! Quadratic terms.
    m = 0
    do n = 1,ndim
       do k = n, ndim
          m = m + 1
          F = F + Coeff_XX(m)*Cell%Face_Centroid(n,face)*Cell%Face_Centroid(k,face)
       end do
    end do

    ! Function return.
    select case(relation)
       case ('<')
          BC_SURFACE = F      <  tolerance
       case ('=')
          BC_SURFACE = ABS(F) <= tolerance
       case ('>')
          BC_SURFACE = F      > tolerance
    end select

  END FUNCTION BC_SURFACE

  SUBROUTINE MATL_INIT (Hits_Vol)
    !=======================================================================
    ! Purpose(s):
    !
    !   Initialize the Matl structure and the velocity part of
    !   the Zone structure
    !
    !=======================================================================
    use cutoffs_module,       only: cutvof
    use interfaces_module,    only: background_body, Body_Mass, Matnum, nbody,   &
                                    Body_Temp
    use matl_module,          only: Matl, SLOT_INCREASE, SLOT_SET
    use mesh_module,          only: Cell
    use parameter_module,     only: mat_slot, mat_slot_new, maxmat, &
                                    mbody, ncells, nmat
    use pgslib_module,        only: PGSLib_GLOBAL_MAXVAL
    use property_data_module, only: background_material
    use restart_variables,    only: restart
    use property_module,      only: density_material

    ! Arguments
    real(r8), dimension(nbody,ncells), intent(INOUT) :: Hits_Vol

    ! Local Variables
    integer :: stat
    logical :: fatal
    integer :: ib, jb, m, s
    integer,  dimension(maxmat) :: Nbodym
    integer,  dimension(mbody,maxmat) :: Body_Num
    logical,  dimension(ncells) :: Insert_mat
    real(r8), dimension(ncells) :: Vofm, Total, Tmp1
    integer,  dimension(ncells) :: Nmtl

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


    if ( restart ) then
       ! Don't do any work
       return
    end if

    ! Set up arrays that count the number of bodies
    ! of a given material and the body number of each.
    Nbodym   = 0
    Body_Num = 0

    Total = 0.0_r8
    do ib = 1,nbody
       m = matnum(ib)
       Nbodym(m) = Nbodym(m) + 1
       Body_Num(Nbodym(m),m) = ib

       where (Hits_Vol(ib,:)/Cell%Volume < cutvof) Hits_Vol(ib,:) = 0.0_r8
       Total = Total + Hits_Vol(ib,:)
    end do

    ! Convert from total volume to volume fraction
    ! and set to 0 or 1 if within cutvof away
    Total = Total / Cell%Volume
    where (ABS(Total) < cutvof) Total = 0.0_r8
    where (ABS(Total - 1.0_r8) < cutvof) Total = 1.0_r8

1   continue

    ! Zero out Matl
    do s = 1,mat_slot
       call SLOT_SET (Matl, s)
    end do
    Nmtl = 0

    ! Count the required number of material slots
    ! Background materials must reside where 0.0 <= Total < 1.0
    where (Total >= 0.0_r8 .and. Total < 1.0_r8)
       Nmtl = 1
       Matl(1)%Cell%Id = background_material
    end where

    ! Set the Matl Id
    SET_MATL_ID: do m = 1,nmat

       Vofm = 0.0_r8
       do jb = 1,Nbodym(m)
          ib = Body_Num(jb,m)
          Vofm = Vofm + Hits_Vol(ib,:) / Cell%Volume
       end do

       Insert_mat = .false.
       where (Vofm > 0.0_r8) Insert_mat = .true.
       do s = 1,mat_slot
          Insert_mat = Insert_mat .and. .not.(Matl(s)%Cell%Id == m)
       end do
       where (Insert_mat) Nmtl = Nmtl + 1
       do s = 1,mat_slot
          where (Nmtl == s .and. Insert_mat)
             Matl(s)%Cell%Id = m
          end where
       end do

    end do SET_MATL_ID

    mat_slot_new = PGSLib_GLOBAL_MAXVAL(Nmtl)

    ! Print out the number of slots needed, the current
    ! number allocated, and enlarge Matl if necessary
#if DEBUG_MAT_SLOT
    write (message, 2) mat_slot_new, mat_slot
2   format(4x,i2,' material slots are needed: ',i2, ' material slots allocated')
    call TLS_info (message)
#endif

    if (mat_slot_new > mat_slot) then
       call SLOT_INCREASE (Matl, mat_slot, mat_slot_new)
       go to 1     ! Go back up to see if mat_slot is sufficient
    end if

    ! Compute volume fractions and body volumes for each material.
    SET_MATL_VOF: do m = 1,nmat

       Vofm = 0.0_r8

       do jb = 1,Nbodym(m)

          ib = Body_Num(jb,m)

          ! Accumulate the volume fractions
          Vofm = Vofm + Hits_Vol(ib,:)/Cell%Volume

          ! Accumulate the body mass
          Body_mass(ib) = SUM (Hits_Vol(ib,:)*density_material(m))

       end do

       ! Put the material quantities into the appropriate slot
       do s = 1,mat_slot

          where (Matl(s)%Cell%Id == m)
             Matl(s)%Cell%Vof = Vofm
          end where

       end do

    end do SET_MATL_VOF

    ! Total the volume fractions over all materials
    Total = 0.0_r8
    do s = 1,mat_slot
       Total = Total + Matl(s)%Cell%Vof
    end do

    ! Fill background with material "background_material", if not already filled
    Tmp1 = MIN(MAX(1.0_r8 - Total, 0.0_r8), 1.0_r8)

    ! Set up background_material material arrays
    Vofm = 0.0_r8
    do s = 1,mat_slot
       where (Matl(s)%Cell%Id == background_material)
          Vofm = Matl(s)%Cell%Vof
       end where
    end do

    where (Tmp1 > cutvof)
       Vofm = Vofm + Tmp1
       Tmp1  = Tmp1*density_material(matnum(background_body))*Cell%Volume
    end where

    ! Put modified background_material materials back into the appropriate slots
    do s = 1,mat_slot
       where (Matl(s)%Cell%Id == background_material)
          Matl(s)%Cell%Vof = Vofm
       end where
    end do

    !! NNC, Feb 2013.  At this point (original commented-out code below) we
    !! check that computed volume fractions (matl) are valid and then go-ahead
    !! and normalize them anyway!  Why the vofs can't be expected to valid at
    !! this point I don't understand and isn't proper.  Needs to be fixed!

!    ! Total the volume fractions
!    Total = 0.0_r8
!    do s = 1,mat_slot
!       Total = Total + Matl(s)%Cell%Vof
!    end do
!
!    ! Check to see if the Sum(Vof) = 1.0
!    call CHECK_VOF (fatal)
!    call check_vof (stat)
!    if (stat /= 0) call TLS_fatal ('computed volume fractions are invalid')
!
!    ! Renormalize Vofm array in case they are under or over-filled.
!    ! Also initialize Rho_Old
!    Total = 1.0_r8/Total
!    do s = 1,mat_slot
!       Matl(s)%Cell%Vof     = Total*Matl(s)%Cell%Vof
!       where (Matl(s)%Cell%Vof == 0.0_r8) Matl(s)%Cell%Id = 0
!    end do

    !! Preserved from the above code -- is it really necessary?
    do s = 1,mat_slot
       where (Matl(s)%Cell%Vof == 0.0_r8) Matl(s)%Cell%Id = 0
    end do

    !! Check that the computed VoFs are valid.
    call check_vof (stat)
    if (stat /= 0) then
      call TLS_info ('  Computed volume fractions are invalid; attempting to normalize.')
      Total = 0.0_r8
      do s = 1,mat_slot
        Total = Total + Matl(s)%Cell%Vof
      end do
      Total = 1.0_r8/Total
      do s = 1,mat_slot
        Matl(s)%Cell%Vof = Total*Matl(s)%Cell%Vof
      end do
      call check_vof (stat)
      if (stat /= 0) then
        call TLS_fatal ('normalized volume fractions are invalid')
      else
        call TLS_info ('  Normalization was successful.')
      end if
    end if

    ! Compute mesh Velocities from body velocities
    call VELOCITY_INIT (Hits_Vol)

  END SUBROUTINE MATL_INIT

  subroutine ZONE_INIT (Hits_Vol)
    !=======================================================================
    ! Purpose(s):
    !
    !   Initialize all cell-average quantities in the Zone structure
    !
    !=======================================================================
    use cutoffs_module,       only: alittle
    use interfaces_module,    only: Body_Temp, Matnum, nbody
    use matl_module,          only: Matl, GATHER_VOF
    use matl_utilities,       only: update_matl
    use mesh_module,          only: Cell

    use parameter_module,     only: mat_slot, ncells, nmat
    use property_module,      only: ENTHALPY_DENSITY_MATERIAL, DENSITY_MATERIAL
    use zone_module,          only: Zone
    use restart_variables,    only: restart
    use tempGrad_module,      only: Body_Temperature, Body_Temp_Grad
    use physics_module, only: heat_transport, heat_species_transport

    ! Arguments
    real(r8), dimension(nbody,ncells), intent(IN) :: Hits_Vol

    ! Local Variables
    integer :: i, ib, iz, m, s
    real(r8), pointer, dimension(:)   :: enth_sum, scratch
    real(r8), pointer, dimension(:,:) :: bzRho
    real(r8), pointer, dimension(:)   :: mass_sum, vof
    real(r8) :: enth_matl, rho, temp
    real(r8), allocatable :: Volume_Fraction(:,:)
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE (Volume_Fraction(0:nmat,ncells), STAT=s)
    call TLS_fatal_if_any (s /= 0, 'ZONE_INIT: error allocating Volume_Fraction')

    ! scratch(1:ncells) is used for various values of that dimension.
    ! Pointers with meaningful names are associated with it prior to use.
    ALLOCATE (scratch(ncells), STAT=s)
    call TLS_fatal_if_any (s /= 0, 'ZONE_INIT: error allocating scratch')

    ! Gather volume fractions of materials using enth_sum as a temporary
    ! Volume Fraction is need later inside of the T_of_H call
    Volume_Fraction = 0.0_r8
    vof => scratch
    do i = 1,nmat
       call GATHER_VOF (i, vof)
       Volume_Fraction(i,:) = vof
    end do
    NULLIFY(vof)

    RESTARTCHECK: if (restart) then

      if (heat_transport .or. heat_species_transport) then
        !FIXME: compute Zone%Enthalpy or is it read from the restart file?
        Zone%Enthalpy_old = Zone%Enthalpy
      else
        Zone%Enthalpy = 0.0
        Zone%Enthalpy_old = 0.0
      end if

    else

      !! NNC, 4 May 2012.  Here we need to compute cell values for temperature,
      !! density, and enthalpy.  The latter two are average quantities based on
      !! the volume fraction of material phases present in a cell -- averaged
      !! material properties.  At this stage the material phase volume fractions
      !! are available.  The difficulty is in assigning temperatures, which are
      !! needed to compute enthalpies (and potentially densitites), as these are
      !! specified by BODY namelist, along with material phases, whose region
      !! does not necessarily conform to the mesh.  The original code assigned
      !! a preliminary temperature (body contributing the most mass to the cell)
      !! which was used to compute average cell enthalpies.  Using these values
      !! as primary, new consistent temperatures and volume fractions were
      !! solved for, and the updated volume fractions pushed back into MATL.
      !! We could conceivably continue doing this last step, but I anticipate
      !! moving initial field values out of BODY namelists and setting them on
      !! the basis of mesh blocks only, making it unnecessary.  Thus I've opted
      !! to drop the last step, and stick with the 'preliminary' temperatures
      !! and computed enthalpies.

      ! NOTE: Some reasonable value must be assigned to Zone%Temp in this
      !       body-loop to ensure that the return values from  T_of_H() are valid.
      !
      !       T_of_H accesses TI_Initial_Guess() which is just Zone%Temp.
      !       If heat_conduction is off (false), T_of_H() simply returns
      !       TI_Initial_Guess(), which is Zone%Temp as calculated herein.
      !
      !       To meet this need, the zonal temperatures are assigned from
      !       the body that contributes the largest mass to the zone, either
      !       the constant value or from the body temperature gradient, if defined.

      allocate (enth_sum(ncells), STAT=s)
      call TLS_fatal_if_any (s /= 0, 'ZONE_INIT: error allocating enth_sum')

      ! Compute a cell-average density
      Zone%Rho = 0.0_r8
      do s = 1,mat_slot
        do i = 1, ncells
          Zone(i)%Rho = Zone(i)%Rho + Matl(s)%Cell(i)%Vof * DENSITY_MATERIAL(Matl(s)%Cell(i)%Id)
        end do
      end do
      Zone%Rho_Old = Zone%Rho

      ! Initialize target arrays to zero
      mass_sum    => scratch
      mass_sum(:)  = 0.0_r8
      enth_sum(:)  = 0.0_r8
      Zone(:)%Temp = 0.0_r8

      ! Accumulate the partial density values for each body in each zone
      ALLOCATE (bzRho(nbody,ncells), STAT=s)
      call TLS_fatal_if_any (s /= 0, 'ZONE_INIT: error allocating bzRho')
      bzRho(:,:) = 0.0_r8
      !FORALL (ib = 1:nbody, iz = 1:ncells, Hits_Vol(ib,iz) .gt. alittle)
      !  bzRho(ib,iz) = DENSITY_MATERIAL(MatNum(ib))*Hits_Vol(ib,iz)
      !END FORALL
      do iz = 1, ncells
        do ib = 1, nbody
          if (Hits_Vol(ib,iz) > alittle) then
            bzRho(ib,iz) = DENSITY_MATERIAL(MatNum(ib))*Hits_Vol(ib,iz)
          end if
        end do
      end do

      TOTAL_ENTHALPY_LOOP: do ib = 1,nbody
        ! The enthalpy/mass of the cell is calculated as the sum of the
        ! full enthalpy of each material divided by the mass of the cell

        m   = MatNum(ib)

        ! Accumulate total mass by zone
        where (Hits_Vol(ib,:) > alittle)
          mass_sum(:) = mass_sum(:) + DENSITY_MATERIAL(MatNum(ib))*Hits_Vol(ib,:)
        end where

        GRADIENT: if (Body_Temp_Grad(ib)%On) then
          ! Apply a temperature gradient, zone by zone

          do iz = 1, ncells
            if (Hits_Vol(ib,iz) .gt. alittle) then
              temp         = Body_Temperature(ib, Body_Temp(ib), Zone(iz), Cell(iz))
              enth_matl    = ENTHALPY_DENSITY_MATERIAL(m,temp)
              enth_sum(iz) = enth_sum(iz) + enth_matl * Hits_Vol(ib,iz)

              ! If the maximum mass contribution is from this body, assign the temp
              if (count(maxloc(bzRho(:,iz)).eq.ib).gt.0) Zone(iz)%Temp = temp

            end if
          end do

        else GRADIENT
          ! Utilize the constant body temperature
          enth_matl = ENTHALPY_DENSITY_MATERIAL(m,Body_Temp(ib))
          enth_sum(:) = enth_sum(:) + enth_matl*Hits_Vol(ib,:)

          do iz = 1,ncells
            ! If the maximum mass contribution is from this body, assign the temp
            if (count(maxloc(bzRho(:,iz)).eq.ib).gt.0) Zone(iz)%Temp = Body_Temp(ib)
          end do

        end if GRADIENT

      end do TOTAL_ENTHALPY_LOOP

      DEALLOCATE (bzRho)

      ! Complete and store Enthaply
      where(mass_sum(:) .gt. alittle)
        enth_sum(:) = enth_sum(:)/mass_sum(:)
      elsewhere
        enth_sum(:) = 0.0_r8
      end where
      Zone(:)%Enthalpy     = enth_sum
      Zone(:)%Enthalpy_Old = Zone%Enthalpy
      Zone(:)%Temp_Old = Zone(:)%Temp
      NULLIFY(mass_sum)

      call update_matl (Volume_Fraction)

      DEALLOCATE (enth_sum)

    end if RESTARTCHECK

    DEALLOCATE (scratch)
    DEALLOCATE (Volume_Fraction)

  end subroutine ZONE_INIT


  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  !! NNC, Feb 2013.  This is a re-implementation of the original VOF_CHECK
  !! subroutine that actually works in parallel.  Additional improvements
  !! to the diagnostic messages have been made.  The routine has been split
  !! into two parts to facilitate testing.  The auxillary routine, which
  !! does all the work, gets its data via arguments and not module use and
  !! thus can be directly tested.  The original code could not be tested,
  !! and obviously wasn't.

  subroutine check_vof (stat)

    use parameter_module, only: nmat, ncells
    use mesh_module, only: unpermute_mesh_vector
    use matl_utilities, only: matl_get_vof
    use property_module, only: get_user_material_id

    integer, intent(out) :: stat

    integer  :: m, mmap(nmat)
    real(r8) :: vfrac(nmat,ncells)

    call matl_get_vof (vfrac) ! we don't want to deal directly with matl

    do m = 1, nmat
      mmap(m) = get_user_material_id(m)
    end do

    call check_vof_aux (vfrac, mmap, unpermute_mesh_vector, stat)

  end subroutine check_vof

  subroutine check_vof_aux (vfrac, mmap, cmap, stat)

    use parallel_communication, only: global_any, global_minval, global_maxval

    real(r8), intent(in)  :: vfrac(:,:)
    integer,  intent(in)  :: mmap(:)    ! user material IDs
    integer,  intent(in)  :: cmap(:)    ! user cell IDs
    integer,  intent(out) :: stat

    integer :: m, j
    real(r8) :: vfmin, vfmax, vfsum, eps
    character(256) :: message(2)

    !! Local data type to keep track of cells with bad vof data.
    integer, parameter :: MAX_REPORT = 5
    type vector
      integer :: n = 0, array(MAX_REPORT), m = 0
    end type
    type(vector) :: list1, list2

    ASSERT(size(vfrac,1) == size(mmap))
    ASSERT(size(vfrac,2) == size(cmap))

    stat = 0

    do m = 1, size(vfrac,1)

      !! Examine the material volume fractions.
      vfmin = 0.0_r8; vfmax = 1.0_r8
      do j = 1, size(vfrac,2)
        if (vfrac(m,j) < 0.0_r8) then
          call append_vector (list1, cmap(j))
          vfmin = min(vfmin, vfrac(m,j))
        else if (vfrac(m,j) > 1.0_r8) then
          call append_vector (list2, cmap(j))
          vfmax = max(vfmax, vfrac(m,j))
        end if
      end do

      10 format(a,i0,a,/,a,es12.5)

      !! Report cells where the volume fraction is negative.
      if (global_any(size_vector(list1) > 0)) then
        stat = -1
        write(message,10) 'material ', mmap(m), ' volume fraction < 0 in cells: ', &
                          'minimum volume fraction: ', global_minval(vfmin)
        call report_errors (message, list1)
      end if

      !! Report cells where the volume fraction exceeds 1.
      if (global_any(size_vector(list2) > 0)) then
        stat = -1
        write(message,10) 'material ', mmap(m), ' volume fraction > 1 in cells: ', &
                          'maximum volume fraction less 1: ', global_maxval(vfmax) - 1.0_r8
        call report_errors (message, list2)
      end if

      call clear_vector (list1)
      call clear_vector (list2)

    end do

    if (stat /= 0) return ! don't bother with the rest of the checks

    !! Examine the volume fraction sums.
    vfmin = 1.0_r8; vfmax = 1.0_r8
    eps = 2*epsilon(1.0_r8) ! I'm uncertain what tolerance we can require.
    do j = 1, size(vfrac,2)
      vfsum = sum(vfrac(:,j))
      if (vfsum < 1.0_r8 - eps) then
        call append_vector (list1, cmap(j))
        vfmin = min(vfmin, vfsum)
      else if (vfsum > 1.0_r8 + eps) then
        call append_vector (list2, cmap(j))
        vfmax = max(vfmax, vfsum)
      end if
    end do

    20 format(a,/,a,es12.5)

    !! Report cells where the volume fraction sum falls short of 1.
    if (global_any(size_vector(list1) > 0)) then
      stat = -1
      write(message,20) 'volume fraction sum < 1 in cells: ', &
                        'minimum volume fraction sum less 1: ', global_minval(vfmin) - 1.0_r8
      call report_errors (message, list1)
    end if

    !! Report cells where the volume fraction sum exceeds 1.
    if (global_any(size_vector(list2) > 0)) then
      stat = -1
      write(message,20) 'volume fraction sum > 1 in cells: ', &
                        'maximum volume fraction sum less 1: ', global_maxval(vfmax) - 1.0_r8
      call report_errors (message, list2)
    end if

  contains

    pure subroutine clear_vector (this)
      type(vector), intent(inout) :: this
      this%n = 0; this%m = 0
    end subroutine

    pure subroutine append_vector (this, value)
      type(vector), intent(inout) :: this
      integer, intent(in) :: value
      if (this%n < size(this%array)) then
        this%n = this%n + 1
        this%array(this%n) = value
      else
        this%m = this%m + 1
      end if
    end subroutine

    pure integer function size_vector (this)
      type(vector), intent(in) :: this
      size_vector = this%n
    end function

    subroutine report_errors (message, list)
      use parallel_communication, only: is_IOP, allocate_collated_array, collate, global_sum
      character(*), intent(inout) :: message(:)
      type(vector), intent(in) :: list
      integer :: n, m
      integer, pointer :: glist(:)
      n = global_sum(list%n)
      m = global_sum(list%m)
      call allocate_collated_array (glist, n)
      call collate (glist, list%array(:list%n))
      !if (is_IOP) message(1) = trim(message(1)) // list_to_string(glist,m)
      if (is_IOP) call append_list_to_string (glist, m, message(1))
      call TLS_error (message)
    end subroutine

    !! This generates a string of the integers in the list up to a maximum,
    !! with a trailing comment on the number omitted, if any.  This uses
    !! some F2003/F2008 features of allocatable, deferred-length character
    !! scalar variable.  May cause problems with some compilers, but not
    !! NAG 5.3.1 or Intel 13.0.  List should have at least 1 element.

    !function list_to_string (list, omitted) result (string)
    !  use string_utilities, only: i_to_c
    !  integer, intent(in) :: list(:), omitted
    !  character(:), allocatable :: string
    !  integer :: n
    !  string = ' ' // i_to_c(list(1))
    !  do n = 2, min(size(list), MAX_REPORT)
    !    string = string // ' ' // i_to_c(list(n))
    !  end do
    !  if (size(list)+omitted > MAX_REPORT) then
    !    string = string // ' [' // i_to_c(size(list)+omitted-MAX_REPORT) // ' more items omitted]'
    !  end if
    !end function

    !! This version does not require the more advanced F2003/F2003 features.

    subroutine append_list_to_string (list, omitted, string)
      use string_utilities, only: i_to_c
      integer, intent(in) :: list(:), omitted
      character(*), intent(inout) :: string
      integer :: n
      do n = 1, min(size(list), MAX_REPORT)
        string = trim(string) // ' ' // i_to_c(list(n))
      end do
      if (size(list)+omitted > MAX_REPORT) then
        string = trim(string) // ' [' // i_to_c(size(list)+omitted-MAX_REPORT) // ' more items omitted]'
      end if
    end subroutine

  end subroutine check_vof_aux

  SUBROUTINE VELOCITY_INIT (Hits_Vol)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute initial cell-centered velocities
    !
    !=======================================================================
    use cutoffs_module,    only: alittle
    use interfaces_module, only: Body_Vel, nbody, Body_Temp, Matnum
    use parameter_module,  only: ndim, ncells
    use zone_module,       only: Zone
    use property_module,   only: density_material

    ! Arguments
    real(r8), dimension(nbody,ncells), intent(IN) :: Hits_Vol

    ! Local Variables
    integer :: m, n
    real(r8), dimension(ncells)      :: Mass, Massc
    real(r8), dimension(ndim,ncells) :: Momentum
    real(r8) :: rhomat

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Loop over bodies, accumulating cell-centered momentum and total mass
    Mass     = 0.0_r8
    Massc    = 0.0_r8
    Momentum = 0.0_r8
    do m = 1,nbody
       rhomat = density_material(Matnum(m))
       do n = 1,ndim
          Momentum(n,:) = Momentum(n,:) + Hits_Vol(m,:)*rhomat*Body_Vel(n,m)
       end do
       Mass = Mass + Hits_Vol(m,:)*rhomat
    end do

    ! Now compute the initial cell-centered, center-of-mass velocities
    do n = 1,ndim

       ! Get the momentum
       Massc = Momentum(n,:)

       ! Compute the center of mass velocity
       where (ABS(Massc) <= alittle)
          Zone%Vc(n) = 0.0_r8
       elsewhere
          Zone%Vc(n) = Massc/(Mass + alittle)
       end where

    end do

  END SUBROUTINE VELOCITY_INIT


  SUBROUTINE BOUNDARYBETWEENMATERIALS (face, material_1, material_2, Mask)
    !=======================================================================
    ! Purpose:
    !
    ! Find all faces that have a cell made of material_1 on one side
    ! and material_2 on the other side.
    !
    !=======================================================================
    use cutoffs_module,    only: cutvof
    use gs_module,         only: EE_GATHER
    use matl_module,       only: GATHER_VOF
    use parameter_module,  only: ncells, nfc

    ! Arguments
    integer, intent(IN) :: face
    integer, intent(IN) :: material_1
    integer, intent(IN) :: material_2
    logical, dimension(ncells), intent(INOUT) :: Mask

    ! Local Variables
    real(r8) :: threshold
    integer  :: c, status
    real(r8), allocatable, dimension(:)   :: Vof_1, Vof_2
    real(r8), allocatable, dimension(:,:) :: Vof_Tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set the threshold.
    threshold = cutvof

    ! Get some working space.
    ALLOCATE(Vof_1(ncells),STAT=status)
    if (status /= 0) call TLS_panic ('BoundaryBetweenMaterials: allocation failed: vof_1')
    ALLOCATE(Vof_2(ncells),STAT=status)
    if (status /= 0) call TLS_panic ('BoundaryBetweenMaterials: allocation failed: vof_2')
    ALLOCATE(Vof_Tmp(nfc,ncells),STAT=status)
    if (status /= 0) call TLS_panic ('BoundaryBetweenMaterials: allocation failed: vof_tmp')

    ! start with a null boundary
    Mask = .false.

    ! get vof of material 1 in each cell
    call GATHER_VOF (material_1, Vof_1)

    ! get vof of material 2 in each cell
    call GATHER_VOF (material_2, Vof_2)

    ! gather vof of material 2 across cell faces
    call EE_GATHER (Vof_Tmp, Vof_2)

    ! look for cells made of (mostly) material 1
    do c = 1, ncells
       if (Vof_1(c) > threshold .and. Vof_Tmp(face,c) > threshold) Mask(c) = .true.
    end do

    ! gather vof of material 1 across cell faces
    call EE_GATHER (Vof_Tmp, Vof_1)

    ! look for cells made of (mostly) material 2
    do c = 1, ncells
       if (Vof_2(c) > threshold .and. Vof_Tmp(face,c) > threshold) Mask(c) = .true.
    end do

    ! release working space
    DEALLOCATE(Vof_Tmp)
    DEALLOCATE(Vof_2)
    DEALLOCATE(Vof_1)

  END SUBROUTINE BOUNDARYBETWEENMATERIALS

  SUBROUTINE EXTERNALMATERIALBOUNDARY (face, material, Mask)
    !=======================================================================
    ! Purpose:
    !
    ! Find all faces that have a cell made of material_1 on one side
    ! and material_2 on the other side.
    !
    !=======================================================================
    use cutoffs_module,    only: cutvof
    use matl_module,       only: GATHER_VOF
    use mesh_module,       only: Mesh
    use parameter_module,  only: ncells

    ! Arguments
    integer, intent(IN) :: face
    integer, intent(IN) :: material
    logical, dimension(ncells), intent(INOUT) :: Mask

    ! Local Variables
    real(r8) :: threshold
    integer  :: c, status
    real(r8), allocatable, dimension(:) :: Vof

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Set the threshold.
    threshold = cutvof

    ! Get some working space.
    ALLOCATE(Vof(ncells),STAT=status)
    if (status /= 0) call TLS_panic ('ExternalMaterialBoundary: allocation failed: vof')

    ! Start with a null boundary.
    Mask = .false.

    ! Get the material volume fraction (vof) in each cell.
    call GATHER_VOF (material, Vof)

    ! Look for cell boundary faces having the specified material.
    do c = 1, ncells
       if (Vof(c) > threshold .and. Mesh(c)%Ngbr_Cell(face) == 0) Mask(c) = .true.
    end do

    ! Release working space.
    DEALLOCATE(Vof)

  END SUBROUTINE EXTERNALMATERIALBOUNDARY

  SUBROUTINE MESHBCFACES(face, set1, Mask, nssets, Ngbr_Mesh_Face_Set )
    !=======================================================================
    ! Purpose:
    !
    ! Find all faces and cells that are in the side sets specified by set1.
    ! At least some mesh generators, e.g. Cubit only list one
    ! cell/face pair for each boundary location.  We need cells and face
    ! ids on both sides of internal boundaries.  Do not include faces on gap
    ! elements.
    !
    !=======================================================================
    use bc_data_module,    only: Mesh_Face_Set
    use parameter_module,  only: ncells, nfc
    use mesh_module,       only: Mesh

    ! Arguments
    integer, intent(IN) :: face, set1, nssets
    logical, dimension(ncells), intent(INOUT) :: Mask
    integer, dimension(nssets,nfc,nfc,ncells), intent(IN) :: Ngbr_Mesh_Face_Set

    ! Local Variables
    integer :: n, nf, icell

    ! Start with a null boundary.
    Mask = .false.

    !nssets = SIZE(Mesh_Face_Set,1)

    ! Check each neighbor of each cell to see if it is in the face set list.
    do n = 1,nssets
       CELL_LOOP: do icell = 1,ncells
          if (Mesh_Face_Set(n,face,icell) == set1) then
             Mask(icell) = .true.
          else
             if (Mesh(icell)%Ngbr_Face(face) > 0) then
                nf = Mesh(icell)%Ngbr_Face(face)
             else
                cycle CELL_LOOP
             end if
             if (Ngbr_Mesh_Face_Set(n,nf,face,icell) == set1) then
                Mask(icell) = .true.
             end if
          end if
       end do CELL_LOOP
    end do

  END SUBROUTINE MESHBCFACES

 !!
 !! COMPUTE_CELL_ENTHALPY
 !!
 !! Computes the cell average enthalpy density given the temperature and the
 !! the volume fractions of the material phases.  This routine is transitional.
 !! The (old) material phases are assigned to the mesh, but the enthalpy is
 !! defined in terms of (new) material systems.  This requires mapping both
 !! volume fractions and materials to the new scheme to evaluate the enthalpy.
 !! In addition, material phases could have been assigned to cells that were
 !! inconsistent with the temperature.  This routine also reapportions the
 !! material phase fractions so they are consistent with the temperature.
 !!

  subroutine compute_cell_enthalpy (T, vof, H)

    use parameter_module, only: nmat
    use phase_property_table
    use material_table
    use material_system
    use material_property
    use material_interop, only: void_material_index, material_to_system, phase_to_material

    real(r8), intent(in)    :: T(:)     ! cell temperature
    real(r8), intent(inout) :: vof(:,:) ! Truchas material phase volume fractions
    real(r8), intent(out)   :: H(:)     ! enthalpy density

    integer :: ncell, m, i, j, k, stat
    integer, allocatable :: matids(:)
    real(r8), allocatable :: vfrac(:,:), pfrac(:)
    real(r8) :: state(1), value
    type(mat_system), pointer :: ms
    type(mat_prop) :: mp
    integer, pointer :: phase_id(:)
    character(128) :: errmsg

    ncell = size(vof,dim=2)

    ASSERT(size(vof,dim=1) == nmat)
    ASSERT(size(T) == ncell)
    ASSERT(size(H) == ncell)

    !! Map material volume fractions (VOF) to material system volume fractions (VFRAC).
    allocate(matids(mt_num_material()), vfrac(ncell,mt_num_material()))
    call mt_get_material_ids (matids)
    vfrac = 0.0_r8
    do m = 1, nmat  ! loop over Truchas material numbers
      if (m == void_material_index) cycle
      do i = size(matids), 1, -1 ! find the destination index
        if (matids(i) == material_to_system(m)) exit
      end do
      INSIST(i > 0)
      vfrac(:,i) = vfrac(:,i) + vof(m,:)
    end do

    ASSERT(ppt_has_property('enthalpy density'))

    H = 0.0_r8

    do i = 1, size(matids) ! loop over material systems
      !! Create the enthalpy density property MP for this material system.
      call mp_create (mp, matids(i), ppt_property_id('enthalpy density'), stat, errmsg)
      if (stat /= 0) call TLS_fatal ('COMPUTE_CELL_ENTHALPY: ' // trim(errmsg))
      !! Accumulate the enthalpy density due to this material system.
      do j = 1, ncell
        if (vfrac(j,i) == 0.0_r8) cycle ! no contribution
        state(1) = T(j)
        call mp_eval (mp, state, value)
        H(j) = H(j) + vfrac(j,i)*value
      end do
      !! Compute the correct phase volume fractions and push back into VOF.
      !! Only relevant for multi-phase material systems.
      ms => mt_get_material(matids(i))
      ASSERT(associated(ms))
      call ms_get_phase_id (ms, phase_id)
      if (size(phase_id) > 1) then  ! multi-phase material system
        allocate(pfrac(size(phase_id)))
        do j = 1, ncell
          if (vfrac(j,i) == 0.0_r8) cycle ! nothing to do
          state(1) = T(j)
          call ms_phase_mixture (ms, state, pfrac)
          do k = 1, size(phase_id)
            m = phase_to_material(phase_id(k))
            vof(m,j) = pfrac(k) * vfrac(j,i)
          end do
        end do
        deallocate(pfrac)
      end if
      deallocate(phase_id)
      call destroy (mp)
    end do

  end subroutine compute_cell_enthalpy

END MODULE INIT_MODULE
