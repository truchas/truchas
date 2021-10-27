!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use material_model_driver,  only: matl_model
  implicit none
  private

  ! Public Procedures
  public :: INITIAL

  ! For testing:
  public :: check_vof_aux

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  subroutine INITIAL (t, dt)
    !=======================================================================
    ! Purpose(s):
    !
    !   Initialize the Matl, Zone, and Alloy (if necessary) derived types.
    !   Initialize all BC data structures and information.
    !   Initialize all probe objects.
    !
    !=======================================================================
    use overwrite_module,       only: OVERWRITE_BC, OVERWRITE_MATL,       &
        OVERWRITE_VEL, OVERWRITE_ZONE
    use parameter_module,       only: nmat
    use legacy_mesh_api,        only: ncells
    use restart_variables,      only: restart
    use restart_driver,         only: restart_matlzone, restart_solid_mechanics, restart_species, restart_ustruc
    use zone_module,            only: Zone
    use solid_mechanics_module, only: SOLID_MECH_INIT
    use diffusion_solver_data,  only: ds_enabled, num_species, &
        ds_sys_type, DS_SPEC_SYS, DS_TEMP_SYS, DS_TEMP_SPEC_SYS
    use diffusion_solver,       only: ds_init, ds_set_initial_state, ds_get_face_temp_view, &
                                      ds_get_temp
    use probes_driver,          only: probes_init
    use ustruc_driver,          only: ustruc_driver_init
    use flow_driver, only: flow_driver_init, flow_enabled, flow_driver_set_initial_state
    use vtrack_driver, only: vtrack_driver_init, vtrack_enabled
    use solid_mechanics_driver, only: solid_mechanics_init
    use physics_module,         only: heat_transport, flow, solid_mechanics, legacy_solid_mechanics
    use toolhead_driver,        only: toolhead_init
    use mesh_manager,           only: unstr_mesh_ptr
    use body_namelist,          only: bodies_params
    use compute_body_volumes_proc
    use material_model_driver,  only: init_material_model
    use unstr_mesh_type
    use output_control, only: output_init

    real(r8), intent(in) :: t, dt

    ! Local Variables
    integer :: m, stat
    real(r8) :: density
    logical :: found
    type(unstr_mesh), pointer :: mesh
    real(r8), allocatable :: phi(:,:), vel_fn(:), hits_vol(:,:)
    real(r8), pointer :: temperature_fc(:) => null()
    character(200) :: errmsg
    character(:), allocatable :: errmsg2
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    !call init_material_model

    mesh => unstr_mesh_ptr('MAIN')

    call property_init

    ! Set up Thermodynamic Reference Temperatures and Enthalpies
    !call THERMO_REFERENCES

    ! vof_initialize returns the volume of each body in each cell
    ! hits_vol is historically used to hold this information
    ! I think this is converted to volume fractions in matl_init.
    ! That's hopelessly confusing and will be corrected later.
    call TLS_info ('')
    call TLS_info ('Computing initial volume fractions ... ')
    call compute_body_volumes (mesh, 3, bodies_params, hits_vol, stat, errmsg2)
    if (stat /= 0) call TLS_fatal(errmsg2)

    ! Either read Zone and Matl from a restart file or initialize them
    if (restart) then

      call restart_matlzone(vel_fn)
      call restart_solid_mechanics

    else

      ! Initialize Matl%Cell and Zone%Vc.
      call MATL_INIT (Hits_Vol)

    end if

    ! Initialize Zone%Rho and Zone%Temp
    call ZONE_INIT (mesh, Hits_Vol)

    if (legacy_solid_mechanics) call BC_INIT()

    if (flow) call flow_driver_init

    ! Allow arbitrary overwriting of the Matl, Zone, and BC types.
    call OVERWRITE_MATL ()
    call OVERWRITE_ZONE ()
    call OVERWRITE_VEL  ()
    call OVERWRITE_BC   ()

    ! Calculate initial stress-strain field and displacements if solid mechanics is active
    if (solid_mechanics) then
      call solid_mechanics_init
    else if (legacy_solid_mechanics) then
      Call SOLID_MECH_INIT
    end if

    call toolhead_init(t)
    call output_init(t)

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
      call ds_init(t)
      select case (ds_sys_type)
      case (DS_SPEC_SYS)
        call ds_set_initial_state (t, dt, conc=phi)
        deallocate(phi)
      case (DS_TEMP_SYS)
        call ds_set_initial_state (t, dt, temp=zone%temp)
        call ds_get_temp (zone%temp)  ! possibly adjusted on void cells
      case (DS_TEMP_SPEC_SYS)
        call ds_set_initial_state (t, dt, temp=zone%temp, conc=phi)
        call ds_get_temp (zone%temp)  ! possibly adjusted on void cells
        deallocate(phi)
      end select
    end if

    ! Initialize the flow solver.
    if (flow) then
      if (ds_enabled) call ds_get_face_temp_view(temperature_fc)
      if (allocated(vel_fn)) then
        call flow_driver_set_initial_state(t, dt, temperature_fc, vel_fn)
      else
        call flow_driver_set_initial_state(t, dt, temperature_fc)
      end if
    end if

    ! Initialize the microstructure modeling driver (if enabled).
    call ustruc_driver_init (t)
    if (restart) call restart_ustruc

    ! Initialize probes.
    call probes_init

  END SUBROUTINE INITIAL

  subroutine property_init

    use physics_module
    use material_utilities

    integer :: stat
    character(:), allocatable :: errmsg

    call required_property_check(matl_model, 'density', stat, errmsg)
    if (stat /= 0) call TLS_fatal ('PROPERTY_INIT: ' // errmsg)

    call constant_property_check(matl_model, 'density', stat, errmsg)
    if (stat /= 0) call TLS_fatal('non-constant density specified for materials: ' // errmsg)

    if (heat_transport) then
      call add_enthalpy_prop(matl_model, stat, errmsg)
      if (stat /= 0) call TLS_fatal('PROPERTY_INIT: ' // errmsg)
    else ! unused but expected for initialization and output
      !FIXME: If the physics don't need this we should not reference it.
      call define_property_default(matl_model, 'specific-enthalpy', default=0.0_r8)
      call add_enthalpy_prop(matl_model, stat, errmsg)
      if (stat /= 0) call TLS_fatal('PROPERTY_INIT: ' // errmsg)
    end if

  end subroutine property_init

!  subroutine flow_property_init
!
!    use material_utilities
!    use fluid_data_module, only: boussinesq_approximation
!    use viscous_data_module, only: inviscid
!
!    integer :: stat
!    character(:), allocatable :: errmsg
!
!    if (boussinesq_approximation) then
!      call define_fluid_property_default(matl_model, 'density-delta', default=0.0_r8)
!    end if
!
!    if (.not.inviscid) then
!      call required_fluid_property_check(matl_model, 'viscosity', stat, errmsg)
!      if (stat /= 0) call TLS_fatal(errmsg)
!    end if
!
!  end subroutine flow_property_init

  subroutine init_conc (hits_vol, phi)

    use legacy_mesh_api, only: ncells, mesh, GAP_ELEMENT_1
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
    use bc_data_module,         only: BC_Type, BC_Name, BC_Variable, nbc_surfaces,          &
                                      Conic_XX, Conic_YY, Conic_ZZ, Conic_XY, Conic_XZ,     &
                                      Conic_YZ, Conic_X, Conic_Y, Conic_Z, Conic_Constant,  &
                                      Conic_Tolerance, Bounding_Box, Surface_Name,          &
                                      Conic_Relation, Surfaces_In_This_BC, Mesh_Surface
    use bc_displacement_init,   only: Initialize_Displacement_BC, Node_Set_BC_Init,         &
                                      append_to_displacement_bc,                            &
                                      Make_Displacement_BC_Atlases, Interface_Surface_Id
    use legacy_mesh_api,        only: ncells, ndim, nfc, nvc, EE_GATHER
    use legacy_mesh_api,        only: Cell, Mesh, DEGENERATE_FACE, mesh_face_set
    use parallel_communication, only: global_count, global_sum
    use physics_module,         only: legacy_solid_mechanics

    ! Local Variables
    integer                              :: f, n, nf, p, &
                                                             c, set1, nssets
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

    ! Thermo-mechanics bcs are handled differently using the "new" bc stuff
    ! We still use the face and bc surface loops to get the masks for each surface
    ! specified in the bc namelists
    ! NB: Code will crash if displacement BCs are defined but solid mechanics is not active.
    if (legacy_solid_mechanics) call Initialize_Displacement_BC

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

          ! Print out how many cells were affected by this BC
          n = global_count(Mask1)
          if (n > 0) then

             Faces(p) = Faces(p) + n
             Tmp2     = MERGE (Cell%Face_Area(f), 0.0_r8, Mask1)
             Area(p)  = Area(p) + global_sum(Tmp2)
             n = global_count(Mask1 .and. Mesh%Ngbr_cell(f) /= 0)

             ! If the BC is internal, warn the user and
             ! set the internal BC flag (BC%Internal)
             if (n > 0) then

                write (message, 20) n, f, p
20              format (i6,' interior #',i1,' faces are affected by BC namelist ',i2,'!')
                !call TLS_warn (message)
                Mask1 = Mask1 .and. Mesh%Ngbr_cell(f) /= 0

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

       end if

    end do BC_SURFACE_CHECK

    if (legacy_solid_mechanics) call Make_Displacement_BC_Atlases

  END SUBROUTINE BC_INIT

  FUNCTION BC_SURFACE (face, Coeff_X, Coeff_XX, constant, tolerance, relation)
    !=======================================================================
    ! Purpose(s):
    !
    !=======================================================================
    use legacy_mesh_api, only: ndim, ncells, Cell

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
    use interfaces_module,    only: matnum, nbody
    use matl_module,          only: matl, slot_increase, slot_set
    use legacy_mesh_api,      only: ncells, Cell, ncells_real
    use parameter_module,     only: mat_slot, mat_slot_new, nmat
    use parallel_communication, only: global_maxval
    use restart_variables,    only: restart

    real(r8), intent(inout) :: Hits_Vol(:,:)

    integer :: b, m, s, j, stat
    integer :: mcount(ncells)
    logical :: bm_mask(nbody,nmat)
    real(r8) :: vofm, total(ncells_real)

    if (restart) return ! nothing to do

    !! Body material mask: BM_MASK(b,m) true if body b is material m.
    bm_mask = .false.
    do b = 1, nbody
      bm_mask(b,matnum(b)) = .true.
    end do

    !! Count the number of materials in each cell and size MATL accordingly
    mcount = 0
    do m = 1, nmat
      do j = 1, ncells
        if (any(hits_vol(:,j) > 0 .and. bm_mask(:,m))) mcount(j) = mcount(j) + 1
      end do
    end do
    mat_slot_new = global_maxval(mcount)
    if (mat_slot_new > mat_slot) call slot_increase(matl, mat_slot, mat_slot_new)

    !! Zero out MATL
    do s = 1, mat_slot
      call slot_set(matl, s)
    end do

    !! Initialize MATL
    mcount = 0
    do m = 1, nmat
      do j = 1, ncells
        vofm = sum(hits_vol(:,j), mask=bm_mask(:,m)) / cell(j)%volume
        if (vofm > 0) then
          mcount(j) = mcount(j) + 1
          s = mcount(j)
          matl(s)%cell(j)%id = m
          matl(s)%cell(j)%vof = vofm
        end if
      end do
    end do

    !! We are guaranteed that HITS_VOL sums to the cell volume within roundoff.
    !! TODO? tweak the vofs to drop small values and adjust to sum exactly to 1

    ! Compute mesh Velocities from body velocities
    call VELOCITY_INIT (Hits_Vol)

  END SUBROUTINE MATL_INIT

  subroutine ZONE_INIT (mesh, Hits_Vol)
    !=======================================================================
    ! Purpose(s):
    !
    !   Initialize all cell-average quantities in the Zone structure
    !
    !=======================================================================
    use interfaces_module,    only: Body_Temp, Matnum, nbody
    use matl_utilities,       only: matl_get_cell_vof
    use zone_module,          only: Zone
    use restart_variables,    only: restart
    use physics_module, only: heat_transport
    use avg_phase_prop_type
    use unstr_mesh_type

    type(unstr_mesh), intent(inout) :: mesh
    real(r8), intent(in) :: hits_vol(:,:)

    integer  :: i, j, pid(matl_model%nphase)
    real(r8) :: vfrac(matl_model%nphase)
    type(avg_phase_prop), allocatable :: func
    character(:), allocatable :: errmsg

    INSIST(size(hits_vol,dim=1) == nbody)

    RESTARTCHECK: if (restart) then

      if (heat_transport) then
        !FIXME: compute Zone%Enthalpy or is it read from the restart file?
        Zone%Enthalpy_old = Zone%Enthalpy
      else
        Zone%Enthalpy = 0.0
        Zone%Enthalpy_old = 0.0
     end if

    else

      !! ZONE%TEMP
      call compute_equil_cell_temp(mesh, matnum(1:nbody), body_temp(1:nbody), hits_vol, zone%temp)
      zone%temp_old = zone%temp

      !! ZONE%RHO
      pid = [(i, i = 1, matl_model%nphase)]
      call matl_model%alloc_avg_phase_prop('density', pid, func, errmsg)
      if (.not.allocated(func)) call TLS_fatal(errmsg)
      do j = 1, mesh%ncell_onP
        call matl_get_cell_vof(j, vfrac)
        call func%compute_value(vfrac, [zone(j)%temp], zone(j)%rho)
      end do
      zone(mesh%ncell_onP+1:)%rho = 0 ! gap elements
      zone%rho_old = zone%rho

      !! ZONE%ENTHALPY
      call matl_model%alloc_avg_phase_prop('enthalpy', pid, func, errmsg)
      if (.not.allocated(func)) call TLS_fatal(errmsg)
      do j = 1, mesh%ncell_onP
        call matl_get_cell_vof(j, vfrac)
        call func%compute_value(vfrac, [zone(j)%temp], zone(j)%enthalpy)
      end do
      zone(mesh%ncell_onP+1:)%enthalpy = 0 ! gap elements
      zone%enthalpy_old = zone%enthalpy

    end if RESTARTCHECK

  end subroutine ZONE_INIT


  subroutine compute_equil_cell_temp(mesh, body_pid, body_temp, body_vol, tcell)

    use scalar_func_containers, only: scalar_func_box
    use equil_temp_type
    use unstr_mesh_type

    type(unstr_mesh), intent(inout) :: mesh
    integer, intent(in) :: body_pid(:)
    type(scalar_func_box), intent(in) :: body_temp(:)
    real(r8), intent(in) :: body_vol(:,:)
    real(r8), intent(out) :: tcell(:)

    integer :: i, j
    integer, allocatable :: nvbi(:)
    type(equil_temp), allocatable :: equil_temp_func
    character(:), allocatable :: errmsg
    real(r8), allocatable :: temps(:)

    ASSERT(size(body_pid) == size(body_temp))
    ASSERT(size(body_pid) == size(body_vol,dim=1))
    ASSERT(size(tcell) == size(body_vol,dim=2))
    ASSERT(size(tcell) >= mesh%ncell_onP)

    !! Default and dummy values
    tcell(:mesh%ncell_onP) = 0   ! default void temperature for true cells
    tcell(mesh%ncell_onP+1:) = 0 ! dummies for possible gap element cells

    !! Array of non-void body indices
    nvbi = pack([(i,i=1,size(body_pid))], body_pid /= matl_model%void_index)
    if (size(nvbi) == 0) return  ! all bodies are void

    call matl_model%alloc_equil_temp(body_pid(nvbi), equil_temp_func, errmsg)
    if (.not.allocated(equil_temp_func)) call TLS_fatal(errmsg)

    allocate(temps(size(nvbi)))
    call mesh%init_cell_centroid
    do j = 1, mesh%ncell_onP
      associate (w => body_vol(nvbi,j))
        if (count(w > 0) == 0) cycle ! entirely void cell
        do i = 1, size(w)
          if (w(i) > 0) temps(i) = body_temp(nvbi(i))%eval(mesh%cell_centroid(:,j))
        end do
        call equil_temp_func%compute(w, temps, tcell(j))
      end associate
    end do

  end subroutine compute_equil_cell_temp

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  !! NNC, Feb 2013.  This is a re-implementation of the original VOF_CHECK
  !! subroutine that actually works in parallel.  Additional improvements
  !! to the diagnostic messages have been made.  The routine has been split
  !! into two parts to facilitate testing.  The auxillary routine, which
  !! does all the work, gets its data via arguments and not module use and
  !! thus can be directly tested.  The original code could not be tested,
  !! and obviously wasn't.

  subroutine check_vof (stat)

    use parameter_module, only: nmat
    use legacy_mesh_api, only: ncells, unpermute_mesh_vector, ncells_real
    use matl_utilities, only: matl_get_vof

    integer, intent(out) :: stat

    integer :: n
    real(r8) :: vfrac(nmat,ncells)

    call matl_get_vof (vfrac) ! we don't want to deal directly with matl

    call check_vof_aux (vfrac(:,:ncells_real), unpermute_mesh_vector(:ncells_real), stat)

  end subroutine check_vof

  subroutine check_vof_aux (vfrac, cmap, stat)

    use parallel_communication, only: global_any, global_minval, global_maxval

    real(r8), intent(in)  :: vfrac(:,:)
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

      10 format(3a,/,a,es12.5)

      !! Report cells where the volume fraction is negative.
      if (global_any(size_vector(list1) > 0)) then
        stat = -1
        write(message,10) 'material "', matl_model%phase_name(m), &
            '" volume fraction < 0 in cells: ', &
            'minimum volume fraction: ', global_minval(vfmin)
        call report_errors (message, list1)
      end if

      !! Report cells where the volume fraction exceeds 1.
      if (global_any(size_vector(list2) > 0)) then
        stat = -1
        write(message,10) 'material "', matl_model%phase_name(m), &
            '" volume fraction > 1 in cells: ', &
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
      use parallel_communication, only: is_IOP, gather, global_sum
      character(*), intent(inout) :: message(:)
      type(vector), intent(in) :: list
      integer :: n, m
      integer, pointer :: glist(:)
      n = global_sum(list%n)
      m = global_sum(list%m)
      allocate(glist(merge(n,0,is_iop)))
      call gather (list%array(:list%n), glist)
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
    use legacy_mesh_api,   only: ndim, ncells
    use zone_module,       only: Zone

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
       if (Matnum(m) == matl_model%void_index) cycle
       rhomat = matl_model%const_phase_prop(Matnum(m), 'density')
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

       Zone%Vc_old(n) = Zone%Vc(n)

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
    use legacy_mesh_api,   only: ncells, nfc, EE_GATHER
    use matl_module,       only: GATHER_VOF

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
    use legacy_mesh_api,   only: ncells, Mesh

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
    use legacy_mesh_api, only: ncells, nfc, Mesh, Mesh_Face_Set

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
    use avg_matl_prop_type
    use material_class

    real(r8), intent(in), contiguous, target :: T(:)     ! cell temperature
    real(r8), intent(inout) :: vof(:,:) ! Truchas material phase volume fractions
    real(r8), intent(out)   :: H(:)     ! enthalpy density

    integer :: ncell, nm, m, p1, p2, j
    real(r8), allocatable :: vfrac(:,:)
    real(r8), pointer :: state(:,:)
    type(avg_matl_prop), allocatable :: avg_H
    character(:), allocatable :: errmsg

    ncell = size(vof,dim=2)

    ASSERT(size(vof,dim=1) == nmat)
    ASSERT(size(T) == ncell)
    ASSERT(size(H) == ncell)

    nm = matl_model%nmatl_real
    allocate(vfrac(ncell,nm))

    !! Generate the material volume fraction array VFRAC from VOF
    do m = 1, nm
      call matl_model%get_matl_phase_index_range(m, p1, p2)
      vfrac(:,m) = sum(vof(p1:p2,:),dim=1)
    end do

    !! Compute the volume fraction weighted average material enthalpy density.
    call matl_model%alloc_avg_matl_prop('enthalpy', avg_H, errmsg)
    ASSERT(allocated(avg_H))
    state(1:1,1:ncell) => T
    call avg_H%compute_value(vfrac, state, H)

    !! Compute the correct phase volume fractions for multi-phase materials
    !! and push them back into VOF.
    do m = 1, nm
      call matl_model%get_matl_phase_index_range(m, p1, p2)
      if (p1 == p2) cycle
      block
        real(r8) :: beta(p2-p1+1)
        do j = 1, ncell
          if (vfrac(j,m) > 0) then
            call matl_model%get_matl_phase_frac(m, T(j), beta)
            vof(p1:p2,j) = vfrac(j,m)*beta
          end if
        end do
      end block
    end do

  end subroutine compute_cell_enthalpy

END MODULE INIT_MODULE
