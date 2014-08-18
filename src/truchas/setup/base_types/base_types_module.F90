MODULE BASE_TYPES_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures to allocate, deallocate, and default base types
  !
  ! Public Interface(s):
  !
  !   * call BASE_TYPES_ALLOCATE ()
  !
  !     Allocate the base types.
  !
  !   * call BASE_TYPES_DEALLOCATE ()
  !
  !     Deallocate the base types.
  !
  !   * call BASE_TYPES_DEFAULT ()
  !
  !     Default the base types, which is defaulted in SLOT_INCREASE.
  !
  !   * call MESH_VERTEX_ALLOCATE (Mesh, Vertex, ncells_input, nnodes_input)
  !
  !     Allocate the mesh types.
  !
  !   * call MESH_VERTEX_DEALLOCATE (Mesh, Vertex)
  !
  !     Deallocate the mesh types.
  !
  !   * call MESH_REALLOCATE (MeshPermute, VertexPermute)
  !
  !     Reallocate the mesh types to a different size.
  !
  ! Contains: BASE_TYPES_ALLOCATE
  !           BASE_TYPES_DEALLOCATE
  !           BASE_TYPES_DEFAULT
  !           MESH_VERTEX_ALLOCATE
  !           MESH_VERTEX_DEALLOCATE
  !           MESH_REALLOCATE
  !           MESH_VERTEX_PERMUTE
  !           MESH_PERMUTE
  !           MESH_RENUMBER
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !            Robert Ferrell (CPCA, Ltd., ferrell@cpca.com)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  ! Public Procedures
  public :: BASE_TYPES_ALLOCATE, BASE_TYPES_DEALLOCATE, BASE_TYPES_DEFAULT, &
            MESH_VERTEX_ALLOCATE, MESH_VERTEX_DEALLOCATE, &
            PERMUTE_MESH, PERMUTE_VERTEX, RENUMBER_CELLS_VERTICES, &
            ANNOUNCE_MESH_SIZES

  !! These procedures are incomplete, in that they don't handle all the
  !! components of the derived type as they should.  They aren't used,
  !! so make them private to prevent their future use. (NNC, 1/13/2005)
  private :: MESH_REALLOCATE, MESH_VERTEX_PERMUTE, MESH_RENUMBER

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE BASE_TYPES_ALLOCATE ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Allocate the Cell, Matl, and Zone base types and Fluxing_Velocity
    !
    !=======================================================================
    use bc_module,           only: BC
    use matl_module,         only: SLOT_INCREASE, Matl
    use mesh_module,         only: Cell
    use parameter_module,    only: nfc, ncells, mat_slot, mat_slot_new, nmat, nprobes
    use zone_module,            only: Zone
    use fluid_data_module,       only: Fluxing_Velocity
    use solid_mechanics_module, only: SOLID_MECHANICS_ALLOCATE
    use solid_mechanics_input,  only: solid_mechanics
    use turbulence_module,      only: TURBULENCE_ALLOCATE

    use probe_module,           only: probes
    use EM_data_proxy,          only: EM_is_on

    ! Arguments

    ! Local Variables
    integer :: memstat
    integer :: ems, sms, smv, smt
    integer :: i
    integer :: nprobevars, nprobescavars, nprobevecvars, nprobetensvars  
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Inform the user of allocation.
    call TLS_info ('')
    call TLS_info ('Allocating base derived types ...', advance=.false.)

    ! Allocate the BC derived type.
    ALLOCATE (BC(ncells), STAT = memstat)
    if (memstat /= 0) call TLS_panic ('BASE_TYPES_ALLOCATE: BC derived type memory allocation error')

    ! BC base type.
    BC%Flag = 0
    BC%Internal = 0

    ! Allocate the Cell derived type.
    ALLOCATE (Cell(ncells), STAT = memstat)
    if (memstat /= 0) call TLS_panic ('BASE_TYPES_ALLOCATE: Cell derived type memory allocation error')

    ! Allocate the Zone derived type.
    ALLOCATE (Zone(ncells), STAT = memstat)
    if (memstat /= 0) call TLS_panic ('BASE_TYPES_ALLOCATE: Zone derived type memory allocation error')

    !Allocate Fluxing Velocity variable
    ALLOCATE (Fluxing_Velocity(nfc,ncells), STAT = memstat)
    if (memstat /= 0) call TLS_panic ('BASE_TYPES_ALLOCATE: Fluxing_Velocity derived type memory allocation error')

    ! Allocate the Matl derived type.
    call SLOT_INCREASE (Matl, mat_slot, mat_slot_new)

    ! Allocate the material property, displacement, strain, and stress arrays
    call SOLID_MECHANICS_ALLOCATE ()

    ! allocate arrays for the turbulence model
    call TURBULENCE_ALLOCATE ()

    ! allocate any probes structures if the PROBE namelist exists in input file

    if (nprobes > 0) then

       ALLOCATE (probes(nprobes), STAT = memstat)
       
       ems  = 0
       if (EM_is_on()) then
          ems = 1
       end if
       sms   = 0
       smv   = 0
       smt   = 0
       if (solid_mechanics) then
          sms = 3
          smv = 1
          smt = 4
       end if

       nprobescavars  = 7 + ems + nmat + sms
       nprobevecvars  = 1 + smv
       nprobetensvars = smt

       nprobevars     = nprobescavars + nprobevecvars + nprobetensvars

       do i=1,nprobes
          ALLOCATE(probes(i)%pid(nprobevars), STAT = memstat)
          ALLOCATE(probes(i)%NameLU(nprobevars), STAT = memstat)
          ALLOCATE(probes(i)%ScalarVarLU(nprobescavars), STAT = memstat)
          ALLOCATE(probes(i)%VectorVarLU(nprobevecvars), STAT = memstat)
          ALLOCATE(probes(i)%TensorVarLU(nprobetensvars), STAT = memstat)
       end do

       call TLS_fatal_if_any ((memstat /= 0), 'BASE_TYPES_ALLOCATE: error allocating probes pointer')

    endif

    ! Set the new arrays to their defaults
    call BASE_TYPES_DEFAULT ()

    ! Inform the user of allocation.
    call TLS_info ('done.')

  END SUBROUTINE BASE_TYPES_ALLOCATE

  SUBROUTINE BASE_TYPES_DEALLOCATE ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Deallocate the base types.
    !
    !=======================================================================
    use bc_module,         only: BC
    use matl_module,       only: SLOT_DECREASE, Matl
    use mesh_module,       only: Cell
    use parameter_module,  only: mat_slot, mat_slot_new, nprobes
    use probe_module,      only: probes
    use zone_module,       only: ZONE
    use fluid_data_module, only: Fluxing_Velocity

    ! Local Variables
    integer :: i
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! deallocate the BC derived type.
    if (ASSOCIATED(BC)) DEALLOCATE (BC)

    ! Deallocate the Cell derived type.
    if (ASSOCIATED(Cell)) DEALLOCATE (Cell)

    ! Deallocate the Zone derived type.
    if (ASSOCIATED(Zone)) DEALLOCATE (Zone)

    ! Deallocate Fluxing_Velocity.
    if (ASSOCIATED(Fluxing_Velocity)) DEALLOCATE (Fluxing_Velocity)

    if (nprobes > 0) then
       do i=1,nprobes
          if (ASSOCIATED(probes(i)%NameLU))      DEALLOCATE(probes(i)%NameLU)
          if (ASSOCIATED(probes(i)%ScalarVarLU)) DEALLOCATE(probes(i)%ScalarVarLU)
          if (ASSOCIATED(probes(i)%VectorVarLU)) DEALLOCATE(probes(i)%VectorVarLU)
          if (ASSOCIATED(probes(i)%TensorVarLU)) DEALLOCATE(probes(i)%TensorVarLU)
       end do
       DEALLOCATE (probes)
    end if

    ! Deallocate the Matl derived type by decreasing to zero slots.
    mat_slot_new = 0
    call SLOT_DECREASE (Matl, mat_slot, mat_slot_new)

  END SUBROUTINE BASE_TYPES_DEALLOCATE

  SUBROUTINE BASE_TYPES_DEFAULT ()
    !=======================================================================
    ! Purpose:
    !
    !   Default the base types.
    !
    !=======================================================================
    use mesh_module,       only: Cell
    use parameter_module,  only: ndim, nfc
    use zone_module,       only: Zone

    ! Local variables
    integer :: f, n
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Cell base type.
    Cell%Volume = 0.0_r8
    do n = 1,ndim
       Cell%Centroid(n) = 0.0_r8
       do f = 1,nfc
          Cell%Face_Normal(n,f) = 0.0_r8
          Cell%Face_Centroid_L(n,f) = 0.0_r8
       end do
    end do
    do f = 1,nfc
       Cell%Face_Area(f) = 0.0_r8
       Cell%Halfwidth(f) = 0.0_r8
    end do

    ! Zone base type.
    Zone%Rho          = 0.0_r8
    Zone%Rho_old      = 0.0_r8
    Zone%Temp         = 0.0_r8
    Zone%Temp_old     = 0.0_r8
    Zone%Enthalpy     = 0.0_r8
    Zone%Enthalpy_old = 0.0_r8
    Zone%P            = 0.0_r8
    do n = 1,ndim
       Zone%Vc(n)     = 0.0_r8
       Zone%Vc_old(n) = 0.0_r8
    end do

  END SUBROUTINE BASE_TYPES_DEFAULT

  SUBROUTINE MESH_VERTEX_ALLOCATE (Mesh, Vertex, ncells_input, nnodes_input)
    !=======================================================================
    ! Purpose(s):
    !
    !   Allocate the mesh types.
    !
    !=======================================================================
    use mesh_module,          only: MESH_CONNECTIVITY, VERTEX_DATA, Vrtx_Bdy
    use parameter_module,     only: ncells, nnodes, ndim

    ! Arguments
    type(MESH_CONNECTIVITY),  dimension(:),  pointer :: Mesh
    type(VERTEX_DATA),        dimension(:),  pointer :: Vertex
    integer, optional, intent(IN)   :: ncells_input
    integer, optional, intent(IN)   :: nnodes_input

    ! Local Variables
    integer :: ncells_local, nnodes_local, n
    logical, save :: first_time = .true.
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Inform the user if this the first time here.
    if (first_time) then
       call TLS_info ('Allocating mesh derived types ... ', advance=.false.)
    end if

    ! If sizes were passed in, use those, otherwise use global ncells and nnodes.
    if (PRESENT(ncells_input)) then
       ncells_local = ncells_input
    else
       ncells_local = ncells
    end if
    if (PRESENT(nnodes_input)) then
       nnodes_local = nnodes_input
    else
       nnodes_local = nnodes
    end if


    ! Allocate the Mesh & Vertex derived type.
    call ALLOCATE_MESH (Mesh, ncells_local)
    call ALLOCATE_VERTEX (Vertex, nnodes_local)

    ! Formerly a side effect of VERTEX_DATA_PRESET() (called from allocate_vertex)
    do n = 1, ndim
      if (associated(Vrtx_Bdy(n)%data)) deallocate(Vrtx_Bdy(n)%data)
    end do

    ! Inform the user of allocation.
    if (first_time) then
       call TLS_info ('done.')
       first_time = .false.
       call ANNOUNCE_MESH_SIZES (Mesh, Vertex)
    end if

  END SUBROUTINE MESH_VERTEX_ALLOCATE

  SUBROUTINE ANNOUNCE_MESH_SIZES(Mesh, Vertex)
    !=======================================================================
    ! Purpose(s):
    !
    !   Announce cell and node sizes on each processor
    !
    !=======================================================================
    use mesh_module,          only: MESH_CONNECTIVITY, VERTEX_DATA
    use parallel_info_module, only: p_info
    use pgslib_module,        only: PGSLib_COLLATE

    ! Arguments
    type(MESH_CONNECTIVITY),  dimension(:),  pointer :: Mesh
    type(VERTEX_DATA),        dimension(:),  pointer :: Vertex

    ! Local Variables
    integer :: pe
    integer, dimension(p_info%nPE) :: Ncells_All, NNodes_All
    character(len=128) :: message
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    call PGSLib_COLLATE (NCells_All, SIZE(Mesh,1))
    call PGSLib_COLLATE (Nnodes_All, SIZE(Vertex,1))

    call TLS_info('', TLS_VERB_NOISY)
    do pe = 1, SIZE(NCells_All,1)
       write(message,'(3(a,i0))') ' On pe ', pe, ' ncells = ', NCells_All(pe), ' nnodes = ', NNodes_All(pe)
       call TLS_info (message, TLS_VERB_NOISY)
    end do

    end subroutine announce_mesh_sizes
       
  SUBROUTINE ALLOCATE_MESH (Mesh, ncells_input)
    !=======================================================================
    ! Purpose(s):
    !
    !   Allocate an instance of the Mesh_Connectivity data structure.
    !   Sets all components to default values
    !
    !=======================================================================
    use mesh_module, only: MESH_CONNECTIVITY

    ! Arguments
    type(MESH_CONNECTIVITY), dimension(:), pointer :: Mesh
    integer, intent(IN) :: ncells_input

    ! Local Variables
    integer :: memstat
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Allocate the Mesh derived type.
    ALLOCATE (Mesh(ncells_input), STAT = memstat)
    if (memstat /= 0) call TLS_panic ('BASE_TYPES_ALLOCATE: Mesh derived type memory allocation error')

  END SUBROUTINE ALLOCATE_MESH

  SUBROUTINE ALLOCATE_VERTEX (Vertex, nnodes_input)
    !=======================================================================
    ! Purpose(s):
    !
    !   Allocate the mesh types.
    !
    !=======================================================================
    use mesh_module, only: VERTEX_DATA

    ! Arguments
    type(VERTEX_DATA), dimension(:), pointer :: Vertex
    integer, optional, intent(IN) :: nnodes_input

    ! Local Variables
    integer :: memstat
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Allocate the Vertex derived type.
    ALLOCATE (Vertex(nnodes_input), STAT = memstat)
    if (memstat /= 0) call TLS_panic ('ALLOCATE_VERTEX: Vertex derived type memory allocation error')

  END SUBROUTINE ALLOCATE_VERTEX

  SUBROUTINE MESH_VERTEX_DEALLOCATE (Mesh, Vertex)
    !=======================================================================
    ! Purpose(s):
    !
    !   Deallocate the Mesh and Vertex base types.
    !
    !=======================================================================
    use mesh_module, only: MESH_CONNECTIVITY, VERTEX_DATA

    ! Arguments
    type(MESH_CONNECTIVITY), dimension(:),  pointer :: Mesh
    type(VERTEX_DATA),       dimension(:),  pointer :: Vertex

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Deallocate the Mesh derived type
    call DEALLOCATE_MESH (Mesh)

    ! Deallocate the Mesh derived type
    call DEALLOCATE_VERTEX (Vertex)

  END SUBROUTINE MESH_VERTEX_DEALLOCATE

  SUBROUTINE DEALLOCATE_MESH (Mesh)
    !=======================================================================
    ! Purpose(s):
    !
    !   Deallocate an instance of the Mesh derived type
    !
    !=======================================================================
    use mesh_module, only: MESH_CONNECTIVITY

    ! Arguments
    type(MESH_CONNECTIVITY), dimension(:),  pointer :: Mesh

    ! Local Variables
    integer :: memstat
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Deallocate the Mesh derived type
    if (ASSOCIATED(Mesh)) DEALLOCATE (Mesh)

  END SUBROUTINE DEALLOCATE_MESH

  SUBROUTINE DEALLOCATE_VERTEX (Vertex)
    !=======================================================================
    ! Purpose(s):
    !
    !   Deallocate the Vertex base types.
    !
    !=======================================================================
    use mesh_module, only: VERTEX_DATA

    ! Arguments
    type(VERTEX_DATA), dimension(:), pointer :: Vertex

    ! Local Variables
    integer :: memstat
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Deallocate the Vertex derived type
    if (ASSOCIATED(Vertex)) DEALLOCATE (Vertex)

  END SUBROUTINE DEALLOCATE_VERTEX

  SUBROUTINE MESH_REALLOCATE (MeshPermute, VertexPermute)
    !=======================================================================
    ! Purpose(s):
    !
    !   Re-Allocate the mesh.
    !   The permutation vectors tell how to permute the old mesh
    !   into the new mesh.  The permutation vectors have the
    !   new data layout.  
    !   The new ncells and nnodes are determined implicitly from
    !   the sizes of MeshPermute and VertexPermute.
    !   The Mesh and Vertex arrays may not have local sizes of
    !   ncells and nnodes on input, but the do on output.
    !   WARNING: ncells and nnodes are changed in this routine.
    !
    !=======================================================================
    use mesh_module,   only: Mesh, Vertex, MESH_CONNECTIVITY, VERTEX_DATA
    use parallel_info_module, only: p_info
    use parameter_module,     only: ncells, nnodes
    use pgslib_module,        only: PGSLib_SUM_PREFIX, PGSLib_PERMUTE, pgslib_collate

    ! Arguments
    integer, dimension(:), intent(IN) :: MeshPermute
    integer, dimension(:), intent(IN) :: VertexPermute

    ! Local Variables
    type(MESH_CONNECTIVITY),  dimension(:), pointer :: Mesh_Local
    type(VERTEX_DATA),        dimension(:), pointer :: Vertex_Local
    integer, dimension(:), pointer :: MeshPermute_OrigLayout, MeshPermute_Index, &
                                      VertexPermute_OrigLayout, VertexPermute_Index
    integer :: pe
    integer, dimension(p_info%nPE) :: Ncells_All, NNodes_All
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Inform the user of re-allocation.
    call TLS_info ('')
    call TLS_info (' Reallocating mesh derived types ...')
    
    ! Make local copies on inputs
    call MESH_VERTEX_ALLOCATE (Mesh_Local, Vertex_Local, SIZE(Mesh,1), SIZE(Vertex,1))
    Mesh_Local   = Mesh
    Vertex_Local = Vertex

    ! Deallocate inputs, reallocate to new sizes
    call MESH_VERTEX_DEALLOCATE (Mesh, Vertex)

    ! New sizes
    ncells = SIZE(MeshPermute,1)
    nnodes = SIZE(VertexPermute,1)

    call MESH_VERTEX_ALLOCATE (Mesh, Vertex, ncells, nnodes)

    ! To do the permuation we need the permute vector to have
    ! the original layout (distribution)
    ALLOCATE (MeshPermute_OrigLayout(SIZE(Mesh_local,1)))
    ALLOCATE (MeshPermute_Index(SIZE(MeshPermute,1)))
    ALLOCATE (VertexPermute_OrigLayout(SIZE(Vertex_local,1)))
    ALLOCATE (VertexPermute_Index(SIZE(VertexPermute,1)))

    ! Use the _OrigLayout vectors as temporaries
    ! We are doing an identity permutation, but the distributions
    ! of the source and dest differ.
    MeshPermute_Index = 1
    MeshPermute_Index = PGSLib_SUM_PREFIX (MeshPermute_Index)
    call PGSLib_PERMUTE (DEST   = MeshPermute_OrigLayout, &
                         SOURCE = MeshPermute,            &
                         INDEX  = MeshPermute_Index)

    VertexPermute_Index = 1
    VertexPermute_Index = PGSLib_SUM_PREFIX (VertexPermute_Index)
    call PGSLib_PERMUTE (DEST   = VertexPermute_OrigLayout, &
                         SOURCE = VertexPermute,            &
                         INDEX = VertexPermute_Index)

    ! Permute the input (which is now stored in the _Local arrays)
    ! into the output (which is now the Mesh and Vertex arrays).
    call MESH_VERTEX_PERMUTE (Mesh_Local, Vertex_Local,  &
                       MeshPermute_OrigLayout,    &
                       VertexPermute_OrigLayout)

    ! Deallocate unneeded arrays and derived types.
    DEALLOCATE (MeshPermute_OrigLayout, VertexPermute_OrigLayout, &
                MeshPermute_Index, VertexPermute_Index) 
    call MESH_VERTEX_DEALLOCATE (Mesh_Local, Vertex_Local)

    ! Renumber mesh so that it points to new cell and vertex locations.
    call MESH_RENUMBER (MeshPermute, VertexPermute)

    ! Inform the user of a successful reallocation.
    call TLS_info (' Mesh derived types reallocated.')

    call PGSLib_COLLATE (NCells_All, SIZE(Mesh,1))
    call PGSLib_COLLATE (Nnodes_All, SIZE(Vertex,1))
    do pe = 1, SIZE(NCells_All,1)
       write(message, 10) PE, NCells_All(pe), NNodes_all(pe)
       call TLS_info (message)
10     format (1x,' On pe ',i8, ' ncells = ',i10, ' nnodes = ',i10)
    end do

  END SUBROUTINE MESH_REALLOCATE

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE PERMUTE_MESH(Mesh, MeshPermutation, SCOPE)
    !=======================================================================
    ! Purpose(s):
    !
    !   Permute the Mesh argument 
    !   This routine permutes all fields of Mesh, so can be used
    !   at any time.
    !   The result has the layout (per-processor array sizes) of MeshPermutation
    !   That may be different from the input layout.
    !
    !   Also, updates UnPermute_Mesh_Vector vector by composing with MeshPermutation
    !   and, updates Permute_Mesh vector
    !
    !=======================================================================
    use mesh_module,   only: MESH_CONNECTIVITY, UnPermute_Mesh_Vector, &
                             Permute_Mesh_Vector,               &
                             UnPermute_Mesh_Initialized
    use parameter_module, only: ncells
    use pgslib_module, only: PGSLib_Deallocate_Trace,     &
         &                   PGSLib_GS_Trace,             &
         &                   PGSLib_Permute,              &
         &                   PGSLib_Redistribute,         &
         &                   PGSLib_Global_ANY,           &
         &                   PGSLib_SUM_PREFIX,           &
         &                   PGSLib_SCOPE, OPERATOR(==),  &
         &                   PGSLib_Local

    ! Arguments
    type (MESH_CONNECTIVITY), dimension(:), POINTER :: Mesh
    integer, dimension(:), TARGET, intent(INOUT) :: MeshPermutation
    type (PGSLib_SCOPE),intent(IN), OPTIONAL :: SCOPE

    ! Local variables
    integer, pointer :: Permuter(:), Offset(:), Id(:)
    type(MESH_CONNECTIVITY), pointer :: Mesh_New(:), Mesh_Old(:)
    type(PGSLib_GS_Trace), pointer :: Mesh_Permute_Trace
    integer :: i
    logical :: LocalPermuter, LocalScope_P

    call TLS_info ('')
    call TLS_info (' Permuting mesh ... ', advance=.false.)

    ! If the scope is local we have to make a global permutation vector 
    ! (sounds silly, but when it comes to renumbering, makes life much easier)
    if (PRESENT(SCOPE)) then
       LocalScope_P = (SCOPE == PGSLib_Local)
    else
       LocalScope_P = .FALSE.
    end if
    
    ! If scope is local, make permutation vector global
    if (LocalScope_P) then
       ALLOCATE(Offset(SIZE(MeshPermutation)))
       Offset = 1
       Offset = PGSLib_SUM_PREFIX(Offset)
       MeshPermutation = MeshPermutation + Offset(1) - 1
       DEALLOCATE(Offset)
    end if
    

    ! We must check that MeshPermutation has the same layout as Mesh
    ! If it doesn't then we have to make a temporary permute vector
    ! which has that layout.
    IF ( PGSLib_Global_ANY ( (/ SIZE(Mesh) /= SIZE(MeshPermutation) /) ) ) then
       ! Permute and mesh have different layouts, so need to first
       ! get the permutation vector to original layout.
       ALLOCATE(Permuter(SIZE(Mesh)))
       
       ! One way to get the distribution changed is to PACK the input
       ! permutation vector into the new one.

       call PGSLib_REDISTRIBUTE (Permuter, SOURCE = MeshPermutation)


       ! We need to remember that we allocated Permuter in this routine
       LocalPermuter = .TRUE.
    ELSE
       Permuter => MeshPermutation
       LocalPermuter = .FALSE.
    END IF

    ! Now we are ready to start the Permuting.  We need to make some space
    call ALLOCATE_MESH (Mesh_New, SIZE(MeshPermutation, 1))
    
    ! Since we use the same pattern repeatedly, save the trace, but
    ! first put it into known state.
    NULLIFY(Mesh_Permute_Trace)
    ! Permute Mesh%Ngbr_Cell
    do i = 1,SIZE(Mesh_New(1)%Ngbr_Cell,1)
       call PGSLib_PERMUTE (DEST   = Mesh_New%Ngbr_Cell(i), &
                            SOURCE = Mesh%Ngbr_Cell(i),     &
                            INDEX  = Permuter,              &
                            TRACE  = Mesh_Permute_Trace)
    end do

    ! Permute Mesh%Ngbr_Cell_Orig
    do i = 1,SIZE(Mesh_New(1)%Ngbr_Cell,1)
       call PGSLib_PERMUTE (DEST   = Mesh_New%Ngbr_Cell_Orig(i), &
                            SOURCE = Mesh%Ngbr_Cell_Orig(i),     &
                            INDEX  = Permuter,                   &
                            TRACE  = Mesh_Permute_Trace)
    end do

    ! Permute Mesh%Ngbr_Cell_PE
    do i = 1,SIZE(Mesh_New(1)%Ngbr_Cell,1)
       call PGSLib_PERMUTE (DEST   = Mesh_New%Ngbr_Cell_PE(i), &
                            SOURCE = Mesh%Ngbr_Cell_PE(i),     &
                            INDEX  = Permuter,                 &
                            TRACE  = Mesh_Permute_Trace)
    end do

    ! Permute Mesh_New%Ngbr_Face.
    do i = 1,SIZE(Mesh_New(1)%Ngbr_Face,1)
       call PGSLib_PERMUTE (DEST   = Mesh_New%Ngbr_Face(i),       &
                            SOURCE = Mesh%Ngbr_Face(i), &
                            INDEX  = Permuter,             &
                            TRACE  = Mesh_Permute_Trace)
    end do

    ! Permute Mesh%Ngbr_Vrtx.
    do i = 1,SIZE(Mesh_New(1)%Ngbr_Vrtx,1)
       call PGSLib_PERMUTE (DEST   = Mesh_New%Ngbr_Vrtx(i), &
                            SOURCE = Mesh%Ngbr_Vrtx(i),     &
                            INDEX  = Permuter,              &
                            TRACE  = Mesh_Permute_Trace)
    end do

    ! Permute Mesh%Ngbr_Vrtx_Orig.
    do i = 1,SIZE(Mesh_New(1)%Ngbr_Vrtx_Orig,1)
       call PGSLib_PERMUTE (DEST   = Mesh_New%Ngbr_Vrtx_Orig(i), &
                            SOURCE = Mesh%Ngbr_Vrtx_Orig(i),     &
                            INDEX  = Permuter,                   &
                            TRACE  = Mesh_Permute_Trace)
    end do

    ! Permute Mesh%Ngbr_Vrtx_PE.
    do i = 1,SIZE(Mesh_New(1)%Ngbr_Vrtx_PE,1)
       call PGSLib_PERMUTE (DEST   = Mesh_New%Ngbr_Vrtx_PE(i), &
                            SOURCE = Mesh%Ngbr_Vrtx_PE(i),     &
                            INDEX  = Permuter,                 &
                            TRACE  = Mesh_Permute_Trace)
    end do

    ! Permute Mesh%cell_shape
    call PGSLib_PERMUTE (DEST   = Mesh_New%cell_shape, &
                         SOURCE = Mesh%cell_shape,     &
                         INDEX  = Permuter,            &
                         TRACE  = Mesh_Permute_Trace)

    ! Permute Mesh%CBlockID
    call PGSLib_PERMUTE (DEST   = Mesh_New%CBlockID, &
                         SOURCE = Mesh%CBlockID,     &
                         INDEX  = Permuter,          &
                         TRACE  = Mesh_Permute_Trace)

    ! Done with the Mesh permutation.
    call PGSLib_DEALLOCATE_TRACE (Mesh_Permute_Trace)

    ! Now point at new mesh
    Mesh_Old => Mesh
    Mesh     => Mesh_New
    
    ! Since we changed the mesh size, need to update ncells
    ncells = SIZE(Mesh,1)

    ! And get rid of the original mesh
    call DEALLOCATE_MESH (Mesh_Old)

    ! Finally, update UnPermute_Mesh_Vector, and Permute_Mesh_Vector

    ! Permute_Mesh_Vector is the total permutation done to date

    ! Need this since Permute_Mesh_Vector needs to be in known state
    if (.NOT. UnPermute_Mesh_Initialized) NULLIFY(Permute_Mesh_Vector)
    if (ASSOCIATED(Permute_Mesh_Vector)) DEALLOCATE(Permute_Mesh_Vector)
    ALLOCATE(Permute_Mesh_Vector(SIZE(Mesh,1)))

    if (UnPermute_Mesh_Initialized) then
       call PGSLib_PERMUTE (DEST   = Permute_Mesh_Vector, &
                           SOURCE = Permuter,            &
                           INDEX  = UnPermute_Mesh_Vector)
       DEALLOCATE(UnPermute_Mesh_Vector)
    else
       UnPermute_Mesh_Initialized = .TRUE.
       call PGSLib_REDISTRIBUTE (DEST   = Permute_Mesh_Vector, &
                                 SOURCE = Permuter)
       NULLIFY(UnPermute_Mesh_Vector)
    end if

    ! done with the permutation vector
    if (LocalPermuter) DEALLOCATE(Permuter)
    
    ! We like to store the inverse of Permute_Mesh_Vector
    ! UnPermute_Mesh_Vector always has the size of the Mesh
    ALLOCATE(UnPermute_Mesh_Vector(SIZE(Mesh)))

    ! To get the inverse, we apply Permute_Mesh_Vector to the idenity
    ALLOCATE(Id(SIZE(Permute_Mesh_Vector)))
    Id = PGSLib_SUM_PREFIX((/ (1, i=1,SIZE(Id)) /) )

    call PGSLib_PERMUTE (DEST   = UnPermute_Mesh_Vector, &
                         SOURCE = Id,             &
                         INDEX  = Permute_Mesh_Vector)

    ! Clean up
    DEALLOCATE(Id)

    call TLS_info ('done.')

  END SUBROUTINE PERMUTE_MESH

  SUBROUTINE PERMUTE_VERTEX(Vertex, VertexPermutation)
    !=======================================================================
    ! Purpose(s):
    !
    !   Permute the Vertex argument 
    !   This routine permutes all fields of Vertex, so can be used
    !   at any time.
    !   The result has the layout (per-processor array sizes) of VertexPermutation
    !   That may be different from the input layout.
    !
    !   Also, updates UnPermute_Vertex_Vector by composing with VertexPermutation,
    !   and updates Permute_Vertex_Vector.
    !
    !=======================================================================
    use mesh_module,   only: Vertex_Data, UnPermute_Vertex_Vector, &
                             Permute_Vertex_Vector, Vrtx_Bdy,      &
                             UnPermute_Vertex_Initialized
    use parameter_module, only: nnodes, ndim
    use pgslib_module, only: PGSLib_Deallocate_Trace,     &
         &                   PGSLib_GS_Trace,             &
         &                   PGSLib_Permute,              &
         &                   PGSLib_Redistribute,         &
         &                   PGSLib_Global_ANY,           &
         &                   PGSLib_SUM_PREFIX

    ! Arguments
    type (VERTEX_DATA), dimension(:), &
                        POINTER       :: Vertex
    integer, dimension(:), &  
                        TARGET,       &
                        intent(IN   ) :: VertexPermutation
    
    ! Local variables
    integer, dimension(:), &
                        POINTER       :: Permuter
    integer, dimension(:), &
                        POINTER       :: Id
    type (VERTEX_DATA), dimension(:), &
                        POINTER       :: Vertex_New
    type (VERTEX_DATA), dimension(:), &
                        POINTER       :: Vertex_Old
    type (PGSLib_GS_Trace), POINTER   :: Vertex_Permute_Trace
    integer :: i, n
    logical :: LocalPermuter

    call TLS_info (' Permuting vertices ... ', advance=.false.)

    ! We must check that VertexPermutation has the same layout as Vertex
    ! If it doesn't then we have to make a temporary permute vector
    ! which has that layout.
    IF ( PGSLib_Global_ANY ( (/ SIZE(Vertex) /= SIZE(VertexPermutation) /) ) ) then
       ! Permute and mesh have different layouts, so need to first
       ! get the permutation vector to original layout.
       ALLOCATE(Permuter(SIZE(Vertex)))
       
       ! One way to get the distribution changed is to PACK the input
       ! permutation vector into the new one.

       call PGSLib_REDISTRIBUTE (Permuter, VertexPermutation)


       ! We need to remember that we allocated Permuter in this routine
       LocalPermuter = .TRUE.
    ELSE
       Permuter => VertexPermutation
       LocalPermuter = .FALSE.
    END IF

    ! Now we are ready to start the Permuting.  We need to make some space
    call ALLOCATE_Vertex(Vertex_New, SIZE(VertexPermutation, 1))

    ! Formerly a side effect of VERTEX_DATA_PRESET() (called from allocate_vertex)
    do n = 1, ndim
      if (associated(Vrtx_Bdy(n)%data)) deallocate(Vrtx_Bdy(n)%data)
    end do
    
    ! Since we use the same pattern repeatedly, save the trace, but
    ! first put it into known state.
    NULLIFY(Vertex_Permute_Trace)

    ! Permute Vertex%Coord
    do i = 1,SIZE(Vertex_New(1)%Coord,1)
       call PGSLib_PERMUTE (DEST   = Vertex_New%Coord(i), &
                            SOURCE = Vertex%Coord(i),     &
                            INDEX  = Permuter,            &
                            TRACE  = Vertex_Permute_Trace)
    end do

    call PGSLib_PERMUTE (DEST   = Vertex_New%RSum_rvol, &
                         SOURCE = Vertex%RSum_rvol,     &
                         INDEX  = Permuter,             &
                         TRACE  = Vertex_Permute_Trace)


    ! Done with the Vertex permutation.
    call PGSLib_DEALLOCATE_TRACE (Vertex_Permute_Trace)

    ! Now point at new vertex
    Vertex_Old => Vertex
    Vertex     => Vertex_New
    
    ! Since we changed the local sizes of Vertex, need to update nnodes
    nnodes = SIZE(Vertex,1)

    ! And get rid of the original vertex
    call DEALLOCATE_VERTEX (Vertex_Old)

    ! Finally, update UnPermute_Vertex_Vector and Permute_Vertex_Vector

    ! Permute_Vertex_Vector is the total permutation done to date

    ! Need to make sure Permute_Vertex_Vector is in a known state
    if (.NOT. UnPermute_Vertex_Initialized) NULLIFY(Permute_Vertex_Vector)
    if (ASSOCIATED(Permute_Vertex_Vector))  DEALLOCATE(Permute_Vertex_Vector)
    ALLOCATE(Permute_Vertex_Vector(SIZE(Vertex)))

    if (UnPermute_Vertex_Initialized) then
       call PGSLib_PERMUTE(DEST   = Permute_Vertex_Vector, &
                           SOURCE = Permuter,              &
                           INDEX  = UnPermute_Vertex_Vector)
       if (ASSOCIATED(UnPermute_Vertex_Vector)) DEALLOCATE(UnPermute_Vertex_Vector)
    else
       UnPermute_Vertex_Initialized = .TRUE.
       call PGSLib_REDISTRIBUTE (DEST   = Permute_Vertex_Vector, &
                                 SOURCE = Permuter)
       NULLIFY(UnPermute_Vertex_Vector)
    end if
    
    ! done with the permutation vector
    if (LocalPermuter) DEALLOCATE(Permuter)
    

    ! We like to store the inverse of Permute_Vertex_Vector
    ! UnPermute_Vertex_Vector always has the size of the Vertex
    ALLOCATE(UnPermute_Vertex_Vector(SIZE(Vertex)))

    ! To get the inverse, we apply Permute_Vertex_Vector to the idenity
    ALLOCATE(Id(SIZE(Permute_Vertex_Vector)))
    Id = PGSLib_SUM_PREFIX((/ (1, i=1,SIZE(Id)) /) )

    call PGSLib_PERMUTE (DEST   = UnPermute_Vertex_Vector, &
                         SOURCE = Id,                      &
                         INDEX  = Permute_Vertex_Vector)

    ! Clean up
    DEALLOCATE(Id)

    call TLS_info ('done.')

  END SUBROUTINE PERMUTE_VERTEX

  SUBROUTINE MESH_VERTEX_PERMUTE (Mesh_Local, Vertex_Local, MeshPermute, VertexPermute)
    !=======================================================================
    ! Purpose(s):
    !
    !   Permute the base Mesh & Vertex types from _Local arrays into final
    !   arrays. (_Local means local to this module.  These arrays are
    !   distributed.)
    !
    !=======================================================================
    use mesh_module,   only: Mesh, Vertex, MESH_CONNECTIVITY, VERTEX_DATA
    use pgslib_module, only: PGSLib_Deallocate_Trace, &
         &                   PGSLib_GS_Trace,             &
         &                   PGSLib_Permute,              &
         &                   PGSLib_Global_SUM

    ! Arguments
    type(MESH_CONNECTIVITY),  dimension(:), intent(IN) :: Mesh_Local
    type(VERTEX_DATA),        dimension(:), intent(IN) :: Vertex_Local
    integer, dimension(SIZE(Mesh_Local,1)),   intent(IN) :: MeshPermute
    integer, dimension(SIZE(Vertex_Local,1)), intent(IN) :: VertexPermute

    ! Local Variables
    integer :: i
    ! Save the trace between permute calls, to improve performance
    type (PGSLib_GS_Trace), POINTER :: Mesh_Trace, VTX_Trace

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
       
    ! Set up to do mesh permutations
    ! First check that we have appropriate source and dest
    if (PGSLib_Global_SUM(SIZE(Mesh,1)) < PGSLib_Global_SUM(SIZE(Mesh_Local,1))) then
       call TLS_panic ('MESH_VERTEX_PERMUTE: Total size of Mesh not large enough to hold Mesh_Local')
    end if

    ! Initialize trace to put it into known state.
    ! Since this is first use of Trace, it must be NULL.  F90 doesn't do that
    ! automatically (although F90 can), so do it manually.
    NULLIFY(Mesh_Trace)
    
    ! Permute Mesh%Ngbr_Cell.
    do i = 1,SIZE(Mesh(1)%Ngbr_Cell,1)
       call PGSLib_PERMUTE (DEST   = Mesh%Ngbr_Cell(i),       &
                            SOURCE = Mesh_Local%Ngbr_Cell(i), &
                            INDEX  = MeshPermute,             &
                            TRACE  = Mesh_Trace)
    end do

    ! Permute Mesh%Ngbr_Face.
    do i = 1,SIZE(Mesh(1)%Ngbr_Face,1)
       call PGSLib_PERMUTE (DEST   = Mesh%Ngbr_Face(i),       &
                            SOURCE = Mesh_Local%Ngbr_Face(i), &
                            INDEX  = MeshPermute,             &
                            TRACE  = Mesh_Trace)
    end do

    ! Permute Mesh%Ngbr_Vrtx.
    do i = 1,SIZE(Mesh(1)%Ngbr_Vrtx,1)
       call PGSLib_PERMUTE (DEST   = Mesh%Ngbr_Vrtx(i),       &
                            SOURCE = Mesh_Local%Ngbr_Vrtx(i), &
                            INDEX  = MeshPermute,             &
                            TRACE  = Mesh_Trace)
    end do

    ! Done with the Mesh permutation.
    call PGSLib_Deallocate_Trace (Mesh_Trace)

    ! Set up to do Vertex permutations
    ! First check that we have appropriate source and dest
    if (PGSLib_Global_SUM(SIZE(Vertex,1)) < PGSLib_Global_SUM(SIZE(Vertex_Local,1))) then
       call TLS_panic ('MESH_VERTEX_PERMUTE: Total size of Vertex not large enough to hold Vertex_Local')
    end if

    ! Initialize trace to put it into known state.
    ! Since this is first use of Trace, it must be NULL.  F90 doesn't do that
    ! automatically (although F90 can), so do it manually.
    NULLIFY(VTX_Trace)

    ! Permute Vertex%Coord.
    do i = 1,SIZE(Vertex(1)%Coord,1)
       call PGSlib_PERMUTE (DEST   = Vertex%Coord(i),       &
                            SOURCE = Vertex_Local%Coord(i), &
                            INDEX  = VertexPermute,         &
                            TRACE  = VTX_Trace)
    end do

    ! Permute Vertex%Rsum_Rvol.
    call PGSLib_PERMUTE (DEST   = Vertex%Rsum_rvol,       &
                         SOURCE = Vertex_Local%Rsum_rvol, &
                         INDEX  = VertexPermute,          &
                         TRACE  = VTX_Trace)

    ! Done with Vertex permutation
    call PGSLib_Deallocate_Trace(VTX_Trace)

  END SUBROUTINE MESH_VERTEX_PERMUTE

  SUBROUTINE MESH_RENUMBER (MeshPermute, VertexPermute)
    !=======================================================================
    ! Purpose(s):
    !
    !   Renumber Mesh%Ngbr_Vrtx and Mesh_Ngbr_Cell to point to new
    !   locations of vertices and cells, based on the permutation vectors.
    !
    !=======================================================================
    use mesh_module,      only: Mesh
    use parameter_module, only: ncells, nvc
    use pgslib_module,    only: PGSLib_GATHER

    ! Arguments
    integer, dimension(ncells), intent(IN) :: MeshPermute
    integer, dimension(:),      intent(IN) :: VertexPermute

    ! Local Variables
    integer, dimension(nvc,ncells) :: Mesh_Vrtx, Mesh_Vrtx_New
    integer :: i

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    call TLS_info (' Renumbering mesh ... ', advance=.false.)
    ! On input Mesh%Ngbr_Vrtx points to the old vertex (node) locations.
    ! VertexPermute tells the new location for each vertex.
    do i = 1, SIZE(Mesh(1)%Ngbr_Vrtx, 1)
       Mesh_Vrtx(i,:) = Mesh%Ngbr_Vrtx(i)
    end do

    ! Mesh_Vrtx_New gets the new vertex number
    call PGSLib_GATHER (DEST   = Mesh_Vrtx_New, &
                        SOURCE = VertexPermute, &
                        INDEX  = Mesh_Vrtx)

    do i = 1, SIZE(Mesh(1)%Ngbr_Vrtx, 1)
       Mesh%Ngbr_Vrtx(i) = Mesh_Vrtx_New(i,:)
    end do
    
    call TLS_info ('done.')

  END SUBROUTINE MESH_RENUMBER

  SUBROUTINE RENUMBER_CELLS_VERTICES (Mesh, VertexPermute)
    !=======================================================================
    ! Purpose(s):
    !
    !   Renumber Mesh%Ngbr_Vrtx and Mesh_Ngbr_Cell to point to new
    !   locations of vertices and cells, based on the permutation vectors.
    !
    !=======================================================================
    use mesh_module,      only: MESH_CONNECTIVITY
    use parameter_module, only: nvc
    use pgslib_module,    only: PGSLib_GATHER

    ! Arguments
    type(MESH_CONNECTIVITY), dimension(:), intent(INOUT) :: Mesh
    integer, dimension(:), TARGET, intent(IN) :: VertexPermute

    ! Local Variables
    integer, dimension(:,:), pointer :: Mesh_Vrtx, Mesh_Vrtx_New
    integer :: i

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    call TLS_info (' Renumbering vertices ... ', advance=.false.)

    ! On input Mesh%Ngbr_Vrtx points to the old vertex (node) locations.
    ! VertexPermute tells the new location for each vertex.

    ALLOCATE(Mesh_Vrtx(nvc, SIZE(Mesh,1)))
    ALLOCATE(Mesh_Vrtx_New(nvc, SIZE(Mesh,1)))

    do i = 1, nvc
       Mesh_Vrtx(i,:) = Mesh%Ngbr_Vrtx(i)
    end do

    ! Mesh_Vrtx_New gets the new vertex number
    call PGSLib_GATHER (DEST   = Mesh_Vrtx_New, &
                        SOURCE = VertexPermute, &
                        INDEX  = Mesh_Vrtx)

    do i = 1, nvc
       Mesh%Ngbr_Vrtx(i) = Mesh_Vrtx_New(i,:)
    end do
    
    DEALLOCATE(Mesh_Vrtx_New, Mesh_Vrtx)
    call TLS_info ('done.')

  END SUBROUTINE RENUMBER_CELLS_VERTICES

END MODULE BASE_TYPES_MODULE
