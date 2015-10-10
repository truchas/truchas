MODULE GS_UTIL
  !=======================================================================
  ! Purpose(s):
  !
  !   Define utilities used to support gather & scatter.
  !
  !=======================================================================
  use gs_info_module, only: EE_TRACE,               &
                            EE_All_Ngbr_Trace,      &
                            EE_Mask_Initialized,    &
                            EL_Nbr_Mask,            &
                            NN_All_Ngbr_Trace
  implicit none
  private

  ! Public procedures
  Public :: GS_INIT_EE_MASK, &
            EE_GS_INIT,      &
            EN_GS_INIT,      &
            NN_GS_INIT

  ! Arrays and variables used only inside this module

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE GS_INIT_EE_MASK()
    !=======================================================================
    ! Purpose(s):
    !
    !=======================================================================
    use truchas_logging_services, only: TLS_fatal_if_any
    use mesh_module,      only: Mesh, DEGENERATE_FACE
    use parameter_module, only: ncells, nfc

    ! Local variables
    integer :: memerror, f

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE (El_Nbr_MASK(nfc, ncells), STAT=memerror)
    call TLS_fatal_if_any (memerror /= 0, 'EE_GATHER_INT: could not allocate El_Nbr_Mask')
    
    do f = 1,nfc
       El_Nbr_Mask(f,:) = Mesh%Ngbr_cell(f) /= 0 .and. &
                          Mesh%Ngbr_cell(f) /= DEGENERATE_FACE
    end do

  END SUBROUTINE GS_INIT_EE_MASK

  SUBROUTINE EE_GS_INIT()
    !=======================================================================
    ! Purpose(s):
    !
    !=======================================================================
    use ArrayAllocate_Module
    use mesh_module,      only: Mesh
    use parameter_module, only: ncells, nfc
    use pgslib_module,    only: PGSLIB_SETUP_TRACE
    use var_vector_module

    ! Local variables
    integer :: f
    integer, dimension(nfc,ncells) :: Mesh_Ngbr_Cell
    integer, dimension(:,:), POINTER :: Mesh_Ngbr_Cell_PE
    integer, dimension(:), POINTER :: Mesh_Ngbr_Cells_All

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    if (.NOT. EE_MASK_Initialized) then
       EE_MASK_Initialized = .true.
       call GS_INIT_EE_MASK()
    end if
    do f = 1,nfc
       Mesh_Ngbr_Cell(f,:) = Mesh%Ngbr_cell(f)
    end do

    ! We want pgslib to return values in Mesh_Ngbr_Cell_PE.  To signal that,
    ! we send it an unallocated pointer.  Since pointers are not NULL by defualt
    ! (in F90, fixed in F95), we must NULLIFY.
    NULLIFY(Mesh_Ngbr_Cell_PE)

    EE_Trace => PGSLib_Setup_Trace(INDEX        = Mesh_Ngbr_Cell,   &
                                   SIZE_OF_DEST = ncells,           &
                                   PE_ARRAY     = Mesh_Ngbr_Cell_PE,&
                                   MASK         = El_Nbr_Mask)
    do f = 1, nfc
       ! Save original Ngbr_Cell
       Mesh%Ngbr_Cell_Orig(f) = Mesh%Ngbr_Cell(f) 
       ! Extract the new index numbers
       Mesh%Ngbr_Cell(f)      = Mesh_Ngbr_Cell(f,:) 
       ! Save the PE number, for global/local stuff
       Mesh%Ngbr_Cell_PE(f)   = MERGE(Mesh_Ngbr_Cell_PE(f,:), -1, EL_Nbr_Mask(f,:))
    end do

    ! This code is for gathering/scattering from/to all neighbor cells
    ! Notice that here we don't need a mask, since the length of the
    ! neighbor list varies.  
    
    ! Preserve the original global cell numbers
    Call CREATE(Mesh%Ngbr_Cells_All_Orig, SIZES = SIZES(Mesh%Ngbr_Cells_All))
    Mesh%Ngbr_Cells_All_Orig = FLATTEN(Mesh%Ngbr_Cells_All)

    ! We need to move the list of neighbors into a single long array
    ! in order to pass it to pgslib for initialization.
    ! That is provided by FLATTENing Mesh%Ngbr_Cells_All
    Mesh_Ngbr_Cells_All => FLATTEN(Mesh%Ngbr_Cells_All)
    
    EE_All_Ngbr_Trace => PGSLib_Setup_Trace(INDEX = Mesh_Ngbr_Cells_All, &
                                            SIZE_OF_DEST = ncells)
    
  END SUBROUTINE EE_GS_INIT

  SUBROUTINE EN_GS_INIT()
    !=======================================================================
    ! Purpose(s):
    !
    !=======================================================================
    use gs_info_module,   only: EN_TRACE
    use mesh_module,      only: Mesh
    use parameter_module, only: ncells, nnodes, nvc
    use pgslib_module,    only: PGSLIB_SETUP_TRACE

    ! Local variables
    integer :: v
    integer, dimension(nvc, ncells) :: Mesh_Ngbr_Vrtx
    integer, dimension(:,:), POINTER :: Mesh_Ngbr_Vrtx_PE

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do v = 1, nvc
       Mesh_Ngbr_Vrtx(v,:) = Mesh%Ngbr_Vrtx(v)
    end do

    ! We want pgslib to return values in Mesh_Ngbr_Vrtx_PE.  To signal that,
    ! we send it an unallocated pointer.  Since pointers are not NULL by defualt
    ! (in F90, fixed in F95), we must NULLIFY.
    NULLIFY(Mesh_Ngbr_Vrtx_PE)

    EN_Trace => PGSLib_Setup_Trace(INDEX        = Mesh_Ngbr_Vrtx,    &
                                   SIZE_OF_DEST = NNodes,            &
                                   PE_ARRAY     = Mesh_Ngbr_Vrtx_PE)
                                
    do v = 1,nvc
       ! Save original vertex number in case we need it
       Mesh%Ngbr_Vrtx_Orig(v) = Mesh%Ngbr_Vrtx(v)
       ! Extract the new vertex number into our data structure
       Mesh%Ngbr_Vrtx(v)      = Mesh_Ngbr_Vrtx(v,:)
       ! Save the PE number, for global/local stuff
       Mesh%Ngbr_Vrtx_PE(v)   = Mesh_Ngbr_Vrtx_PE(v,:)
    end do
  
  END SUBROUTINE EN_GS_INIT

  SUBROUTINE NN_GS_INIT()
    !=======================================================================
    ! Purpose(s): Initialize the trace for node-node gather & scatter
    !
    !=======================================================================
    use ArrayAllocate_Module
    use mesh_module,      only: Vertex_Ngbr_All, Vertex_Ngbr_All_Orig
    use parameter_module, only: nnodes
    use var_vector_module
    use pgslib_module,    only: PGSLIB_SETUP_TRACE

    ! Local variables
    integer, dimension(:), POINTER :: Vertex_Ngbr_List

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! This code is for gathering/scattering from/to all neighbor vertices
    ! Notice that here we don't need a mask, since the length of the
    ! neighbor list varies.  
    
    ! Preserve the original global cell numbers
    ALLOCATE(Vertex_Ngbr_All_Orig(SIZE(Vertex_Ngbr_All)))
    CALL CREATE(Vertex_Ngbr_All_Orig, SIZES = SIZES(Vertex_Ngbr_All))
    Vertex_Ngbr_All_Orig = FLATTEN(Vertex_Ngbr_All)

    ! We need to move the list of neighbors into a single long array
    ! in order to pass it to pgslib for initialization.
    Vertex_Ngbr_List => FLATTEN(Vertex_Ngbr_All)
    
    NN_All_Ngbr_Trace => PGSLib_Setup_Trace(INDEX = Vertex_Ngbr_List, &
                                            SIZE_OF_DEST = nnodes)

  END SUBROUTINE NN_GS_INIT

END MODULE GS_UTIL
