MODULE MESH_GEN_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures for generating the mesh.
  !
  !   Public Interface:
  !
  !     * call MESH_GEN ()
  !
  !        Read or create the mesh; then compute the mesh connectivity.
  !
  !     * call FLAG_FACE_NEIGHBORS ()
  !
  !        Sets the Nbrg_Cells_Face field of Mesh to indicated face neighbors.
  !
  !	
  ! Contains: CONNECTIVITY
  !           RCM
  !           MESH_AXIS
  !           MESH_GEN
  !           FLAG_FACE_NEIGHBORS
  !
  ! Author(s): Robert Ferrell (ferrell.cpca.com)
  !            Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  public :: MESH_GEN, FLAG_FACE_NEIGHBORS

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Private parameters
  integer, parameter :: NOT_LOCAL_INDEX = -1

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE FLAG_FACE_NEIGHBORS ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Fill Ngbr_Cells_All_Face field in the Mesh data structure
    !   This identifies which neighbor cells share a node
    !   with each of the nfc faces.
    !   For each neighbor of a reference cell we will compare all of its nvc
    !   vertices with the faces (as defined by the vertices) of the reference
    !   cell.  If the neighbor cell has a vertex in that face, then it is
    !   a face neighbor (neighbor of that face) and gets flagged as such.
    !   To avoid killing ourselves with too much memory, we can gather one 
    !   of the nvc vertices at a time.
    !
    !=======================================================================
    use gs_module,        only: EE_GATHER
    use mesh_module,      only: Mesh, Face_Vrtx, Initialize_Face_Bit_Mask,   &
                                Set_Face_Neighbor, Clear_Face_Neighbor,      &
                                DEGENERATE_FACE
    use parameter_module, only: ncells, nfc, nvf, nvc
    use var_vector_module

    ! Local Variables
    type(int_var_vector), pointer, dimension(:) :: Ngbr_Vertices  ! Holds the gathered vertices
    integer, pointer, dimension(:) :: Ngbr_Data      ! Points to neighbors of Ref_Cell
    integer, pointer, dimension(:) :: Ngbr_Face_Bits ! Points to neighbors of Ref_Cell
    integer, dimension(nvf) :: Reference_Face ! Current face of interest
    integer :: ref_cell, ngbr_vtx, f, v, ngbr

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    call TLS_info ('')
    call TLS_info (' Finding face neighbors ... ', advance=.false.)

    ! We must initialize the face bit mask array
    call INITIALIZE_FACE_BIT_MASK ()
    
    ! Allocate the storage to hold the face neighbor flags, then clear the flags
    call CREATE (Mesh%Ngbr_Cells_Face, SIZES=SIZES(Mesh%Ngbr_Cells_All))

    CLEAR_ALL: do Ref_Cell = 1, ncells
       Ngbr_Face_Bits => FLATTEN(Mesh(Ref_Cell)%Ngbr_Cells_Face)
       ! The next line initializes ngbr_face_bits, which was pointed
       ! to allocated but uninitialized storage in the create call
       ! above.  This avoids an intermittent complaint about
       ! CLEAR_FACE_NEIGHBOR using an uninitialized variable.  This
       ! line could probably replace the contents of the entire loop
       Ngbr_Face_Bits = 0
       do ngbr = 1, SIZE(Ngbr_Face_Bits)
          do f = 1, nfc
             call CLEAR_FACE_NEIGHBOR (Ngbr_Face_Bits(ngbr), f)
          end do
       end do
    end do CLEAR_ALL

    ! We need space to hold the neighbors vertices
    ALLOCATE(Ngbr_Vertices(ncells))
    call CREATE (Ngbr_Vertices, SIZES=SIZES(Mesh%Ngbr_Cells_All))

    ! Loop over all nvc vertices.
    do ngbr_vtx = 1, nvc
       
       call EE_GATHER (Ngbr_Vertices, Mesh%Ngbr_Vrtx_Orig(ngbr_vtx))

       CELLS_LOOP: do ref_cell = 1, ncells

          FACE_LOOP: do f = 1, nfc

             if (Mesh(ref_cell)%Ngbr_Cell_Orig(f) == DEGENERATE_FACE) cycle FACE_LOOP

             do v = 1, nvf
                Reference_Face(v) = Mesh(ref_cell)%Ngbr_Vrtx_Orig(Face_Vrtx(f,v))
             end do

             Ngbr_Data      => FLATTEN(Ngbr_Vertices(ref_cell))
             Ngbr_Face_Bits => FLATTEN(Mesh(ref_cell)%Ngbr_Cells_Face)
             
             NGBR_LOOP: do ngbr = 1, SIZE(Ngbr_Data)
             
                if (Touches_Face( Reference_Face, Ngbr_Data(ngbr) )) then
                   call SET_FACE_NEIGHBOR (Ngbr_Face_Bits(ngbr), f)
                end if
             
             end do NGBR_LOOP
          
          end do FACE_LOOP

       end do CELLS_LOOP

    end do
                
    ! Deallocate the Ngbr_Vertices temporary
    ! First dismiss the var_vector components, then the array itself
    Call DESTROY(Ngbr_Vertices)
    DEALLOCATE(Ngbr_Vertices)

    call TLS_info ('done.')

  END SUBROUTINE FLAG_FACE_NEIGHBORS
      
  SUBROUTINE MESH_GEN ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Read or create the mesh; then compute the mesh connectivity.
    !
    !=======================================================================
    use ArrayAllocate_Module,   only: ARRAYCREATE, ARRAYDESTROY
    use base_types_B_module,    only: PERMUTE_MESH, &
                                      PERMUTE_VERTEX, RENUMBER_CELLS_VERTICES, &
                                      ANNOUNCE_MESH_SIZES
    use bc_data_module,         only: Mesh_Face_Set, Mesh_Face_Set_Tot
    use mesh_input_module,      only: MESH_READ, mesh_file, &
                                      coordinate_scale_factor, &
                                      use_RCM, MESH_READ_SIDE_SETS
    use mesh_partition_module,  only: MESH_PARTITIONS
    use mesh_module,            only: Mesh, Vertex,             &
                                      MESH_CONNECTIVITY,        &
                                      VERTEX_DATA, Vrtx_Bdy,    &
                                      Permute_Mesh_Vector,      &
                                      Permute_Vertex_Vector,    &
                                      OPERATOR(.DISTRIBUTE.),   &
                                      Vertex_Ngbr_All,          &
                                      mesh_has_cblockid_data
    use mesh_utilities,         only: NODE_CONNECTIVITY, SET_DEGENERATE_FACES
    use parallel_info_module,   only: p_info
    use parameter_module,       only: ncells, ncells_tot, nnodes, &
                                      nnodes_tot, ndim, nfc, nssets
    use pgslib_module,          only: PGSLib_DIST,            &
                                      PGSLIB_PERMUTE,         &
                                      PGSLib_LOCAL,           &
                                      PGSLib_BCAST
    use restart_variables,      only: restart
    use restart_driver,         only: restart_mesh, restart_side_sets
    use partitioner_data,       only: PARTITIONER_INIT
    use var_vector_module
    use input_utilities, only: NULL_C

    ! Local Variables
    logical :: fatal, read_mesh, read_face_sets
    integer :: memerror, n, m
    type(VERTEX_DATA), allocatable, dimension(:) :: Vertex_Tot
    type(MESH_CONNECTIVITY), allocatable, dimension(:) :: Mesh_Tot
    integer, pointer, dimension(:) :: MeshPermute, VertexPermute, Node_BC_Tmp, &
                                      RCM_Permute, Face_Set_Tmp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Nullify pointers.
    NULLIFY (MeshPermute)
    NULLIFY (VertexPermute)
    NULLIFY (RCM_Permute)
    NULLIFY (Mesh_Face_Set_Tot)


    ! make sure the user partitioning info is right
    call PARTITIONER_INIT ()

    ! Set various read mesh flags.
    read_mesh         = (mesh_file /= NULL_C) .or. restart
    read_face_sets    = .false.

    ! Assign vertices by reading or creating mesh 
    ASSIGN_VERTEX: if (read_mesh) then

       ! Read existing mesh.
       if (restart) then

          call restart_mesh ()
          call set_degenerate_faces(Mesh)
          call mesh_read_side_sets () ! defines Mesh_Face_Set_Tot
          read_face_sets = ASSOCIATED(Mesh_Face_Set_Tot)
          call restart_side_sets () ! skip side set data in older restart files

       else

          ! Read mesh; read into temporaries of total size,
          ! then distribute mesh across pe's.
          if (p_info%IOP) then
             ALLOCATE (Mesh_tot(ncells_tot), Vertex_tot(nnodes_tot), STAT = memerror)
          else
             ALLOCATE (Mesh_tot(0), Vertex_tot(0), STAT = memerror)
          end if
          fatal = (memerror /= 0)
          call TLS_fatal_if_any (fatal, 'MESH_GEN: could not allocate space for mesh temporary arrays')

          !! NNC, Dec 2013.  The Vertex_Data type is now default initialized
          !! making the following assignment unnecessary.  However the rhs
          !! function reference had the side effect of resetting the Vrtx_Bdy
          !! structure to its initial state (unbelievable!)  That is now done here.
          !! Vertex_tot = VERTEX_DATA_PRESET()
          do n = 1, ndim
            if (associated(Vrtx_Bdy(n)%data)) deallocate(Vrtx_Bdy(n)%data)
          end do

          ! Set flag if mesh file has cell 'material' data: IDEAS, PATRAN, TOPAZ, EXODUSII.
          mesh_has_cblockid_data = .true.

          ! Read the mesh.
          if (p_info%IOP) then

             ! I/O processor reads the mesh into total arrays.
             call MESH_READ (Mesh_Tot, Vertex_Tot)

             ! If face set(s) were read from mesh
             read_face_sets    = ASSOCIATED(Mesh_Face_Set_Tot)
          end if
 
         ! Distribute the Mesh just read in.
          Mesh   = .DISTRIBUTE. Mesh_Tot
          Vertex = .DISTRIBUTE. Vertex_Tot
          DEALLOCATE (Mesh_Tot, Vertex_Tot)

          call set_degenerate_faces(Mesh)

       end if

    else

       ! Create mesh; create into temporaries of total size,
       ! then distribute mesh across pe's.
       if (p_info%IOP) then
          ALLOCATE (Mesh_Tot(ncells_tot), Vertex_Tot(nnodes_tot), STAT = memerror)
       else
          ALLOCATE (Mesh_Tot(0),          Vertex_Tot(0),          STAT = memerror)
       end if
       fatal = (memerror /= 0)
       call TLS_fatal_if_any (fatal, 'MESH_GEN: could not allocate space for mesh temporary arrays')

       if (p_info%IOP) then
          call CREATE_MESH (Mesh_Tot, Vertex_Tot)
       end if
       !! NNC, Dec 2013.  The Vertex_Data type is now default initialized
       !! making the following assignment unnecessary.  However the rhs
       !! function reference had the side effect of resetting the Vrtx_Bdy
       !! structure to its initial state (unbelievable!)  That is now done here.
       !! The assignment was originally done in CREATE_MESH, meaning it was
       !! only done on the IO process -- probably not correct, but it didn't
       !! matter because in this execution path the data structure had never
       !! been allocated.
       !! Vertex_tot = VERTEX_DATA_PRESET() ! moved from inside create_mesh
       do n = 1, ndim
         if (associated(Vrtx_Bdy(n)%data)) deallocate(Vrtx_Bdy(n)%data)
       end do

       Mesh   = .DISTRIBUTE. Mesh_Tot
       Vertex = .DISTRIBUTE. Vertex_Tot
       DEALLOCATE (Mesh_Tot, Vertex_Tot)

    end if ASSIGN_VERTEX

    ! Broadcast nodal BC input flags in case they have changed.
    if (.not. p_info%UseGlobalServices) then
       call PGSLib_BCAST (read_face_sets)
       call PGSLib_BCAST (nssets)
    end if

    ! Scale the coordinates equally; default scale factor is unity.
    ! Only do this if this is NOT a restart.  If this is a restart,
    ! then the scaled coordinates were dumped, so don't rescale.
    if (.NOT. Restart) then
       do n = 1,ndim
          Vertex%Coord(n) = coordinate_scale_factor*Vertex%Coord(n)
       end do
    end if

    ! Find the Mesh partitions, and the permutation vectors which
    ! take current mesh to partitioned mesh.
    call MESH_PARTITIONS (MeshPermute, VertexPermute)

    ! Reallocate mesh so that it will be properly partitioned.
    ! Also renumber Mesh%Ngbr_Vrtx so that it points to correct
    ! nodes after reallocation.

    call PERMUTE_MESH (Mesh, MeshPermute)
    call PERMUTE_VERTEX (Vertex, VertexPermute)
    call RENUMBER_CELLS_VERTICES (Mesh, VertexPermute)

    call ANNOUNCE_MESH_SIZES(Mesh, Vertex)
    ! Done with MeshPermute and VertexPermutation
    DEALLOCATE (MeshPermute, VertexPermute)

    if (use_RCM) then

       ! Need All_Ngbr_CONNECTIVITY for RCM
       call ALL_NGBR_CONNECTIVITY (Mesh)

       ! Renumber the mesh via reverse Cuthill-McKee renumbering

       ! Need to allocate RCM_Permute on the fly, since size of each
       ! domain changes in mesh_reallocate.
       ALLOCATE(RCM_Permute(SIZE(Mesh)))
       RCM_Permute = COMPUTE_RCM (Mesh)

       ! Permute mesh according to RCM permutation.  Since we haven't
       ! done anything with node numbering (yet) don't need to permute nodes.
       ! Scope is local because this is a per-domain re-ordering.

       call PERMUTE_MESH (Mesh, RCM_Permute, SCOPE=PGSLib_Local)

       ! do not have to renumber because we didn't permute nodes
   
       ! Done with RCM_Permute
       DEALLOCATE(RCM_Permute)
   
    end if

    ! don't have to renumber because we didn't permute the nodes
    ! Do have to rebuild cell connectivity, since that has changed.

    ! Compute mesh connectivity.  (This work should be combined with All_Ngbr_CONNECTIVITY)
    call CONNECTIVITY (Mesh)

    ! Need to recompute, since go clobbered inside Permute_Mesh
    call ALL_NGBR_CONNECTIVITY (Mesh)

    ! Need node connectivity
    ! Have to allocate this here, rather than at the top, since nnodes changes during 
    ! this routine.

    ALLOCATE(Vertex_Ngbr_All(nnodes))
    Vertex_Ngbr_All = NODE_CONNECTIVITY()

    ! Distribute face set data if available
    if (read_face_sets) then
       if (.not. p_info%IOP) ALLOCATE (Mesh_Face_Set_Tot(nssets,nfc,0))
       ALLOCATE (Mesh_Face_Set(nssets,nfc,ncells), STAT = memerror)
       if (memerror /= 0) call TLS_panic ('Mesh_Gen: could not Allocate Mesh_Face_Set')
       do n = 1,nssets
          do m = 1,nfc
             call PGSLib_DIST (Mesh_Face_Set(n,m,:), Mesh_Face_Set_Tot(n,m,:))
          end do
       end do

       DEALLOCATE (Mesh_Face_Set_Tot)
       ALLOCATE (Face_Set_Tmp(ncells), STAT = memerror)
       if (memerror /= 0) call TLS_panic ('Mesh_Gen: could not Allocate Face_Set_Tmp')
       do n = 1,nssets
          do m = 1,nfc
             Face_Set_Tmp = Mesh_Face_Set(n,m,:)
             call PGSLib_PERMUTE (DEST = Mesh_Face_Set(n,m,:), SOURCE = Face_Set_Tmp, &
                                  INDEX  = Permute_Mesh_Vector)
          end do
       end do
       DEALLOCATE (Face_Set_Tmp)
    end if

  END SUBROUTINE MESH_GEN

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE CREATE_MESH (Mesh_Tot, Vertex_Tot)
    !======================================================================
    ! Purpose:
    !
    !   Create a mesh, based on input parameters
    !
    !   Creates mesh only on IO processor.  The mesh must then be
    !   distributed across other processors (outside of this routine)
    !======================================================================
    use mesh_input_module,      only: Coord, Heps,                &
                                      Fuzz, Ncell, Nseg, Ratio,   &
                                      Coord_label
    use mesh_module,            only: MESH_CONNECTIVITY,          &
                                      VERTEX_DATA
    use parameter_module,       only: ncells_tot, nnodes_tot,     &
                                      Nx_tot, Mx_tot,             &
                                      ndim, nvc, nvf
    use random_module,          only: GENERATE_RANDOM
    
    ! Arguments
    type(VERTEX_DATA),       intent(INOUT), dimension(:) :: Vertex_Tot
    type(MESH_CONNECTIVITY), intent(INOUT), dimension(:) :: Mesh_Tot

    integer, dimension(0:nvc) :: Nvtx
    real(r8), pointer, dimension(:,:) :: Xv_Tot
    real(r8), dimension(ndim,nnodes_tot) :: Del_Tot
    real(r8), dimension(nnodes_tot)      :: Rnum_Tot
    real(r8) :: SHIFT, PI, xi, eta, zeta, dx, dy, dz
    logical :: orthog_mesh
    integer :: i, iend, j, jend, k, kend, n, nc, nv, maxdim, skip, v, vertex_start

    ! Test that the mesh and vertex structures are the proper size

    if (SIZE(Mesh_Tot,1) /= ncells_tot) then
       call TLS_panic ('CREATE_MESH: Mesh_tot not size ncells_tot')
    end if
    
    if (SIZE(Vertex_Tot,1) /= nnodes_tot) then
       call TLS_panic ('CREATE_MESH: Vertex_tot not size nnodes_tot')
    end if
    
    PI  = 4.0*atan(1.0)

    ! Scale Heps so not too large
    Heps = 0.5*Heps

    ! Create a logical mesh; allocate arrays
    maxdim = MAXVAL(Mx_Tot)
    ALLOCATE(Xv_Tot(ndim, maxdim))
    Xv_Tot = 0.

          ! Create vertex axes; assign orthogonal mesh flag
          orthog_mesh = .true.
          iend = 1; jend = 1; kend = 1
          do n = 1,ndim
             orthog_mesh = orthog_mesh .and. Fuzz(n) == 0.
             call MESH_AXIS (Coord_label(n), Nseg(n), Ncell(n,:), Coord(n,:), Ratio(n,:), &
                Mx_tot(n), maxdim, Xv_Tot(n,:))
             select case(n)
             case (1)
                iend = Mx_tot(n)
             case (2)
                jend = Mx_tot(n)
             case (3)
                kend = Mx_tot(n)
             end select
          end do

          dx = Xv_Tot(1,2) - Xv_Tot(1,1)
          dy = Xv_Tot(2,2) - Xv_Tot(2,1)
          dz = Xv_Tot(3,2) - Xv_Tot(3,1)

          ! Assign vertex coordinates
          nv = 1
          do k = 1, kend
             do j = 1, jend
                do i = 1, iend

                   xi   = (i-1)*dx
                   eta  = (j-1)*dy
                   zeta = (k-1)*dz

                   SHIFT = Heps*SIN( 2*PI*xi/(Xv_Tot(1,iend)-Xv_Tot(1,1)) )           &
                               *SIN( 2*PI*eta/(Xv_Tot(2,jend)-Xv_Tot(2,1)) )          &
                               *SIN( 2*PI*zeta/(Xv_Tot(3,kend)-Xv_Tot(3,1)) )

                   ! Assign vertex coordinates
                   do n = 1,ndim
                      select case(n)
                      case (1)
                         Vertex_Tot(nv)%Coord(n) = Xv_Tot(n,i) + SHIFT
                      case (2)
                         Vertex_Tot(nv)%Coord(n) = Xv_Tot(n,j) + SHIFT
                      case (3)
                         Vertex_Tot(nv)%Coord(n) = Xv_Tot(n,k) + SHIFT
                      end select
                   end do

                   ! Increment vertex number
                   nv = nv + 1

                end do
             end do
          end do

          ! Non-Orthogonal mesh (Fuzz interior vertices only!)
          RANDOM_MESH: if (.not.(orthog_mesh)) then
             ! Make sure kend is correct for 2-D case
             if (ndim == 2 .and. kend == 1) kend = 3

             ! Random number [-1.0, +1.0]
             call GENERATE_RANDOM(-1.0_r8,1.0_r8,nnodes_tot,Rnum_Tot)

             ! Initialize Min-Delta arrays
             Del_Tot = 0.

             ! Min-Cell width about each vertex
             vertex_start = 2
             do n = 1,ndim-1
                select case(n)
                case (1)
                   vertex_start = vertex_start + Mx_tot(n)
                case (2)
                   vertex_start = vertex_start + Mx_tot(n)*Mx_tot(n-1)
                end select
             end do
             do k = 2, kend - 1
                do j = 2, jend - 1
                   do i = 2, iend - 1

                      ! Compute vertex number
                      nv = vertex_start
                      do n = 1,ndim
                         select case(n)
                         case (1)
                            nv = nv + i - 2
                         case (2)
                            nv = nv + (j - 2)*Mx_tot(n-1)
                         case (3)
                            nv = nv + (k - 2)*Mx_tot(n-1)*Mx_tot(n-2)
                         end select
                      end do

                      ! Compute min-cell width
                      do n = 1,ndim
                         select case(n)
                         case (1)
                            Del_Tot(n,nv) = MIN((Xv_Tot(n,i+1)-Xv_Tot(n,i)), (Xv_Tot(n,i)-Xv_Tot(n,i-1)))
                         case (2)
                            Del_Tot(n,nv) = MIN((Xv_Tot(n,j+1)-Xv_Tot(n,j)), (Xv_Tot(n,j)-Xv_Tot(n,j-1)))
                         case (3)
                            Del_Tot(n,nv) = MIN((Xv_Tot(n,k+1)-Xv_Tot(n,k)), (Xv_Tot(n,k)-Xv_Tot(n,k-1)))
                         end select
                      end do

                   end do
                end do
             end do

             ! Limit min-delta arrays
             ! (The limiting factor must be less than one-half to avoid an
             ! overlap of the fuzzed vertices. A factor much less than
             ! one-half is not necessary but it is recommended; this further
             ! limits the vertex position adjustments.)
             do n = 1,ndim
                Del_Tot(n,:) = (1./2.)*Del_Tot(n,:)
             end do

             ! Fuzz interior vertex coordinates
             do k = 2, kend - 1
                do j = 2, jend - 1
                   do i = 2, iend - 1

                      ! Compute vertex number
                      nv = vertex_start
                      do n = 1,ndim
                         select case(n)
                         case (1)
                            nv = nv + i - 2
                         case (2)
                            nv = nv + (j - 2)*Mx_tot(n-1)
                         case (3)
                            nv = nv + (k - 2)*Mx_tot(n-1)*Mx_tot(n-2)
                         end select
                      end do

                      ! Fuzz interior vertices
                      do n = 1,ndim
                         Vertex_Tot(nv)%Coord(n) = Vertex_Tot(nv)%Coord(n) + &
                            Fuzz(n)*Rnum_Tot(nv)*Del_Tot(n,nv)
                      end do

                   end do
                end do
             end do

          end if RANDOM_MESH

          ! Destroy working arrays
          DEALLOCATE(Xv_Tot)

          ! Assign ending do loop indices
          iend = 1; jend = 1; kend = 1
          do n = 1,ndim
             select case(n)
             case (1)
                iend = Nx_tot(n)
             case (2)
                jend = Nx_tot(n)
             case (3)
                kend = Nx_tot(n)
             end select
          end do

          ! Get the vertices
          nc = 1
          do k = 1, kend
             do j = 1, jend
                do i = 1, iend

                   ! Compute cell vertices
                   do n = 1,ndim
                      select case(n)
                      case (1)
                         Nvtx(nvf) = i
                      case (2)
                         Nvtx(nvf) = Nvtx(nvf) + (j - 1)*Mx_tot(n-1)
                      case (3)
                         Nvtx(nvf) = Nvtx(nvf) + (k - 1)*Mx_tot(n-1)*Mx_tot(n-2)
                      end select
                   end do
                   Nvtx(0) = Nvtx(nvf)
                   do v = 1,nvf-1
                      select case(v)
                      case (1)
                         Nvtx(v) = Nvtx(v-1) + 1
                      case (2)
                         Nvtx(v) = Nvtx(v-1) + Mx_tot(1)
                      case (3)
                         Nvtx(v) = Nvtx(v-1) - 1
                      end select
                   end do
                   do n = 1,ndim
                      select case(n)
                      case (1)
                         skip = 1
                      case (2)
                         skip = skip*Mx_tot(n-1)
                      case (3)
                         skip = skip*Mx_tot(n-1)
                      end select
                   end do
                   do v = 1,nvc/2
                      Nvtx(v+nvf) = Nvtx(v) + skip
                   end do

                   ! Assign cell vertices
                   do v = 1,nvc
                      Mesh_Tot(nc)%Ngbr_Vrtx(v) = Nvtx(v)
                   end do

                   ! Increment cell number
                   nc = nc + 1

                end do
             end do
          end do

  END SUBROUTINE CREATE_MESH

  SUBROUTINE ALL_NGBR_CONNECTIVITY (Mesh)
    !=======================================================================
    ! Purpose(s):
    !
    !   Initialize the neighbor arrays Mesh%Ngbr_Cell_All, which lists all 
    !   neighbor cells.  Cells are neighbors if they share at least one vertex.
    !   Notice that this ismore general than Ngbr_Cell, which lists only
    !   cells which share a face.
    !
    !======================================================================

    use ArrayAllocate_Module, only: ArrayCreate, ArrayDestroy
    use mesh_module,          only: MESH_CONNECTIVITY
    use pgslib_module,        only: PGSLib_Deallocate_Trace, &
                                    PGSLib_Gather,           &
                                    PGSLIB_GLOBAL_EOSHIFT,   &
                                    PGSLib_GRADE_UP,         &
                                    PGSLib_GRADE_UP_LOCAL,   &
                                    PGSLib_GS_Trace,         &
                                    PGSLIB_PARITY_PREFIX,    &
                                    PGSLib_Permute,          &
                                    PGSLIB_SUM_PREFIX,       &
                                    PGSLib_SUM_SUFFIX,       &
                                    PGSLib_SCATTER_SUM    
    use parameter_module,     only: ncells, nvc, nnodes
    use var_vector_module

    ! Arguments
    type(MESH_CONNECTIVITY), dimension(ncells), intent(INOUT) :: Mesh

    ! Local variables
    integer :: c, v, i, l, p, n, ngbr_count
    integer :: Total_CV_Pairs, Total_Neighbors

    ! Local arrays
    integer, POINTER, dimension(:,:) :: CV_Pairs_ALL
    integer, POINTER, dimension(:)   :: Vertex_Rank
    integer, POINTER, dimension(:)   :: Global_Cell_Number
    integer, POINTER, dimension(:)   :: CV_Pairs_Temp
    logical, POINTER, dimension(:)   :: CV_Pairs_MASK
    logical, POINTER, dimension(:)   :: CV_Pairs_Seg
    integer, POINTER, dimension(:)   :: SegmentLength
    integer, POINTER, dimension(:)   :: SegmentOffset
    integer, POINTER, dimension(:)   :: Offset_V
    integer, POINTER, dimension(:)   :: Length_V
    integer, POINTER, dimension(:,:) :: Offset_C
    integer, POINTER, dimension(:,:) :: Length_C
    integer, POINTER, dimension(:,:) :: Ngbr_Vrtx
    integer, POINTER, dimension(:)   :: C_To_Seg_Index
    integer, POINTER, dimension(:)   :: C_Ngbr_List
    integer, POINTER, dimension(:)   :: C_to_C_Segments

    integer,    dimension(ncells)    :: Length_Too_Many
    integer,    dimension(ncells)    :: Length_Unique
    type(int_var_vector), dimension(ncells)    :: C_Ngbrs_Too_Many
    type(int_var_vector), dimension(ncells)    :: C_Ngbrs_Rank
    type(int_var_vector), dimension(ncells)    :: C_Ngbrs_Temp
    type(log_var_vector), dimension(ncells)    :: C_Ngbrs_Mask
    integer,    POINTER, dimension(:):: Too_Many_Temp, temp_ngbr_too_many
    integer,    POINTER, dimension(:):: Temp_Ngbr_List
    logical,    POINTER, dimension(:):: temp_ngbr_mask

    integer, PARAMETER :: Cell_Slot = 1
    integer, PARAMETER :: Vertex_Slot = 2

    ! Communication buffers
    type (PGSLib_GS_Trace), POINTER                     :: CV_Pairs_TRACE
    type (PGSLib_GS_Trace), POINTER                     :: S_To_V_Trace 
    type (PGSLib_GS_Trace), POINTER                     :: V_To_C_Trace

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Each hex has nvc vertices around it.  We will pair vertices and cells.
    ! Then we will group vertices.  Cells in the same group (segment) are neighbors.

    call TLS_info (' Finding all cell neighbors ... ', advance=.false.)

    ! The first list is all the (cell,vertex) pairs
    Total_CV_Pairs = ncells*nvc

    call ARRAYCREATE (CV_Pairs_All, 1, 2, 1, Total_CV_Pairs, "CV_Pairs__All in All_Ngbr_CONNECTIVITY")

    ! Initialize with (cell,vertex) pairs
    ! Need global cell number
    ALLOCATE(Global_Cell_Number(NCells))
    Global_Cell_Number = PGSLib_SUM_PREFIX((/ (1,c=1,ncells) /))
    p = 1
    do c = 1, ncells
       do v = 1, nvc
          CV_Pairs_All(Cell_Slot,  p) = Global_Cell_Number(c)
          CV_Pairs_All(Vertex_Slot,p) = Mesh(c)%Ngbr_Vrtx(v)
          p = p + 1
       end do
    end do
    DEALLOCATE(Global_Cell_Number)
    
    ! Permute CV_Pairs_All to create segments of contstant vertex number
    call ARRAYCREATE (Vertex_Rank, 1, Total_CV_Pairs, "Vertex_Rank in CONNECTIVITY")
    call ARRAYCREATE (CV_Pairs_Temp, 1, Total_CV_Pairs, "CV_Pairs_Temp in CONNECTIVITY")
    Vertex_Rank = PGSLib_GRADE_UP(CV_Pairs_All(Vertex_Slot,:))
    
    NULLIFY(CV_Pairs_TRACE)
    call PGSLib_PERMUTE (DEST   = CV_Pairs_Temp,     &
                         SOURCE = CV_Pairs_All(Cell_Slot,:), &
                         INDEX  = Vertex_Rank,        &
                         TRACE  = CV_Pairs_TRACE)
    CV_Pairs_All(Cell_Slot,:) = CV_Pairs_Temp

    call PGSLib_PERMUTE (DEST   = CV_Pairs_Temp,     &
                         SOURCE = CV_Pairs_All(Vertex_Slot,:), &
                         INDEX  = Vertex_Rank,        &
                         TRACE  = CV_Pairs_TRACE)
    CV_Pairs_All(Vertex_Slot,:) = CV_Pairs_Temp

    call PGSLib_DEALLOCATE_TRACE (CV_Pairs_TRACE)

    ! Now we need to find the segments of neighbors
    call ARRAYCREATE (CV_Pairs_Mask, 1, Total_CV_Pairs, "CV_Pairs_Mask in CONNECTIVITY")
    call ARRAYCREATE (CV_Pairs_Seg,  1, Total_CV_Pairs, "CV_Pairs_Seg in CONNECTIVITY")
    
    CV_Pairs_Mask = CV_Pairs_All(Vertex_Slot,:) /= &
                    PGSLib_Global_EOSHIFT(CV_Pairs_All(Vertex_Slot,:), &
                                          SHIFT = -1,        &
                                          BOUNDARY = -1)
    CV_Pairs_Seg = PGSLib_PARITY_PREFIX(CV_Pairs_Mask)
    
    ! Enumerate the segments.  The offset and length of each segment
    ! will get mapped to the nodes, then to the cells.  Notice that
    ! there is a 1<->1 mapping between nodes and segments (by construction
    ! of a segment).

    call ARRAYCREATE (SegmentLength, 1, Total_CV_Pairs, "SegmentLength in CONNECTIVITY")
    call ARRAYCREATE (SegmentOffset, 1, Total_CV_Pairs, "SegmentOffset in CONNECTIVITY")

    ! This puts the length of each segment into the head slot of each segment
    SegmentLength = 1
    SegmentLength = PGSLib_SUM_SUFFIX(SegmentLength, SEGMENT = CV_Pairs_Seg)

    ! This puts the offset of the head of each segment into the head slot of each segment
    SegmentOffset = 1
    SegmentOffset = PGSLib_SUM_PREFIX(SegmentOffset)
       
    ! Scatter to the vertices
    call ARRAYCREATE (Offset_V, 1, nnodes, "Offset_V in CONNECTIVITY")
    call ARRAYCREATE (Length_V, 1, nnodes, "Length_V in CONNECTIVITY")
    ! Make sure this is empty so it will get initialized.
    NULLIFY(S_To_V_TRACE)
    Offset_V = 0
    Length_V = 0
    call PGSLib_SCATTER_SUM (DEST   = Offset_V,          &
                             SOURCE = SegmentOffset,     &
                             INDEX  = CV_Pairs_All(Vertex_Slot,:), &
                             TRACE  = S_To_V_TRACE,      &
                             MASK   = CV_Pairs_MASK)

    call PGSLib_SCATTER_SUM (DEST   = Length_V,          &
                             SOURCE = SegmentLength,     &
                             INDEX  = CV_Pairs_All(Vertex_Slot,:), &
                             TRACE  = S_To_V_TRACE,      &
                             MASK   = CV_Pairs_MASK)
    call PGSLib_DEALLOCATE_TRACE (S_To_V_TRACE)


    ! We are done with a bunch of the temporaries
    call ARRAYDESTROY (Vertex_Rank)
    call ARRAYDESTROY (CV_Pairs_Temp)
    call ARRAYDESTROY (CV_Pairs_Mask)
    call ARRAYDESTROY (CV_Pairs_Seg )
    call ARRAYDESTROY (SegmentOffset)
    call ARRAYDESTROY (SegmentLength)

    ! Move the segment offsets and lengths from the vertices to the cells
    ! This is a (vector) gather.  (Vector because each cell gathers a single
    ! item from each vertex, which makes a total of nvc items gathered
    ! to each cell.)
    call ARRAYCREATE (Offset_C, 1, nvc, 1, ncells, "Offset_C in CONNECTIVITY")
    call ARRAYCREATE (Length_C, 1, nvc, 1, ncells, "Length_C in CONNECTIVITY")
    call ARRAYCREATE (Ngbr_Vrtx,1, nvc, 1, ncells, "Ngbr_Vrtx in CONNECTIVITY")

    ! Load up Ngbr_Vrtx to be in nice array format
    do c = 1, ncells
       do v = 1, nvc
          Ngbr_Vrtx(v,c) = Mesh(c)%Ngbr_Vrtx(v)
       end do
    end do
    
    NULLIFY(V_To_C_TRACE)
    call PGSLib_GATHER (DEST   = Offset_C,  &
                        SOURCE = Offset_V,  &
                        INDEX  = Ngbr_Vrtx, &
                        TRACE  = V_To_C_TRACE)
    call PGSLib_GATHER (DEST   = Length_C,  &
                        SOURCE = Length_V,  &
                        INDEX  = Ngbr_Vrtx, &
                        TRACE  = V_To_C_TRACE)

    ! Cleup up some more temporaries
    call ARRAYDESTROY (Offset_V)
    call ARRAYDESTROY (Length_V)
    call ARRAYDESTROY (Ngbr_Vrtx)
    call PGSLib_DEALLOCATE_TRACE (V_To_C_TRACE)


    ! Before we gather the segments (ie the neighbor lists) we need to make
    ! space for them.  The data goes into a ragged array, which is in Mesh.
    ! We count (on each PE) the number of neighbors total.  This is how
    ! much total ragged array we need.
    Total_Neighbors = SUM(Length_C)
    
    ! We will work with a flattened array
    call ARRAYCREATE (C_To_Seg_Index, 1, Total_Neighbors, "C_To_Seg_Index in CONNECTIVITY")
    call ARRAYCREATE (C_Ngbr_List, 1, Total_Neighbors, "C_Ngbr_List in CONNECTIVITY")
    call ARRAYCREATE (C_to_C_Segments, 1, ncells, "C_to_C_Segments in CONNECTIVITY")

    ! Fill this with the pointers to CV_Pairs_All.  This is where we need
    ! Offset_C and Length_C
    i = 1
    do c = 1, ncells
       ! Point to the start of each new segment of neighbors
       C_to_C_Segments(c) = i
       do v = 1, nvc
          do l = 1, Length_C(v,c)
             C_To_Seg_Index(i) = Offset_C(v,c) + l - 1
             i = i + 1
          end do
       end do
    end do
    
    
    ! Gather the list of neighbors
    call PGSLib_GATHER (DEST   = C_Ngbr_List,       &
                        SOURCE = CV_Pairs_All(Cell_Slot,:), &
                        INDEX  = C_To_Seg_Index)
    
    

    ! Done with CV_PAirs_All, so dump it
    call ARRAYDESTROY (CV_Pairs_All)
    ! We need to eliminate cells which got listed more than once.
    ! That is easiest doen in cell-centric data structure.
    ! so, move C_Ngbr_List into a ragged array.
    ! Each vetor in the ragged array has size Length_C(i).
    Length_Too_Many = SUM(Length_C, DIM=1)

    ! Done with a few temporaries, so we can ditch them before using more memory
    call ARRAYDESTROY (Offset_C)
    call ARRAYDESTROY (Length_C)

    call CREATE (C_Ngbrs_Too_Many, SIZES = Length_Too_Many)
    ! Transfer the data from a flat array into the ragged array
    C_Ngbrs_Too_Many =   C_Ngbr_List

    ! Get rid of a bunch of temporaries
    call ARRAYDESTROY (C_To_Seg_Index)
    call ARRAYDESTROY (C_Ngbr_List)
    call ARRAYDESTROY (C_to_C_Segments)

    ! We need to maks out all the duplicates in C_Ngbrs_Too_Many
    ! First, construct the mask, with T for each unique neighbor
    call CREATE (C_Ngbrs_Mask, SIZES = Length_Too_Many)
    call CREATE (C_Ngbrs_Rank, SIZES = Length_Too_Many)
    call CREATE (C_Ngbrs_Temp, SIZES = Length_Too_Many)

    do c = 1, ncells
       ! First find the rank of each neighbor
       C_Ngbrs_Rank(c) = PGSLib_GRADE_UP_Local(FLATTEN(C_Ngbrs_Too_Many(c)))

       ! Sort each varying vector, to get all same neighbor adjacent
       ! Need a temporary, since we are permuting and overwrting
       C_Ngbrs_Temp(c) = FLATTEN(C_Ngbrs_Too_Many(c))
       Too_Many_Temp => FLATTEN(C_Ngbrs_Temp(c))
       Too_Many_Temp(FLATTEN(C_Ngbrs_Rank(c))) = FLATTEN(C_Ngbrs_Too_Many(c))
       C_Ngbrs_Too_Many(c) = FLATTEN(C_Ngbrs_Temp(c))
       
       C_Ngbrs_Mask(c) = EOSHIFT(FLATTEN(C_Ngbrs_Too_Many(c)), &
                                 SHIFT = -1,                   &
                                 BOUNDARY = -1               ) &
                         /= FLATTEN(C_Ngbrs_Too_Many(c))
       Length_Unique(c) = COUNT(FLATTEN(C_Ngbrs_Mask(c)))
    end do

    ! We now know the total number of unique cell neighbors, so we can create
    ! the space to store them.
    call CREATE (Mesh%Ngbr_Cells_All, SIZES = Length_Unique)

    ! Finally, we can PACK out the duplicates.
    do c = 1, ncells
       temp_ngbr_list     => FLATTEN(Mesh(c)%Ngbr_Cells_All)
       temp_ngbr_too_many => FLATTEN(C_Ngbrs_Too_Many(c)   )
       temp_ngbr_mask     => FLATTEN(C_Ngbrs_Mask(c))
       ngbr_count = 1
       ! Pack out duplicates, explicit looping for easier debugging
!       Mesh(c)%Ngbr_Cells_All = PACK(FLATTEN(C_Ngbrs_Too_Many(c)),   &
!                                     MASK = FLATTEN(C_Ngbrs_Mask(c)),&
!                                     VECTOR=FLATTEN(Mesh(c)%Ngbr_Cells_All))

       do n = 1, SIZES(C_Ngbrs_Too_Many(c))
          if (temp_ngbr_mask(n)) then
             temp_ngbr_list(ngbr_count) = temp_ngbr_too_many(n)
             ngbr_count = ngbr_count + 1
          end if
       end do
       if ((ngbr_count - 1) /= SIZES(Mesh(c)%Ngbr_Cells_ALL)) then
          call TLS_panic ('ALL_NGBR_CONNECTIVITY: wrong number of neigbhors found')
       end if
          
#ifdef DEBUG_MESH_GEN
       Temp_Ngbr_List => FLATTEN(Mesh(c)%Ngbr_Cells_All)
       if (SIZES(Mesh(c)%Ngbr_Cells_All) > 0) then
          do n = 1, SIZES(Mesh(c)%Ngbr_Cells_All)
             write(25,*) c, SIZES(Mesh(c)%Ngbr_Cells_All), Temp_Ngbr_List(n)
          end do
       else
          write(25,*)    c, SIZES(Mesh(c)%Ngbr_Cells_All), -1
       end if
#endif       
    end do

#ifdef DEBUG_MESH_GEN
    close(25)
#endif

    ! Clean up the remaining local arrays
    call DESTROY (C_Ngbrs_Too_Many)
    call DESTROY (C_Ngbrs_Mask)
    call DESTROY (C_Ngbrs_Rank)
    call DESTROY (C_Ngbrs_Temp)

    call TLS_info ('done.')
    
  END SUBROUTINE ALL_NGBR_CONNECTIVITY

  SUBROUTINE CONNECTIVITY (Mesh)
    !=======================================================================
    ! Purpose(s):
    !
    !   Initialize the neighbor arrays Ngbr_Cell and Ngbr_Face, which are
    !   components of the Mesh derived type. Given the list of vertices
    !   associated with each cell (Mesh%Ngbr_Vrtx), find the cell numbers
    !   opposite each face (Mesh%Ngbr_Cell) of our logical cell (a square
    !   in 2-D and a cube in 3-D. Also find the face numbers of the neighbor
    !   cells sharing faces with each cell (Mesh%Ngbr_Face).
    !
    !   This task is carried out in the following fashion. Given a list of
    !   cells and the vertices surrounding each cell, we:
    !    1) Construct the list of faces around each cell, noting that a
    !       face is uniquely identified by the vertices at its corners.
    !    2) Some nodes may be degenerate, so indicate those.
    !    3) Sort the nodes of each face into descending order.
    !    4) Sort the faces, using a segmented sort.  Ultimately, same
    !       faces will become adjacent in the Faces list.
    !    5) Match faces with face up or down in list, using global_cshift.
    !    6) Send connecivity info back to Mesh data structure.
    !
    !   Discussion: The package that generated the problem has set up a 1-D
    !               array of cells and vertices. Each cell has a list of
    !               vertices that form that cell. This code must now
    !               determine the linear index of each of the cells that it
    !               shares a face with. It must also determine which of the
    !               possible nfc faces it shares with each cell. This 'orients'
    !               each cell with its neighbor cells. This is parallelized
    !               as well (using parallel sorts), authored by Robert C.
    !               Ferrell (ferrell@cpca.com) of CPCA, Ltd. on 5/29/97.
    !
    !=======================================================================
    use mesh_module,          only: Face_Vrtx, MESH_CONNECTIVITY,       &
                                    DEGENERATE_FACE, degenerate_points, &
                                    degenerate_lines, triangle_faces,   &
                                    quad_faces, DEGENERATE_NODE, EXTERNAL_FACE
    use pgslib_module,        only: PGSLib_GS_Trace,         &
                                    PGSLIB_SUM_PREFIX,       &
                                    PGSLIB_PARITY_PREFIX,    &
                                    PGSLIB_GRADE_UP,         &
                                    PGSLIB_PERMUTE,          &
                                    PGSLib_Deallocate_Trace, &
                                    PGSLIB_GLOBAL_EOSHIFT,   &
                                    PGSLIB_Setup_Trace,      &
                                    PGSLIB_SCATTER_BUFFER,   &
                                    PGSLIB_GLOBAL_COUNT,     &
                                    PGSLib_GLOBAL_SUM,       &
                                    PGSLib_Size_Of_Dup,      &
                                    PGSLib_Size_Of_Sup,      &
                                    PGSLib_Dup_Index

    use parameter_module,     only: boundary_faces, ncells, nfc, boundary_faces_tot, nvf

    ! Arguments
    type(MESH_CONNECTIVITY), dimension(ncells), intent(INOUT) :: Mesh

    ! Local variables
    integer :: allfaces
    
    ! The array Faces is the main work array of this routine.  It
    ! holds the list of faces, gets sorted and is used to 
    ! identify neighboring cells.
    ! The entries in Faces are listed here
    integer, dimension(nvf, SIZE(Mesh,1)*SIZE(Mesh(1)%Ngbr_Face,1)) :: Faces
    logical, dimension(     SIZE(Mesh,1)*SIZE(Mesh(1)%Ngbr_Face,1)) :: FaceSegment
    logical, dimension(     SIZE(Mesh,1)*SIZE(Mesh(1)%Ngbr_Face,1)) :: FaceBoundary
    logical, dimension(     SIZE(Mesh,1)*SIZE(Mesh(1)%Ngbr_Face,1)) :: FaceMask
    integer, dimension(     SIZE(Mesh,1)*SIZE(Mesh(1)%Ngbr_Face,1)) :: FaceRank
    ! HomeCell is the (face, cell) pair for the face
    ! AdjCell  is the (face, cell) pair for the adjacent face
    integer, dimension(2, SIZE(Mesh,1)*SIZE(Mesh(1)%Ngbr_Face,1)) :: HomeCell
    integer, dimension(2, SIZE(Mesh,1)*SIZE(Mesh(1)%Ngbr_Face,1)) :: AdjCell
    integer, dimension(   SIZE(Mesh,1)*SIZE(Mesh(1)%Ngbr_Face,1)) :: FTemp

    ! These are scratch arrays that are needed a bit.
    integer, dimension(SIZE(Mesh,1))      :: MeshGlobalNumber, ITemp

    ! Communication buffers
    type (PGSLib_GS_Trace), POINTER                     :: GS_Trace
    integer, dimension(:,:),   pointer :: AFSup, AFDup
    integer, dimension(:,:),   pointer :: ACSup, ACDup

    ! Misc integers
    integer :: tf1, f, v, v2, FaceIndex, c
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Each hex has nvf faces, some of which may be degenerate.
    allfaces = SIZE(Mesh,1)*SIZE(Mesh(1)%Ngbr_Face,1)

    write (message, 5) PGSLib_GLOBAL_SUM(allfaces)
5   format (1x,'Establishing mesh connectivity for', i10,' faces ...')
    call TLS_info ('')
    call TLS_info (message)
    call TLS_info ('')

    ! Construct the list of faces surrounding each cell
    ITemp = 1
    MeshGlobalNumber = PGSLib_SUM_PREFIX(ITemp)
    FaceIndex = 1
    do c = 1, SIZE(Mesh)
       do f = 1, SIZE(Mesh(1)%Ngbr_Face,1)
          do v = 1, nvf
             Faces(v, FaceIndex) = Mesh(c)%Ngbr_Vrtx(Face_Vrtx(f,v))
          end do
          HomeCell(1,FaceIndex) = f
          HomeCell(2,FaceIndex) = MeshGlobalNumber(c)
          FaceIndex = FaceIndex + 1
       end do
    end do
    
    ! If any node occurs twice in a face, then it is a degenerate node.
    ! Find degenerate nodes with a double loop
    do v = 1, nvf
       do v2 = v-1, 1, -1
          where (Faces(v,:) == Faces(v2,:)) Faces(v,:) = DEGENERATE_NODE
       end do
    end do
    
    ! Sort the nodes of each face into descending so that
    ! the nodes will uniquely identify the face.  Since nvc == 4 usually,
    ! n^2 sort is fine for this.
    do f = 1, Allfaces
       do v = 1, nvf
          tf1 = Faces(nvf,f)
          do v2 = nvf-1, v, -1
             if (tf1 > Faces(v2,f)) then ! swap 
                Faces(v2+1,f) = Faces(v2,f)
                Faces(v2,f)   = tf1
             else                        ! leave in place, 
                tf1 = Faces(v2,f)                
             end if
          end do
       end do
    end do

    ! Now sort the global face list to find matches.  
    FaceBoundary = .false.
    do v = 1, nvf
       write (message, '(a,i1,a)') '   Sorting face vertex ',v,' ...'
       call TLS_info (message)

       ! Segments are constructed from segment boundaries
       FaceSegment = PGSLib_PARITY_PREFIX (FaceBoundary)

       ! Rank the vertex 
       FaceRank = PGSLib_GRADE_UP (Faces(v,:), SEGMENT=FaceSegment)
       ! Permute the faces
       ! We need to use the same pattern many times, so save the setup
       ! in a trace.  Trace must be NULL for that to work.
       NULLIFY(GS_Trace)
       do v2 = 1, nvf
          FTemp = Faces(v2,:)
          call PGSLib_PERMUTE (DEST   = Faces(v2,:), &
                               SOURCE = FTemp,       &
                               INDEX  = FaceRank,    &
                               TRACE  = GS_Trace)
       end do
       ! Permute the (face,cell) pairs
       FTemp = HomeCell(1,:) 
       call PGSLib_PERMUTE (DEST = HomeCell(1,:), SOURCE = FTemp, INDEX = FaceRank, TRACE=GS_Trace)
       FTemp = HomeCell(2,:)
       call PGSLib_PERMUTE (DEST = HomeCell(2,:), SOURCE = FTemp, INDEX = FaceRank, TRACE=GS_Trace)
       ! Don't need to permute AdjCell because it is still empty.
       call PGSLib_Deallocate_Trace(GS_Trace)

       ! Segment boundaries for next iteration are the just sorted vertices
       FaceBoundary = FaceBoundary .or. &
            (Faces(v,:) /= PGSLib_GLOBAL_EOSHIFT (Faces(v,:), SHIFT=-1, BOUNDARY=-1))
    end do
    
    ! Find matching faces
    ! Initialize AdjCell to assume we won't find a neighbor.
    AdjCell = EXTERNAL_FACE
    ! Check for match from upward
    FaceMask = .true.
    do v = 1, nvf
       FaceMask = FaceMask .and. &
                  (Faces(v,:) == PGSLib_GLOBAL_EOSHIFT(Faces(v,:), SHIFT=-1, BOUNDARY=-1))
    end do
    ! For all matches, we have found the adjacent face
    AdjCell(1,:) = MERGE (PGSLib_GLOBAL_EOSHIFT (HomeCell(1,:), SHIFT=-1, BOUNDARY=-1), &
                          AdjCell(1,:), FaceMask)
    AdjCell(2,:) = MERGE (PGSLib_GLOBAL_EOSHIFT (HomeCell(2,:), SHIFT=-1, BOUNDARY=-1), &
                          AdjCell(2,:), FaceMask)
    ! Check for match from downward
    FaceMask = .true.
    do v = 1, nvf
       FaceMask = FaceMask .and. &
                  (Faces(v,:) == PGSLib_GLOBAL_EOSHIFT (Faces(v,:), SHIFT=1, BOUNDARY=-1))
    end do
    ! For all matches, we have found the adjacent face
    AdjCell(1,:) = MERGE (PGSLib_GLOBAL_EOSHIFT (HomeCell(1,:), SHIFT=1, BOUNDARY=-1), &
                          AdjCell(1,:), FaceMask)
    AdjCell(2,:) = MERGE (PGSLib_GLOBAL_EOSHIFT (HomeCell(2,:), SHIFT=1, BOUNDARY=-1), &
                          AdjCell(2,:), FaceMask)
             
    ! Degenerate faces get identified
    where (COUNT((Faces == DEGENERATE_NODE),DIM=1) >= 2)
       AdjCell(1,:) = DEGENERATE_FACE
       AdjCell(2,:) = DEGENERATE_FACE
    end where
    
    ! Finally, have to send AdjCell back to mesh structure.
    ! The destination cell is in HomeCell(2,:)
    GS_Trace => PGSLib_Setup_Trace (INDEX        = HomeCell(2,:), &
                                    SIZE_OF_DEST = SIZE(Mesh,1)   ) 

    ALLOCATE (AFSup(SIZE(Mesh(1)%Ngbr_Face,1),PGSLib_Size_OF_Sup(GS_Trace)))
    ALLOCATE (ACSup(SIZE(Mesh(1)%Ngbr_Cell,1),PGSLib_Size_OF_Sup(GS_Trace)))

    ! Initialize
    AFSup = -1
    ACSup = -1
    do f = 1, allfaces
       if (HomeCell(2,f) > 0) then
          Mesh(HomeCell(2,f))%Ngbr_Face(HomeCell(1,f)) = AdjCell(1,f)
          Mesh(HomeCell(2,f))%Ngbr_Cell(HomeCell(1,f)) = AdjCell(2,f)
       else
          AFSup(HomeCell(1,f), - HomeCell(2,f)) = AdjCell(1,f)
          ACSup(HomeCell(1,f), - HomeCell(2,f)) = AdjCell(2,f)
       end if
    end do
    
    ALLOCATE (AFDup(SIZE(Mesh(1)%Ngbr_Face,1), PGSLib_Size_Of_Dup(GS_Trace)))
    ALLOCATE (ACDup(SIZE(Mesh(1)%Ngbr_Cell,1), PGSLib_Size_Of_Dup(GS_Trace)))

    AFDup = -1
    ACDup = -1

    AFDup = PGSLib_Scatter_Buffer(AFSup, GS_Trace)
    ACDup = PGSLib_Scatter_Buffer(ACSup, GS_Trace)

    ! Set the connectivity in the Mesh type.
    do c = 1, PGSLib_Size_Of_Dup(GS_Trace)
       do f = 1, SIZE(AFDup,1)
          if (AFDup(f,c) >= 0) then
             Mesh(PGSLib_Dup_Index(GS_Trace, c))%Ngbr_Face(f) = AFDup(f,c)
          end if
       end do
       do f = 1, SIZE(ACDup,1)
          if (ACDup(f,c) >= 0) then
             Mesh(PGSLib_Dup_Index(GS_Trace, c))%Ngbr_Cell(f) = ACDup(f,c)
             Mesh(GS_Trace%Duplicate_Indices(c))%Ngbr_Cell(f) = ACDup(f,c)
          end if
       end do
    end do

    ! Clean up the temporary stuff.
    call PGSLib_Deallocate_Trace (GS_Trace)
    DEALLOCATE (AFDup, ACDup, AFSup, ACSup)

    ! Count the number of degenerate faces.
    degenerate_points = PGSLib_GLOBAL_COUNT (COUNT((Faces == DEGENERATE_NODE),DIM=1) == 3)
    degenerate_lines  = PGSLib_GLOBAL_COUNT (COUNT((Faces == DEGENERATE_NODE),DIM=1) == 2)
    triangle_faces    = PGSLib_GLOBAL_COUNT (COUNT((Faces == DEGENERATE_NODE),DIM=1) == 1)
    quad_faces        = PGSLib_GLOBAL_COUNT (COUNT((Faces == DEGENERATE_NODE),DIM=1) == 0)

#ifdef DEBUG_MESH_GEN
    call TLS_info ('Diagnostics from CONNECTIVITY')
    write(message,'(a,i0)') '  Points    ', degenerate_points
    call TLS_info (message)
    write(message,'(a,i0)') '  Lines     ', degenerate_lines
    call TLS_info (message)
    write(message,'(a,i0)') '  Tri_Faces ', triangle_faces
    call TLS_info (message)
    write(message,'(a,i0)') '  Qud_Faces ', quad_faces
    call TLS_info (message)
#endif
    
    ! Count the number of external boundary faces.
    boundary_faces = 0
    do f = 1,nfc
       boundary_faces = boundary_faces + COUNT(Mesh%Ngbr_Cell(f) == 0)
    end do

    boundary_faces_tot = PGSLib_GLOBAL_SUM(boundary_faces)
    write (message, 6) boundary_faces_tot
6   format (3x,'There are ',i10,' external boundary faces.')
    call TLS_info (message)
    call TLS_info (' Mesh connectivity established.')
    call TLS_info ('')

  END SUBROUTINE CONNECTIVITY

  SUBROUTINE MESH_AXIS (axis, nsegs, Ncell, Coord, Ratio, nvx, mvx, Vertex)
    !=======================================================================
    ! Purpose(s):
    !
    !   Generate an axis using mesh segments for an orthogonal grid.
    !
    !=======================================================================
    ! Arguments
    character, intent(IN) :: axis
    integer, intent(IN)  :: nvx
    integer, intent(IN)  :: nsegs
    integer, intent(IN)  :: mvx
    integer,  dimension(nsegs),   intent(IN)  :: Ncell
    real(r8), dimension(nsegs+1), intent(IN)  :: Coord
    real(r8), dimension(nsegs),   intent(IN)  :: Ratio
    real(r8), dimension(mvx),     intent(OUT) :: Vertex

    ! Local Variables
    integer :: i, j, m, n, ncl
    real(r8) :: dbl_coord, expf, expq, h1, h2, rn, width
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Write segment information
    write (message, 1) axis
1   format (1x,a,'-direction mesh:')
    call TLS_info ('')
    call TLS_info (message)

    write (message, 2) axis, axis, axis
2   format (3x,'Seg',6x,a,7x,'Cells',3x,'Ratio', 5x,'First d',a,6x,'Last  d',a)
    call TLS_info ('')
    call TLS_info (message)

    ! Fatal i nsegs is negative
    if (nsegs <= 0) then
       write (message, 5) axis, nsegs
5      format('MESH_AXIS: ',a,'-axis number mesh segments =',i9,' is <= 0!')
       call TLS_panic (message)
    end if

    ! Number of cells
    ncl = 0

    ! First coordinate
    Vertex(1) = Coord(1)

    ! Loop over segments
    MESH_SEGMENTS: do i = 1, nsegs

       ! Width of next mesh segment
       width = Coord(i+1) - Coord(i)

       ! Fatal if width is negative
       if (width <= 0.0_r8) then
          write (message, 6) axis, i, width
6         format('MESH_AXIS: ',a,'-axis mesh width(',i2,') =',1pe12.5,' is <= 0!')
          call TLS_panic (message)
       end if

       ! Fatal if Ncells is negative
       if (Ncell(i) <= 0) then
          write (message, 7) axis, i, Ncell(i)
7         format ('MESH_AXIS: ',a,'-axis ncell(',i2,') =',i10,' <= 0!')
          call TLS_panic (message)
       end if

       EXPANSION_FACTOR: if (Ratio(i) > 0.0_r8 .and. Ratio(i) /= 1.0_r8) then

          expf = Ratio(i)                ! Expansion factor
          rn = expf**Ncell(i)            ! Set r**n
          h1 = width*(expf-1.0_r8)/(rn-1.0_r8) ! First cell size
          expq = 1.0_r8/expf                ! Inverse expansion factor
          h2 = h1*rn*expq                ! Last cell size

       else if (Ratio(i) == 0.0_r8 .or. Ratio(i) == 1.0_r8) then

          h1 = width/real(Ncell(i))      ! Constant mesh
          h2 = h1
          expf = 1.0_r8                     ! No expansion
          expq = 1.0_r8

       else

          ! Fatal if Ratio is negative
          write (message, 8) axis, i, Ratio(i)
8         format ('MESH_AXIS: ',a,'-axis Ratio(',i2,') =',1pe12.5,' < 0  or  > 1!')
          call TLS_panic (message)

       end if EXPANSION_FACTOR

       m   = ncl + 2        ! Seg first interior vertex index
       ncl = ncl + Ncell(i) ! Seg last interior vertex index
       n   = (ncl + m) / 2  ! Seg midpoint index

       ! Write generation output
       write (message, 3) i, Coord(i), Ncell(i), Ratio(i), h1, h2
3      format (3x,i2,1pe14.6,i5,0pf11.7,2(1pe14.6))
       call TLS_info (message)

       dbl_coord = Coord(i)           ! Starting segment coordinate

       do j = m, n                    ! Generate first half going foward
          dbl_coord  = dbl_coord + h1
          Vertex(j)  = dbl_coord
          h1 = h1*expf
       end do

       dbl_coord     = Coord(i+1)     ! Force seg end coord to be exact
       Vertex(ncl+1) = Coord(i+1)

       do j = ncl, n+1, -1            ! Generate last half going backwards
          dbl_coord = dbl_coord - h2
          Vertex(j) = dbl_coord
          h2 = h2*expq
       end do

    end do MESH_SEGMENTS

    ! Write coordinate
    write(message,'(5x,1pe14.6," ----")') Coord(nsegs+1)
    call TLS_info (message)
    write(message,'(19x,i5)') ncl
    call TLS_info (message)

  END SUBROUTINE MESH_AXIS

  !-----------------------------------------------------------------------------

  FUNCTION COMPUTE_RCM (Mesh) RESULT(RCM2)
     !--------------------------------------------------------------------------
     ! Purpose: returns a permutation vector RCM such that
     !          Mesh_New_Order(RMC(i)) = Mesh_Old_Order(i)
     !
     !          reordering algorithm is Reverse Cuthill-Mckee
     !          "Iterative Methods for Sparse Linear Systems,"
     !          Yousef Saad, Pg 77
     !--------------------------------------------------------------------------
    use mesh_module,       only: MESH_CONNECTIVITY
    use parameter_module,  only: ncells
    use var_vector_module, only: SIZES, FLATTEN

    ! arguments
    type (MESH_CONNECTIVITY), dimension(ncells), intent(IN) :: Mesh
    integer,       dimension(ncells)             :: RCM

    ! local variables
    integer :: i, j, cellid, rcount, scount
    integer :: rmax, next, ni, num_neighbors
    integer, pointer, dimension(:) :: Ngbr_List
    integer, dimension(ncells) :: M
    integer, dimension(ncells) :: RCM2 

    ! these are working vectors which need be allowed to grow later
    integer, dimension(ncells) :: R
    integer, dimension(ncells) :: S
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize all vectors to -1

    RCM = -1
    M   = -1
    R   = -1
    S   = -1

    ! Setup initial level set information

    cellid = Next_Cell(M) ! starting node. Can be any node.
    R(1) = cellid;        ! level set R begins with one starting node 
                          ! or cellid[] to more general later
    rcount = 1           ! number of items currently in R
    M(cellid) = 1         ! mark first node as discoverd
    RCM(ncells) = cellid  ! first addition to permutation vector is starting node
    next = ncells - 1     ! next index for permutation vector RCM
    rmax = 0              ! for max number of members in R and S

    ni = Make_Local(-1, ncells) ! Initialize Make_Local function
       
    do

        if (next < 1) exit

        ! Reset next level set S
        scount = 0
        S      = -1

        ! SortR() add sort later

        ! Now loop through current level set R
        R_LOOP: do j = 1,rcount

            Ngbr_List => FLATTEN (Mesh(R(j))%Ngbr_Cells_All)
            ! number of neighbors for R(j)
            num_neighbors = SIZES(Mesh(R(j))%Ngbr_Cells_All)

            GET_NEIGHBORS: if (SIZES(Mesh(R(j))%Ngbr_Cells_All) > 0) then

                NEIGHBOR_LOOP: do i=1,num_neighbors

                    ! get node number ni for the ith neighbor
                    ni = Make_Local(Ngbr_List(i))
                    if (ni == NOT_LOCAL_INDEX) cycle NEIGHBOR_LOOP ! If not a local index, ignore
                    if (M(ni) > 0) cycle NEIGHBOR_LOOP
                    ! add unmarked neighbor to next level set
                    scount = scount + 1
                    S(scount) = ni
                    ! Mark the neighbor as discovered
                    M(ni) = 1
                    ! add this neighbor to the permutation vector
                    RCM(next) = ni
                    ! increment counter for permutation vector P
                    next = next - 1

                end do NEIGHBOR_LOOP

            end if GET_NEIGHBORS

        end do R_LOOP
        
        ! Now copy next level set S into current level set R
        rcount = scount
        R = S

        ! Need to find a new starting cell for non-connected regions
        if (rcount == 0) then

           cellid = Next_Cell(M) ! Starting node for this region. Can be any node.
           R(1) = cellid;        ! level set R begins with one starting node 
                                 ! or cellid[] to more general later
           rcount = 1            ! number of items currently in R
           M(cellid) = 1         ! mark first node as discoverd
           RCM(next) = cellid    ! first addition to permutation vector is starting node
           next = next - 1       ! next index for permutation vector RCM

        end if

    end do

    ! create inverse of permutation vector
    do j=1,ncells
        RCM2(RCM(j)) = j
    end do

  END FUNCTION COMPUTE_RCM

  !-----------------------------------------------------------------------------

  FUNCTION NEXT_CELL (M)
    !=======================================================================
    ! Purpose(s):
    !
    !=======================================================================

    ! Arguments
    integer, dimension(:) :: M

    ! Local Variables
    integer :: n, Next_Cell
    integer, save :: next_cell_ptr = 1

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Looking for next unmarked cell
    do n =  next_cell_ptr, SIZE(M)
       ! is cell not marked, if not this is the next cell
       if (M(n) < 0) then
          ! pointer for next time
          next_cell_ptr = n + 1
          ! set next_cell
          next_cell     = n
          ! break out of loop having found next_cell
          exit
       end if

    end do

  END FUNCTION NEXT_CELL

  FUNCTION MAKE_LOCAL (global_index, ncells)
    !=======================================================================
    ! Purpose(s):
    !   Determine if a index is local to this processor
    !   If it is, return the local index
    !   If it isn't, return NOT_LOCAL_INDEX
    !   This routine will get replaced
    !   This routine assumes that the mesh is ncells on each processor
    !   and the indices form a compact ordered set running from 1 to ncells_tot.
    !   The first call must supply ncells.
    !
    !=======================================================================
    use pgslib_module, only: PGSLib_SUM_PREFIX

    ! Arguments
    integer :: global_index
    integer, optional :: ncells

    ! Function Return
    integer :: MAKE_LOCAL

    ! Local variables
    integer, save                  :: lower, upper
    logical, save                  :: initialized = .FALSE.
    integer, pointer, dimension(:) :: Offset
    
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! If not initialized, set lower and upper
    if (.not. initialized) then
       ALLOCATE (Offset(ncells))
       Offset = 1
       Offset = PGSLib_SUM_PREFIX(Offset)
       Lower = Offset(1)
       Upper = Offset(ncells)
       initialized = .TRUE.
       DEALLOCATE (Offset)
    end if

    if ((lower <= global_index) .AND. (global_index <= upper)) then
       MAKE_LOCAL = global_index - lower + 1
    else
       MAKE_LOCAL = NOT_LOCAL_INDEX
    end if

  END FUNCTION MAKE_LOCAL
    
  FUNCTION TOUCHES_FACE (Face, vertex)
    !=======================================================================
    ! Purpose(s):  Tests if vertex is in Face.  Returns .TRUE.
    !              if so, .FALSE. otherwise.
    !
    !=======================================================================
    use parameter_module, only: nvf

    ! Arguments
    integer, dimension(nvf), intent(IN) :: Face
    integer, intent(IN) :: vertex

    ! Function return
    logical :: TOUCHES_FACE

    ! Local variables
    integer :: v

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    TOUCHES_FACE = .false.
    do v = 1, nvf
       if (vertex == Face(v)) then
          TOUCHES_FACE = .true.
          exit
       end if
    end do

  END FUNCTION TOUCHES_FACE

END MODULE MESH_GEN_MODULE
