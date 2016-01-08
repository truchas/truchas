!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE MESH_INPUT_MODULE
  !=======================================================================
  ! Purpose:
  !
  !   Define procedures for the input of mesh parameters
  !
  !   Public Interface:
  !
  !     * call MESH_INPUT ()
  !
  !       sets defaults, reads, checks, and broadcasts input variables
  !       in the MESH namelist
  !
  !     * call MESH_READ (Mesh, Vertex)
  !
  !       reads various format mesh files
  !
  !     * call MESH_SIZES ()
  !
  !       set the problem mesh dimensions
  !
  ! Contains: MESH_INPUT
  !           MESH_READ
  !           MESH_READ_SIZE
  !           MESH_SIZES
  !           MESH_CHECK
  !           MESH_DEFAULT
  !           MESH_INPUT_PARALLEL
  !
  ! Author(s): Jerry S. Brock, LANL (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL (dbk@lanl.gov)
  !            Bryan Lally, LANL (lally@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only: mseg, ndim, mbody
  use truchas_logging_services
  use parallel_communication
  implicit none
  private

  public :: MESH_INPUT, MESH_READ, MESH_SIZES, MESH_READ_SIDE_SETS

  ! Magic values used to detect variables not initialized by input
  character, parameter :: NULL_C = char(0)
  integer, parameter :: NULL_I = huge(1)
  real(r8), parameter :: NULL_R = huge(1.0_r8)

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! MESH namelist input variables
  character(120), public, save :: mesh_file
  character(120), public, save :: mesh_file_format
  real(r8), public, save :: coordinate_scale_factor
  integer,  public, save, dimension(ndim,mseg) :: Ncell
  real(r8), public, save, dimension(ndim) :: Fuzz
  real(r8), public, save :: Heps
  real(r8), public, save, dimension(ndim,mseg) :: Ratio
  real(r8), public, save, dimension(ndim,mseg+1) :: Coord
  integer,  public, save, dimension(mbody) :: gap_element_blocks
  integer,  public, save, dimension(127) :: interface_side_sets
  integer,  public, save :: exodus_block_modulus

  ! Derived MESH namelist quantities
  integer, public, save, dimension(ndim) :: Nseg
  character, public, save, dimension(3) :: Coord_label = (/ 'X','Y','Z' /)

  ! RCM renumbering flag
  logical, public, save :: use_RCM
  
  ! Logical unit the mesh file is opened on
  integer, save :: msh_lun

CONTAINS

  SUBROUTINE MESH_INPUT (lun)
    !=======================================================================
    ! Purpose:
    !
    !   Read the MESH namelist
    !
    !=======================================================================
    use input_utilities,        only: seek_to_namelist, NULL_C
    use string_utilities,       only: i_to_c
    use truchas_env,            only: input_dir
    use parallel_info_module,   only: p_info
    use pgslib_module,          only: pgslib_bcast
    use parameter_module,       only: mseg, ndim, Nx_tot
    use restart_variables,      only: restart
    use parameter_module,       only: ncells_tot, nnodes_tot
    use mesh_gen_data,          only: partitions_total, set_generated_mesh

    integer, intent(in) :: lun

    ! local variables
    logical :: fatal, found
    integer :: ios
    integer :: n
    integer :: total_cells
    character(256) :: message

    ! mesh namelist specification
    Namelist /MESH/ Ncell, Coord, Ratio, Fuzz, Heps, mesh_file, mesh_file_format, &
                    coordinate_scale_factor, use_RCM, &
                    gap_element_blocks, interface_side_sets, exodus_block_modulus

    call TLS_info ('')
    call TLS_info ('Reading MESH Namelist ...')

    ! Set defaults.
    call MESH_DEFAULT ()

    ! Input is done in IO PE, results are broadcast to other PEs.
    if (p_info%IOP) then
      rewind lun
      call seek_to_namelist (lun, 'MESH', found, iostat=ios)
    end if
    
    call pgslib_bcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading input file')

    call pgslib_bcast (found)
    if (.not.found) call TLS_fatal ('MESH namelist not found')

    ! Read the namelist.
    if (p_info%IOP) read(lun, nml=mesh, iostat=ios)
    call pgslib_bcast (ios)
    if (ios /= 0) call TLS_fatal ('error reading MESH namelist; iostat=' // i_to_c(ios))

    ! Broadcast total mesh size variables.
    call MESH_INPUT_PARALLEL ()

    ! all is OK, let's set the partition numbers
    partitions_total = p_info%npe

    ! Read in the mesh from the mesh file if this is not a restart.
    ! Set flag for indicating whether we can cartesian partition
    ! The default .FALSE.
    call set_Generated_Mesh(.FALSE.)

    if (mesh_file /= NULL_C) then
       mesh_file = adjustl(trim(mesh_file))
       if (mesh_file(1:1) /= '/' .and. mesh_file(1:1) /= '~') then
          mesh_file = trim(input_dir) // trim(mesh_file)
       end if
       if (mesh_file_format == NULL_C) mesh_file_format = 'ExodusII'
    end if
    READ_THE_MESH: if (mesh_file /= NULL_C .and. .not.(restart)) then
       ! First read the sizes and broadcast them.
       call MESH_READ_SIZE ()

       if (exodus_block_modulus < 0) &
           call TLS_fatal ('negative value specified for Exodus_Block_Modulus')

       if (nnodes_tot <= 0 .or. ncells_tot <= 0) &
           call TLS_fatal ('number of nodes and/or cells in ' // trim(mesh_file) // ' is <= 0')

       write(message,'(9x,a)') 'Opened mesh file ' // trim(mesh_file)
       call TLS_info (message)
       write(message,'(11x,3(a,:,i0))') 'containing ', nnodes_tot, ' nodes and ', ncells_tot, ' cells'
       call TLS_info (message)

    end if READ_THE_MESH

    ! If we made the mesh, check it and announce.
    GENERATE_MESH: if (mesh_file == NULL_C) then

       ! If we are generating a mesh, okay to cartesian partition
       call set_Generated_Mesh( .TRUE. )

       ! we built the mesh - Check for user input errors
       total_cells = 1
       do n = 1,ndim
          fatal = .false.
          call MESH_CHECK(Coord_label(n), mseg, Coord(n,:), Fuzz(n), Heps, &
                          Ratio(n,:), Ncell(n,:), fatal, Nx_tot(n), Nseg(n))
          total_cells = total_cells*Nx_tot(n)
          if (fatal) call TLS_fatal ('terminating execution due to previous input errors')
       end do

       if (ndim == 2) then
          write (message, 15) (Nx_tot(n), n = 1,ndim), total_cells
15        format (9x,'Mesh Size: ',i3,' x',i3,' = ',i7,' cells')
          if (.not.restart) call TLS_info (message)
       else if (ndim == 3) then
          write (message, 20) (Nx_tot(n), n = 1,ndim), total_cells
20        format (9x,'Mesh Size: ',i3,' x',i3,' x',i3,' = ',i7,' cells')
          if (.not.restart) call TLS_info (message)
       end if

    end if GENERATE_MESH

  END SUBROUTINE MESH_INPUT

  SUBROUTINE MESH_READ (Mesh, Vertex)
    !=======================================================================
    ! Purpose:
    !
    !   Read vertex (node) coordinates and element <=> node connectivity
    !   from the mesh file
    !
    !=======================================================================
    use bc_data_module,       only: Mesh_Face_Set_Tot
    use mesh_module,          only: MESH_CONNECTIVITY, VERTEX_DATA
    use parameter_module,     only: ncells_tot, ndim, nnodes_tot, nfc, nssets
    use exodus_mesh_type
    use exodus_mesh_io, only: read_exodus_mesh
    use string_utilities, only: i_to_c

    ! Arguments
    type(MESH_CONNECTIVITY), dimension(ncells_tot), intent(INOUT) :: Mesh
    type(VERTEX_DATA),       dimension(nnodes_tot), intent(INOUT) :: Vertex

    ! Local Variables
    logical :: hybrid_mesh, tet, prism, pyramid
    integer :: i, j, lc, n, mem_stat, status
    integer :: nc_count, ntmp, new_id
    integer, pointer  :: ctemp(:), ftemp(:)
    type(exodus_mesh), target :: exo_mesh
    character(:), allocatable :: errmsg, msg

    call TLS_info ('')
    call TLS_info ('Reading ExodusII mesh file ' // trim(mesh_file) // ' ...')

    ! Read the Exodus II mesh file.
    call read_exodus_mesh (mesh_file, exo_mesh, status, errmsg)
    if (status /= 0) then
      call TLS_panic ('READ_MESH: error reading ExodusII mesh file: ' // errmsg)
    end if

    ! Copy the node coordinates.
    do n = 1, nnodes_tot
      Vertex(n)%Coord = exo_mesh%coord(:,n)
    end do

    ! Copy the element connectivity, organized by blocks.
    nc_count = 0
    do n = 1, exo_mesh%num_eblk
       ! Overwrite the block ID with its value modulo EXODUS_BLOCK_MODULUS.
       if (exodus_block_modulus > 0) then
          associate (id => exo_mesh%eblk(n)%id)
             new_id = modulo(exo_mesh%eblk(n)%id, exodus_block_modulus)
             if (new_id /= id) then
                msg = ' Element block ' // i_to_c(id) // ' merged with block ' // i_to_c(new_id)
                call TLS_info (msg, TLS_VERB_NORMAL)
                id = new_id
             end if
          end associate
       end if
       lc = size(exo_mesh%eblk(n)%connect,dim=1)
       do i = 1, exo_mesh%eblk(n)%num_elem
          nc_count = nc_count + 1
          Mesh(nc_count)%Ngbr_vrtx(1:lc) = exo_mesh%eblk(n)%connect(:,i)
          Mesh(nc_count)%CBlockID = exo_mesh%eblk(n)%ID
          ! Fix node order for prisms (not sure if this matters except for gap elements
          ! with triangular faces.  Swap vertices 2 and 4, 5 and 3.  This assumes that 
          ! all elements in a block have the same number of non-zero nodes
          if (lc == 6) then
             ntmp = Mesh(nc_count)%Ngbr_vrtx(2)
             Mesh(nc_count)%Ngbr_vrtx(2) = Mesh(nc_count)%Ngbr_vrtx(4)
             Mesh(nc_count)%Ngbr_vrtx(4) = ntmp
             ntmp = Mesh(nc_count)%Ngbr_vrtx(3)
             Mesh(nc_count)%Ngbr_vrtx(3) = Mesh(nc_count)%Ngbr_vrtx(5)
             Mesh(nc_count)%Ngbr_vrtx(5) = ntmp
          end if
       end do
    end do

    ! Allocate the side set arrays.
    nssets = exo_mesh%num_sset
    if (nssets > 0) then
       ALLOCATE(Mesh_Face_Set_Tot(nssets,nfc,ncells_tot), STAT=mem_stat)
       if (mem_stat /= 0) call TLS_panic ('MESH_READ: could not allocate Mesh_Face_Set_Tot')
       Mesh_Face_Set_Tot = 0
    end if

    ! Copy the side set data, if it is present.
    if (nssets > 0) then
       do n = 1, nssets
          ctemp => exo_mesh%sset(n)%elem
          ftemp => exo_mesh%sset(n)%face
          do i = 1, exo_mesh%sset(n)%num_side
             ! Convert face numbers from Exodus (Cubit?) to Truchas scheme
             ! Find type of cell that has the boundary/interface face

             tet = .false.; prism = .false.; pyramid = .false.

             j = ctemp(i)
             if (Mesh(j)%Ngbr_Vrtx(7) <= 0 .and. Mesh(j)%Ngbr_Vrtx(8) <= 0) then
                if (Mesh(j)%Ngbr_Vrtx(6) <= 0) then
                   tet = Mesh(j)%Ngbr_Vrtx(5) <= 0
                   pyramid = .not. tet
                else
                   prism = .true.
                end if
             end if

             if (tet) then
                select case(ftemp(i))
                case(1)
                   ftemp(i) = 4
                case(2)
                   ftemp(i) = 1
!                      case(3)
!                         ftemp(i) = 3
                case(4)
                   ftemp(i) = 5
                end select
             else if (prism) then
                select case(ftemp(i))
                case(1)
                   ftemp(i) = 5
                case(2)
                   ftemp(i) = 1
                case(3)
                   ftemp(i) = 2
                case(4)
                   ftemp(i) = 3
                case(5)
                   ftemp(i) = 4
                end select
             else if (pyramid) then
                select case (ftemp(i))
                case(1)
                  ftemp(i) = 2
                case(2)
                  ftemp(i) = 4
                case(3)
                  ftemp(i) = 1
                case(4)
                  ftemp(i) = 3
                end select
             else
                select case(ftemp(i))
                case(1)
                   ftemp(i) = 2
                case(2)
                   ftemp(i) = 4
                case(3)
                   ftemp(i) = 1
                case(4)
                   ftemp(i) = 3
                end select
             end if

             Mesh_Face_Set_Tot(n,ftemp(i),ctemp(i)) = exo_mesh%sset(n)%ID
          end do

       end do
    end if

    ! Done; announce.
    call TLS_info (' Closed ExodusII mesh file ' // trim(mesh_file))

    ! Check for nonhex elements if we're in 3-D. if they're present,
    ! then set degenerate nodes.
    if (ndim == 3) then

       ! Is this a hybrid mesh?
       hybrid_mesh = ANY(Mesh%Ngbr_Vrtx(7) <= 0 .and. Mesh%Ngbr_Vrtx(8) <= 0)

       if (hybrid_mesh) then

          ELEMENT_CHECK: do i = 1,ncells_tot
             ! Is this a non-hex cell?
             tet = .false.; prism = .false.; pyramid = .false.
             if (Mesh(i)%Ngbr_Vrtx(7) <= 0 .and. Mesh(i)%Ngbr_Vrtx(8) <= 0) then
                if (Mesh(i)%Ngbr_Vrtx(6) <= 0) then
                   tet = Mesh(i)%Ngbr_Vrtx(5) <= 0
                   pyramid = .not. tet
                else
                   prism = .true.
                end if
             end if

             ! if this cell is not a hex, reassign the vertices and set degenerate faces.
             if (tet) then
                ! Tet: Vertices 5-8 are the 4th vertex read in, 4 is
                !      the 3rd, 3 is the 2nd, and 2 is the 1st. Faces
                !      2 and 6 are degenerate (zero area).
                Mesh(i)%Ngbr_Vrtx(5:8) = Mesh(i)%Ngbr_Vrtx(4)
                Mesh(i)%Ngbr_Vrtx(4)    = Mesh(i)%Ngbr_Vrtx(3)
                Mesh(i)%Ngbr_Vrtx(3)    = Mesh(i)%Ngbr_Vrtx(2)
                Mesh(i)%Ngbr_Vrtx(2)    = Mesh(i)%Ngbr_Vrtx(1)
             else if (pyramid) then
                ! Pyramid: Vertices 6-8 are equivalent to the 5th vertex
                !          read in. Face 6 is degenerate.
                Mesh(i)%Ngbr_Vrtx(6:8) = Mesh(i)%Ngbr_Vrtx(5)
             else if (prism) then
                ! Prism: Vertex 7 is equivalent to the 6th vertex
                !        read in and vertex 8 is equivalent to the
                !        5th vertex read in. Face 6 is degenerate.
                Mesh(i)%Ngbr_Vrtx(8) = Mesh(i)%Ngbr_Vrtx(5)
                Mesh(i)%Ngbr_Vrtx(7) = Mesh(i)%Ngbr_Vrtx(6)
             end if
          end do ELEMENT_CHECK

       end if
    end if

  END SUBROUTINE MESH_READ

  SUBROUTINE MESH_SIZES ()
    !=======================================================================
    ! Purpose:
    !
    !   Set the relevant problem mesh dimensions
    !
    !=======================================================================
    use mesh_decomposition_module, only: MESH_BLOCK_DECOMPOSITION, &
                                         MESH_DOMAIN_SIZES
    use parameter_module,          only: boundary_faces, boundary_faces_tot, &
                                         Mx, Mx_tot, ncells, ncells_tot,     &
                                         nnodes, nnodes_tot, Nx, Nx_tot,     &
                                         ndim
    use parallel_info_module,      only: p_info
    use pgslib_module,             only: PGSLib_Global_ALL
    use restart_variables,         only: restart, restart_ncells, restart_nnodes

    ! Local Variables
    logical :: read_from_file, all_zero
    integer :: n

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! File read flag.
    all_zero = ALL(NX_tot == 0)
    read_from_file = restart .or. PGSLib_Global_ALL(all_zero)

    ! Set max dimensions for a single block mesh
    do n = 1,ndim
       Mx_tot(n) = Nx_tot(n) + 1
    end do

    ! Set other mesh parameters
    MESH_PARAMETERS: if (read_from_file) then

       ! if this is a restart, get the total mesh sizes.
       if (restart) then
          ncells_tot = restart_ncells
          nnodes_tot = restart_nnodes
       end if

       ! if mesh is input from a file, then we need to know
       ! the sizes of each of the domains.  This is determined by
       ! a combination of: default choices, user input, mesh file format
       call MESH_DOMAIN_SIZES (msh_lun, ncells_tot, nnodes_tot, ncells, nnodes)

       do n = 1,ndim
          Mx(n) = Nx(n) + 1
       end do

    else

       ! Since we are generating mesh on the fly we can decompose it on the fly too.
       ! We need to find out what the mesh decomposition is.
       ! This routine works only for the simple brick mesh that we generate.
       ! The decomposition is BLOCK on the third dimension.
       call MESH_BLOCK_DECOMPOSITION (Nx_tot, Nx)

       ! The distribution of the vertices follows the distribution of the cells.
       !.not.ce that only the last pe has an additional plane of boundary vertices.
       if (p_info%IsParallel) then
          if (p_info%thisPE < p_info%nPE)  then
             ! No additional plane of vertices for internal boundaries
             do n = 1,ndim-1
                Mx(n) = Nx(n) + 1
             end do
             Mx(ndim) = Nx(ndim)
          else
             ! Last PE gets an extra plane for external boundary
             do n = 1,ndim
                Mx(n) = Nx(n) + 1
             end do
          end if
       else
          ! This is a serial run.
          do n = 1,ndim
             Mx(n) = Nx(n) + 1
          end do
       end if

       ! Total number of cells and vertices.
       ncells = 1; ncells_tot = 1; nnodes = 1; nnodes_tot = 1
       do n = 1,ndim
          ncells     = ncells*Nx(n)
          ncells_tot = ncells_tot*Nx_tot(n)
          nnodes     = nnodes*Mx(n)
          nnodes_tot = nnodes_tot*Mx_tot(n)
       end do

       ! Total number of boundary faces.
       boundary_faces = 1; boundary_faces_tot = 1
       do n = 1,ndim
          boundary_faces     = boundary_faces*(Mx(n) + 1)
          boundary_faces_tot = boundary_faces_tot*(Mx_tot(n) + 1)
       end do
       boundary_faces     = boundary_faces - ncells
       boundary_faces_tot = boundary_faces_tot - ncells_tot

    end if MESH_PARAMETERS

  END SUBROUTINE MESH_SIZES

  SUBROUTINE MESH_CHECK (axis, mseg, Coord, fuzz, Heps, Ratio, Ncell, fatal, n, nseg)
    !=======================================================================
    ! Purpose:
    !
    !   Check MESH namelist parameters for orthogonal mesh generation.
    !   return n, the total number of cells for this axis.
    !=======================================================================
    ! Arguments
    character, intent(IN) :: axis
    integer, intent(IN) :: mseg
    real(r8), dimension(mseg+1), intent(IN) :: Coord
    real(r8), intent(IN) :: fuzz
    real(r8), intent(IN) :: Heps
    real(r8), dimension(mseg), intent(INOUT) :: Ratio
    integer, dimension(mseg), intent(INOUT) :: Ncell
    logical, intent(INOUT) :: fatal
    integer, intent(OUT) :: n, nseg

    ! Local Variables
    integer :: i, msegm1
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Basic constants
    n      = 0
    nseg   = 0
    msegm1 = mseg - 1

    ! Check mesh segments
    SEGMENT_LOOP: do i = 1, msegm1

       ! Check number of cells for this segment.
       CELL_NUMBER_CHECK: if (Ncell(i) < 0) then

          ! Negative number of cells; no good.
          write (message,10) axis, i, Ncell(i)
10        format(a,'-axis ncell(',i3,') =',i10,' is negative!')
          call TLS_error (message)
          fatal = .true.

       ! Positive number of cells. . .
       else if (Ncell(i) > 0) then

          n    = n + Ncell(i)
          nseg = nseg + 1

          ! Non-contiguous segment.
          if (nseg /= i) then
             write (message,15) axis, i
15           format(a,'-axis mesh segment ',i0,' is not contiguous with other mesh segments!')
             call TLS_error (message)
             fatal = .true.
          end if

          ! This seems to be a valid segment, now check to be sure the
          ! user input all the coordinates.
          if (Coord(i) == NULL_R .or. Coord(i+1) == NULL_R) then
             call TLS_error ('COORD data for ' // axis // '-axis not fully specified')
             fatal = .true.
          endif

          ! Non-monotonic coordinates.
          if (Coord(i) >= Coord(i+1)) then
             write (message,20) axis,i,Coord(i),i+1,Coord(i+1)
20           format(a,'-axis coord(',i0,')=',1pe13.5,' Equal or greater than coord(',i0,')=',1pe13.5)
             call TLS_error (message)
             fatal = .true.
          end if

          ! Zero expansion gives uniform zoning.
          if (Ratio(i) == 0.0_r8) Ratio(i) = 1.0_r8

          ! Outside allowed bounds for Ratio.
          if (Ratio(i) < 0.8 .or. Ratio(i) > 1.2) then
             write (message,25) axis, i!, Ratio(i)
25           format(a,'-axis ratio(',i3,') < 0.8 or > 1.2 !')
             call TLS_warn (message)
          end if

       end if CELL_NUMBER_CHECK

    end do SEGMENT_LOOP

    ! Check mseg, Fuzz, and Heps.
    ! Fuzz is greater than one; not allowed.
    if (Heps > 0 .and. fuzz > 0) then
       write (message,30) Heps,fuzz
30     format('Heps = ',1pe12.5,' > 0 and Fuzz = ',1pe12.5,' > 0. ', &
              'Only one value (fuzz or Heps) can be nonzero. ')
       call TLS_error (message)
       fatal = .true.
    end if

    if (fuzz > 0.95) then
       write (message,35) axis, fuzz
35     format(a,'-axis randomization MAY be too large! Fuzz = ',1pe12.5, &
              '.  Fuzz greater than 0.95 may give negative volumes!')
       call TLS_warn (message)
       fatal = .false.
    end if

    if (nseg == msegm1 .and. Ncell(mseg) > 0) then
       write (message,40) axis, msegm1
40     format('number of ',a,'-axis mesh segments exceeds internal limit MSEG=',i0)
       call TLS_error (message)
       fatal = .true.
    end if

    ! Fuzz is greater than one; not allowed.
    if (Heps > 1.0_r8) then
       call TLS_error ('Heps must be <= 1')
       fatal = .true.
    end if

    ! Fuzz is less than zero; not allowed.
    if (Heps < 0.0_r8) then
       call TLS_error ('Heps must be >= 0')
       fatal = .true.
    end if

    ! Fuzz is greater than one; not allowed.
    if (fuzz > 1.0_r8) then
       write (message,55) axis, fuzz
55     format(a,'-axis randomization is too large: Fuzz = ',1pe12.5, &
              '.  Fuzz must be less than one.')
       call TLS_error (message)
       fatal = .true.
    end if


    ! Fuzz is negative; not allowed.
    if (fuzz < 0.0_r8) then
       write (message,60) axis, fuzz
60     format(a,'-axis randomization is negative: Fuzz = ',1pe12.5, &
              '.  Fuzz should be positive.')
       call TLS_error (message)
       fatal = .true.
    end if

    ! Make sure coordinate scale factor is positive.
    if (coordinate_scale_factor <= 0.0_r8) then
       write (message,65) Coordinate_Scale_Factor
65     format(3(1pe12.5,1x),'is not a valid coordinate scale factor!')
       call TLS_error (message)
       fatal = .true.
    end if

  END SUBROUTINE MESH_CHECK

  SUBROUTINE MESH_DEFAULT ()
    !=======================================================================
    ! Purpose:
    !
    !   Default MESH namelist
    !=======================================================================
    use input_utilities,  only: NULL_C
    use parameter_module, only: boundary_faces, ncells,   &
                                nfaces, nfc, nnodes, nvc, &
                                nvf, nfv

    ! local variables
    integer :: f, i, j
    integer :: v
    integer :: vertex_number

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Mesh array sizes
    boundary_faces = 1 ! Number of boundary faces
    ncells         = 1 ! Number of cells
    nfaces         = 1 ! Number of unique faces
    nnodes         = 1 ! Number of nodes

    ! Number of cells for each mesh segment.
    Ncell = 0

    ! Fuzzing parameters
    Fuzz = 0.0_r8
    
    Heps = 0.0_r8
    Ratio = 0.0_r8

    ! Mesh file name
    mesh_file = NULL_C

    ! Mesh file format
    mesh_file_format = NULL_C

    ! Vertex coordinate scale factor.
    coordinate_scale_factor = 1.0_r8

    ! RCM off by default
    use_RCM = .false.

    ! No gap elements or internal interfaces
    gap_element_blocks = 0
    interface_side_sets = 0

    ! Modulus for determining congruent element block IDs.
    exodus_block_modulus = 10000

   ! Initialize the coordinates for in-line generation to NULL_R.
   ! This permits error-trapping above to catch situations where the 
   ! user failed to specify a coordinate.
   do  j = 1, mseg+1
       do i = 1, ndim
          Coord(i,j) = NULL_R
       end do
    end do

  END SUBROUTINE MESH_DEFAULT

  SUBROUTINE MESH_INPUT_PARALLEL ()
    !======================================================================
    ! Purpose:
    !
    !   broadcast and distribute MESH namelist values from IO_ROOT_PE to
    !   all PE's in the system
    !======================================================================
    use parallel_info_module, only: p_info
    use pgslib_module,        only: PGSLIB_BCAST

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Broadcast Data
    BROADCAST_VARIABLES: if (.not. p_info%UseGlobalServices) then
       call PGSLib_BCast (Coord)
       call PGSLib_BCast (Fuzz)
       call PGSLib_BCast (Ncell)
       call PGSLib_BCast (Ratio)
       call PGSLib_BCast (mesh_file)
       call PGSLib_BCast (mesh_file_format)
       call PGSLib_BCast (coordinate_scale_factor)
       call PGSLib_BCast (use_RCM)
       call PGSLib_BCast (gap_element_blocks)
       call PGSLib_BCast (interface_side_sets)
       call PGSLib_BCast (exodus_block_modulus)
    end if BROADCAST_VARIABLES

  END SUBROUTINE MESH_INPUT_PARALLEL

  subroutine mesh_read_size ()

    use parameter_module, only: ncells_tot, nnodes_tot
    use exodus_truchas_hack, only: read_exodus_mesh_size
    use string_utilities, only: raise_case

    logical :: file_exists
    integer :: stat

    !! Verify the mesh file exists.
    if (is_IOP) inquire(file=trim(mesh_file),exist=file_exists)
    call broadcast (file_exists)
    if (.not.file_exists) call TLS_fatal ('mesh file "' // trim(mesh_file) // '" not found')

    !! The default and only supported format is ExodusII.
    if (mesh_file_format == NULL_C) mesh_file_format = 'EXODUSII'
    mesh_file_format = raise_case(mesh_file_format)
    if (mesh_file_format /= 'EXODUSII') then
      call TLS_fatal ('Unknown MESH_FILE_FORMAT value: "' // trim(mesh_file_format) // '"')
    end if
    
    !! Read the number of nodes and cells from the file and broadcast.
    if (is_IOP) call read_exodus_mesh_size (mesh_file, nnodes_tot, ncells_tot, stat)
    call broadcast (stat)
    if (stat /= 0) call TLS_fatal ('Error reading ExodusII file "' // trim(mesh_file) // '"')
    call broadcast (ncells_tot)
    call broadcast (nnodes_tot)

  end subroutine mesh_read_size
  
  SUBROUTINE MESH_READ_SIDE_SETS ()
  
    use parameter_module, only: nssets
    use parallel_info_module, only: p_info
    use bc_data_module, only: Mesh_Face_Set_Tot
    use exodus_truchas_hack, only: read_exodus_side_sets
    use pgslib_module,only: PGSLib_BCast
    use string_utilities, only: raise_case
    
    logical :: fatal
    integer :: status
    character(128) :: Fatal_Error_String
  
    if (raise_case(mesh_file_format) == "EXODUSII") then
    
      fatal = .false.
    
      if (p_info%IOP) then
        if (associated(mesh_face_set_tot)) then
          Fatal_Error_String = 'mesh_face_set_tot already defined'
          fatal = .true.
        end if
        call read_exodus_side_sets (mesh_file, mesh_face_set_tot, status)
        if (status /= 0) then
          Fatal_Error_String = 'error reading side sets from ExodusII file ' // trim(mesh_file)
          fatal = .true.
        else
          nssets = 0
          if (associated(mesh_face_set_tot)) nssets = size(mesh_face_set_tot,dim=1)
        end if
      end if
      
      call TLS_fatal_if_any (fatal, Fatal_Error_String)
      
      !! Broadcast the number of side sets;
      !! the side set data is distributed elsewhere.
      call PGSLib_BCast (nssets)
      
    end if
  
  END SUBROUTINE MESH_READ_SIDE_SETS
    
END MODULE MESH_INPUT_MODULE
