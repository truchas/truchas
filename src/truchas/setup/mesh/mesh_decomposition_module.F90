!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE MESH_DECOMPOSITION_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Determine mesh decomposition. This module contains two routines.
  !
  !   Mesh_Domain_Decomposition is called if the mesh is read in from a file
  !     in that case, the decomposition is simply dividing the data into
  !     segments.  The size of each segment depends on other input data 
  !     (eg if the mesh has been domain decomposed, use one domain
  !      per processor unless user instructs otherwise).
  !
  !   Mesh_Block_Decomposition is called if the mesh is generated at run-time.
  !     in this case the mesh is a rectangular block, and this routine
  !     does a block decomposition of the block to the processors.
  !
  !   NOTE: These routines only determine array extents, they do not
  !     actually distribute the data. That is done in distribute_module.
  !
  ! Contains: BLOCK_PARAMETERS
  !           MESH_BLOCK_DECOMPOSITION
  !           MESH_DOMAIN_SIZES
  !
  ! Author(s): Robert Ferrell (CPCA, Ltd., ferrell@cpca.com)
  !
  !=======================================================================
  use truchas_logging_services
  implicit none
  private

  public :: MESH_BLOCK_DECOMPOSITION, MESH_DOMAIN_SIZES

CONTAINS

  SUBROUTINE BLOCK_PARAMETERS (max_block, min_block, critPE, N_Tot)
    !=======================================================================
    ! Purpose(s):
    !   For a given N_Tot, determine block distribution for N_tot.
    !   PE's <= critPE have block size max_block.
    !   PE's > critPT have blcok size = min_block.
    !   The algorithm used here has (max_block - min_block) = 1.
    !=======================================================================
    use parallel_info_module, only: p_info

    ! Argument List
    integer, intent(IN)  :: N_Tot
    integer, intent(OUT) :: critPE, max_block, min_block

    ! Local Variables
    integer :: nPE

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Number of Processors
    nPE = p_info%nPE

    ! First find the average number of data items in a block,
    ! then truncate to an integer, to get min_block.
    min_block = FLOOR((REAL(N_tot)/REAL(nPE)))

    ! Max_block is one larger than min_block.
    max_block = min_block + 1
    
    ! critPE is the cutoff PE
    ! (Counting of PE's is one based at F90 level)
    critPE = (N_Tot - min_block*nPE) 

  END SUBROUTINE BLOCK_PARAMETERS

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE MESH_BLOCK_DECOMPOSITION (Nx_tot, Nx)
    !=======================================================================
    ! Purpose(s):
    !   Determine decomposition of a block mesh onto nPEs.
    !   Input:  nx_tot, ny_tot, nz_tot  - these are global mesh sizes.
    !   Output: nx, ny, nz              - these are local mesh sizes.
    !           (nx and ny are copied from nx_tot and ny_tot)
    !           (nz_tot is distributed as evenly as possible over nPEs)
    !=======================================================================
    use parallel_info_module, only: p_info
    use parameter_module,     only: ndim
    use pgslib_module,        only: PGSLib_BCAST, PGSLib_COLLATE

    ! Argument List
    integer, dimension(ndim), intent(IN)  :: Nx_tot
    integer, dimension(ndim), intent(OUT) :: Nx

    ! Local Variables
    integer :: nPE, max_block, min_block, critPE, pe, n
    integer, allocatable, dimension(:,:) :: Nx_all
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Number of Processors
    nPE = p_info%nPE

    LOCAL_MESH_SIZE: if (p_info%IsParallel) then

       call BLOCK_PARAMETERS (max_block, min_block, critPE, Nx_tot(ndim))
       if (p_info%thisPE <= critPE) then
          Nx(ndim) = max_block
       else
          Nx(ndim) = min_block
       end if

    else

       Nx(ndim) = Nx_tot(ndim)

    end if LOCAL_MESH_SIZE

    do n = 1,ndim-1
       Nx(n) = Nx_tot(n)
    end do

    ! allocate temporary arrays to hold collated data
    ALLOCATE (Nx_all(ndim,nPE))

    do n = 1,ndim
       call PGSLib_COLLATE (Nx_all(n,:), Nx(n))
       call PGSLib_BCAST (Nx_all(n,:))
    end do

    do pe = 1,nPE
       if (ndim == 2) then
          write (message, 10) pe, (Nx_all(n,pe), n = 1,ndim)
10        format (1x,' On pe ',i8, ' nx = ',i8, ' ny = ',i8)
       else if (ndim == 3) then
          write (message, 11) pe, (Nx_all(n,pe), n = 1,ndim)
11        format (1x,' On pe ',i8, ' nx = ',i8, ' ny = ',i8, ' nz = ',i8)
       end if
    end do

    if (ndim == 2) then
       write (message,15) (Nx(n), n = 1,ndim)
15     format (1x,'nx = ',i8,', ny = ',i8)
    else if (ndim == 3) then
       write (message,16) (Nx(n), n = 1,ndim)
16     format (1x,'nx = ',i8,', ny = ',i8,', nz = ',i8)
    end if

    if (ndim == 2) then
       write (message,20) (Nx_tot(n), n = 1,ndim)
20     format (1x,'nx_tot = ',i8,', ny_tot = ',i8)
    else if (ndim == 3) then
       write (message,21) (Nx_tot(n), n = 1,ndim)
21     format (1x,'nx_tot = ',i8,', ny_tot = ',i8,', nz_tot = ',i8)
    end if

    ! Deallocate temporaries
    DEALLOCATE (Nx_all)

  END SUBROUTINE MESH_BLOCK_DECOMPOSITION

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE MESH_DOMAIN_SIZES (msh_lun, ncells_tot, nnodes_tot, ncells, nnodes)
    !=======================================================================
    ! Purpose(s):
    !   Figure out how the mesh is decomposed, and determine how to assign 
    !   domains to processors. Determine size of domain on each processor.
    !=======================================================================
    use mesh_distribute_module, only: number_domains, Ncells_List, Nnodes_List
    use parallel_info_module,   only: p_info
    use pgslib_module,          only: PGSLib_BCAST

    ! Argument List
    integer, intent(in) :: msh_lun
    integer, intent(IN)  :: ncells_tot
    integer, intent(IN)  :: nnodes_tot
    integer, intent(OUT) :: ncells
    integer, intent(OUT) :: nnodes

    ! Local Variables
    logical :: fatal
    integer :: mesh_file_number_domains, nPE, max_block, &
                                min_block, critPE, memerror, dom
    integer, pointer, dimension(:) :: Mesh_File_Nnodes_List => null()
    integer, pointer, dimension(:) :: Mesh_File_Ncells_List => null()
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! If this routine is called,  then the mesh is being read in from
    ! a file.  In that case, there are three ways we might determine
    ! the extents of the domains.
    !   1. If the user did not provide any information about domains,
    !      and the mesh does not have any domain information, then
    !      we simply divide ncells_tot and nnodes_tot evenly across all PEs.
    !   2. If the user did not provide any information about domains
    !      but the mesh has domain information, then we read the mesh
    !      domain information and determine how to map to the PEs.
    !   3. If the user did provide domain information, we make sure
    !      it is consistent with the mesh, and use as much user
    !      information as possible.

    ! Number of Processors
    nPE = p_info%nPE
    fatal = .false.

    ! Check if the file is not of type "decomposed"
    mesh_file_number_domains = -1

    ! Check if the user provided any information about domain sizes, either
    ! through input deck or through mesh file.
    ! If no specific information about domain sizes, then use block 
    ! decomposition.
    ! Also, if mesh file has wrong number of domains, use block decomposition
    BLOCK_DECOMPOSITION: if ( ((.NOT. ANY(Ncells_List >= 0))   .OR. & ! No user domain sizes
         &                     (.NOT. ANY(Nnodes_List >= 0)) ) .OR. &
         &                    ( mesh_file_number_domains < 0)  .OR. & ! No mesh file domain sizes
         &                    ((mesh_file_number_domains > 0) .AND. & ! Wrong number of domains in file
         &                     (mesh_file_number_domains/= nPE)) ) then


       ! If mesh was not of type "decomposed", then we do simple decomposition
       ! Also, if mesh was decomposed, but number of domains /= number of processors,
       ! the we still do a simple decomposition, effectively ignoring the 
       ! provided mesh decomposition.
       if ((number_domains > 0) .and.(number_domains /= nPE)) then
          write(message,7) number_domains, nPE
7         format('User requested ',i0,' domains, but using ',i0)
          call TLS_warn (message)
       end if

       if ((mesh_file_number_domains >= 0) .and. (mesh_file_number_domains /= nPE)) then
          write(message,10) mesh_file_number_domains, number_domains
10        format('Mesh file decomposed into ',i0,' domains, but using ',i0)
          call TLS_warn (message)
       end if

       ! Do a simple block decomposition of the linear mesh data structures

       if (p_info%IsParallel) then

          ! True parallel run
          call BLOCK_PARAMETERS (max_block, min_block, critPE, ncells_tot)
          if (p_info%thisPE <= critPE) then
             ncells = max_block
          else
             ncells = min_block
          end if
          call BLOCK_PARAMETERS (max_block, min_block, critPE, nnodes_tot)
          if (p_info%thisPE <= critPE) then
             nnodes = max_block
          else
             nnodes = min_block
          end if

       else

          ! Serial mode
          ncells = ncells_tot
          nnodes = nnodes_tot

       end if

    end if BLOCK_DECOMPOSITION

    ! If user specified domain sizes, then use those domain sizes,
    ! No matter what was in the mesh file

    USER_DECOMPOSITION: if ((ANY(Ncells_List >= 0)) .and. & ! User domain sizes
         &                  (ANY(Nnodes_List >= 0))) then

       ! Use domain information provided by user
       ! For now, insist that number_domains == nPE
       fatal = .false.
       if (number_domains /= nPE) then
          fatal = .true.
          write (message, 12) number_domains, nPE
12        format ('user requested ',i8,' domains, but ran on ',i8,' processors')
       end if
       call TLS_fatal_if_any (fatal, 'MESH_DOMAIN_DECOMPOSITION: ' // message)

       if ((mesh_file_number_domains >= 0) .and. &
           (mesh_file_number_domains /= Number_Domains)) then
          write (message, 15) mesh_file_number_domains, number_domains
15        format (9x,'WARNING: Mesh file decomposed into ',i8,' domains, but input deck specified ',i8)
          mesh_file_number_domains = -1 !Reset, so will fall through on next IF check.
       end if

       if (p_info%IsParallel) then
          nnodes = Nnodes_List(p_info%thisPE)
          ncells = Ncells_List(p_info%thisPE)
       end if

    end if USER_DECOMPOSITION

    ! If mesh file specified correct number of domains use those domain sizes
    MESH_FILE_DECOMPOSITION: if (mesh_file_number_domains == nPE) then

       ! If mesh was decomposed into domains, and that information is in
       ! the mesh file (of type "decomposed") then use thta information
       ! to specify domain sizes
       number_domains = mesh_file_number_domains

       ! We are assuming one domain per PE (see IF clause above)

       if (p_info%IsParallel) then

          ! Parallel mode
          ncells = Mesh_File_Ncells_List(p_info%thisPE)
          nnodes = Mesh_File_Nnodes_List(p_info%thisPE)

       else

          ! Serial mode
          ncells = ncells_tot
          nnodes = nnodes_tot

          ! Nullify the pointers to be safe
          NULLIFY (Mesh_File_Ncells_List, Mesh_File_Nnodes_List) 

       end if

    end if MESH_FILE_DECOMPOSITION

    ! Deallocate temporaries
    if (ASSOCIATED(Mesh_File_Ncells_List) .and. ASSOCIATED(Mesh_File_Nnodes_List)) &
        DEALLOCATE (Mesh_File_Ncells_List, Mesh_File_Nnodes_List)

  END SUBROUTINE MESH_DOMAIN_SIZES

END MODULE MESH_DECOMPOSITION_MODULE
