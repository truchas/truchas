Module GRID_MAPPING_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !    Provide routines for mapping function fields between different
  !    grids which occupy the same physical space.
  !
  !    Public Interface(s):
  !
  !      compute_int_volumes
  !      map_cell_field
  !      grid_vol_fracs
  !      grid_vols
  !      write_int_volumes
  !      read_int_volumes
  !      right_int_volumes
  !      destroy_grid_int_vols
  !      destroy_gm_mesh
  !
  !      type gm_mesh
  !      type grid_int_vols
  !
  ! Author(s): Andrew Kuprat (kuprat@lanl.gov)
  !
  !=======================================================================
  USE hpsort, only: hpsortimp,hpsortim
  USE gm_mesh_type, only : gm_mesh, destroy_gm_mesh
  USE grid_mapping_utils, only: write_msg
  implicit none

  PRIVATE

  PUBLIC :: gm_mesh, destroy_gm_mesh

  PUBLIC :: compute_int_volumes
  PUBLIC :: map_cell_field
  PUBLIC :: grid_vol_fracs
  PUBLIC :: grid_vols
  PUBLIC :: write_int_volumes
  PUBLIC :: read_int_volumes
  PUBLIC :: right_int_volumes
  PUBLIC :: destroy_grid_int_vols

  character(len=80), dimension(3) :: message
  integer, parameter :: dp = KIND(1.0d0)

  type :: gm_mesh_checksum
     integer :: nnod = 0
     integer :: nelt = 0
     integer :: nodes_per_elt = 0
     real(dp) :: pos_node_0 = 0.d0
     real(dp) :: pos_node_1 = 0.d0
     real(dp) :: node_elt_0 = 0.d0
     real(dp) :: node_elt_1 = 0.d0
     real(dp) :: block_elt_0 = 0.d0
     real(dp) :: block_elt_1 = 0.d0
  end type gm_mesh_checksum

  type :: sparse_mtx_by_rows
     integer                        :: nrows=0, ncols=0, numentries=0
     integer, dimension(:), pointer :: irowst => null()
     integer, dimension(:), pointer :: jcn => null()
     real(dp),    dimension(:), pointer :: val => null()
  end type sparse_mtx_by_rows

  type, public :: grid_int_vols
     private
     type(sparse_mtx_by_rows)        :: BA
     real(dp), dimension(:), pointer :: Vol_eltA => null()
     real(dp), dimension(:), pointer :: Vol_eltB => null()
     real(dp), dimension(:), pointer :: VolinBmesh_eltA => null()
     real(dp), dimension(:), pointer :: VolinAmesh_eltB => null()
     real(dp), dimension(:), pointer :: VolBlockinBmesh_eltA => null()
     real(dp), dimension(:), pointer :: VolBlockinAmesh_eltB => null()
     integer                         :: nodes_per_elt_Amesh = 0
     integer                         :: nodes_per_elt_Bmesh = 0
     type(gm_mesh_checksum)          :: checksum_Amesh
     type(gm_mesh_checksum)          :: checksum_Bmesh
  end type grid_int_vols

  integer, parameter :: BOUNDARY=-1
  integer, parameter :: DEGENERATE_FACE=-2
  integer, parameter :: BEST=4
  integer, parameter :: FAT=5

  integer :: NUMWALKPOPS=0

  character(LEN=256) :: outstring

  integer, dimension(4,6), save :: node_hexface = reshape ( (/ &
       3, 7, 8, 4, &
       5, 6, 2, 1, &
       1, 4, 8, 5, &
       6, 7, 3, 2, &
       1, 2, 3, 4, &
       6, 5, 8, 7  &
       /), (/ 4, 6/) )  ! rt-hand rule would give INWARD normals.

  integer, dimension(3,4), save :: node_tetface = reshape ( (/ &
       1, 2, 3, &
       4, 2, 1, &
       4, 3, 2, &
       4, 1, 3  &
       /), (/ 3, 4/) )  ! rt-hand rule would give INWARD normals.

CONTAINS

  subroutine compute_int_volumes(mesh_a,mesh_b,int_vols,maxwarn,ier)

    type(gm_mesh), intent(in) :: mesh_a, mesh_b
    type(grid_int_vols), intent(inout) :: int_vols
    integer, optional, intent(in) :: maxwarn
    integer, optional, intent(out) :: ier

    integer, dimension(:,:), pointer :: volij,tempvolij => null()
    real(dp), dimension(:), pointer :: vol,tempvol => null()
    real(dp), dimension(:), pointer :: volint => null()
    integer, dimension(:), pointer :: eltb_found => null()
    integer, dimension(:), pointer :: iprm, iprmeltb => null()
    real(dp), dimension(:,:,:), pointer :: inplane_eltfaceb => null()
    logical, dimension(:), pointer :: isgap => null()
    integer, dimension(:,:), pointer :: adjelt_elta => null()
    integer, dimension(:,:), pointer :: adjface_elta => null()
    integer, dimension(:,:), pointer :: adjelt_eltb => null()
    integer, dimension(:,:), pointer :: adjface_eltb => null()
    logical, dimension(:), pointer :: in_mesh_b_block => null()
    real(dp), dimension(:,:), pointer :: pos_elta => null()
    integer :: currelta,&
         curreltb,matind,ind,i,k,num_int,matsize,irow,&
         currblockA,curreltb_old,&
         minblocka,maxblocka,minblockb,maxblockb,minblock,&
         maxblock,alloc_error,numwarn,maxwarn1
    real(dp) :: ascend,eps_dist,eps_vol,ref_vol
    real(dp), parameter :: eps_vol_frac=1.d-8
    real(dp), parameter :: eps_dist_frac=1.d-8
    real(dp), parameter :: pi=3.1415926535897932d0

    real(dp), dimension(3) :: xq
    real(dp) :: maxdist, dist_out

! Arrays that keep their state between calls to subroutines
    logical, dimension(:), pointer :: nextelt_logical => null()
    integer, dimension(:), pointer :: nextelt_integer => null()
    logical, dimension(:), pointer :: walk_logical => null()
    integer, dimension(:), pointer :: walk_integer1 => null()
    integer, dimension(:), pointer :: walk_integer2 => null()
    logical, dimension(:), pointer :: getvols_logical => null()
    integer, dimension(:), pointer :: getvols_integer1 => null()
    integer, dimension(:), pointer :: getvols_integer2 => null()

    if (present(ier)) then
       ier=0
    endif

    if (present(maxwarn)) then
       maxwarn1=maxwarn
    else
       maxwarn1=0
    endif

    call destroy_grid_int_vols(int_vols)

    if (size(mesh_a%node_elt,1).eq.4) then
       int_vols%nodes_per_elt_Amesh=4
    elseif (size(mesh_a%node_elt,1).eq.8) then
       int_vols%nodes_per_elt_Amesh=8
    else
       write(message,'(a/a/)') &
            'FATAL: Mesh A is neither hex or tet mesh!', &
            'Abort from routine COMPUTE_INT_VOLUMES.'
       call write_msg(message(:3))
       call error_handler
       return
    endif

    if (size(mesh_b%node_elt,1).eq.4) then
       int_vols%nodes_per_elt_Bmesh=4
    elseif (size(mesh_b%node_elt,1).eq.8) then
       int_vols%nodes_per_elt_Bmesh=8
    else
       write(message,'(a/a/)') &
            'FATAL: Mesh B is neither hex or tet mesh!', &
            'Abort from routine COMPUTE_INT_VOLUMES.'
       call write_msg(message(:3))
       call error_handler
       return
    endif

! Compute mesh checksums
    call get_checksum(mesh_a,int_vols%checksum_Amesh)
    call get_checksum(mesh_b,int_vols%checksum_Bmesh)
    
    allocate (iprmeltb(mesh_b%nelt))
    allocate (volint(mesh_b%nelt))
    allocate (eltb_found(mesh_b%nelt))

    allocate(int_vols%Vol_eltA(mesh_a%nelt))
    allocate(int_vols%Vol_eltB(mesh_b%nelt))
    allocate(int_vols%VolinBmesh_eltA(mesh_a%nelt))
    allocate(int_vols%VolinAmesh_eltB(mesh_b%nelt))
    allocate(int_vols%VolBlockinBmesh_eltA(mesh_a%nelt))
    allocate(int_vols%VolBlockinAmesh_eltB(mesh_b%nelt))

! Initialize mesh-element overlap volumes to zero.    
    int_vols%VolinBmesh_eltA=0.d0
    int_vols%VolinAmesh_eltB=0.d0
    int_vols%VolBlockinBmesh_eltA=0.d0
    int_vols%VolBlockinAmesh_eltB=0.d0

    call get_adjelt_relation(mesh_a%nelt,mesh_a%node_elt,adjelt_elta,adjface_elta)
    call get_adjelt_relation(mesh_b%nelt,mesh_b%node_elt,adjelt_eltb,adjface_eltb)

    call get_vols(mesh_a%nelt,mesh_a%node_elt,adjface_elta,mesh_a%pos_node,int_vols%Vol_eltA)
    call get_vols(mesh_b%nelt,mesh_b%node_elt,adjface_eltb,mesh_b%pos_node,int_vols%Vol_eltB)

    allocate(isgap(mesh_b%nelt))
    call get_isgap(mesh_b%nelt,mesh_b%node_elt,mesh_b%pos_node,int_vols%Vol_eltB, &
         eps_vol_frac,isgap)

    call get_inplane_eltface(mesh_b%nelt,mesh_b%node_elt,mesh_b%pos_node,isgap,adjface_eltb, &
         eps_dist_frac,inplane_eltfaceb)

    allocate (volij(2,mesh_a%nelt)) ! MESH_A%NELT is a 'chunksize' ; array will grow if necessary
    allocate (vol(mesh_a%nelt))     ! MESH_A%NELT is a 'chunksize' ; array will grow if necessary

    minblocka=minval(mesh_a%block_elt)
    maxblocka=maxval(mesh_a%block_elt)
    minblockb=minval(mesh_b%block_elt)
    maxblockb=maxval(mesh_b%block_elt)
    minblock=min(minblocka,minblockb)
    maxblock=max(maxblocka,maxblockb)
    allocate(in_mesh_b_block(minblock:maxblock),STAT=alloc_error)
    if (alloc_error.ne.0) then
       write(message,'(a/a/)') &
            'FATAL: Could not allocate element block range.', &
            'Abort from routine COMPUTE_INT_VOLUMES.'
       call write_msg(message(:3))
       call error_handler
       return
    endif
    in_mesh_b_block=.false.
    do i=1,mesh_b%nelt
       in_mesh_b_block(mesh_b%block_elt(i))=.true.
    enddo

    allocate(pos_elta(3,int_vols%nodes_per_elt_Amesh))

    currelta=0
    call get_nextelt(currelta,adjelt_elta,mesh_a%nelt, &
         nextelt_logical,nextelt_integer)
    matind=0
    curreltb=0
    numwarn=0
    do while (currelta.ne.0)

       if (in_mesh_b_block(mesh_a%block_elt(currelta))) then

          currblockA=mesh_a%block_elt(currelta)

          do k=1,int_vols%nodes_per_elt_Amesh
             pos_elta(:,k)=mesh_a%pos_node(:,mesh_a%node_elt(k,currelta))
          enddo

          call get_centroid_and_distance(pos_elta,xq,maxdist)

          ref_vol=4.d0/3.d0*pi*(maxdist**3)
          eps_vol=eps_vol_frac*ref_vol

! If CURRELTA is a gap element, then skip it.
          if (int_vols%Vol_eltA(currelta).le.eps_vol) then
             call get_nextelt(currelta,adjelt_elta,mesh_a%nelt, &
                  nextelt_logical,nextelt_integer)
             cycle
          endif

          eps_dist=eps_dist_frac*maxdist

          curreltb_old=curreltb
          call walk_mesh_pt(curreltb,xq,inplane_eltfaceb, &
               adjelt_eltb,mesh_b%nelt,isgap,dist_out,&
               walk_logical,walk_integer1,walk_integer2)

          if (dist_out.gt.maxdist) then
             if (maxwarn1.eq.0) then
                write(message,'(a,i10/a/)') &
                     'FATAL: WALK_MESH_PT fails for element #',currelta, &
                     'Abort from routine COMPUTE_INT_VOLUMES.'
                call write_msg(message(:3))
                call error_handler
                return
             endif
             numwarn=numwarn+1
             write(message,'(a,i10/a/)') &
                  'Warning: WALK_MESH_PT fails for element #',currelta, &
                  'Warning from routine COMPUTE_INT_VOLUMES.'
             call write_msg(message(:3))
             if (numwarn.le.maxwarn1) then
                if (curreltb_old.eq.0) then
                   curreltb=1
                else
                   curreltb=curreltb_old
                endif
                call get_nextelt(currelta,adjelt_elta,mesh_a%nelt, &
                     nextelt_logical,nextelt_integer)
                cycle
             else        
                write(message,'(a/a/)') &
                     'FATAL: MAXWARN exceeded.', &
                     'Abort from routine COMPUTE_INT_VOLUMES.'
                call write_msg(message(:3))
                call error_handler
                return
             endif
          endif

          curreltb_old=curreltb
          call get_vols_around_elta(curreltb,currelta, &
               adjelt_eltb,mesh_b%nelt,volint,eltb_found,num_int, &
               pos_elta,adjface_elta(:,currelta),inplane_eltfaceb, &
               int_vols%Vol_eltA(currelta),int_vols%Vol_eltB, &
               eps_dist,maxdist,isgap, BEST, &
               getvols_logical,getvols_integer1,getvols_integer2)
          ! restore CURRELTB to elt at centroid of CURRELTA
          curreltb=curreltb_old
          if (num_int.eq.0) then
             call get_vols_around_elta(curreltb,currelta, &
                  adjelt_eltb,mesh_b%nelt,volint,eltb_found,num_int, &
                  pos_elta,adjface_elta(:,currelta),inplane_eltfaceb, &
                  int_vols%Vol_eltA(currelta),int_vols%Vol_eltB, &
                  eps_dist,maxdist,isgap, FAT, &
                  getvols_logical,getvols_integer1,getvols_integer2)
             ! restore CURRELTB to elt at centroid of CURRELTA
             curreltb=curreltb_old
             if (num_int.eq.0) then
                if (maxwarn1.eq.0) then
                   write(message,'(a,i10/a/)') &
                        'FATAL: GET_VOLS_AROUND_ELTA fails for element #',currelta, &
                        'Abort from routine COMPUTE_INT_VOLUMES.'
                   call write_msg(message(:3))
                   call error_handler
                   return
                endif
                numwarn=numwarn+1
                write(message,'(a,i10/a/)') &
                     'Warning: GET_VOLS_AROUND_ELTA fails for element #',currelta, &
                     'Warning from routine COMPUTE_INT_VOLUMES.'
                call write_msg(message(:3))
                if (numwarn.le.maxwarn1) then
                   call get_nextelt(currelta,adjelt_elta,mesh_a%nelt, &
                        nextelt_logical,nextelt_integer)
                   cycle
                else        
                   write(message,'(a/a/)') &
                        'FATAL: MAXWARN exceeded.', &
                        'Abort from routine COMPUTE_INT_VOLUMES.'
                   call write_msg(message(:3))
                   call error_handler
                   return
                endif
             endif
          endif

          ! Add MESH B intersection data for ELTA into the unsorted row-column data
          ! structure for the volume intersection matrix

          do ind=1,num_int
             matind=matind+1

             if (matind.gt.size(volij,2)) then
                tempvolij => volij
                allocate(volij(2,size(tempvolij,2)+mesh_a%nelt))
                volij(1:2,1:size(tempvolij,2))=tempvolij
                deallocate(tempvolij)

                tempvol => vol
                allocate(vol(size(tempvol)+mesh_a%nelt))
                vol(1:size(tempvol))=tempvol
                deallocate(tempvol)
             endif

             volij(1,matind)=eltb_found(ind)
             volij(2,matind)=currelta

             if (mesh_b%block_elt(eltb_found(ind)).eq.currblockA) then
                int_vols%VolBlockinBmesh_eltA(currelta)=int_vols%VolBlockinBmesh_eltA(currelta)+volint(ind)
                int_vols%VolBlockinAmesh_eltB(eltb_found(ind))=int_vols%VolBlockinAmesh_eltB(eltb_found(ind))+volint(ind)
                vol(matind)=volint(ind)
             else
! We signify intersection volumes where the element blocks do not match
! by negating them.
                vol(matind)=-volint(ind)
             endif
             int_vols%VolinBmesh_eltA(currelta)=int_vols%VolinBmesh_eltA(currelta)+volint(ind)
             int_vols%VolinAmesh_eltB(eltb_found(ind))=int_vols%VolinAmesh_eltB(eltb_found(ind))+volint(ind)

          enddo

       endif

       call get_nextelt(currelta,adjelt_elta,mesh_a%nelt, &
            nextelt_logical,nextelt_integer)

    enddo

    matsize=matind

    ! Sort unsorted row-column data structure for BAint volumes

    allocate(iprm(matsize))
    ascend=1.0d0
    do i=1,matsize
       iprm(i)=i
    enddo
    call hpsortimp(matsize,2,2,volij,ascend,iprm)

    ! Create packed-row format matrix

    allocate(int_vols%BA%val(matsize))
    allocate(int_vols%BA%jcn(matsize))
    allocate(int_vols%BA%irowst(mesh_b%nelt+1))
    int_vols%BA%nrows=mesh_b%nelt
    int_vols%BA%ncols=mesh_a%nelt
    int_vols%BA%numentries=matsize

    matind=1
    do irow=1,int_vols%BA%nrows
       int_vols%BA%irowst(irow)=matind
       if (matind.gt.matsize) then
          cycle
       endif

       do while(volij(1,iprm(matind)).eq.irow)
          int_vols%BA%jcn(matind)=volij(2,iprm(matind))
          int_vols%BA%val(matind)=vol(iprm(matind))
          matind=matind+1
          if (matind.gt.matsize) then
             exit
          endif
       enddo
    enddo
    int_vols%BA%irowst(int_vols%BA%nrows+1)=matind          

    call cleanup

  contains

    subroutine cleanup
      
      if (associated(adjelt_elta)) deallocate(adjelt_elta)
      if (associated(adjface_elta)) deallocate(adjface_elta)
      if (associated(adjelt_eltb)) deallocate(adjelt_eltb)
      if (associated(adjface_eltb)) deallocate(adjface_eltb)
      if (associated(isgap)) deallocate(isgap)
      if (associated(inplane_eltfaceb)) deallocate(inplane_eltfaceb)
      if (associated(volij)) deallocate(volij)
      if (associated(vol)) deallocate(vol)
      if (associated(iprm)) deallocate(iprm)
      if (associated(volint)) deallocate(volint)
      if (associated(eltb_found)) deallocate(eltb_found)
      if (associated(iprmeltb)) deallocate(iprmeltb)
      if (associated(in_mesh_b_block)) deallocate(in_mesh_b_block)
      if (associated(pos_elta)) deallocate(pos_elta)
      if (associated(nextelt_logical)) deallocate(nextelt_logical)
      if (associated(nextelt_integer)) deallocate(nextelt_integer)
      if (associated(walk_logical)) deallocate(walk_logical)
      if (associated(walk_integer1)) deallocate(walk_integer1)      
      if (associated(walk_integer2)) deallocate(walk_integer2)      
      if (associated(getvols_logical)) deallocate(getvols_logical)
      if (associated(getvols_integer1)) deallocate(getvols_integer1)      
      if (associated(getvols_integer2)) deallocate(getvols_integer2)      

    end subroutine cleanup

    subroutine error_handler

      call destroy_grid_int_vols(int_vols)
      call cleanup
      if (.not.present(ier)) stop
      ier=-1

    end subroutine error_handler

  end subroutine compute_int_volumes

  subroutine get_adjelt_relation(nelt,node_elt,adjelt_elt,adjface_elt)

    integer, intent(in) :: nelt
    integer, dimension(:,:), intent(in) :: node_elt
    integer, dimension(:,:), pointer :: adjelt_elt
    integer, dimension(:,:), pointer :: adjface_elt

    if (size(node_elt,1).eq.4) then
       !   Tet-mesh
       allocate (adjelt_elt(4,nelt))
       allocate (adjface_elt(4,nelt))
       call get_adjtet_relation(nelt,node_elt,adjelt_elt,adjface_elt)
    else
       !   Hex-mesh
       allocate (adjelt_elt(6,nelt))
       allocate (adjface_elt(6,nelt))
       call get_adjhex_relation(nelt,node_elt,adjelt_elt,adjface_elt)
    endif

  end subroutine get_adjelt_relation

  subroutine get_adjtet_relation(ntet,node_tet,adjtet_tet,adjface_tet)

    integer, intent(in) :: ntet
    integer, dimension(:,:), intent(in) :: node_tet
    integer, dimension(:,:), intent(out) :: adjtet_tet
    integer, dimension(:,:), intent(out) :: adjface_tet
    integer :: i,j,k,nkey,facei,teti,faceip1,tetip1
    integer :: temp(3)
    integer, dimension(:,:), pointer :: ikey => null()

    allocate(ikey(5,ntet*4))
    nkey=0
    do i=1,ntet
       do j=1,4
          do k=1,3
             temp(k)=node_tet(node_tetface(k,j),i)
          enddo
          nkey=nkey+1
          ikey(1,nkey)=minval(temp)
          ikey(3,nkey)=maxval(temp)
          ikey(2,nkey)=sum(temp)-ikey(1,nkey)-ikey(3,nkey)
          ikey(4,nkey)=i
          ikey(5,nkey)=j
       enddo
    enddo

    call hpsortim(nkey,3,5,ikey)

    i=1
    do while(i.le.nkey)
       if (i.eq.nkey) then
          teti=ikey(4,i)
          facei=ikey(5,i)
          adjtet_tet(facei,teti)=BOUNDARY
          adjface_tet(facei,teti)=BOUNDARY
          exit
       endif
       if (count(ikey(1:3,i)==ikey(1:3,i+1))==3) then
          teti=ikey(4,i)
          facei=ikey(5,i)
          tetip1=ikey(4,i+1)
          faceip1=ikey(5,i+1)
          adjtet_tet(facei,teti)=tetip1
          adjface_tet(facei,teti)=faceip1
          adjtet_tet(faceip1,tetip1)=teti
          adjface_tet(faceip1,tetip1)=facei
          i=i+2
       else
          teti=ikey(4,i)
          facei=ikey(5,i)
          adjtet_tet(facei,teti)=BOUNDARY
          adjface_tet(facei,teti)=BOUNDARY
          i=i+1
       endif
    enddo

    deallocate (ikey)

  end subroutine get_adjtet_relation

  subroutine get_adjhex_relation(nhex,node_hex,adjhex_hex,adjface_hex)

    integer, intent(in) :: nhex
    integer, dimension(:,:), intent(in) :: node_hex
    integer, dimension(:,:), intent(out) :: adjhex_hex
    integer, dimension(:,:), intent(out) :: adjface_hex
    integer :: nkey,i,j,k,facei,hexi,faceip1,hexip1
    integer :: minvl, maxvl
    integer, dimension(1):: minpos
    integer, dimension(1):: maxpos
    integer, dimension(4):: temp
    logical, dimension(4):: mask
    integer, dimension(:,:), pointer :: ikey => null()

    allocate(ikey(6,nhex*6))

    adjhex_hex=BOUNDARY
    adjface_hex=BOUNDARY

    nkey=0
    do i=1,nhex
       do j=1,6
          do k=1,4
             temp(k)=node_hex(node_hexface(k,j),i)
          enddo
          nkey=nkey+1
          minpos=minloc(temp)
          minvl=minval(temp)
          maxpos=maxloc(temp)
          maxvl=maxval(temp)

! Put IKEY(1:4,nkey) in ascending order
          ikey(1,nkey)=minvl
          ikey(4,nkey)=maxvl
          mask=.true.
          mask(minpos)=.false.
          mask(maxpos)=.false.
          ikey(2,nkey)=minval(temp,mask=mask)
          ikey(3,nkey)=maxval(temp,mask=mask)
          ikey(5,nkey)=i
          ikey(6,nkey)=j

!   Check for degenerate face which is when there are topologically less
!   than 3 distinct vertices on a face.  If there are three distinct
!   vertices, this is not a degenerate face.
          if (ikey(1,nkey).eq.ikey(2,nkey)) then
             if (ikey(2,nkey).eq.ikey(3,nkey)) then
                adjhex_hex(j,i)=DEGENERATE_FACE
                adjface_hex(j,i)=DEGENERATE_FACE
             elseif (ikey(3,nkey).eq.ikey(4,nkey)) then
                adjhex_hex(j,i)=DEGENERATE_FACE
                adjface_hex(j,i)=DEGENERATE_FACE
             endif
          elseif ((ikey(2,nkey).eq.ikey(3,nkey)).and. &
             (ikey(3,nkey).eq.ikey(4,nkey))) then
                adjhex_hex(j,i)=DEGENERATE_FACE
                adjface_hex(j,i)=DEGENERATE_FACE
          endif
       enddo
    enddo

! Multikey sort puts IKEY(*,1:NKEY) in lexicographic ascending order based
! on the multikeys IKEY(1:4,*).
    call hpsortim(nkey,4,6,ikey)

! Match up pairs of faces
    i=1
    do while(i.lt.nkey)
       if (count(ikey(1:4,i)==ikey(1:4,i+1))==4) then
          hexi=ikey(5,i)
          facei=ikey(6,i)
          hexip1=ikey(5,i+1)
          faceip1=ikey(6,i+1)
          adjhex_hex(facei,hexi)=hexip1
          adjface_hex(facei,hexi)=faceip1
          adjhex_hex(faceip1,hexip1)=hexi
          adjface_hex(faceip1,hexip1)=facei
          i=i+2
       else
          i=i+1
       endif
    enddo

    deallocate(ikey)

  end subroutine get_adjhex_relation

  subroutine get_inplane_eltface(nelt,node_elt,pos_mesh,isgap,adjface_elt, &
       eps_dist_frac,inplane_eltface)

    integer, intent(in) :: nelt
    integer, dimension(:,:), intent(in) :: node_elt
    real(dp), dimension(:,:), intent(in) :: pos_mesh 
    logical, dimension(:), intent(in) :: isgap
    integer, dimension(:,:), intent(in) :: adjface_elt
    real(dp), intent(in) :: eps_dist_frac
    real(dp), pointer, dimension(:,:,:) :: inplane_eltface

    real(dp), dimension(3) :: x1,x2,x3,x4,x12,x23,x34,xq,areavec,normvec
    real(dp), dimension(:,:), pointer :: pos_elt
    real(dp) :: norm,maxdist,eps_dist
    integer :: i,j,k


    if (size(node_elt,1).eq.4) then
       !   Tet-mesh
       allocate(inplane_eltface(5,4,nelt))
       allocate(pos_elt(3,4))
       do i=1,nelt
! Inward pointing planes only defined for non-gap elements

          if (.not.isgap(i)) then
             do k=1,4
                pos_elt(:,k)=pos_mesh(:,node_elt(k,i))
             enddo
             call get_centroid_and_distance(pos_elt,xq,maxdist)
             eps_dist=eps_dist_frac*maxdist
             
             do j=1,4
                x1= pos_elt(:,node_tetface(1,j))
                x2= pos_elt(:,node_tetface(2,j))
                x3= pos_elt(:,node_tetface(3,j))
                areavec=areavector_tri(x1,x2,x3)
                norm=sqrt(areavec(1)**2+areavec(2)**2+areavec(3)**2)
                normvec=areavec/norm
                inplane_eltface(1:3,j,i)=normvec
                inplane_eltface(BEST,j,i)=dot_product(normvec,x1)
                inplane_eltface(FAT,j,i)=inplane_eltface(BEST,j,i)-eps_dist
             enddo
          endif

       enddo
       deallocate(pos_elt)
    else
       !   Hex-mesh
       allocate(inplane_eltface(5,6,nelt))
       allocate(pos_elt(3,8))
       do i=1,nelt
! Inward pointing planes only defined for non-gap elements

          if (.not.isgap(i)) then
             do k=1,8
                pos_elt(:,k)=pos_mesh(:,node_elt(k,i))
             enddo
             call get_centroid_and_distance(pos_elt,xq,maxdist)
             eps_dist=eps_dist_frac*maxdist

             do j=1,6
                if (adjface_elt(j,i).ne.DEGENERATE_FACE) then
                   x1=pos_elt(:,node_hexface(1,j))
                   x2=pos_elt(:,node_hexface(2,j))
                   x3=pos_elt(:,node_hexface(3,j))
                   x4=pos_elt(:,node_hexface(4,j))
                   x12=(x1+x2)*0.5d0
                   x23=(x2+x3)*0.5d0
                   x34=(x3+x4)*0.5d0
                   
                   areavec=areavector_tri(x12,x23,x34)
                   norm=sqrt(areavec(1)**2+areavec(2)**2+areavec(3)**2)
                   normvec=areavec/norm
                   inplane_eltface(1:3,j,i)=normvec
                   inplane_eltface(BEST,j,i)=dot_product(normvec,x12)
                   inplane_eltface(FAT,j,i)=min(dot_product(normvec,x1), &
                        dot_product(normvec,x2),dot_product(normvec,x3), &
                        dot_product(normvec,x4)) - eps_dist
                endif
             enddo
          endif
       enddo
       deallocate(pos_elt)

    endif

  end subroutine get_inplane_eltface

  subroutine get_centroid_and_distance(pos_elt,xq,dist)
    real(dp), dimension(:,:), intent(in) :: pos_elt
    real(dp), dimension(:), intent(out) :: xq
    real(dp), intent(out) :: dist

    integer :: k
    real(dp) :: distk

    xq=sum(pos_elt,2)/size(pos_elt,2)
    dist=0.d0
    do k=1,size(pos_elt,2)
       distk=((xq(1)-pos_elt(1,k))**2 + &
            (xq(2)-pos_elt(2,k))**2 + &
            (xq(3)-pos_elt(3,k))**2)
       dist=max(dist,distk)
    enddo
    dist=sqrt(dist)

  end subroutine get_centroid_and_distance

  subroutine get_inplane_hex(x,adjface_elt,inplane_hexface)

    real(dp), dimension(:,:), intent(in) :: x
    integer, dimension(:), intent(in) :: adjface_elt
    real(dp), dimension(:,:), intent(out) :: inplane_hexface

    real(dp), dimension(3) :: x12,x23,x34,areavec,normvec
    real(dp) :: norm
    integer :: j

    do j=1,6
       if (adjface_elt(j).ne.DEGENERATE_FACE) then
          x12=( x(:,node_hexface(1,j)) + x(:,node_hexface(2,j)) ) * 0.5d0
          x23=( x(:,node_hexface(2,j)) + x(:,node_hexface(3,j)) ) * 0.5d0
          x34=( x(:,node_hexface(3,j)) + x(:,node_hexface(4,j)) ) * 0.5d0
          areavec=areavector_tri(x12,x23,x34)
          norm=sqrt(areavec(1)**2+areavec(2)**2+areavec(3)**2)
          normvec=areavec/norm
          inplane_hexface(1:3,j)=normvec
          inplane_hexface(4,j)=dot_product(normvec,x12)
       endif
    enddo

  end subroutine get_inplane_hex

  subroutine get_vols(nelt,node_elt,adjface_elt,pos_mesh,vols)

    integer, intent(in) :: nelt
    integer, dimension(:,:), intent(in) :: node_elt
    integer, dimension(:,:), intent(in) :: adjface_elt
    real(dp), dimension(:,:), intent(in) :: pos_mesh 
    real(dp), dimension(:), intent(out) :: vols

    real(dp), dimension(3) :: x1,x2,x3,x4,xref,x12,x23,x34,areavec
    integer :: i,j

    if (size(node_elt,1).eq.4) then
       !   Tet-mesh
       do i=1,nelt
          x1= pos_mesh(:,node_elt(1,i))
          x2= pos_mesh(:,node_elt(2,i))
          x3= pos_mesh(:,node_elt(3,i))
          areavec=areavector_tri(x1,x2,x3)
          x4= pos_mesh(:,node_elt(4,i))
          vols(i)=dot_product(areavec,x4-x1)/3.d0
       enddo

    else
       !   Hex-mesh
       do i=1,nelt
          vols(i)=0.d0
          xref=pos_mesh(:,node_elt(1,i))
          do j=1,6
             if (adjface_elt(j,i).ne.DEGENERATE_FACE) then
                x12=( pos_mesh(:,node_elt(node_hexface(1,j),i)) + &
                     pos_mesh(:,node_elt(node_hexface(2,j),i)) ) * 0.5
                x23=( pos_mesh(:,node_elt(node_hexface(2,j),i)) + &
                     pos_mesh(:,node_elt(node_hexface(3,j),i)) ) * 0.5
                x34=( pos_mesh(:,node_elt(node_hexface(3,j),i)) + &
                     pos_mesh(:,node_elt(node_hexface(4,j),i)) ) * 0.5
                areavec=areavector_tri(x12,x23,x34)
                vols(i)=vols(i)+dot_product(areavec,xref-x12)*4.d0/3.d0
             endif
          enddo
       enddo

    endif

  end subroutine get_vols

  subroutine get_isgap(nelt,node_elt,pos_mesh,vols,eps_vol_frac,isgap)

    integer, intent(in) :: nelt
    integer, dimension(:,:), intent(in) :: node_elt
    real(dp), dimension(:,:), intent(in) :: pos_mesh 
    real(dp), dimension(:), intent(in) :: vols
    real(dp), intent(in) :: eps_vol_frac
    logical, dimension(:), intent(out) :: isgap

    real(dp), dimension(3) :: xq
    real(dp), dimension(:,:), pointer :: pos_elt
    real(dp) :: maxdist, ref_vol, eps_vol
    integer :: i,k,nodes_per_elt
    real(dp), parameter :: pi=3.1415926535897932d0

    nodes_per_elt=size(node_elt,1)
    allocate(pos_elt(3,nodes_per_elt))
    do i=1,nelt
       do k=1,nodes_per_elt
          pos_elt(:,k)=pos_mesh(:,node_elt(k,i))
       enddo

       call get_centroid_and_distance(pos_elt,xq,maxdist)
 
       ref_vol=4.d0/3.d0*pi*(maxdist**3)
       eps_vol=eps_vol_frac*ref_vol
       if (vols(i).le.eps_vol) then
          isgap(i)=.true.
       else
          isgap(i)=.false.
       endif
    enddo
    deallocate(pos_elt)
    
  end subroutine get_isgap

  function areavector_tri(x1,x2,x3)

    real(dp), dimension(3), intent(in) :: x1,x2,x3
    real(dp), dimension(3) :: areavector_tri

    real(dp), dimension(3) :: a,b

    a=x2-x1
    b=x3-x1

    areavector_tri(1)=a(2)*b(3)-a(3)*b(2)
    areavector_tri(2)=a(3)*b(1)-a(1)*b(3)
    areavector_tri(3)=a(1)*b(2)-a(2)*b(1)

    areavector_tri=areavector_tri*0.5d0

  end function areavector_tri

  subroutine get_nextelt(currelt,adjelt_elt,nelt,onstack,s)

    integer, intent(inout) :: currelt
    integer, dimension(:,:), intent(in) :: adjelt_elt
    integer, intent(in) :: nelt
    logical, dimension(:), pointer :: onstack
    integer, dimension(:), pointer :: s

    integer :: i,oppelt
    integer, save :: top, next_pos_seed

    ! CURRELT = 0 happens upon first invocation ; allocate and initialize
    if (currelt.eq.0) then
       allocate(onstack(nelt))
       onstack=.false.
       allocate(s(nelt))
       top=0
       next_pos_seed=1
    endif
    ! TOP = 0 happens upon first invocation or when we have just exhaustively
    ! searched a face-connected component of mesh.
    ! Cycle sequentially through the mesh looking for an unvisited element.
    ! If one is found, this will start a new face-connected component that we
    ! will exhaustively search.
    if (top.eq.0) then
       do i=next_pos_seed,nelt
          if (.not.onstack(i)) then
             currelt=i
             call push(s,top,currelt)
             onstack(currelt)=.true.
             next_pos_seed=i+1
             exit
          endif
       enddo
       ! We have exhaustively searched all components of the mesh
       if (top.eq.0) then
          currelt=0
          return
       endif
    endif
    ! Pop next element off mesh ; this element will be returned
    call pop(s,top,currelt)
    ! Place all unvisited neighbors on the stack
    do i=1,size(adjelt_elt,1)
       oppelt=adjelt_elt(i,currelt)
       if ((oppelt.ne.BOUNDARY).and.(oppelt.ne.DEGENERATE_FACE)) then
          if (.not.onstack(oppelt)) then
             call push(s,top,oppelt)
             onstack(oppelt)=.true.
          endif
       endif
    enddo

  end subroutine get_nextelt

  function tripleproduct(x1,x2,x3)

    real(dp), dimension(3), intent(in) :: x1,x2,x3
    real(dp)               :: tripleproduct

    real(dp), dimension(3) :: x1_cross_x2

    x1_cross_x2(1)=x1(2)*x2(3)-x1(3)*x2(2)
    x1_cross_x2(2)=x1(3)*x2(1)-x1(1)*x2(3)
    x1_cross_x2(3)=x1(1)*x2(2)-x1(2)*x2(1)

    tripleproduct=dot_product(x1_cross_x2,x3)

  end function tripleproduct

  subroutine walk_mesh_pt(currelt,x,inplane_eltface, &
       adjelt_elt,nelt,isgap,dist_out,onstack,onstacklist,eltstack)

    integer, intent(inout) :: currelt
    integer, intent(in) :: nelt 
    real(dp), dimension(3), intent(in) :: x
    real(dp), dimension(:,:,:), intent(in) :: inplane_eltface
    integer, dimension(:,:), intent(in) :: adjelt_elt
    logical, dimension(:), intent(in) :: isgap
    real(dp), intent(out) :: dist_out
    logical, dimension(:), pointer :: onstack
    integer, dimension(:), pointer :: onstacklist
    integer, dimension(:), pointer :: eltstack

    integer :: lenlist,top,h_top,i,h_opp,next_poss_seed_elt
    integer :: best_choice
    real(dp) :: dist_out_elt, dist_out_i, dist_out_htop
!
!  This algorithm differs somewhat from Algorithm 2 in the 
!  Physics_Algorithms manual (which
!  refers to the previous version of this algorithm).  Now when we
!  walk from element to element on the mesh, we keep track of 
!  DIST_OUT which is a measure of how much the best-element-so-far
!  failed in containing the query point X.  If no elements
!  in the mesh contain X, we return the one with lowest DIST_OUT value.
!  Although this causes the whole mesh to be searched, this is hopefully
!  rare, and often it will lead to a good intersection calculation
!  for mesh elements near a boundary where the meshes do not match up well.
!  That is, the element returned, although not containing X, 
!  will provide an adequate seed for subroutine GET_VOLS_AROUND_ELTA
!  and good intersection volumes will be generated.
!
    best_choice=0

    ! If CURRELT is zero, allocate work arrays and use '1' as the 
    ! first elt for starting the search for the elt containing the
    ! point X.
    if (currelt.eq.0) then
       allocate(onstack(nelt))
       allocate(onstacklist(nelt))
       allocate(eltstack(nelt))
       onstack=.false.
       currelt=1
       next_poss_seed_elt=2
    else
       next_poss_seed_elt=1
    endif

    lenlist=0

    ! Put CURRELT on stack
    do while (currelt.ne.0)
       top=1
       eltstack(top)=currelt
       onstack(currelt)=.true.
       lenlist=lenlist+1
       onstacklist(lenlist)=currelt

       ! Pop next elt off stack
       do while (top>0)
          currelt=eltstack(top)
          top=top-1
          NUMWALKPOPS=NUMWALKPOPS+1
          ! If CURRELT is a gap element, we do not try to locate X
          ! within this element, but put its neighbors on the stack.
          if (isgap(currelt)) then
             do i=1,size(adjelt_elt,1)
                h_opp=adjelt_elt(i,currelt)
                if ((h_opp.ne.BOUNDARY).and.(h_opp.ne.DEGENERATE_FACE)) then
                   if (.not.onstack(h_opp)) then
                      top=top+1
                      eltstack(top)=h_opp
                      onstack(h_opp)=.true.
                      lenlist=lenlist+1
                      onstacklist(lenlist)=h_opp
                   endif
                endif
             enddo
             cycle
          endif

! We loop over the faces and see how much 'outside' of each face X is.
! If X is not outside of any face, then X is inside CURRELT.
! When we look at each face, we put the opposite element on the stack
! (if that element exists and if it has not been put on the stack
! before).  However, as we go through the faces, if X is currently
! outside of face I by a distance DIST_OUT_I and that distance is 
! greater than DIST_OUT_HTOP, and the opposite element
! has not been put on the stack before, we call that element H_TOP.
! (In this case, DIST_OUT_I becomes the new DIST_OUT_HTOP.)
! We will delay putting H_TOP on the stack until the end, so that
! H_TOP is popped off first and we proceed immediately in a favorable
! direction in the mesh.  In the case that a better candidate for 
! H_TOP is encountered, we place the previous H_TOP candidate onto
! the stack.
          h_top=0
          do i=1,size(adjelt_elt,1)
             h_opp=adjelt_elt(i,currelt)
             if (h_opp.ne.DEGENERATE_FACE) then
                dist_out_i=inplane_eltface(FAT,i,currelt)- &
                     dot_product(x,inplane_eltface(1:3,i,currelt))
                if (i.eq.1) then
                   dist_out_elt=dist_out_i
                else
                   dist_out_elt=max(dist_out_elt,dist_out_i)
                endif
                if (h_opp.ne.BOUNDARY) then
                   if (.not.onstack(h_opp)) then
! we will either put H_OPP on the stack
! or make it become the new H_TOP.
                      if (h_top.eq.0) then
! H_TOP unassigned, so H_OPP becomes our first H_TOP...
                         h_top=h_opp
                         dist_out_htop=dist_out_i
                      elseif (dist_out_i.gt.dist_out_htop) then
! If DIST_OUT_I exceeds DIST_OUT_HTOP then H_TOP is updated...
! The rationale is that if X is a greater distance normal to the outside
! of this face than any other face, we should traverse this face
! instead.
                         if (.not.onstack(h_top)) then
! H_TOP is being replaced, so put it on the stack
                            top=top+1
                            eltstack(top)=h_top
                            onstack(h_top)=.true.
                            lenlist=lenlist+1
                            onstacklist(lenlist)=h_top
                         endif
                         h_top=h_opp
                         dist_out_htop=dist_out_i
                      else
! Just put H_OPP on the stack
                         top=top+1
                         eltstack(top)=h_opp
                         onstack(h_opp)=.true.
                         lenlist=lenlist+1
                         onstacklist(lenlist)=h_opp
                      endif
                   endif
                endif
             endif
          enddo

          if (best_choice.eq.0) then
             best_choice=currelt
             dist_out=dist_out_elt
          elseif (dist_out_elt.lt.dist_out) then
             best_choice=currelt
             dist_out=dist_out_elt
          endif

          if (dist_out.le.0) then
             ! The point X has been located within CURRELT.
             ! Clear ONSTACK array for next time.
             do i=1,lenlist
                onstack(onstacklist(i))=.false.
             enddo
             return
          else
             ! Current elt not viable.  Put the adjacent elt H_TOP onto the stack, so 
             ! that it will be processed next.  H_TOP is our candidate for the best
             ! element to walk to now.
             if (h_top.ne.0) then
                if (.not.onstack(h_top)) then
                   top=top+1
                   eltstack(top)=h_top
                   onstack(h_top)=.true.
                   lenlist=lenlist+1
                   onstacklist(lenlist)=h_top
                endif
             endif
          endif
       enddo

       ! Couldn't find elt in this connected component of elts.
       ! Look for an unvisited elt and use it as a seed for searching
       ! another connected component of elts.
       currelt=0
       do i=next_poss_seed_elt,nelt
          if (.not.onstack(i)) then
             currelt=i
             next_poss_seed_elt=i+1
             exit
          endif
       enddo

    enddo

    ! If we are here, then we failed to locate X within an elt, but
    ! we return best_choice

    currelt=best_choice
    onstack=.false.

  end subroutine walk_mesh_pt

  subroutine get_vols_around_elta(currelt,elta, &
       adjelt_elt,nelt,volint,elt,num_int, &
       x_elta,adjface_elta,inplane_eltface,volelta,volb, &
       eps_dist,maxdist,isgap, &
       b_elt_type,onstack,onstacklist,eltstack)

    use overlap_module, only: vector_type,poly3D_type, &
         Load_Hex, Load_Tet, RK, &
         plane_type, plane_poly_int3d, volm_Poly3D

    integer, intent(inout) :: currelt
    integer, intent(out) :: num_int
    integer, intent(in) :: elta
    integer, intent(in) :: nelt
    integer, dimension(:,:), intent(in) :: adjelt_elt
    real(dp), dimension(:), intent(inout) :: volint ! CANNOT USE INTENT(OUT)!!
    integer, dimension(:), intent(inout) :: elt     ! CANNOT USE INTENT(OUT)!!
    real(dp), dimension(:,:), intent(in) :: x_elta
    integer, dimension(:), intent(in) :: adjface_elta
    real(dp), dimension(:,:,:), intent(in) :: inplane_eltface
    real(dp), intent(in) :: volelta,eps_dist,maxdist
    real(dp), dimension(:), intent(in) :: volb
    logical, dimension(:), intent(in) :: isgap
    integer, intent(in) :: b_elt_type
    logical, dimension(:), pointer :: onstack
    integer, dimension(:), pointer :: onstacklist
    integer, dimension(:), pointer :: eltstack

    real(RK) :: Volint_poly
    type(vector_type) :: t1,t2,t3,t4
    type(vector_type) :: h1,h2,h3,h4,h5,h6,h7,h8
    type(poly3D_type) :: poly_elta
    type(poly3D_type) :: tmp, int_poly
    type(vector_type) :: nrml
    type(plane_type) :: plane
    real(RK) :: dist

    integer :: lenlist,top,i,h_opp,j
    real(dp)            :: distance,dmin,dmax
    logical, dimension(6) :: lplanecuts  ! 6 is maximum number of faces on elements of mesh B
    logical :: lint_nonempty
    real(dp), dimension(4,6) :: inplane_hexface ! specific for hexes
    real(dp) :: max_x,min_x,max_y,min_y,max_z,min_z
    logical :: planar, elta_inside_seedelt

    num_int=0

    ! Allocate work arrays (and initialize ONSTACK) on first invocation.
    if (.not.associated(onstack)) then
       allocate(onstack(nelt))
       allocate(onstacklist(nelt))
       allocate(eltstack(nelt))
       onstack=.false.
    endif

    if (size(x_elta,2).eq.4) then
       t1%x=x_elta(:,1)
       t2%x=x_elta(:,2)
       t3%x=x_elta(:,3)
       t4%x=x_elta(:,4)
       poly_elta=Load_Tet(t1,t2,t3,t4)
    else
       ! ELTA is a hex.  Determine if it is planar.  If so, use Load_Hex with vertex data.
       ! Otherwise, use Load_Hex to create an enclosing hex, and then intersect this
       ! hex with the six planarized faces of the nonplanar hex.

       call get_inplane_hex(x_elta,adjface_elta,inplane_hexface)
       planar=.true.
       CHKPLANAR: do j=1,6
          if (adjface_elta(j).ne.DEGENERATE_FACE) then
             do i=1,4
                distance=abs(dot_product(x_elta(:,node_hexface(i,j)),inplane_hexface(1:3,j)) &
                     -inplane_hexface(4,j))
                if (distance.gt.eps_dist*0.1d0) then
                   planar=.false.
                   exit CHKPLANAR
                endif
             enddo
          endif
       enddo CHKPLANAR

       if (planar.and.(count(adjface_elta.eq.DEGENERATE_FACE).eq.0)) then

          h1%x=x_elta(:,1)
          h2%x=x_elta(:,2)
          h3%x=x_elta(:,3)
          h4%x=x_elta(:,4)
          h5%x=x_elta(:,5)
          h6%x=x_elta(:,6)
          h7%x=x_elta(:,7)
          h8%x=x_elta(:,8)
          poly_elta=Load_Hex(h1,h2,h3,h4,h5,h6,h7,h8)

       else

! Make a Cartesian hex that is guaranteed to be bigger than ELT_A
          max_x=maxval(x_elta(1,:))+10.*maxdist
          min_x=minval(x_elta(1,:))-10.*maxdist
          max_y=maxval(x_elta(2,:))+10.*maxdist
          min_y=minval(x_elta(2,:))-10.*maxdist
          max_z=maxval(x_elta(3,:))+10.*maxdist
          min_z=minval(x_elta(3,:))-10.*maxdist
          h1%x=(/min_x,min_y,min_z/)
          h2%x=(/max_x,min_y,min_z/)
          h3%x=(/max_x,max_y,min_z/)
          h4%x=(/min_x,max_y,min_z/)
          h5%x=(/min_x,min_y,max_z/)
          h6%x=(/max_x,min_y,max_z/)
          h7%x=(/max_x,max_y,max_z/)
          h8%x=(/min_x,max_y,max_z/)
          poly_elta=Load_Hex(h1,h2,h3,h4,h5,h6,h7,h8)

! Intersect the bigger hex with the averaged face planes of ELT_A to 
! get a 'planarized' version of ELT_A.
          do i=1,6
             if (adjface_elta(i).ne.DEGENERATE_FACE) then
                nrml=vector_type(-1.d0*inplane_hexface(1:3,i))
                dist=-1.d0*inplane_hexface(4,i)
                plane=plane_type(nrml,dist)
                call plane_poly_int3d(plane,poly_elta,tmp)
                poly_elta=tmp
                if (poly_elta%numCrnr<4) then
                   write(outstring,*) 
                   write(message,'(a,i10,a/a/)') &
                        'Warning: planarization of hex ',elta,' annihilates it!', &
                        'Warning from routine GET_VOLS_AROUND_ELTA.'
                   call write_msg(message(:3))
                   exit
                endif
             endif
          enddo

       endif
    endif

! See if ELTA fits entirely within the FAT version of the seed element
! (initial value of CURRELT) of the B mesh.  If so, record the volume and exit.

    elta_inside_seedelt= .true.
    CHKINSIDE: do i=1,size(inplane_eltface,2) ! loop over faces of CURRELT
       if (adjelt_elt(i,currelt).ne.DEGENERATE_FACE) then
          do j=1,size(x_elta,2)  ! loop over nodes of ELTA
             if (dot_product(x_elta(:,j),inplane_eltface(1:3,i,currelt)).lt. &
                inplane_eltface(FAT,i,currelt)) then
                elta_inside_seedelt=.false.
                exit CHKINSIDE
             endif
          enddo
       endif
    enddo CHKINSIDE
    
    if (elta_inside_seedelt) then
       num_int=1
       elt(num_int)=currelt
       volint(num_int)=volelta
       return
    endif

    lenlist=0

    ! Put CURRELT on stack
    top=1
    eltstack(top)=currelt
    onstack(currelt)=.true.
    lenlist=lenlist+1
    onstacklist(lenlist)=currelt

    ! Pop next elt off stack
    do while (top>0)
       currelt=eltstack(top)
       top=top-1

       ! If CURRELT is a gap element, we avoid intersecting it
       ! with ELTA, but put its unvisited neighbors on the stack.
       if (isgap(currelt)) then
          do i=1,size(adjelt_elt,1)
             h_opp=adjelt_elt(i,currelt)
             if ((h_opp.ne.BOUNDARY).and.(h_opp.ne.DEGENERATE_FACE)) then
                if (.not.onstack(h_opp)) then
                   top=top+1
                   eltstack(top)=h_opp
                   onstack(h_opp)=.true.
                   lenlist=lenlist+1
                   onstacklist(lenlist)=h_opp
                endif
             endif
          enddo
          cycle
       endif

       ! Get volume of intersection between ELTA and CURRELT.

       lplanecuts=.false.

       int_poly=poly_elta
       lint_nonempty=.true.
       do i=1,size(inplane_eltface,2)  ! Loop over faces of CURRELT
          if (adjelt_elt(i,currelt).ne.DEGENERATE_FACE) then
! If the nodes of ELTA all lie on the wrong side of the planarized face,
! the intersection with CURRELT is null.  If the nodes of ELTA all lie
! on the right side of the planarized face, then the intersection of this
! face with the current intersection polyhedron can be skipped.
             dmin=1.d50+10.d0*eps_dist ! i.e. +infinity
             dmax=-dmin             ! i.e. -infinity
             do j=1,size(x_elta,2)        ! Loop over nodes of ELTA
                distance=dot_product(x_elta(:,j),inplane_eltface(1:3,i,currelt)) &
                     - inplane_eltface(b_elt_type,i,currelt)
                dmin=min(dmin,distance)
                dmax=max(dmax,distance)
             enddo
          endif
          if (dmax.le.eps_dist) then
             lint_nonempty=.false.
             exit
          elseif (dmin.ge.-eps_dist) then
             cycle
          endif
          lplanecuts(i)=.true.
          nrml=vector_type(-1.d0*inplane_eltface(1:3,i,currelt))
          dist=-1.d0*inplane_eltface(b_elt_type,i,currelt)
          plane=plane_type(nrml,dist)
          call plane_poly_int3d(plane,int_poly,tmp)
          if (tmp%numCrnr<4) then
             lint_nonempty=.false.
             exit
          else
             int_poly=tmp
          endif
       enddo
       if (lint_nonempty) then
          volint_poly=volm_Poly3D(int_poly)
!  The following two lines of code impose the most basic restriction on the 
!  intersection volume:  The intersection cannot have a negative volume,
!  and cannot exceed the volume of either element.
          volint_poly=max(volint_poly,0.d0)
          volint_poly=min(volint_poly,volelta,volb(currelt))
       else
          volint_poly=0.d0
       endif

       if (volint_poly.gt.0.d0) then
          num_int=num_int+1
          volint(num_int)=volint_poly
          elt(num_int)=currelt
          ! Place adjacent elements that cut elta onto stack
          do i=1,size(adjelt_elt,1)
             if (lplanecuts(i)) then
                h_opp=adjelt_elt(i,currelt)
                if ((h_opp.ne.BOUNDARY).and.(h_opp.ne.DEGENERATE_FACE)) then
                   if (.not.onstack(h_opp)) then
                      top=top+1
                      eltstack(top)=h_opp
                      onstack(h_opp)=.true.
                      lenlist=lenlist+1
                      onstacklist(lenlist)=h_opp
                   endif
                endif
             endif
          enddo
       endif
    enddo
    do i=1,lenlist
       onstack(onstacklist(i))=.false.
    enddo

  end subroutine get_vols_around_elta

  subroutine push(s,top,idat)
    integer, dimension(:), intent(inout) :: s
    integer, intent(inout) :: top
    integer, intent(in) :: idat

    if (top.ge.size(s)) then
! This is not supposed to happen under any circumstances,
! so we do not provide a way of continuing execution.
       write(message,'(a/a/)') &
            'FATAL: Stack overflow.', &
            'Abort from routine PUSH in GRID_MAPPING_MODULE.'
       call write_msg(message(:3))
       stop
    endif
    top=top+1
    s(top)=idat

  end subroutine push

  subroutine pop(s,top,idat)
    integer, dimension(:), intent(inout) :: s
    integer, intent(inout) :: top
    integer, intent(out) :: idat

    if (top.eq.0) then
! This is not supposed to happen under any circumstances,
! so we do not provide a way of continuing execution.
       write(message,'(a/a/)') &
            'FATAL: Stack empty.', &
            'Abort from routine POP in GRID_MAPPING_MODULE.'
       call write_msg(message(:3))
       stop
    else
       idat=s(top)
       top=top-1
    endif

  end subroutine pop

  subroutine MAP_CELL_FIELD(src,dest,int_vols, &                ! necessary
       reverse_order,exactly_conservative,preserve_constants, & ! optional
       strict,defval,ier)                                       ! optional 

    real(dp), dimension(:), intent(in) :: src
    real(dp), dimension(:), intent(out) :: dest
    type(grid_int_vols), intent(in) :: int_vols
    logical, intent(in), optional :: reverse_order, exactly_conservative, &
         preserve_constants,strict
    real(dp), intent(in), optional :: defval
    integer, intent(out), optional :: ier

    logical  :: used_default, reverse_order1, exactly_conservative1, preserve_constants1, strict1
    real(dp) :: defval1
    integer :: i,j,k,count

    if (present(ier)) ier = 0

    if (present(reverse_order)) then
       reverse_order1=reverse_order
    else
       reverse_order1=.false.
    endif

    if (present(exactly_conservative)) then
       exactly_conservative1=exactly_conservative
    else
       exactly_conservative1=.false.
    endif

    if (present(preserve_constants)) then
       preserve_constants1=preserve_constants
    else
       preserve_constants1=.false.
    endif

    if (exactly_conservative1.and.preserve_constants1) then
       write(message,'(a/a/)') &
            'FATAL: EXACTLY_CONSERVATIVE and PRESERVE_CONSTANTS both true', &
            'Abort from routine MAP_CELL_FIELD'
       call write_msg(message(:3))
       if (present(ier)) then
          ier=-1
          return
       else
          stop
       endif
    endif

    if (present(strict)) then
       strict1=strict
    else
       strict1=.true.
    endif

    if (present(defval)) then
       defval1=defval
    else
       defval1=0.d0
    endif

    used_default=.false.

    if (.not.reverse_order1) then
       if (.not.strict1) then
          do i=1,int_vols%BA%nrows
             if (int_vols%BA%irowst(i)<int_vols%BA%irowst(i+1)) then

                dest(i)=0.0d0

                do k=int_vols%BA%irowst(i),int_vols%BA%irowst(i+1)-1
                   if (exactly_conservative1) then
                      dest(i)=dest(i)+abs(int_vols%BA%val(k))*src(int_vols%BA%jcn(k))* &
                           int_vols%Vol_eltA(int_vols%BA%jcn(k))/int_vols%VolinBmesh_eltA(int_vols%BA%jcn(k))
                   else
                      dest(i)=dest(i)+abs(int_vols%BA%val(k))*src(int_vols%BA%jcn(k))
                   endif
                enddo
                if (preserve_constants1) then
                   dest(i)=dest(i)/int_vols%VolinAmesh_eltB(i)
                else
                   dest(i)=dest(i)/int_vols%Vol_eltB(i)
                endif

             else
                dest(i)=defval1
                used_default = .true.
             endif
          enddo
       else
          do i=1,int_vols%BA%nrows
             count=0
             dest(i)=0.0d0

             do k=int_vols%BA%irowst(i),int_vols%BA%irowst(i+1)-1
                if (int_vols%BA%val(k).gt.0.d0) then
                   count=count+1
                   if (exactly_conservative1) then
                      dest(i)=dest(i)+int_vols%BA%val(k)*src(int_vols%BA%jcn(k))* &
                           int_vols%Vol_eltA(int_vols%BA%jcn(k))/int_vols%VolBlockinBmesh_eltA(int_vols%BA%jcn(k))
                   else
                      dest(i)=dest(i)+int_vols%BA%val(k)*src(int_vols%BA%jcn(k))
                   endif
                endif
             enddo
             if (count.gt.0) then
                if (preserve_constants1) then
                   dest(i)=dest(i)/int_vols%VolBlockinAmesh_eltB(i)
                else
                   dest(i)=dest(i)/int_vols%Vol_eltB(i)
                endif
             else
                dest(i)=defval1
                used_default = .true.
             endif
          enddo
       endif
    else
       do j=1,int_vols%BA%ncols
          dest(j)=0.d0
       enddo

       if (.not.strict1) then
          do i=1,int_vols%BA%nrows
             do k=int_vols%BA%irowst(i),int_vols%BA%irowst(i+1)-1
                if (exactly_conservative1) then
                   dest(int_vols%BA%jcn(k))=dest(int_vols%BA%jcn(k))+abs(int_vols%BA%val(k))*src(i)* &
                        int_vols%Vol_eltB(i)/int_vols%VolinAmesh_eltB(i)
                else
                   dest(int_vols%BA%jcn(k))=dest(int_vols%BA%jcn(k))+abs(int_vols%BA%val(k))*src(i)
                endif
             enddo
          enddo
          do j=1,int_vols%BA%ncols
             if (int_vols%VolinBmesh_eltA(j).eq.0.d0) then
                dest(j)=defval1
                used_default=.true.
             else
                if (preserve_constants1) then
                   dest(j)=dest(j)/int_vols%VolinBmesh_eltA(j)
                else
                   dest(j)=dest(j)/int_vols%Vol_eltA(j)
                endif
             endif
          enddo
       else
          do i=1,int_vols%BA%nrows
             do k=int_vols%BA%irowst(i),int_vols%BA%irowst(i+1)-1
                if (int_vols%BA%val(k).gt.0.d0) then
                   if (exactly_conservative1) then
                      dest(int_vols%BA%jcn(k))=dest(int_vols%BA%jcn(k))+int_vols%BA%val(k)*src(i)* &
                           int_vols%Vol_eltB(i)/int_vols%VolBlockinAmesh_eltB(i)
                   else
                      dest(int_vols%BA%jcn(k))=dest(int_vols%BA%jcn(k))+int_vols%BA%val(k)*src(i)
                   endif
                endif
             enddo
          enddo
          do j=1,int_vols%BA%ncols
             if (int_vols%VolBlockinBmesh_eltA(j).eq.0.d0) then
                dest(j)=defval1
                used_default=.true.
             else
                if (preserve_constants1) then
                   dest(j)=dest(j)/int_vols%VolBlockinBmesh_eltA(j)
                else
                   dest(j)=dest(j)/int_vols%Vol_eltA(j)
                endif
             endif
          enddo
       endif
    endif

    if (used_default .and. (.not.present(defval))) then
       write(message,'(a/a/)') &
            'Warning: default value needed and none provided; used zero.', &
            'Warning from routine MAP_CELL_FIELD'
       call write_msg(message(:3))
    endif

  end subroutine map_cell_field

  subroutine grid_vol_fracs(int_vols,vf_a,vf_b,strict,ier)

    type(grid_int_vols), intent(in) :: int_vols
    real(dp), dimension(:), intent(out) :: vf_a,vf_b
    logical, intent(in), optional :: strict
    integer, intent(out), optional :: ier

    logical strict1

    if (present(ier)) ier=0

    if (size(vf_a).ne.int_vols%BA%ncols) then 
       write(message,'(a/a/)') 'FATAL:  Volume fraction array VF_A has incorrect length.', &
            'Abort from routine GRID_VOL_FRACS'
       call write_msg(message(:3))
       if (.not.present(ier)) stop
       ier=-1
       return
    endif
    if (size(vf_b).ne.int_vols%BA%nrows) then 
       write(message,'(a/a/)') 'FATAL:  Volume fraction array VF_B has incorrect length.', &
            'Abort from routine GRID_VOL_FRACS'
       call write_msg(message(:3))
       if (.not.present(ier)) stop
       ier=-1
       return
    endif

    if (present(strict)) then
       strict1=strict
    else
       strict1=.true.
    endif

    if (strict1) then
       if (.not.associated(int_vols%VolBlockinBmesh_eltA).or. &
            .not.associated(int_vols%VolBlockinAmesh_eltB)) then
          write(message,'(a/a/)') 'FATAL:  Fractional volume arrays not defined.', &
               'Abort from routine GRID_VOL_FRACS'
          call write_msg(message(:3))
          if (.not.present(ier)) stop
          ier=-1
          return
       endif

       where (int_vols%Vol_eltA > 0.d0) 
          vf_a=int_vols%VolBlockinBmesh_eltA/int_vols%Vol_eltA
       elsewhere
          vf_a=0.d0
       endwhere

       where (int_vols%Vol_eltB > 0.d0)
          vf_b=int_vols%VolBlockinAmesh_eltB/int_vols%Vol_eltB
       elsewhere
          vf_b=0.d0
       endwhere

    else
       if (.not.associated(int_vols%VolinBmesh_eltA).or. &
            .not.associated(int_vols%VolinAmesh_eltB)) then
          write(message,'(a/a/)') 'FATAL:  Fractional volume arrays not defined.', &
               'Abort from routine GRID_VOL_FRACS'
          call write_msg(message(:3))
          if (.not.present(ier)) stop
          ier=-1
          return
       endif
       
       where (int_vols%Vol_eltA > 0.d0) 
          vf_a=int_vols%VolinBmesh_eltA/int_vols%Vol_eltA
       elsewhere
          vf_a=0.d0
       endwhere

       where (int_vols%Vol_eltB > 0.d0)
          vf_b=int_vols%VolinAmesh_eltB/int_vols%Vol_eltB
       elsewhere
          vf_b=0.d0
       endwhere

    endif

  end subroutine grid_vol_fracs

  subroutine grid_vols(int_vols,vol_a,vol_b,ier)

    type(grid_int_vols), intent(in) :: int_vols
    real(dp), dimension(:), intent(out) :: vol_a,vol_b
    integer, intent(out), optional :: ier

    if (present(ier)) ier=0

    if (size(vol_a).ne.int_vols%BA%ncols) then 
       write(message,'(a/a/)') 'FATAL:  Element volume array VOL_A has incorrect length.', &
            'Abort from routine GRID_VOLS'
       call write_msg(message(:3))
       if (.not.present(ier)) stop
       ier=-1
       return
    endif
    if (size(vol_b).ne.int_vols%BA%nrows) then 
       write(message,'(a/a/)') 'FATAL:  Element volume array VOL_B has incorrect length.', &
            'Abort from routine GRID_VOLS'
       call write_msg(message(:3))
       if (.not.present(ier)) stop
       ier=-1
       return
    endif

    if (.not.associated(int_vols%Vol_eltA).or. &
         .not.associated(int_vols%Vol_eltB)) then
       write(message,'(a/a/)') 'FATAL: Element volume arrays not defined.', &
            'Abort from routine GRID_VOLS'
       call write_msg(message(:3))
       if (.not.present(ier)) stop
       ier=-1
       return
    endif

    Vol_A=int_vols%Vol_eltA
    Vol_B=int_vols%Vol_eltB

  end subroutine grid_vols

  subroutine write_int_volumes(int_vols, iunit)

    type(grid_int_vols), intent(in) :: int_vols
    integer, intent(in) :: iunit

    integer :: i

    write(iunit) int_vols%BA%nrows, int_vols%BA%ncols, int_vols%BA%numentries
    write(iunit) (int_vols%BA%irowst(i),i=1,int_vols%BA%nrows+1)
    write(iunit) (int_vols%BA%jcn(i),i=1,int_vols%BA%numentries)
    write(iunit) (int_vols%BA%val(i),i=1,int_vols%BA%numentries)
    write(iunit) (int_vols%Vol_eltA(i), i=1,int_vols%BA%ncols)
    write(iunit) (int_vols%Vol_eltB(i), i=1,int_vols%BA%nrows)
    write(iunit) (int_vols%VolinBmesh_eltA(i), i=1,int_vols%BA%ncols)
    write(iunit) (int_vols%VolinAmesh_eltB(i), i=1,int_vols%BA%nrows)
    write(iunit) (int_vols%VolBlockinBmesh_eltA(i), i=1,int_vols%BA%ncols)
    write(iunit) (int_vols%VolBlockinAmesh_eltB(i), i=1,int_vols%BA%nrows)
    write(iunit) int_vols%nodes_per_elt_Amesh,int_vols%nodes_per_elt_Bmesh
    write(iunit) int_vols%checksum_Amesh
    write(iunit) int_vols%checksum_Bmesh

  end subroutine write_int_volumes

  subroutine read_int_volumes(int_vols, iunit, ier)

    type(grid_int_vols), intent(inout) :: int_vols
    integer, intent(in) :: iunit
    integer, optional, intent(out) :: ier

    integer :: i,ierr

    if (present(ier)) ier=0

    call destroy_grid_int_vols(int_vols)

    read(iunit,IOSTAT=ierr) int_vols%BA%nrows, int_vols%BA%ncols, int_vols%BA%numentries
    if (ierr.eq.0) allocate(int_vols%BA%irowst(int_vols%BA%nrows+1),STAT=ierr)
    if (ierr.eq.0) read(iunit,IOSTAT=ierr) (int_vols%BA%irowst(i),i=1,int_vols%BA%nrows+1)
    if (ierr.eq.0) allocate(int_vols%BA%jcn(int_vols%BA%numentries),STAT=ierr)
    if (ierr.eq.0) read(iunit,IOSTAT=ierr) (int_vols%BA%jcn(i),i=1,int_vols%BA%numentries)
    if (ierr.eq.0) allocate(int_vols%BA%val(int_vols%BA%numentries),STAT=ierr)
    if (ierr.eq.0) read(iunit,IOSTAT=ierr) (int_vols%BA%val(i),i=1,int_vols%BA%numentries)
    if (ierr.eq.0) allocate(int_vols%Vol_eltA(int_vols%BA%ncols),STAT=ierr)
    if (ierr.eq.0) read(iunit,IOSTAT=ierr) (int_vols%Vol_eltA(i),i=1,int_vols%BA%ncols)
    if (ierr.eq.0) allocate(int_vols%Vol_eltB(int_vols%BA%nrows),STAT=ierr)
    if (ierr.eq.0) read(iunit,IOSTAT=ierr) (int_vols%Vol_eltB(i),i=1,int_vols%BA%nrows)
    if (ierr.eq.0) allocate(int_vols%VolinBmesh_eltA(int_vols%BA%ncols),STAT=ierr)
    if (ierr.eq.0) read(iunit,IOSTAT=ierr) (int_vols%VolinBmesh_eltA(i),i=1,int_vols%BA%ncols)
    if (ierr.eq.0) allocate(int_vols%VolinAmesh_eltB(int_vols%BA%nrows),STAT=ierr)
    if (ierr.eq.0) read(iunit,IOSTAT=ierr) (int_vols%VolinAmesh_eltB(i),i=1,int_vols%BA%nrows)
    if (ierr.eq.0) allocate(int_vols%VolBlockinBmesh_eltA(int_vols%BA%ncols),STAT=ierr)
    if (ierr.eq.0) read(iunit,IOSTAT=ierr) (int_vols%VolBlockinBmesh_eltA(i),i=1,int_vols%BA%ncols)
    if (ierr.eq.0) allocate(int_vols%VolBlockinAmesh_eltB(int_vols%BA%nrows),STAT=ierr)
    if (ierr.eq.0) read(iunit,IOSTAT=ierr) (int_vols%VolBlockinAmesh_eltB(i),i=1,int_vols%BA%nrows)
    if (ierr.eq.0) read(iunit,IOSTAT=ierr) int_vols%nodes_per_elt_Amesh,int_vols%nodes_per_elt_Bmesh
    if (ierr.eq.0) read(iunit,IOSTAT=ierr) int_vols%checksum_Amesh
    if (ierr.eq.0) read(iunit,IOSTAT=ierr) int_vols%checksum_Bmesh

    if (ierr.ne.0) then
       write(message,'(a/a/)') 'FATAL: Problem reading intersection volumes.', &
            'Problem in routine READ_INT_VOLUMES'
       call write_msg(message(:3))
       if (.not.present(ier)) stop
       ier=-1
       return
    endif       

  end subroutine read_int_volumes

  function right_int_volumes(mesh_a, mesh_b, int_vols)

    logical :: right_int_volumes
    type(gm_mesh), intent(in) :: mesh_a, mesh_b 
    type(grid_int_vols), intent(in) :: int_vols

    type(gm_mesh_checksum) :: checksum
    
    right_int_volumes=.false.

    call get_checksum(mesh_a,checksum)
    if (.not.checksums_equal(int_vols%checksum_Amesh,checksum)) return
    call get_checksum(mesh_b,checksum)
    if (.not.checksums_equal(int_vols%checksum_Bmesh,checksum)) return

    right_int_volumes=.true.

  end function right_int_volumes

  subroutine get_checksum(mesh,checksum)
  
    type(gm_mesh), intent(in) :: mesh
    type(gm_mesh_checksum), intent(out) :: checksum

    real(dp) :: rind
    integer :: i,j

    checksum%nnod=mesh%nnod
    checksum%nelt=mesh%nelt
    checksum%nodes_per_elt=size(mesh%node_elt,1)

    checksum%pos_node_0=0.d0
    checksum%pos_node_1=0.d0
    rind=0.d0
    do i=1,checksum%nnod
       do j=1,3
          rind=rind+1.d0
          checksum%pos_node_0=checksum%pos_node_0+real(mesh%pos_node(j,i),kind=dp)
          checksum%pos_node_1=checksum%pos_node_1+rind*real(mesh%pos_node(j,i),kind=dp)
       enddo
    enddo

    checksum%node_elt_0=0.d0
    checksum%node_elt_1=0.d0
    rind=0.d0
    do i=1,checksum%nelt
       do j=1,checksum%nodes_per_elt
          rind=rind+1.d0
          checksum%node_elt_0=checksum%node_elt_0+real(mesh%node_elt(j,i),kind=dp)
          checksum%node_elt_1=checksum%node_elt_1+rind*real(mesh%node_elt(j,i),kind=dp)
       enddo
    enddo

    checksum%block_elt_0=0.d0
    checksum%block_elt_1=0.d0
    rind=0.d0
    do i=1,checksum%nelt
       rind=rind+1.d0
       checksum%block_elt_0=checksum%block_elt_0+real(mesh%block_elt(i),kind=dp)
       checksum%block_elt_1=checksum%block_elt_1+rind*real(mesh%block_elt(i),kind=dp)
    enddo

  end subroutine get_checksum

  function checksums_equal(checksum1,checksum2)

  logical :: checksums_equal
  type(gm_mesh_checksum) :: checksum1,checksum2

  checksums_equal=.false.

  if (checksum1%nnod.ne.checksum2%nnod) return
  if (checksum1%nelt.ne.checksum2%nelt) return
  if (checksum1%nodes_per_elt.ne.checksum2%nodes_per_elt) return
  if (checksum1%pos_node_0.ne.checksum2%pos_node_0) return
  if (checksum1%pos_node_1.ne.checksum2%pos_node_1) return
  if (checksum1%node_elt_0.ne.checksum2%node_elt_0) return
  if (checksum1%node_elt_1.ne.checksum2%node_elt_1) return
  if (checksum1%block_elt_0.ne.checksum2%block_elt_0) return
  if (checksum1%block_elt_1.ne.checksum2%block_elt_1) return

  checksums_equal=.true.

  end function checksums_equal

  subroutine destroy_grid_int_vols(int_vols)

    type(grid_int_vols), intent(inout) :: int_vols

    if (associated(int_vols%BA%irowst)) deallocate(int_vols%BA%irowst)
    if (associated(int_vols%BA%jcn)) deallocate(int_vols%BA%jcn)
    if (associated(int_vols%BA%val)) deallocate(int_vols%BA%val)
    int_vols%BA%nrows=0
    int_vols%BA%ncols=0
    int_vols%BA%numentries=0
    if (associated(int_vols%Vol_eltA)) deallocate(int_vols%Vol_eltA)
    if (associated(int_vols%Vol_eltB)) deallocate(int_vols%Vol_eltB)
    if (associated(int_vols%VolinBmesh_eltA)) deallocate(int_vols%VolinBmesh_eltA)
    if (associated(int_vols%VolinAmesh_eltB)) deallocate(int_vols%VolinAmesh_eltB)
    if (associated(int_vols%VolBlockinBmesh_eltA)) deallocate(int_vols%VolBlockinBmesh_eltA)
    if (associated(int_vols%VolBlockinAmesh_eltB)) deallocate(int_vols%VolBlockinAmesh_eltB)
    int_vols%nodes_per_elt_Amesh=0
    int_vols%nodes_per_elt_Bmesh=0

    int_vols%checksum_Amesh%nnod=0
    int_vols%checksum_Amesh%nelt=0
    int_vols%checksum_Amesh%nodes_per_elt=0
    int_vols%checksum_Amesh%pos_node_0=0.d0
    int_vols%checksum_Amesh%pos_node_1=0.d0
    int_vols%checksum_Amesh%node_elt_0=0.d0
    int_vols%checksum_Amesh%node_elt_1=0.d0
    int_vols%checksum_Amesh%block_elt_0=0.d0
    int_vols%checksum_Amesh%block_elt_1=0.d0

    int_vols%checksum_Bmesh%nnod=0
    int_vols%checksum_Bmesh%nelt=0
    int_vols%checksum_Bmesh%nodes_per_elt=0
    int_vols%checksum_Bmesh%pos_node_0=0.d0
    int_vols%checksum_Bmesh%pos_node_1=0.d0
    int_vols%checksum_Bmesh%node_elt_0=0.d0
    int_vols%checksum_Bmesh%node_elt_1=0.d0
    int_vols%checksum_Bmesh%block_elt_0=0.d0
    int_vols%checksum_Bmesh%block_elt_1=0.d0

  end subroutine destroy_grid_int_vols

end module GRID_MAPPING_MODULE

