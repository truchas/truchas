!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module grid_mapping_utils
  !=======================================================================
  ! Purpose(s):
  !
  !    Provide miscellaneous routines used in conjunction with grid
  !    to grid mapping.
  !
  !    Public Interface(s):
  !
  !      write_msg
  !      writegmv_gm_mesh
  !      read_ascii_gm_mesh
  !      write_ascii_gm_mesh
  !
  ! Author(s): Andrew Kuprat (kuprat@lanl.gov)
  !
  !=======================================================================

  use gm_mesh_type, only: gm_mesh

  implicit none

  PRIVATE

  PUBLIC :: write_msg
  PUBLIC :: writegmv_gm_mesh
  PUBLIC :: read_ascii_gm_mesh
  PUBLIC :: write_ascii_gm_mesh

contains

  subroutine write_msg(message)
!
! Subroutine called by GRID_MAPPING_MODULE for warning/error output.
! Replace or revise this subroutine if a simple 'print' is 
! inadequate for this purpose.
!
    character(len=*),dimension(:),intent(in) :: message
    print '(a)', message
  end subroutine write_msg

  subroutine writegmv_gm_mesh(iunit,mesh1,mesh2, &
       field,fieldname,fields,fieldnames,matnames1,matnames2, &
       probtime,cycleno,ier)
!
! This routine takes one or two GM_MESH derived type data and writes
! out a GMV graphics file.  If both MESH1 and MESH2 are present, then
! the two meshes are superimposed, so that GMV will reveal how good
! a match the two meshes are to each other.  (When mapping between
! grids, they are supposed to be a reasonable match so that, e.g., 
! all the centroids of the elements of MESH1 reside within the 
! space taken up by MESH2.)  
!
! If only MESH1 is present, then one of FIELD or FIELDS can be present.
! If FIELD is present, the GMV file contains the cell field with data
! given by FIELD and the name of this FIELD should be specified by
! the user using FIELDNAME.  If FIELDS is present, multiple fields are
! present in the GMV file, and the corresponding names should be
! provided by the user in FIELDNAMES.
! 
! Optionally the user can specify material names for MESH1 using
! MATNAMES1 (and for MESH2 using MATNAMES2 if present).  The user
! can also optionally provide PROBTIME and CYCLENO.
!
! If IER is present, the subroutine will return IER=-1 upon error;
! otherwise it will stop.
!
    integer, parameter :: dp = kind(1.0d0)

    integer, intent(in) :: iunit
    type(gm_mesh), intent(in) :: mesh1
    type(gm_mesh), optional, intent(in) :: mesh2
    real(dp), optional, dimension(:), intent(in) :: field
    real(dp), optional, dimension(:,:), intent(in) :: fields
    character(LEN=*), optional, intent(in) :: fieldname
    character(LEN=*), optional, dimension(:), intent(in) :: fieldnames
    character(LEN=*), optional, dimension(:), intent(in) :: matnames1, matnames2
    real(dp), optional, intent(in) :: probtime
    integer, optional, intent(in) :: cycleno
    integer, optional, intent(out) :: ier

    integer i,j,nodes_per_elt,datatype,nmat,nmat1,nmat2
    character(len=72), dimension(3) :: message
    
    if (present(ier)) then
       ier=0
    endif

    if (present(field).and.present(fields)) then
       write(message,'(a/a/)') 'FATAL:  FIELD and FIELDS cannot both be present.', &
            'Abort from routine WRITEGM_GM_MESH.'
       call write_msg(message(:3))
       if (present(ier)) then
          ier=-1
          return
       else
          stop
       endif
    endif
    if (present(mesh2)) then
       if (present(field).or.present(fields)) then
          write(message,'(a/a/)') &
               'FATAL:  When MESH2 is present, neither FIELD or FIELDS can be present.', &
               'Abort from routine WRITE_GM_MESH.'
          call write_msg(message(:3))
          if (present(ier)) then
             ier=-1
             return
          else
             stop
          endif
       endif
    endif

    write(iunit,'(a)') 'gmvinput ascii'
    write(iunit,'(a)',advance='no') 'nodev '
    if (present(mesh2)) then
       write(iunit,*) mesh1%nnod+mesh2%nnod
    else
       write(iunit,*) mesh1%nnod
    endif
    do i=1,mesh1%nnod
       write(iunit,*) (mesh1%pos_node(j,i),j=1,3)
    enddo
    if (present(mesh2)) then
       do i=1,mesh2%nnod
          write(iunit,*) (mesh2%pos_node(j,i),j=1,3)
       enddo
    endif
    write(iunit,'(a)',advance='no') 'cells '
    if (present(mesh2)) then
       write(iunit,*) mesh1%nelt+mesh2%nelt
    else
       write(iunit,*) mesh1%nelt
    endif
    nodes_per_elt=size(mesh1%node_elt,1)
    do i=1,mesh1%nelt
       if (nodes_per_elt.eq.4) then
          write(iunit,'(a)',advance='no') 'tet '
       else
          write(iunit,'(a)',advance='no') 'hex '
       endif
       write(iunit,*) nodes_per_elt,(mesh1%node_elt(j,i),j=1,nodes_per_elt)
    enddo
    if (present(mesh2)) then
       nodes_per_elt=size(mesh2%node_elt,1)
       do i=1,mesh2%nelt
          if (nodes_per_elt.eq.4) then
             write(iunit,'(a)',advance='no') 'tet '
          else
             write(iunit,'(a)',advance='no') 'hex '
          endif
          write(iunit,*) nodes_per_elt,(mesh1%nnod+mesh2%node_elt(j,i),j=1,nodes_per_elt)
       enddo
    endif

    write(iunit,'(a)',advance='no') 'material '

    if (present(mesh2)) then
       if (minval(mesh1%block_elt).lt.1) then
          write(message,'(a/a/)') &
               'FATAL:  Nonpositive block number specified.', &
               'Abort from routine WRITE_GM_MESH.'
          call write_msg(message(:3))
          if (present(ier)) then
             ier=-1
             return
          else
             stop
          endif
       endif

       nmat1=maxval(mesh1%block_elt)
       if (present(matnames1)) then
          if (size(matnames1).lt.nmat1) then
             write(message,'(a/a/)') &
               'FATAL: Max blk no. exceeds no. of mat names', &
               'Abort from routine WRITE_GM_MESH.'
             call write_msg(message(:3))
             if (present(ier)) then
                ier=-1
                return
             else
                stop
             endif
          else
             nmat1=max(nmat1,size(matnames1))
          endif
       endif

       if (minval(mesh2%block_elt).lt.1) then
          write(message,'(a/a/)') &
               'FATAL:  Nonpositive block number specified.', &
               'Abort from routine WRITE_GM_MESH.'
          call write_msg(message(:3))
          if (present(ier)) then
             ier=-1
             return
          else
             stop
          endif
       endif

       nmat2=maxval(mesh2%block_elt)
       if (present(matnames2)) then
          if (size(matnames2).lt.nmat2) then
             write(message,'(a/a/)') &
               'FATAL: Max blk no. exceeds no. of mat names', &
               'Abort from routine WRITE_GM_MESH.'
             call write_msg(message(:3))
             if (present(ier)) then
                ier=-1
                return
             else
                stop
             endif
          else
             nmat2=max(nmat2,size(matnames2))
          endif
       endif

       datatype=0 ! 0=cells; 1=nodes
       write(iunit,*) nmat1+nmat2,datatype

       if (present(matnames1)) then
          do i=1,nmat1
             write(iunit,'(a)') trim(matnames1(i))
          enddo
       else
          do i=1,nmat1
             if (i.le.9) then
                write(iunit,'(a,i1)') 'fmat_',i
             else
                write(iunit,'(a,i2)') 'fmat_',i
             endif
          enddo
       endif

       if (present(matnames2)) then
          do i=1,nmat2
             write(iunit,'(a)') trim(matnames2(i))
          enddo
       else
          do i=1,nmat2
             if (i.le.9) then
                write(iunit,'(a,i1)') 'smat_',i
             else
                write(iunit,'(a,i2)') 'smat_',i
             endif
          enddo
       endif
       do i=1,mesh1%nelt
          write(iunit,*) mesh1%block_elt(i)
       enddo
       do i=1,mesh2%nelt
          write(iunit,*) nmat1+mesh2%block_elt(i)
       enddo
    else
       if (minval(mesh1%block_elt).lt.1) then
          write(message,'(a/a/)') &
               'FATAL:  Nonpositive block number specified.', &
               'Abort from routine WRITE_GM_MESH.'
          call write_msg(message(:3))
          if (present(ier)) then
             ier=-1
             return
          else
             stop
          endif
       endif

       nmat=maxval(mesh1%block_elt)
       if (present(matnames1)) then
          if (size(matnames1).lt.nmat) then
             write(message,'(a/a/)') &
               'FATAL: Max blk no. exceeds no. of mat names', &
               'Abort from routine WRITE_GM_MESH.'
             call write_msg(message(:3))
             if (present(ier)) then
                ier=-1
                return
             else
                stop
             endif
          else
             nmat=max(nmat,size(matnames1))
          endif
       endif

       datatype=0 ! 0=cells; 1=nodes
       write(iunit,*) nmat,datatype

       if (present(matnames1)) then
          do i=1,nmat
             write(iunit,'(a)') trim(matnames1(i))
          enddo
       else
          do i=1,nmat
             if (i.le.9) then
                write(iunit,'(a,i1)') 'mat_',i
             else
                write(iunit,'(a,i2)') 'mat_',i
             endif
          enddo
       endif

       do i=1,mesh1%nelt
          write(iunit,*) mesh1%block_elt(i)
       enddo
    endif

    if (present(field)) then
       write(iunit,'(a)') 'variable'
       datatype=0 ! 0=cells; 1=nodes; 2=faces
       if (present(fieldname)) then
          write(iunit,'(a)',advance='no') trim(fieldname)//' '
          write(iunit,*) datatype
       else
          write(iunit,'(a)',advance='no') 'cfield '
          write(iunit,*) datatype
       endif
       do i=1,mesh1%nelt
          write(iunit,*) field(i)
       enddo
       write(iunit,'(a)') 'endvars'
    elseif (present(fields)) then
       write(iunit,'(a)') 'variable'
       datatype=0 ! 0=cells; 1=nodes; 2=faces
       if (present(fieldnames)) then
          if (size(fields,2).ne.size(fieldnames)) then
             write(message,'(a/a/)') &
               'FATAL: No. of field names doesn''t equal no. of fields', &
               'Abort from routine WRITE_GM_MESH.'
             call write_msg(message(:3))
             if (present(ier)) then
                ier=-1
                return
             else
                stop
             endif
          endif
          do i=1,size(fields,2)
             write(iunit,'(a,a,i1)') trim(fieldnames(i)),' ',datatype
             do j=1,mesh1%nelt
                write(iunit,*) fields(j,i)
             enddo
          enddo
          write(iunit,'(a)') 'endvars'
       else
          do i=1,size(fields,2)
             if (i.le.9) then
                write(iunit,'(a,i1,a,i1)') 'cfld_',i,' ',datatype
             else
                write(iunit,'(a,i2,a,i1)') 'cfld_',i,' ',datatype
             endif
             do j=1,mesh1%nelt
                write(iunit,*) fields(j,i)
             enddo
          enddo
          write(iunit,'(a)') 'endvars'
       endif
    endif
       
    if (present(probtime)) then
       write(iunit,*) 'probtime ',probtime
    endif
    if (present(cycleno)) then
       write(iunit,*) 'cycleno ',cycleno
    endif

    write(iunit,'(a)') 'endgmv'

  end subroutine writegmv_gm_mesh

  subroutine read_ascii_gm_mesh(iunit,mesh)
!
! Simple subroutine for reading ascii mesh data to fill
! an object of type GM_MESH.
!
    type(gm_mesh), intent(out) :: mesh
    integer, intent(in) :: iunit

    integer :: i,j,nodes_per_elt

    read(iunit,*) mesh%nelt,mesh%nnod,nodes_per_elt

    allocate(mesh%node_elt(nodes_per_elt,mesh%nelt))
    do i=1,mesh%nelt
       read(iunit,*) (mesh%node_elt(j,i),j=1,nodes_per_elt)
    enddo

    allocate(mesh%pos_node(3,mesh%nnod))
    do i=1,mesh%nnod
       read(iunit,*) (mesh%pos_node(j,i),j=1,3)
    enddo

    ! Assume block elt info is all 1.
    allocate(mesh%block_elt(mesh%nelt))
    do i=1,mesh%nelt
       mesh%block_elt(i)=1
    enddo

  end subroutine read_ascii_gm_mesh

  subroutine write_ascii_gm_mesh(iunit,mesh)
!
! Simple subroutine for writing an ascii data file corresponding
! to an object of type GM_MESH.
!
    type(gm_mesh), intent(in) :: mesh
    integer, intent(in) :: iunit

    integer :: i,j,nodes_per_elt

    nodes_per_elt=size(mesh%node_elt,1)

    write(iunit,*) mesh%nelt,mesh%nnod,nodes_per_elt

    do i=1,mesh%nelt
       write(iunit,*) (mesh%node_elt(j,i),j=1,nodes_per_elt)
    enddo

    do i=1,mesh%nnod
       write(iunit,*) (mesh%pos_node(j,i),j=1,3)
    enddo

    ! Right now we don't write out block_elt

  end subroutine write_ascii_gm_mesh

end module grid_mapping_utils
