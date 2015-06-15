!!
!! UNSTR_MESH_TYPE
!!
!! This module provides a derived type that encapsulates the data describing a
!! distributed unstructured mesh that is comprised either entirely of hex cells
!! or entirely of tet cells.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Revised May 2015
!!

#include "f90_assert.fpp"

module unstr_mesh_type

  use kinds, only: r8
  use base_mesh_class
  use parallel_communication
  use index_partitioning
  use bitfield_type
  implicit none
  private

  type, extends(base_mesh), public :: unstr_mesh
    integer, allocatable :: xcnode(:), cnode(:) ! cell nodes
    integer, allocatable :: xcface(:), cface(:) ! cell faces
    integer, allocatable :: xfnode(:), fnode(:) ! face nodes
  contains
!    procedure :: get_face_set_ids
!    procedure :: get_cell_set_bitmask
!    procedure :: get_face_set_bitmask
!    procedure :: get_global_cnode_array
!    procedure :: get_global_cface_array
!    procedure :: get_global_cblock_array
!    procedure :: get_global_x_array
!    procedure :: get_global_volume_array
!    procedure :: write_profile
    procedure :: compute_geometry
    procedure :: write_profile
    final :: unstr_mesh_delete
  end type unstr_mesh

contains

  !! Final subroutine for UNSTR_MESH objects.
  subroutine unstr_mesh_delete (this)
    type(unstr_mesh), intent(inout) :: this
    call destroy (this%node_ip)
    call destroy (this%face_ip)
    call destroy (this%cell_ip)
  end subroutine unstr_mesh_delete
  
  !! Compute the geometric data components from the node coordinates.
  subroutine compute_geometry (this)
    use cell_geometry, only: cell_volume, face_normal, vector_length
    class(unstr_mesh), intent(inout) :: this
    integer :: j
    ASSERT(allocated(this%volume))
    ASSERT(allocated(this%normal))
    ASSERT(allocated(this%area))
    do j = 1, this%ncell
      associate (cell_nodes => this%cnode(this%xcnode(j):this%xcnode(j+1)-1))
        this%volume(j) = cell_volume(this%x(:,cell_nodes))
      end associate
    end do
    do j = 1, this%nface
      associate(face_nodes => this%fnode(this%xfnode(j):this%xfnode(j+1)-1))
        this%normal(:,j) = face_normal(this%x(:,face_nodes))
        this%area(j) = vector_length(this%normal(:,j))
      end associate
    end do
  end subroutine compute_geometry

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! WRITE_DIST_MESH_PROFILE
 !!
 !! Writes to the tty and output file a profile of the distributed mesh:
 !! numbers of noded, edges, faces, and cells assigned to each processor;
 !! numbers of  on-process and off-process objects.
 !!

  subroutine write_profile (this)

    class(unstr_mesh), intent(in) :: this

    integer :: n
    character(len=80) :: line
    integer, dimension(nPE) :: nnode_vec, nface_vec, ncell_vec
    integer, dimension(2,nPE) :: nvec, fvec, cvec

    call collate (nnode_vec, this%nnode)
    call collate (nface_vec, this%nface)
    call collate (ncell_vec, this%ncell)

    call broadcast (nnode_vec)
    call broadcast (nface_vec)
    call broadcast (ncell_vec)

    call wline ('  UNSTR_MESH Profile:')
    write(line,fmt='(4x,a3,a,4a9)') 'PE', '|', 'nnode', 'nface', 'ncell'
    call wline (line)
    call wline ('    ---+'//repeat('-',27))
    do n = 1, nPE
      write(line,fmt='(4x,i3,a,3i9)') n, '|', nnode_vec(n),  nface_vec(n), ncell_vec(n)
      call wline (line)
    end do

    if (defined(this%node_ip)) then
      call collate (nvec(1,:), this%node_ip%offP_size())
      call collate (nvec(2,:), this%node_ip%onP_size())
      call broadcast (nvec)
    else
      nvec = 0
    end if

    if (defined(this%face_ip)) then
      call collate (fvec(1,:), this%face_ip%offP_size())
      call collate (fvec(2,:), this%face_ip%onP_size())
      call broadcast (fvec)
    else
      fvec = 0
    end if

    if (defined(this%cell_ip)) then
      call collate (cvec(1,:), this%cell_ip%offP_size())
      call collate (cvec(2,:), this%cell_ip%onP_size())
      call broadcast (cvec)
    else
      cvec = 0
    end if

    call wline ('  Mesh Communication Profile:')
    write(line,fmt='(4x,3x,1x,a11,2a16)')  'Nodes', 'Faces', 'Cells'
    call wline (line)
    write(line,fmt='(4x,a3,a,3a16)') 'PE', '|', ('off-PE   on-PE', n=1,3)
    call wline (line)
    call wline ('    ---+'//repeat('-',48))
    do n = 1, nPE
      write(line,fmt='(4x,i3,a,3(i7,i9))') n, '|', nvec(:,n), fvec(:,n), cvec(:,n)
      call wline (line)
    end do

  contains

    subroutine wline (line)
      use truchas_logging_services
      character(len=*), intent(in) :: line
      call TLS_info (line)
    end subroutine wline

  end subroutine write_profile

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!
! !! GET_GLOBAL_CNODE_ARRAY
! !! GET_GLOBAL_CFACE_ARRAY
! !!
! !! These routines return the so-named global index array that is associated
! !! with a distributed mesh.  The collated global index array is returned on
! !! the IO processor and a 0-sized array on all others.  The returned array
! !! pointer is allocated by these procedures.
! !!
! !! N.B.: Any storage the pointer may have been associated with at entry
! !! is _not_ deallocated before allocating the pointer anew.
! !!
!
!  subroutine get_global_cnode_array (this, cnode)
!    class(dist_mesh), intent(in) :: this
!    integer, allocatable, intent(out) :: cnode(:,:)
!    ASSERT(associated(this%cnode))
!    ASSERT(defined(this%cell_ip))
!    ASSERT(defined(this%node_ip))
!    ASSERT(size(this%cnode,2) == this%cell_ip%local_size())
!    ASSERT(minval(this%cnode) >= 1)
!    ASSERT(maxval(this%cnode) <= this%node_ip%local_size())
!    allocate(cnode(size(this%cnode,1),merge(this%cell_ip%global_size(),0,is_IOP)))
!    call collate (cnode, this%node_ip%global_index(this%cnode(:,:this%ncell_onP)))
!    ASSERT(minval(cnode) >= 1 .and. maxval(cnode) <= this%node_ip%global_size())
!  end subroutine get_global_cnode_array
!
!  subroutine get_global_cface_array (this, cface)
!    class(dist_mesh), intent(in) :: this
!    integer, allocatable, intent(out) :: cface(:,:)
!    ASSERT(associated(this%cface))
!    ASSERT(defined(this%cell_ip))
!    ASSERT(defined(this%face_ip))
!    ASSERT(size(this%cface,2) == this%cell_ip%local_size())
!    ASSERT(minval(this%cface) >= 1)
!    ASSERT(maxval(this%cface) <= this%face_ip%local_size())
!    allocate(cface(size(this%cface,1),merge(this%cell_ip%global_size(),0,is_IOP)))
!    call collate (cface, this%face_ip%global_index(this%cface(:,:this%ncell_onP)))
!    ASSERT(minval(cface) >= 1 .and. maxval(cface) <= this%face_ip%global_size())
!  end subroutine get_global_cface_array
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!
! !! GET_GLOBAL_CBLOCK_ARRAY
! !!
! !! This routine returns the global cell block ID array that is associated
! !! with a distributed mesh.  The collated global array is returned on the IO
! !! processor and a 0-sized array on all others.  The returned array pointer
! !! is allocated by this procedure, and any storage the pointer may have been
! !! associated with at entry is _not_ deallocated before allocating the pointer
! !! anew.
! !!
!
!  subroutine get_global_cblock_array (this, cblock)
!    class(dist_mesh), intent(in) :: this
!    integer, allocatable, intent(out) :: cblock(:)
!    ASSERT(allocated(this%cblock))
!    ASSERT(defined(this%cell_ip))
!    allocate(cblock(merge(this%cell_ip%global_size(),0,is_IOP)))
!    call collate (cblock, this%cblock(:this%ncell_onP))
!  end subroutine get_global_cblock_array

end module unstr_mesh_type
