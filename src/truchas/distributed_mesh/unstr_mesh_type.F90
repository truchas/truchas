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

    integer, allocatable :: cfpar(:)  ! relative cell face orientation (bit mask)
    
    real(r8), allocatable :: normal(:,:)

    !! Mesh interface links.
    integer :: nlink = 0, nlink_onP = 0
    integer, pointer :: lface(:,:) => null()  ! pointer due to localize_index_array
    integer, allocatable :: link_set_id(:)    ! user-assigned ID for each link block
    type(bitfield), allocatable :: link_set_mask(:)  ! link block index
    type(ip_desc) :: link_ip
  contains
    procedure :: get_global_cnode_array
!    procedure :: get_global_cface_array
    procedure :: compute_geometry
    procedure :: write_profile
    final :: unstr_mesh_delete
  end type unstr_mesh

contains

  !! Final subroutine for UNSTR_MESH objects.
  subroutine unstr_mesh_delete (this)
    type(unstr_mesh), intent(inout) :: this
    if (associated(this%lface)) deallocate(this%lface)
    call destroy (this%node_ip)
    call destroy (this%face_ip)
    call destroy (this%cell_ip)
    call destroy (this%link_ip)
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

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! GET_GLOBAL_CNODE_ARRAY
 !! GET_GLOBAL_CFACE_ARRAY
 !!
 !! These routines return the so-named global index array that is associated
 !! with the mesh.  The collated global index array is returned on the IO
 !! processor and a 0-sized array on all others.  The returned arrays are
 !! allocated by these procedures.
 !!

  subroutine get_global_cnode_array (this, xcnode, cnode)
    class(unstr_mesh), intent(in) :: this
    integer, allocatable, intent(out) :: xcnode(:), cnode(:)
    associate (xcnode_onP => this%xcnode(:this%ncell_onP+1), &
                cnode_onP => this%cnode(:this%xcnode(this%ncell_onP+1)-1))
      call get_global_ragged_array (xcnode_onP, this%node_ip%global_index(cnode_onP), xcnode, cnode)
    end associate
    ASSERT(minval(cnode) >= 1 .and. maxval(cnode) <= this%node_ip%global_size())
  end subroutine get_global_cnode_array

  subroutine get_global_ragged_array (xarray_l, array_l, xarray, array)
    integer, intent(in) :: xarray_l(:), array_l(:)
    integer, allocatable, intent(out) :: xarray(:), array(:)
    integer :: offset
    ASSERT(size(xarray_l) >= 1)
    ASSERT(size(array_l) == xarray_l(size(xarray_l))-1)
    ASSERT(global_all(xarray_l(2:) - xarray_l(:size(xarray_l)-1) >= 4))
    ASSERT(global_all(xarray_l(2:) - xarray_l(:size(xarray_l)-1) <= 8))
#ifdef NAG_COMPILER_WORKAROUND
    !! NAGFOR 6.0(1052)/GCC 4.8.3 produce bad code under optimization (-O2)
    !! with the one-line allocation in parallel, getting the wrong value for
    !! the global_sum expression, apparently it gets the value of its argument
    !! instead.  Looking at the intermediate C code suggests that tha problem
    !! may not be here but in global_sum and perhaps a race condition?
    offset = global_sum(size(xarray_l)-1)
    allocate(xarray(1+merge(offset,0,is_IOP)))
#else
    allocate(xarray(1+merge(global_sum(size(xarray_l)-1),0,is_IOP)))
#endif
    offset = excl_prefix_sum(size(array_l))
    xarray(1) = 1
    call collate (xarray(2:), xarray_l(2:)+offset)
#ifdef NAG_COMPILER_WORKAROUND
    !! Same comments as above.
    offset = global_sum(size(array_l))
    allocate(array(merge(offset,0,is_IOP)))
#else
    allocate(array(merge(global_sum(size(array_l)),0,is_IOP)))
#endif
    call collate (array, array_l)
    if (is_IOP) then
      ASSERT(size(xarray) >= 1)
      ASSERT(xarray(1) == 1)
      ASSERT(all(xarray(2:) - xarray(:size(xarray)-1) >= 4))
      ASSERT(all(xarray(2:) - xarray(:size(xarray)-1) <= 8))
      ASSERT(size(array) == xarray(size(xarray))-1)
    end if
  contains
    integer function excl_prefix_sum (n) result (psum)
      use parallel_communication, only: nPE, is_IOP, collate, distribute
      integer, intent(in) :: n
      integer :: j
      integer, allocatable :: array(:)
      allocate(array(merge(nPE,0,is_IOP)))
      call collate (array, n)
      if (is_IOP) then
        do j = 2, nPE
          array(j) = array(j) + array(j-1)
        end do
      end if
      call distribute (psum, array)

      psum = psum - n
    end function
  end subroutine get_global_ragged_array

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

end module unstr_mesh_type
