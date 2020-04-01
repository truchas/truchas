!!
!! VF_MATRIX_PATCH_TYPE
!!
!! A concrete implementation for the abstract base class VF_MATRIX.
!! This implementation operates on patch-based view factor data.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 18 July 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#include "f90_assert.fpp"

module vf_matrix_patch_type

  use kinds, only: r8
  use vf_matrix_class, only: vf_matrix
  use rad_encl_type
  use rad_encl_file_type
  use parallel_communication
  implicit none

  type, extends(vf_matrix), public :: vf_matrix_patch
    integer :: npatch       ! number of patches (number of rows) on this process
    integer :: npatch_tot   ! total number of patches (number of columns)
    integer, allocatable :: f2p_map(:)    ! maps local faces to local patches
    integer, allocatable :: f2p_map_g(:)  ! maps global faces to global patches
                                          !   only available on IO process
    real(r8), allocatable :: w(:)         ! w(j) is the fraction of total patch
                                          !  area occupied by face j
    contains
      procedure, public  :: init => init_VFP_PATCH
      procedure, public  :: partition_ER_faces => partition_ER_faces_VFP_PATCH
      procedure, public  :: load_view_factors => load_view_factors_VFP_PATCH
      procedure, public  :: phi_x => phi_x_VFP_PATCH
      procedure, private :: get_face_weights
  end type

contains


  !! Initialize VF_MATRIX_PATCH and define a distribution for the matrix rows.
  subroutine init_VFP_PATCH (this, file)

    use,intrinsic :: iso_c_binding, only: c_size_t

    class(vf_matrix_patch), intent(out) :: this
    type(rad_encl_file), intent(in) :: file

    integer :: nface_tot, npatch_tot
    integer(c_size_t) :: nnonz
    integer :: n

    if (is_IOP) call file%get_vf_dims (nface_tot, npatch_tot, nnonz)

    call broadcast(nface_tot)
    call broadcast(npatch_tot)
    this%nface_tot = nface_tot
    this%npatch_tot = npatch_tot

    !! Partition matrix rows in blocks
    this%npatch = this%npatch_tot/nPE
    if (this_PE <= modulo(this%npatch_tot,nPE)) this%npatch = 1 + this%npatch

    !! Get global face to patch map
    n = merge(this%nface_tot, 0, is_IOP)
    allocate(this%f2p_map_g(n))
    if (is_IOP) call file%get_f2p_map(this%f2p_map_g)

  end subroutine init_VFP_PATCH


  !! Partitions the ER mesh according to the distribution of matrix rows
  !! defined by THIS%INIT.
  subroutine partition_ER_faces_VFP_PATCH (this, color)

    class(vf_matrix_patch), intent(inout) :: this
    integer, allocatable, intent(out) :: color(:)

    integer, allocatable :: color_p(:)    ! global patch coloring
    integer, allocatable :: color_p_l(:)  ! local patch coloring
    integer :: n

    !! Block coloring of the enclosure patches.
    n = merge(this%npatch_tot, 0, is_IOP)
    allocate(color_p(n), color_p_l(this%npatch))
    color_p_l = this_PE
    call collate (color_p, color_p_l)
    deallocate(color_p_l)

    n = merge(this%nface_tot, 0, is_IOP)
    allocate (color(n))
    if (is_IOP) color = color_p(this%f2p_map_g)

  end subroutine partition_ER_faces_VFP_PATCH


  !! Read and distribute the VF matrix data according to the distribution
  !! defined in THIS%INIT.
  subroutine load_view_factors_VFP_PATCH (this, file, encl)

    class(vf_matrix_patch), intent(inout) :: this
    type(rad_encl_file), intent(in) :: file
    type(rad_encl), intent(in) :: encl

    real, allocatable :: amb_vf_face(:)
    integer, allocatable :: f2p_map_col(:)  ! collated map from local faces to global patches
    integer, allocatable :: fidx_col(:)     ! collated map from local faces to global faces
    integer :: n, bsize(nPE), patch_offset

    call this%distribute_vf_rows(file, this%npatch, this%npatch_tot)

    this%nface = encl%nface_onP

    !! Get face to patch map
    n = merge(this%nface_tot, 0, is_IOP)
    allocate(f2p_map_col(n), fidx_col(n))

    call collate (fidx_col, encl%face_map)
    if (is_IOP) f2p_map_col = this%f2p_map_g(fidx_col)

    allocate(this%f2p_map(this%nface))
    call distribute (this%f2p_map, f2p_map_col)

    !! Convert global patch indices to local
    call collate (bsize, this%npatch)
    call broadcast (bsize)
    patch_offset = sum(bsize(1:this_PE-1))
    this%f2p_map = this%f2p_map - patch_offset
    ASSERT ( all(1 <= this%f2p_map) )
    ASSERT ( all(this%f2p_map <= this%npatch) )

    !! Compute face weights
    this%w = this%get_face_weights(encl)

    !! Grow patch-based amb_vf into face-based version
    if (allocated(this%amb_vf)) then
      ASSERT( size(this%amb_vf) == this%npatch )
      allocate(amb_vf_face(this%nface))
      amb_vf_face = this%amb_vf( this%f2p_map )
      call move_alloc(amb_vf_face, this%amb_vf)
    end if

  end subroutine load_view_factors_VFP_PATCH


  !! Computes the global product PHI*x, where x is a local vector.
  subroutine phi_x_VFP_PATCH (this, lhs, x)

    class(vf_matrix_patch), intent(in) :: this
    real(r8), intent(out) :: lhs(:)
    real(r8), intent(in) :: x(:)

    real(r8) :: xp(this%npatch)             ! Local patch-based x
    real(r8) :: global_xp(this%npatch_tot)  ! Global patch-based x
    real(r8) :: lhsp(this%npatch)           ! Local patch-based result
    integer :: i, j

    ASSERT(size(lhs) == this%nface)
    ASSERT(size(x) == this%nface)

    !! Compute intermediate patch-based result
    xp = 0.0_r8
    do i = 1, this%nface
      xp(this%f2p_map(i)) = xp(this%f2p_map(i)) + this%w(i)*x(i)
    end do

    call collate (global_xp, xp)
    call broadcast (global_xp)

    !! Compute patch-based mat-vec
    lhsp = 0.0_r8
    do i = 1, this%npatch
      do j = this%ia(i), this%ia(i+1)-1
        lhsp(i) = lhsp(i) + this%vf(j) * global_xp(this%ja(j))
      end do
    end do

    !! Compute face-based result
    lhs = lhsp( this%f2p_map )

  end subroutine phi_x_VFP_PATCH


  !! Computes the fraction of total patch area occupied by each face.
  function get_face_weights(this, encl) result(ret)

    use cell_geometry, only: face_normal, vector_length

    class(vf_matrix_patch), intent(in) :: this
    type(rad_encl), intent(in) :: encl
    real(r8), allocatable :: ret(:)

    integer :: i
    real(r8), allocatable :: f_area(:), p_area(:)

    !! Get face and patch areas
    allocate(f_area(this%nface),p_area(this%npatch))
    p_area = 0.0_r8
    do i = 1, this%nface
      associate(face_nodes => encl%fnode(encl%xface(i):encl%xface(i+1)-1))
        f_area(i) = vector_length( face_normal(encl%coord(:,face_nodes)) )
        p_area( this%f2p_map(i) ) = p_area( this%f2p_map(i) ) + f_area(i)
      end associate
    end do

    !! Compute weights
    allocate(ret(this%nface))
    do i = 1, this%nface
      ret(i) = f_area(i) / p_area( this%f2p_map(i) )
    end do
    ASSERT( all(0.0_r8 < ret) )
    ASSERT( all(ret <= 1.0_r8) )

    deallocate(f_area, p_area)

  end function get_face_weights

end module vf_matrix_patch_type
