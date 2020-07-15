!TODO: finish documentation
!! MFD_2D_DIFF_MATRIX_TYPE
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! May 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module mfd_2d_diff_matrix_type

  use kinds, only: r8
  use unstr_2d_mesh_type
  use mfd_2d_disc_type
  use pcsr_matrix_type
  use parallel_communication
  use index_partitioning
  implicit none
  private

  type, public :: mfd_2d_diff_matrix
    type(mfd_2d_disc),   pointer :: disc => null()  ! reference only -- do not own
    type(unstr_2d_mesh), pointer :: mesh => null()  ! reference only -- do not own
    !TODO: include a11, a22
    ! real(r8), allocatable :: a11(:)     ! the cell-cell submatrix
    ! real(r8), allocatable :: a12(:,:)   ! the cell-face submatrix
    type(pcsr_matrix)     :: a22        ! the face-face submatrix
    integer, allocatable  :: dir_faces(:)
  contains
    procedure, private :: init_disc
    procedure, private :: init_mold
    generic   :: init => init_disc, init_mold
    procedure :: compute
    procedure :: set_dir_faces
  end type mfd_2d_diff_matrix

contains

  subroutine init_disc (this, disc)

    class(mfd_2d_diff_matrix), intent(out) :: this
    type(mfd_2d_disc), intent(in), target :: disc

    integer :: j
    type(pcsr_graph), pointer :: g
    type(ip_desc), pointer :: row_ip

    this%disc => disc
    this%mesh => disc%mesh

    !TODO: include a11, a12
    ! allocate(this%a11(this%mesh%ncell))
    ! allocate(this%a12_val(size(this%mesh%cface)))

    !! Create a CSR matrix graph for the A22 submatrix.
    allocate(g)
    row_ip => this%mesh%face_ip
    call g%init (row_ip)
      do j = 1, this%mesh%ncell
        associate (cface => this%mesh%cface(this%mesh%cstart(j):this%mesh%cstart(j+1)-1))
          call g%add_clique (cface)
        end associate
      end do
    do j = 1, this%mesh%nlink
      call g%add_clique (this%mesh%lface(:,j))
    end do
    call g%add_complete

    !! Create the A22 submatrix.
    call this%a22%init (g, take_graph=.true.)

  end subroutine init_disc


  subroutine init_mold (this, mold)

    class(mfd_2d_diff_matrix), intent(out) :: this
    class(mfd_2d_diff_matrix), intent(in)  :: mold

    type(pcsr_graph), pointer :: g

    this%disc => mold%disc
    this%mesh => mold%mesh
    !TODO: include a11, a22
    ! allocate(this%a11(size(mold%a11)))
    ! allocate(this%a12_val(size(mold%a12_val)))
    g => mold%a22%graph_ptr()
    call this%a22%init (g, take_graph=.false.)

  end subroutine init_mold


  !! Compute the diffusion matrix given the cell-based array of diffusion
  !! coefficients.  Any existing dirichlet faces are dropped
  subroutine compute (this, coef)

    class(mfd_2d_diff_matrix), intent(inout) :: this
    real(r8), intent(in) :: coef(:)

    integer :: j, l, ir, ic
    real(r8), allocatable :: minv(:)

    ASSERT(size(coef) == this%mesh%ncell)

    call this%a22%set_all (0.0_r8)

    do j = 1, this%mesh%ncell
      !TODO: include a11, a22
      ! associate(a12 => this%a12_val(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
      !   if (coef(j) == 0.0_r8) then
      !     this%a11(j) = 0.0_r8
      !     a12 = 0.0_r8
      !     cycle
      !   end if
      !   minv = coef(j) * this%disc%minv(this%disc%xminv(j):this%disc%xminv(j+1)-1)
      !   !! Fill the A11 and A12 submatrices
      !   allocate(w(size(a12)))
      !   call upm_col_sum (minv, w)
      !   this%a11(j) = sum(w)
      !   a12 = -w
      !   deallocate(w)
      ! end associate
      !! Assemble the A22 CSR submatrix.
      minv = coef(j) * this%disc%minv(this%disc%xminv(j):this%disc%xminv(j+1)-1)
      associate(index => this%mesh%cface(this%mesh%cstart(j):this%mesh%cstart(j+1)-1))
        l = 1
        do ic = 1, size(index)
          do ir = 1, ic-1
            call this%a22%add_to (index(ir), index(ic), minv(l))
            call this%a22%add_to (index(ic), index(ir), minv(l))
            l = l + 1
          end do
          call this%a22%add_to (index(ic), index(ic), minv(l))
          l = l + 1
        end do
      end associate
    end do

    if (allocated(this%dir_faces)) deallocate(this%dir_faces)

  end subroutine compute


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DM_SET_DIR_FACES
 !!
 !! This fix-up subroutine sets selected faces to be Dirichlet faces.  The
 !! implication is that the unknowns associated with those faces are not
 !! actually unknowns and should not be included in the diffusion matrix.
 !! This is handled by projecting out the corresponding face rows and columns
 !! and setting unit diagonal values for those faces (effectively replacing
 !! equations for those unknows with dummy equations that are decoupled from
 !! all other unknowns).  The face-face submatrix is directly modified.  The
 !! modified cell-face submatrix is kept in factored form as the product of
 !! the original cell-face submatrix and the projection matrix described by
 !! the list of Dirichlet faces.  The array pointer DIR_FACES is a list of
 !! local indices of the Dirichlet faces.  The diffusion matrix maintains
 !! an internal reference to this array, and so it should not be modified.
 !!

  subroutine set_dir_faces (this, dir_faces)

    class(mfd_2d_diff_matrix), intent(inout) :: this
    integer, intent(in) :: dir_faces(:)

    integer :: j, n
    integer, allocatable :: tmp(:)

    ASSERT(minval(dir_faces) >= 1)
    ASSERT(maxval(dir_faces) <= this%mesh%nface)

    if (allocated(this%dir_faces)) then
      if (size(dir_faces) > 0) then
        n = size(this%dir_faces) + size(dir_faces)
        allocate(tmp(n))
        n = size(this%dir_faces)
        tmp(:n) = this%dir_faces
        tmp(n+1:) = dir_faces
        deallocate(this%dir_faces)
        call move_alloc (tmp, this%dir_faces)
      end if
    else
      allocate(this%dir_faces(size(dir_faces)))
      this%dir_faces = dir_faces
    end if

    do j = 1, size(dir_faces)
      n = dir_faces(j)
      call this%a22%project_out (n)
      call this%a22%set (n, n, 1.0_r8)
    end do

  end subroutine set_dir_faces

end module mfd_2d_diff_matrix_type
