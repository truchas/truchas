#include "f90_assert.fpp"

module diffusion_matrix

  use kinds, only: r8
  use dist_mesh_type
  use mfd_disc_type
  use pcsr_matrix_type
  use parallel_communication
  use index_partitioning
  implicit none
  private

  type, public :: dist_diff_matrix
    type(mfd_disc),  pointer :: disc => null()  ! reference only -- do not own
    type(dist_mesh), pointer :: mesh => null()  ! reference only -- do not own
    real(r8), allocatable :: a11(:)     ! the cell-cell submatrix
    real(r8), allocatable :: a12(:,:)   ! the cell-face submatrix
    type(pcsr_matrix)     :: a22        ! the face-face submatrix
    integer, allocatable  :: dir_faces(:)
  contains
    procedure, private :: init_disc
    procedure, private :: init_mold
    generic   :: init => init_disc, init_mold
    procedure :: compute
    procedure :: set_dir_faces
    procedure :: incr_cell_diag
    procedure :: incr_face_diag
    procedure :: incr_interface_flux
    procedure :: incr_interface_flux2
    procedure :: compute_face_schur_matrix
  end type dist_diff_matrix
  
contains

  subroutine init_disc (this, disc)

    class(dist_diff_matrix), intent(out) :: this
    type(mfd_disc), intent(in), target :: disc

    integer :: j
    type(pcsr_graph), pointer :: g
    type(ip_desc), pointer :: row_ip

    this%disc => disc
    this%mesh => disc%mesh
    allocate(this%a11(this%mesh%ncell))
    allocate(this%a12(size(this%mesh%cface,1),this%mesh%ncell))

    !! Create a CSR matrix graph for the A22 submatrix.
    allocate(g)
    row_ip => this%mesh%face_ip
    call g%init (row_ip)
    do j = 1, this%mesh%ncell
      call g%add_clique (this%mesh%cface(:,j))
    end do
    do j = 1, this%mesh%nlink
      call g%add_clique (this%mesh%lface(:,j))
    end do
    call g%add_complete

    !! Create the A22 submatrix.
    call this%a22%init (g, take_graph=.true.)

  end subroutine init_disc
  
  subroutine init_mold (this, mold)
  
    class(dist_diff_matrix), intent(out) :: this
    class(dist_diff_matrix), intent(in)  :: mold
    
    type(pcsr_graph), pointer :: g
    
    this%disc => mold%disc
    this%mesh => mold%mesh
    allocate(this%a11(size(mold%a11)))
    allocate(this%a12(size(mold%a12,1),size(mold%a12,2)))
    g => mold%a22%graph_ptr()
    call this%a22%init (g, take_graph=.false.)
    
  end subroutine init_mold

  subroutine compute (this, d)

    use upper_packed_matrix, only: upm_col_sum
  
    class(dist_diff_matrix), intent(inout) :: this
    real(r8), intent(in) :: d(:)
    
    integer :: j, l, ir, ic
    integer, pointer :: index(:)
    real(r8) :: w(size(this%mesh%cface,dim=1))
    real(r8) :: minv(size(w)*(size(w)+1)/2)
    !type(mfd_hex) :: tmp

    ASSERT(size(d) == this%mesh%ncell)
    
    call this%a22%set_all (0.0_r8)
    
    do j = 1, this%mesh%ncell
      if (d(j) == 0.0_r8) then
        this%a11(j) = 0.0_r8
        this%a12(:,j) = 0.0_r8
        cycle
      end if
      !call init_mfd_hex (tmp, this%mesh%x(:,this%mesh%cnode(:,j)))
      !call compute_flux_matrix (tmp, d(j), minv, invert=.true.)
      minv = d(j) * this%disc%minv(:,j)
      !! Fill the A11 and A12 submatrices
      call upm_col_sum (minv, w)
      this%a11(j) = sum(w)
      this%a12(:,j) = -w
      !! Assemble the A22 CSR submatrix.
      index => this%mesh%cface(:,j)
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

    class(dist_diff_matrix), intent(inout) :: this
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
  
  !! This subroutine increments the (entire) diagonal cell-cell diffusion
  !! submatrix with the specified values.  The intended use is to modify
  !! the base diffusion matrix to account for the discretization of a time
  !! derivative term.
  
  subroutine incr_cell_diag (this, values)
    class(dist_diff_matrix), intent(inout) :: this
    real(r8), intent(in) :: values(:)
    ASSERT(size(values) == size(this%a11))
    this%a11 = this%a11 + values
  end subroutine incr_cell_diag

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DM_INCR_FACE_DIAG
 !!
 !! This subroutine increments selected diagonal elements of the face-face
 !! diffusion submatrix.  The array INDICES is a list of local face indices
 !! and VALUES the corresponding array of values to add to the diagonal.
 !! The intended use is to modify the base diffusion matrix to account for
 !! the discretization of certain types of flux boundary conditions.
 !!

  subroutine incr_face_diag (this, indices, values)

    class(dist_diff_matrix), intent(inout) :: this
    integer,  intent(in) :: indices(:)
    real(r8), intent(in) :: values(:)

    integer :: j

    ASSERT(size(indices) == size(values))
    ASSERT(minval(indices) >= 1)
    ASSERT(maxval(indices) <= this%mesh%nface)

    do j = 1, size(indices)
      call this%a22%add_to (indices(j), indices(j), values(j))
    end do

  end subroutine incr_face_diag
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DM_INCR_INTERFACE_FLUX
 !!
 !! This subroutine increments certain elements of the face-face diffusion
 !! submatrix in a specific manner to account for the contribution from the
 !! internal interface flux conditions.  For a matching pair of faces, I and J,
 !! spanning the iterface, the linearization of the condition adds to the
 !! {I,J} submatrix the 2x2 matrix [A, -A; -A, A], for some scalar A.
 !!

  subroutine incr_interface_flux (this, faces, values)
  
    class(dist_diff_matrix), intent(inout) :: this
    integer,  intent(in) :: faces(:,:)
    real(r8), intent(in) :: values(:)
    
    integer :: j, n1, n2
    
    ASSERT(size(faces,1) == 2)
    ASSERT(size(faces,2) == size(values))
    ASSERT(minval(faces) >= 1)
    ASSERT(maxval(faces) <= this%mesh%nface)
  
    do j = 1, size(faces,2)
      n1 = faces(1,j)
      n2 = faces(2,j)
      call this%a22%add_to (n1, n1,  values(j))
      call this%a22%add_to (n1, n2, -values(j))
      call this%a22%add_to (n2, n1, -values(j))
      call this%a22%add_to (n2, n2,  values(j))
    end do
  
  end subroutine incr_interface_flux
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DM_INCR_INTERFACE_FLUX2
 !!
 !! This subroutine increments certain elements of the face-face diffusion
 !! submatrix in a specific manner to account for the contribution from the
 !! internal gap radiation conditions.  For a matching pair of faces, I and J,
 !! spanning the iterface, the linearization of the condition adds to the
 !! {I,J} submatrix the 2x2 matrix [A, -B; -A, B], for some scalars A and B.
 !!

  subroutine incr_interface_flux2 (this, faces, values)
  
    class(dist_diff_matrix), intent(inout) :: this
    integer,  intent(in) :: faces(:,:)
    real(r8), intent(in) :: values(:,:)
    
    integer :: j, n1, n2
    
    ASSERT(size(faces,1) == 2)
    ASSERT(size(values,1) == 2)
    ASSERT(size(faces,2) == size(values,2))
    ASSERT(minval(faces) >= 1)
    ASSERT(maxval(faces) <= this%mesh%nface)
  
    do j = 1, size(faces,2)
      n1 = faces(1,j)
      n2 = faces(2,j)
      call this%a22%add_to (n1, n1,  values(1,j))
      call this%a22%add_to (n1, n2, -values(2,j))
      call this%a22%add_to (n2, n1, -values(1,j))
      call this%a22%add_to (n2, n2,  values(2,j))
    end do
  
  end subroutine incr_interface_flux2
  
  subroutine compute_face_schur_matrix (this, Sff)
  
    class(dist_diff_matrix), intent(in) :: this
    type(pcsr_matrix), intent(inout) :: Sff
    
    integer :: j, n, ir, ic
    integer, pointer :: indices(:)
    real(r8) :: value
    
    ASSERT(associated(this%a22%graph, Sff%graph))
    
    Sff%values = this%a22%values
    do j = 1, this%mesh%ncell
      indices => this%mesh%cface(:,j)
      do ir = 1, size(indices)
        do ic = 1, size(indices)
          value = -this%a12(ir,j)*this%a12(ic,j)/this%a11(j)
          call Sff%add_to (indices(ir), indices(ic), value)
        end do
      end do
    end do
    
    !! Apply the Dirichlet projections.
    do j = 1, size(this%dir_faces)
      n = this%dir_faces(j)
      call Sff%project_out (n)
      call Sff%set (n, n, 1.0_r8)
    end do

  end subroutine compute_face_schur_matrix

end module diffusion_matrix
