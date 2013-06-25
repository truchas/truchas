#include "f90_assert.fpp"

module diffusion_matrix

  use kinds
  use parallel_csr_matrix
  use distributed_mesh, only: dist_mesh
  use mfd_disc_type
  use parallel_communication
  use index_partitioning
  implicit none
  private

  public :: DM_create, DM_destroy, DM_compute
  public :: DM_set_dir_faces, DM_incr_cell_diag, DM_incr_face_diag
  public :: DM_incr_interface_flux, DM_incr_interface_flux2
  public :: DM_compute_face_schur_matrix
  public :: DM_get_Dff

  type, public :: dist_diff_matrix
    type(mfd_disc),  pointer :: disc => null()
    type(dist_mesh), pointer :: mesh => null()
    real(r8), pointer :: a11(:) => null(), a12(:,:) => null()
    type(pcsr_matrix), pointer :: a22 => null()
    integer, pointer :: dir_faces(:) => null()
  end type dist_diff_matrix
  
  interface DM_create
    module procedure DM_create, DM_create_mold
  end interface

contains

  subroutine DM_create (this, disc)

    type(dist_diff_matrix), intent(out) :: this
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
    call pcsr_graph_create (g, row_ip)
    do j = 1, this%mesh%ncell
      call pcsr_graph_add_clique (g, this%mesh%cface(:,j))
    end do
    do j = 1, this%mesh%nlink
      call pcsr_graph_add_clique (g, this%mesh%lface(:,j))
    end do
    call pcsr_graph_fill_complete (g)

    !! Create the A22 submatrix.
    allocate(this%a22)
    call pcsr_matrix_create (this%a22, g, take_graph=.true.)

  end subroutine DM_create
  
  subroutine DM_create_mold (this, mold)
  
    type(dist_diff_matrix), intent(out) :: this
    type(dist_diff_matrix), intent(in)  :: mold
    
    type(pcsr_graph), pointer :: g
    
    this%disc => mold%disc
    this%mesh => mold%mesh
    allocate(this%a11(size(mold%a11)))
    allocate(this%a12(size(mold%a12,1),size(mold%a12,2)))
    allocate(this%a22)
    g => pcsr_matrix_graph(mold%a22)
    call pcsr_matrix_create(this%a22, g, take_graph=.false.)
    
  end subroutine DM_create_mold
  
  subroutine DM_destroy (this)

    type(dist_diff_matrix), intent(inout) :: this

    type(dist_diff_matrix) :: default

    if (associated(this%a11)) deallocate(this%a11)
    if (associated(this%a12)) deallocate(this%a12)
    if (associated(this%a22)) then
      call pcsr_matrix_delete (this%a22)
      deallocate(this%a22)
    end if

    this = default  ! assign default initialization values

  end subroutine DM_destroy
  
  subroutine DM_get_Dff (this, Dff)
    type(dist_diff_matrix), intent(in) :: this
    type(pcsr_matrix), pointer :: Dff
    Dff => this%a22
  end subroutine DM_get_Dff

  subroutine DM_compute (this, d)

    use upper_packed_matrix, only: upm_col_sum
  
    type(dist_diff_matrix), intent(inout) :: this
    real(r8), intent(in) :: d(:)
    
    integer :: j, l, ir, ic
    integer, pointer :: index(:)
    real(r8) :: w(size(this%mesh%cface,dim=1))
    real(r8) :: minv(size(w)*(size(w)+1)/2)
    !type(mfd_hex) :: tmp

    ASSERT(size(d) == this%mesh%ncell)
    
    call pcsr_matrix_set_all_values (this%a22, 0.0_r8)
    
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
          call pcsr_matrix_add_to_value (this%a22, index(ir), index(ic), minv(l))
          call pcsr_matrix_add_to_value (this%a22, index(ic), index(ir), minv(l))
          l = l + 1
        end do
        call pcsr_matrix_add_to_value (this%a22, index(ic), index(ic), minv(l))
        l = l + 1
      end do
    end do
    
    if (associated(this%dir_faces)) deallocate(this%dir_faces)

  end subroutine DM_compute

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

  subroutine DM_set_dir_faces (this, dir_faces)

    type(dist_diff_matrix), intent(inout) :: this
    integer, intent(in) :: dir_faces(:)

    integer :: j, n
    integer, pointer :: tmp(:)
    
    ASSERT(minval(dir_faces) >= 1)
    ASSERT(maxval(dir_faces) <= this%mesh%nface)

    if (associated(this%dir_faces)) then
      if (size(dir_faces) > 0) then
        n = size(this%dir_faces) + size(dir_faces)
        allocate(tmp(n))
        n = size(this%dir_faces)
        tmp(:n) = this%dir_faces
        tmp(n+1:) = dir_faces
        deallocate(this%dir_faces)
        this%dir_faces => tmp
      end if
    else
      allocate(this%dir_faces(size(dir_faces)))
      this%dir_faces = dir_faces
    end if

    do j = 1, size(dir_faces)
      n = dir_faces(j)
      call pcsr_matrix_project_out_index (this%a22, n)
      call pcsr_matrix_set_value (this%a22, n, n, 1.0_r8)
    end do

  end subroutine DM_set_dir_faces
  
  !! This subroutine increments the (entire) diagonal cell-cell diffusion
  !! submatrix with the specified values.  The intended use is to modify
  !! the base diffusion matrix to account for the discretization of a time
  !! derivative term.
  
  subroutine DM_incr_cell_diag (this, values)
    type(dist_diff_matrix), intent(inout) :: this
    real(r8), intent(in) :: values(:)
    ASSERT(size(values) == size(this%a11))
    this%a11 = this%a11 + values
  end subroutine DM_incr_cell_diag

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

  subroutine DM_incr_face_diag (this, indices, values)

    type(dist_diff_matrix), intent(inout) :: this
    integer,  intent(in) :: indices(:)
    real(r8), intent(in) :: values(:)

    integer :: j

    ASSERT(size(indices) == size(values))
    ASSERT(minval(indices) >= 1)
    ASSERT(maxval(indices) <= this%mesh%nface)

    do j = 1, size(indices)
      call pcsr_matrix_add_to_value (this%a22, indices(j), indices(j), values(j))
    end do

  end subroutine DM_incr_face_diag
  
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

  subroutine DM_incr_interface_flux (this, faces, values)
  
    type(dist_diff_matrix), intent(inout) :: this
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
      call pcsr_matrix_add_to_value (this%a22, n1, n1,  values(j))
      call pcsr_matrix_add_to_value (this%a22, n1, n2, -values(j))
      call pcsr_matrix_add_to_value (this%a22, n2, n1, -values(j))
      call pcsr_matrix_add_to_value (this%a22, n2, n2,  values(j))
    end do
  
  end subroutine DM_incr_interface_flux
  
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

  subroutine DM_incr_interface_flux2 (this, faces, values)
  
    type(dist_diff_matrix), intent(inout) :: this
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
      call pcsr_matrix_add_to_value (this%a22, n1, n1,  values(1,j))
      call pcsr_matrix_add_to_value (this%a22, n1, n2, -values(2,j))
      call pcsr_matrix_add_to_value (this%a22, n2, n1, -values(1,j))
      call pcsr_matrix_add_to_value (this%a22, n2, n2,  values(2,j))
    end do
  
  end subroutine DM_incr_interface_flux2
  
  subroutine DM_compute_face_schur_matrix (this, Sff)
  
    type(dist_diff_matrix), intent(in) :: this
    type(pcsr_matrix), intent(inout) :: Sff
    
    integer :: j, n, ir, ic
    integer, pointer :: indices(:)
    real(r8) :: value
    
    ASSERT(associated(this%a22%graph, Sff%graph))
    
    Sff%data = this%a22%data
    do j = 1, this%mesh%ncell
      indices => this%mesh%cface(:,j)
      do ir = 1, size(indices)
        do ic = 1, size(indices)
          value = -this%a12(ir,j)*this%a12(ic,j)/this%a11(j)
          call pcsr_matrix_add_to_value (Sff, indices(ir), indices(ic), value)
        end do
      end do
    end do
    
    !! Apply the Dirichlet projections.
    do j = 1, size(this%dir_faces)
      n = this%dir_faces(j)
      call pcsr_matrix_project_out_index (Sff, n)
      call pcsr_matrix_set_value (Sff, n, n, 1.0_r8)
    end do

  end subroutine DM_compute_face_schur_matrix

end module diffusion_matrix
