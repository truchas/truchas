#include "f90_assert.fpp"

module FHT_precon_type

  use kinds
  use FHT_model_type
  use base_mesh_class
  use index_partitioning
  use diffusion_matrix
  use diff_precon_type
  use data_layout_type
  use rad_problem_type
  use boundary_data
  use interface_data
  use property_mesh_function
  use timing_tree
  implicit none
  private
  
  type, public :: FHT_precon
    type(FHT_model), pointer :: model => null()
    class(base_mesh), pointer :: mesh => null()
    integer, pointer :: vfr_precon_coupling(:) => null()
    type(dist_diff_matrix), pointer :: matrix => null()
    type(diff_precon) :: precon
  end type FHT_precon
  !TODO! Consider having the diff_precon objects take ownership of the dist_diff_matrix
  !TODO! objects so that the FHT_precon object doesn't need to hold them directly
  
  public :: FHT_precon_init, FHT_precon_delete
  public :: FHT_precon_compute, FHT_precon_apply
  
  type, public :: FHT_precon_params
    type(diff_precon_params) :: HC_precon_params
    character(len=16), pointer :: vfr_precon_coupling(:) => null()
  end type
  public :: diff_precon_params, ssor_precon_params, boomer_amg_precon_params
  
  !! Methods of coupling heat conduction and radiosity preconditioning.
  integer, parameter :: VFR_JAC = 1 ! Jacobi (radiosity and conduction decoupled)
  integer, parameter :: VFR_FGS = 2 ! Forward Gauss-Seidel (radiosity, then conduction)
  integer, parameter :: VFR_BGS = 3 ! Backward Gauss-Seidel (conduction, then radiosity)
  integer, parameter :: VFR_FAC = 4 ! Factorization (radiosity, conduction, radiosity)
  
contains

  subroutine FHT_precon_init (this, model, params)
  
    type(FHT_precon), intent(out) :: this
    type(FHT_model), intent(in), target :: model
    type(FHT_precon_params), intent(in) :: params
    
    integer :: j, n

    this%model => model
    this%mesh => model%mesh
    
    !! Create the preconditioner for the heat equation.
    allocate(this%matrix)
    call this%matrix%init (model%disc)
    call diff_precon_init (this%precon, this%matrix, params%HC_precon_params)
    
    !! Initialize the heat equation/view factor radiation
    if (associated(model%vf_rad_prob)) then
      n = size(model%vf_rad_prob)
      INSIST(associated(params%vfr_precon_coupling))
      INSIST(size(params%vfr_precon_coupling) == n)
      allocate(this%vfr_precon_coupling(n))
      do j = 1, n
        select case (params%vfr_precon_coupling(j))
        case ('JACOBI')
          this%vfr_precon_coupling(j) = VFR_JAC
        case ('FORWARD GS')
          this%vfr_precon_coupling(j) = VFR_FGS
        case ('BACKWARD GS')
          this%vfr_precon_coupling(j) = VFR_BGS
        case ('FACTORIZATION')
          this%vfr_precon_coupling(j) = VFR_FAC
        case default
          INSIST(.false.)
        end select
      end do
    end if

  end subroutine FHT_precon_init
  
  subroutine FHT_precon_delete (this)
    type(FHT_precon), intent(inout) :: this
    if (associated(this%matrix)) deallocate(this%matrix)
    call diff_precon_delete (this%precon)
    if (associated(this%vfr_precon_coupling)) deallocate(this%vfr_precon_coupling)
  end subroutine FHT_precon_delete

  subroutine FHT_precon_apply (this, t, u, f)
  
    type(FHT_precon), intent(inout) :: this
    real(r8), intent(in) :: t, u(:)
    real(r8), intent(inout) :: f(:)
    target :: u, f
  
    integer :: index, j, n
    real(r8), pointer :: f2(:), fq(:), Tface(:)
    real(r8) :: f1x(this%mesh%ncell), f2x(this%mesh%nface)
    real(r8), allocatable :: z(:)
    
    ASSERT(size(u) == FHT_model_size(this%model))
    ASSERT(size(f) == FHT_model_size(this%model))
    
    !! Glossary:
    !! Pointers into segments of the F array argument:
    !!  f1  -- cell temperature segment
    !!  f2  -- face temperature segment
    !!  fq  -- a radiosity segment
    !! Local arrays:
    !!  f1x -- off-process extended copy of f1
    !!  f2x -- off-process extended copy of f2
    
    call start_timer ('FHT precon apply')

    !! Precondition the radiosity components:
    !! Factorization and forward Gauss-Seidel coupling methods.
    if (associated(this%model%vf_rad_prob)) then
      call start_timer ('VF rad precon')
      call FHT_model_get_face_temp_view (this%model, f, f2)
      do index = 1, size(this%model%vf_rad_prob)
        if (this%vfr_precon_coupling(index) == VFR_FGS .or. &
            this%vfr_precon_coupling(index) == VFR_FAC) then
          call FHT_model_get_radiosity_view (this%model, index, f, fq)
          allocate(z(size(fq)))
          z = fq
          call this%model%vf_rad_prob(index)%precon (t, z)
          if (this%vfr_precon_coupling(index) == VFR_FGS) fq = z
          !! Update the heat equation face residual.
          call this%model%vf_rad_prob(index)%precon_matvec1 (t, z)
          do j = 1, size(z)
            n = this%model%vf_rad_prob(index)%faces(j)
            f2(n) = f2(n) + this%mesh%area(n) * z(j)
          end do
          deallocate(z)
        end if
      end do
      call stop_timer ('VF rad precon')
    end if

    !! Heat equation cell residual.
    call FHT_model_get_cell_temp_copy (this%model, f, f1x)
    call gather_boundary (this%mesh%cell_ip, f1x)

    !! Heat equation face residual (with radiosity residuals optionally eliminated).
    call FHT_model_get_face_temp_copy (this%model, f, f2x)
    call gather_boundary (this%mesh%face_ip, f2x)

    !! Precondition the heat equation.
    call diff_precon_apply (this%precon, f1x, f2x)
    call FHT_model_set_cell_temp (this%model, f1x, f)
    call FHT_model_set_face_temp (this%model, f2x, f)

    !! Precondition the radiosity components:
    !! Factorization, backward Gauss-Seidel and Jacobi coupling methods.
    if (associated(this%model%vf_rad_prob)) then
      call start_timer ('VF rad precon')
      call FHT_model_get_face_temp_view (this%model, f, f2)
      call FHT_model_get_face_temp_view (this%model, u, Tface)
      do index = 1, size(this%model%vf_rad_prob)
        if (this%vfr_precon_coupling(index) == VFR_JAC .or. &
            this%vfr_precon_coupling(index) == VFR_BGS .or. &
            this%vfr_precon_coupling(index) == VFR_FAC) then
          !! Update the radiosity residual components.
          call FHT_model_get_radiosity_view (this%model, index, f, fq)
          if (this%vfr_precon_coupling(index) /= VFR_JAC) then
            allocate(z(size(fq)))
            call this%model%vf_rad_prob(index)%rhs_deriv (t, Tface(this%model%vf_rad_prob(index)%faces), z)
            fq = fq + z * f2(this%model%vf_rad_prob(index)%faces)
            deallocate(z)
          end if
          call this%model%vf_rad_prob(index)%precon (t, fq)
        end if
      end do
      call stop_timer ('VF rad precon')
    end if

    call stop_timer ('FHT precon apply')
   
  end subroutine FHT_precon_apply
  
  subroutine FHT_precon_compute (this, t, u, h)
  
    type(FHT_precon), intent(inout) :: this
    real(r8), intent(in) :: t, u(:), h
    target :: u

    integer :: index, n, j, n1, n2
    real(r8) :: fdinc, term
    real(r8) :: state(this%mesh%ncell,1), D(this%mesh%ncell), A(this%mesh%ncell), Tface(this%mesh%nface)
    real(r8), allocatable :: values(:), values2(:,:)
    type(dist_diff_matrix), pointer :: dm
    integer, pointer :: faces(:) => null()
    integer, allocatable :: more_dir_faces(:)
    
    ASSERT(size(u) == FHT_model_size(this%model))
    ASSERT(h > 0.0_r8)
    
    call start_timer ('FHT precon compute')
    
    !! Initialize the STATE array.
    call FHT_model_get_cell_temp_copy (this%model, u, state(:,1))
    call gather_boundary (this%mesh%cell_ip, state(:,1))
    
    call FHT_model_get_face_temp_copy (this%model, u, Tface)
    call gather_boundary (this%mesh%face_ip, Tface)

    !! Finite difference approximant to dH/dT.
    fdinc = sqrt(epsilon(1.0d0)) !TODO! fix this naive choice of FD increment.
    call pmf_eval (this%model%H_of_T, state+fdinc, A)
    call pmf_eval (this%model%H_of_T, state, D)  ! D used as temp
    A = this%mesh%volume * ((A - D)/fdinc) / h
    
    call pmf_eval (this%model%conductivity, state, D)
    
    !! Correct data on void cells.
    if (associated(this%model%void_cell)) then
      where (this%model%void_cell)
        D = 0.0_r8
        A = 1.0_r8
      end where
    end if

    !! Jacobian of the basic heat equation that ignores nonlinearities
    !! in the conductivity.  This has the H/T relation eliminated.
    dm => diff_precon_matrix(this%precon)
    call dm%compute (D)
    call dm%incr_cell_diag (A)
    
    !! Dirichlet boundary condition fixups.
    call bd_data_eval (this%model%bc_dir, t)
    call dm%set_dir_faces (this%model%bc_dir%faces)
    
    !! External HTC boundary condition contribution.
    call bd_data_eval (this%model%bc_htc, t)
    allocate(values(size(this%model%bc_htc%faces)))
    values = this%mesh%area(this%model%bc_htc%faces) * this%model%bc_htc%values(1,:)
    call dm%incr_face_diag (this%model%bc_htc%faces, values)
    deallocate(values)

    !! Simple radiation boundary condition contribution.
    call bd_data_eval (this%model%bc_rad, t)
    faces => this%model%bc_rad%faces
    allocate(values(size(faces)))
    values = 4 * this%model%sbconst * this%mesh%area(faces) * this%model%bc_rad%values(1,:) * &
        (Tface(faces) - this%model%abszero)**3
    call dm%incr_face_diag (faces, values)
    deallocate(values)

    !! Internal HTC interface condition contribution.
    call if_data_eval (this%model%ic_htc, t)
    allocate(values(size(this%model%ic_htc%faces,dim=2)))
    values = this%mesh%area(this%model%ic_htc%faces(1,:)) * this%model%ic_htc%values(1,:)
    if (associated(this%model%void_face)) then
      do j = 1, size(values)
        if (any(this%model%void_face(this%model%ic_htc%faces(:,j)))) values(j) = 0.0_r8
      end do
    end if
    call dm%incr_interface_flux (this%model%ic_htc%faces, values)
    deallocate(values)

    !! Internal gap radiation condition contribution.
    call if_data_eval (this%model%ic_rad, t)
    allocate(values2(2,size(this%model%ic_rad%faces,dim=2)))
    do j = 1, size(values2,2)
      n1 = this%model%ic_rad%faces(1,j)
      n2 = this%model%ic_rad%faces(2,j)
      term = 4 * this%model%sbconst * this%mesh%area(n1) * this%model%ic_rad%values(1,j)
      values2(1,j) = term * (Tface(n1) - this%model%abszero)**3
      values2(2,j) = term * (Tface(n2) - this%model%abszero)**3
    end do
    if (associated(this%model%void_face)) then
      do j = 1, size(values2,2)
        if (any(this%model%void_face(this%model%ic_rad%faces(:,j)))) values2(:,j) = 0.0_r8
      end do
    end if
    call dm%incr_interface_flux2 (this%model%ic_rad%faces, values2)
    deallocate(values2)

    if (associated(this%model%void_face)) then
      n = count(this%model%void_face)
      if (n > 0) then
        allocate(more_dir_faces(n))
        n = 0
        do j = 1, this%mesh%nface
          if (this%model%void_face(j)) then
            n = n + 1
            more_dir_faces(n) = j
          end if
        end do
        call dm%set_dir_faces (more_dir_faces)
        deallocate(more_dir_faces)
      end if
    end if
                                       
    !! Enclosure radiation contributions to the preconditioner.
    !! TODO: what about factorization coupling?  Is this still correct?
    if (associated(this%model%vf_rad_prob)) then
      do index = 1, size(this%model%vf_rad_prob)
        faces => this%model%vf_rad_prob(index)%faces
        allocate(values(size(faces)))
        call this%model%vf_rad_prob(index)%rhs_deriv (t, Tface(faces), values)
        call dm%incr_face_diag (faces, this%mesh%area(faces) * values)
        deallocate(values)
      end do
    end if

    !! The matrix is now complete; re-compute the preconditioner.
    call diff_precon_compute (this%precon)
      
    call stop_timer ('FHT precon compute')
    
  end subroutine FHT_precon_compute

end module FHT_precon_type
