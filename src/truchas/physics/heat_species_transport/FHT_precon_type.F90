!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module FHT_precon_type

  use kinds
  use FHT_model_type
  use unstr_mesh_type
  use mfd_diff_matrix_type
  use mfd_diff_precon_type
  use data_layout_type
  use rad_problem_type
  use truchas_timers
  implicit none
  private
  
  type, public :: FHT_precon
    type(FHT_model), pointer :: model => null()
    type(unstr_mesh), pointer :: mesh => null()
    integer, allocatable :: vfr_precon_coupling(:)
    type(mfd_diff_precon) :: precon
  end type FHT_precon
  !TODO! Consider having the MFD_DIFF_PRECON objects take ownership of the MFD_DIFF_MATRIX
  !TODO! objects so that the FHT_PRECON object doesn't need to hold them directly
  
  public :: FHT_precon_init
  public :: FHT_precon_compute, FHT_precon_apply
  
  !! Methods of coupling heat conduction and radiosity preconditioning.
  integer, parameter :: VFR_JAC = 1 ! Jacobi (radiosity and conduction decoupled)
  integer, parameter :: VFR_FGS = 2 ! Forward Gauss-Seidel (radiosity, then conduction)
  integer, parameter :: VFR_BGS = 3 ! Backward Gauss-Seidel (conduction, then radiosity)
  integer, parameter :: VFR_FAC = 4 ! Factorization (radiosity, conduction, radiosity)
  
contains

  subroutine FHT_precon_init (this, model, params)
  
    use parameter_list_type
    use truchas_logging_services

    type(FHT_precon), intent(out) :: this
    type(FHT_model), intent(in), target :: model
    type(parameter_list) :: params
    
    integer :: j, n, stat
    type(parameter_list), pointer :: plist
    character(:), allocatable :: string_array(:), errmsg
    type(mfd_diff_matrix), allocatable :: matrix

    this%model => model
    this%mesh => model%mesh
    
    !! Create the preconditioner for the heat equation.
    allocate(matrix)
    call matrix%init(model%disc)
    call this%precon%init(matrix, params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('FHT_PRECON_INIT: '//errmsg)
    
    !! Initialize the heat equation/view factor radiation
    if (associated(model%vf_rad_prob)) then
      n = size(model%vf_rad_prob)
      allocate(this%vfr_precon_coupling(n))
      call params%get('vfr-precon-coupling', string_array)
      INSIST(size(string_array) == n)
      do j = 1, n
        select case (string_array(j))
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
            if (this%model%vf_rad_prob(index)%fmask(j)) then
              n = this%model%vf_rad_prob(index)%faces(j)
              f2(n) = f2(n) + this%mesh%area(n) * z(j)
            end if
          end do
          deallocate(z)
        end if
      end do
      call stop_timer ('VF rad precon')
    end if

    !! Heat equation cell residual.
    call FHT_model_get_cell_temp_copy (this%model, f, f1x)
    call this%mesh%cell_imap%gather_offp(f1x)

    !! Heat equation face residual (with radiosity residuals optionally eliminated).
    call FHT_model_get_face_temp_copy (this%model, f, f2x)
    call this%mesh%face_imap%gather_offp(f2x)

    !! Precondition the heat equation.
    call this%precon%apply(f1x, f2x)
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
    type(mfd_diff_matrix), pointer :: dm
    integer, pointer :: faces(:) => null()
    integer, allocatable :: more_dir_faces(:)
    
    ASSERT(size(u) == FHT_model_size(this%model))
    ASSERT(h > 0.0_r8)
    
    call start_timer ('FHT precon compute')
    
    !! Initialize the STATE array.
    call FHT_model_get_cell_temp_copy (this%model, u, state(:,1))
    call this%mesh%cell_imap%gather_offp(state(:,1))
    
    call FHT_model_get_face_temp_copy (this%model, u, Tface)
    call this%mesh%face_imap%gather_offp(Tface)

    !! Finite difference approximant to dH/dT.
    fdinc = sqrt(epsilon(1.0d0)) !TODO! fix this naive choice of FD increment.
    !TODO: Hey there has been a compute_deriv -- why was it not used?
    call this%model%H_of_T%compute_value(state+fdinc, A)
    call this%model%H_of_T%compute_value(state, D)  ! D used as temp
    A = this%mesh%volume * ((A - D)/fdinc) / h
    
    call this%model%conductivity%compute_value(state, D)
    
    !! Correct data on void cells.
    if (associated(this%model%void_cell)) then
      where (this%model%void_cell)
        D = 0.0_r8
        A = 1.0_r8
      end where
    end if

    !! Jacobian of the basic heat equation that ignores nonlinearities
    !! in the conductivity.  This has the H/T relation eliminated.
    dm => this%precon%matrix_ref()
    call dm%compute (D)
    call dm%incr_cell_diag (A)
    
    !! Dirichlet boundary condition fixups.
    if (allocated(this%model%bc_dir)) then
      call this%model%bc_dir%compute(t)
      call dm%set_dir_faces(this%model%bc_dir%index)
    end if
    
    !! External HTC boundary condition contribution.
    if (allocated(this%model%bc_htc)) then
      call this%model%bc_htc%compute_deriv(t, Tface)
      call dm%incr_face_diag(this%model%bc_htc%index, this%model%bc_htc%deriv)
    end if

    !! Simple radiation boundary condition contribution.
    if (allocated(this%model%bc_rad)) then
      call this%model%bc_rad%compute_deriv(t, Tface)
      call dm%incr_face_diag(this%model%bc_rad%index, this%model%bc_rad%deriv)
    end if

    !! Internal HTC interface condition contribution.
    if (allocated(this%model%ic_htc)) then
      call this%model%ic_htc%compute_deriv(t, Tface)
      associate (index => this%model%ic_htc%index, &
                 deriv => this%model%ic_htc%deriv)
        if (associated(this%model%void_face)) then
          do j = 1, size(index,2) !FIXME? possibly bad form to modify %deriv
            if (any(this%model%void_face(index(:,j)))) deriv(:,j) = 0.0_r8
          end do
        end if
        call dm%incr_interface_flux3(index, deriv)
      end associate
    end if

    !! Internal gap radiation condition contribution.
    if (allocated(this%model%ic_rad)) then
      call this%model%ic_rad%compute_deriv(t, Tface)
      associate (index => this%model%ic_rad%index, &
                 deriv => this%model%ic_rad%deriv)
        if (associated(this%model%void_face)) then
          do j = 1, size(index,2) !FIXME? possibly bad form to modify %deriv
            if (any(this%model%void_face(index(:,j)))) deriv(:,j) = 0.0_r8
          end do
        end if
        call dm%incr_interface_flux3(index, deriv)
      end associate
    end if

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
        where (.not.this%model%vf_rad_prob(index)%fmask) values = 0
        call dm%incr_face_diag (faces, this%mesh%area(faces) * values)
        deallocate(values)
      end do
    end if

    !! The matrix is now complete; re-compute the preconditioner.
    call this%precon%compute
      
    call stop_timer ('FHT precon compute')
    
  end subroutine FHT_precon_compute

end module FHT_precon_type
