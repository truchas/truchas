!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module ht_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ht_vector_type
  use ht_model_type
  use rad_problem_type
  use unstr_mesh_type
  use mfd_diff_precon_type
  use mfd_diff_matrix_type
  use truchas_timers
  implicit none
  private

  type, public :: ht_precon
    type(ht_model), pointer :: model => null()
    type(unstr_mesh), pointer :: mesh  => null()
    real(r8) :: dt
    real(r8), allocatable :: dHdT(:)
    integer, allocatable :: vfr_precon_coupling(:)
    type(mfd_diff_precon) :: pc
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type ht_precon

  !! Methods of coupling heat conduction and radiosity preconditioning.
  integer, parameter :: VFR_JAC = 1 ! Jacobi (radiosity and conduction decoupled)
  integer, parameter :: VFR_FGS = 2 ! Forward Gauss-Seidel (radiosity, then conduction)
  integer, parameter :: VFR_BGS = 3 ! Backward Gauss-Seidel (conduction, then radiosity)
  integer, parameter :: VFR_FAC = 4 ! Factorization (radiosity, conduction, radiosity)

contains

  subroutine init(this, model, params, stat, errmsg)

    use parameter_list_type

    class(ht_precon), intent(out) :: this
    type(ht_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(mfd_diff_matrix), allocatable :: matrix

    this%model => model
    this%mesh  => model%mesh

    allocate(this%dHdT(this%mesh%ncell))
    allocate(matrix)
    call matrix%init(model%disc)
    call this%pc%init(matrix, params, stat, errmsg)

    if (associated(model%vf_rad_prob)) then
      block
        integer :: n, j
        character(:), allocatable :: array(:)
        n = size(model%vf_rad_prob)
        allocate(this%vfr_precon_coupling(n))
        call params%get('vfr-precon-coupling', array)
        INSIST(size(array) == n)
        do j = 1, n
          select case (array(j))
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
      end block
    end if

  end subroutine init


  subroutine compute(this, t, u, dt)

    class(ht_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    type(ht_vector), intent(inout) :: u
    target :: u

    real(r8), pointer :: state(:,:) => null()
    integer, allocatable :: more_dir_faces(:)

    integer :: index, j, n, n1, n2
    real(r8) :: D(this%mesh%ncell), A(this%mesh%ncell), term
    real(r8), allocatable :: values(:), values2(:,:)
    type(mfd_diff_matrix), pointer :: dm

    ASSERT(dt > 0.0_r8)

    call start_timer('ht-precon-compute')

    !! Generate list of void faces; these are treated like Dirichlet BC.
    if (associated(this%model%void_face)) then
      n = count(this%model%void_face)
      allocate(more_dir_faces(n))
      n = 0
      do j = 1, this%mesh%nface
        if (this%model%void_face(j)) then
          n = n + 1
          more_dir_faces(n) = j
        end if
      end do
    else
      allocate(more_dir_faces(0))
    end if

    !TODO: The existing prop_mesh_func%compute_value function expects a rank-2
    !      state array. This is a workaround until prop_mesh_func is redesigned.
    state(1:this%mesh%ncell,1:1) => u%tc

    call this%mesh%cell_imap%gather_offp(u%tc)
    call this%mesh%face_imap%gather_offp(u%tf)

    this%dt = dt

    !! Hardwired assumption that T is the first component of state -- FIXME!
    call this%model%H_of_T%compute_deriv(state, 1, this%dHdT)
    call this%model%conductivity%compute_value(state, D)
    A = this%mesh%volume * this%dHdT / dt

    !! Correct data on void cells.
    if (associated(this%model%void_cell)) then
      where (this%model%void_cell)
        this%dHdT = 0.0_r8
        D = 0.0_r8
        A = 1.0_r8
      end where
    end if

    !! Jacobian of the basic heat equation that ignores nonlinearities
    !! in the conductivity.  This has the H/T relation eliminated.
    dm => this%pc%matrix_ref()
    call dm%compute(D)
    call dm%incr_cell_diag(A)

    !! Dirichlet boundary condition fixups.
    if (allocated(this%model%bc_dir)) then
      call this%model%bc_dir%compute(t)
      call dm%set_dir_faces(this%model%bc_dir%index)
    end if

    !! External HTC boundary condition contribution.
    if (allocated(this%model%bc_htc)) then
      call this%model%bc_htc%compute_deriv(t, u%tf)
      associate (index => this%model%bc_htc%index, &
                 deriv => this%model%bc_htc%deriv)
        call dm%incr_face_diag(index, deriv)
      end associate
    end if

    !! Simple radiation boundary condition contribution.
    if (allocated(this%model%bc_rad)) then
      call this%model%bc_rad%compute(t, u%tf)
      associate (index => this%model%bc_rad%index, &
                 deriv => this%model%bc_rad%deriv)
        call dm%incr_face_diag(index, deriv)
      end associate
    end if

    !! Experimental evaporation heat flux contribution.
    if (allocated(this%model%evap_flux)) then
      call this%model%evap_flux%compute_deriv(t, u%tf)
      associate (index => this%model%evap_flux%index, &
                 deriv => this%model%evap_flux%deriv)
        call dm%incr_face_diag(index, this%mesh%area(index)*deriv)
      end associate
    endif

    !! Internal HTC interface condition contribution.
    if (allocated(this%model%ic_htc)) then
      call this%model%ic_htc%compute(t, u%tf)
      associate (index => this%model%ic_htc%index, &
                 deriv => this%model%ic_htc%deriv)
        if (associated(this%model%void_face)) then
          do j = 1, size(index,2) !FIXME? Bad form to modify deriv?
            if (any(this%model%void_face(index(:,j)))) deriv(:,j) = 0.0_r8
          end do
        end if
        call dm%incr_interface_flux3(index, deriv) !TODO: rename these methods
      end associate
    end if

    !! Internal gap radiation condition contribution.
    if (allocated(this%model%ic_rad)) then
      call this%model%ic_rad%compute(t, u%tf)
      associate (index => this%model%ic_rad%index, &
                 deriv => this%model%ic_rad%deriv)
        if (associated(this%model%void_face)) then
          do j = 1, size(index,2) !FIXME? Bad form to modify deriv?
            if (any(this%model%void_face(index(:,j)))) deriv(:,j) = 0.0_r8
          end do
        end if
        call dm%incr_interface_flux3(index, deriv) !TODO: rename these methods
      end associate
    end if

    !! Dirichlet fix-ups for void faces.
    call dm%set_dir_faces(more_dir_faces)

    !! Enclosure radiation contributions to the preconditioner.
    !! TODO: what about factorization coupling?  Is this still correct?
    if (associated(this%model%vf_rad_prob)) then
      do index = 1, size(this%model%vf_rad_prob)
        associate (faces => this%model%vf_rad_prob(index)%faces)
          allocate(values(size(faces)))
          call this%model%vf_rad_prob(index)%rhs_deriv(t, u%tf(faces), values)
          where (.not.this%model%vf_rad_prob(index)%fmask) values = 0
          call dm%incr_face_diag(faces, this%mesh%area(faces) * values)
          deallocate(values)
        end associate
      end do
    end if

    !! The matrix is now complete; re-compute the preconditioner.
    call this%pc%compute

    call stop_timer('ht-precon-compute')

  end subroutine compute


  subroutine apply(this, t, u, f)

    class(ht_precon), intent(in) :: this
    real(r8), intent(in) :: t
    type(ht_vector), intent(inout) :: u   ! data is intent(in)
    type(ht_vector), intent(inout) :: f   ! data is intent(inout)

    integer :: index, j, n
    real(r8), allocatable :: z(:)

    call start_timer('ht-precon-apply')

    !! Precondition the radiosity components:
    !! Factorization and forward Gauss-Seidel coupling methods.
    if (associated(this%model%vf_rad_prob)) then
      call start_timer('vf-rad-precon')
      !call HTSD_model_get_face_temp_view (this%model, f, f2)
      do index = 1, size(this%model%vf_rad_prob)
        if (this%vfr_precon_coupling(index) == VFR_FGS .or. &
            this%vfr_precon_coupling(index) == VFR_FAC) then
          z = f%encl(index)%qrad
          call this%model%vf_rad_prob(index)%precon(t, z)
          if (this%vfr_precon_coupling(index) == VFR_FGS) f%encl(index)%qrad = z
          !! Update the heat equation face residual.
          call this%model%vf_rad_prob(index)%precon_matvec1(t, z)
          do j = 1, size(z)
            if (this%model%vf_rad_prob(index)%fmask(j)) then
              n = this%model%vf_rad_prob(index)%faces(j)
              f%tf(n) = f%tf(n) + this%mesh%area(n) * z(j)
            end if
          end do
          deallocate(z)
        end if
      end do
      call stop_timer('vf-rad-precon')
    end if

    !! Heat equation cell residual with the H/T relation residual eliminated.
    !call HTSD_model_get_cell_heat_view (this%model, f, f0)
    !call HTSD_model_get_cell_temp_view (this%model, f, f1)
    !f1x(:this%mesh%ncell_onP) = f1 - (this%mesh%volume(:this%mesh%ncell_onP)/this%dt)*f0
    !f%tc = f%tc - (this%mesh%volume/dt)*f%hc

    if (associated(this%model%void_cell)) then
      where (.not.this%model%void_cell) f%tc = f%tc - (this%mesh%volume/this%dt)*f%hc
    else
      f%tc = f%tc - (this%mesh%volume/this%dt)*f%hc
    end if
    call this%mesh%cell_imap%gather_offp(f%tc)

    !! Heat equation face residual (with radiosity residuals optionally eliminated).
    !call HTSD_model_get_face_temp_copy (this%model, f, f2x)
    !this%mesh%face_imap%gather_offp(f2x)
    call this%mesh%face_imap%gather_offp(f%tf)

    !! Precondition the heat equation.
    !call diff_precon_apply (this%hcprecon, f1x, f2x)
    call this%pc%apply(f%tc, f%tf)
    !call HTSD_model_set_cell_temp (this%model, f1x, f)
    !call HTSD_model_set_face_temp (this%model, f2x, f)

    !! Backsubstitute to obtain the preconditioned H/T-relation residual.
    !f0 = f0 + this%dHdT(:this%mesh%ncell_onP)*f1
    f%hc = f%hc + this%dHdT*f%tc

    !! Precondition the radiosity components:
    !! Factorization, backward Gauss-Seidel and Jacobi coupling methods.
    if (associated(this%model%vf_rad_prob)) then
      call start_timer('vf-rad-precon')
      !call HTSD_model_get_face_temp_view (this%model, f, f2)
      !call HTSD_model_get_face_temp_view (this%model, u, Tface)
      do index = 1, size(this%model%vf_rad_prob)
        if (this%vfr_precon_coupling(index) == VFR_JAC .or. &
            this%vfr_precon_coupling(index) == VFR_BGS .or. &
            this%vfr_precon_coupling(index) == VFR_FAC) then
          associate (fq => f%encl(index)%qrad, faces => this%model%vf_rad_prob(index)%faces)
            !! Update the radiosity residual components.
            !TODO: call HTSD_model_get_radiosity_view (this%model, index, f, fq)
            if (this%vfr_precon_coupling(index) /= VFR_JAC) then
              allocate(z(size(fq)))
              call this%model%vf_rad_prob(index)%rhs_deriv(t, u%tf(faces), z)
              fq = fq + z * f%tf(faces)
              deallocate(z)
            end if
            call this%model%vf_rad_prob(index)%precon(t, fq)
          end associate
        end if
      end do
      call stop_timer('vf-rad-precon')
    end if

    call stop_timer('ht-precon-apply')

  end subroutine apply

end module ht_precon_type
