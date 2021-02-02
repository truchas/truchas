!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module HTSD_precon_type

  use kinds, only: r8
  use HTSD_model_type
  use rad_problem_type
  use unstr_mesh_type
  use mfd_diff_precon_type
  use mfd_diff_matrix_type
  use index_partitioning
  use truchas_timers
  implicit none
  private

  type, public :: HTSD_precon
    type(HTSD_model), pointer :: model => null()
    type(unstr_mesh), pointer :: mesh  => null()
    !! Heat transfer preconditioning data
    real(r8) :: dt
    real(r8), pointer :: dHdT(:) => null()
    type(mfd_diff_precon), pointer :: hcprecon => null()
    integer, pointer :: vfr_precon_coupling(:) => null()
    !! Species diffusion preconditioning data
    type(mfd_diff_precon), pointer :: sdprecon(:) => null()
    logical :: have_soret_coupling = .false.
  end type HTSD_precon

  public :: HTSD_precon_init
  public :: HTSD_precon_delete
  public :: HTSD_precon_compute
  public :: HTSD_precon_apply

  !! Methods of coupling heat conduction and radiosity preconditioning.
  integer, parameter :: VFR_JAC = 1 ! Jacobi (radiosity and conduction decoupled)
  integer, parameter :: VFR_FGS = 2 ! Forward Gauss-Seidel (radiosity, then conduction)
  integer, parameter :: VFR_BGS = 3 ! Backward Gauss-Seidel (conduction, then radiosity)
  integer, parameter :: VFR_FAC = 4 ! Factorization (radiosity, conduction, radiosity)

contains

  subroutine HTSD_precon_init (this, model, params)

    use parameter_list_type
    use truchas_logging_services

    type(HTSD_precon), intent(out) :: this
    type(HTSD_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params

    integer :: n, j, stat
    type(mfd_diff_matrix), allocatable, target :: matrix
    type(mfd_diff_matrix), pointer :: mold_matrix => null()
    character(:), allocatable :: string_array(:), errmsg

    this%model => model
    this%mesh  => model%mesh

    if (associated(model%ht)) then
      !! Heat density/temperature relation derivative.
      allocate(this%dHdT(this%mesh%ncell))
      !! Create the preconditioner for the heat equation.
      allocate(matrix)
      call matrix%init(model%disc)
      allocate(this%hcprecon)
      call this%hcprecon%init(matrix, params, stat, errmsg)
      if (stat /= 0) call TLS_fatal('HTSD_PRECON_INIT: ' // errmsg)
      mold_matrix => this%hcprecon%matrix_ref()
      !! Initialize the heat equation/view factor radiation
      if (associated(model%ht%vf_rad_prob)) then
        n = size(model%ht%vf_rad_prob)
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
    end if

    if (associated(model%sd)) then
      allocate(this%sdprecon(model%num_comp))
      do n = 1, this%model%num_comp
        allocate(matrix)
        if (associated(mold_matrix)) then
          call matrix%init(mold=mold_matrix)
        else
          call matrix%init(model%disc)
        end if
        call this%sdprecon(n)%init(matrix, params, stat, errmsg)
        if (stat /= 0) call TLS_fatal('HTSD_PRECON_INIT: ' // errmsg)
        if (.not.associated(mold_matrix)) mold_matrix => this%sdprecon(n)%matrix_ref()
        if (associated(model%sd(n)%soret)) this%have_soret_coupling = .true.
      end do
      this%have_soret_coupling = (associated(model%ht) .and. this%have_soret_coupling)
    end if

  end subroutine HTSD_precon_init

  subroutine HTSD_precon_delete (this)

    type(HTSD_precon), intent(inout) :: this

    integer :: n

    if (associated(this%dHdT)) deallocate(this%dHdT)
    if (associated(this%vfr_precon_coupling)) deallocate(this%vfr_precon_coupling)
    if (associated(this%hcprecon)) deallocate(this%hcprecon)
    if (associated(this%sdprecon)) deallocate(this%sdprecon)

  end subroutine HTSD_precon_delete

  subroutine HTSD_precon_compute (this, t, u, dt, errc)

    type(HTSD_precon), intent(inout) :: this
    real(r8), intent(in) :: t, u(:), dt
    integer, intent(out) :: errc
    target :: u

    integer :: n, j
    real(r8), pointer :: state(:,:) => null()
    integer, allocatable :: more_dir_faces(:)

    ASSERT(size(u) == HTSD_model_size(this%model))
    ASSERT(dt > 0.0_r8)

    call start_timer ('HTSD precon compute')

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

    state => HTSD_model_new_state_array(this%model, u)

    !! Compute the HT preconditioner.
    if (associated(this%model%ht)) call HT_precon_compute

    !! Compute the SD preconditioner.
    if (associated(this%model%sd)) then
      do n = 1, size(this%sdprecon)
        call SD_precon_compute (n)
      end do
    end if

    deallocate(state)

    errc = 0

    call stop_timer ('HTSD precon compute')

  contains

    subroutine HT_precon_compute ()

      integer :: index, j, n1, n2
      real(r8) :: D(this%mesh%ncell), A(this%mesh%ncell), Tface(this%mesh%nface), term
      real(r8), allocatable :: values(:), values2(:,:)
      type(mfd_diff_matrix), pointer :: dm
      integer, pointer :: faces(:) => null()

      call HTSD_model_get_face_temp_copy (this%model, u, Tface)
      call gather_boundary (this%mesh%face_ip, Tface)

      !! The time step size.
      this%dt = dt

      !! Hardwired assumption that T is the first component of state -- FIXME!
      call this%model%ht%H_of_T%compute_deriv(state, 1, this%dHdT)
      call this%model%ht%conductivity%compute_value(state, D)
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
      dm => this%hcprecon%matrix_ref()
      call dm%compute (D)
      call dm%incr_cell_diag (A)

      !! Dirichlet boundary condition fixups.
      if (allocated(this%model%ht%bc_dir)) then
        call this%model%ht%bc_dir%compute(t)
        call dm%set_dir_faces(this%model%ht%bc_dir%index)
      end if

      !! External HTC boundary condition contribution.
      if (allocated(this%model%ht%bc_htc)) then
        call this%model%ht%bc_htc%compute_deriv(t, Tface)
        associate (index => this%model%ht%bc_htc%index, &
                   deriv => this%model%ht%bc_htc%deriv)
          call dm%incr_face_diag(index, deriv)
        end associate
      end if

      !! Simple radiation boundary condition contribution.
      if (allocated(this%model%ht%bc_rad)) then
        call this%model%ht%bc_rad%compute(t, Tface)
        associate (index => this%model%ht%bc_rad%index, &
                   deriv => this%model%ht%bc_rad%deriv)
          call dm%incr_face_diag(index, deriv)
        end associate
      end if

      !! Experimental evaporation heat flux contribution.
      if (allocated(this%model%ht%evap_flux)) then
        call this%model%ht%evap_flux%compute_deriv(t, Tface)
        associate (index => this%model%ht%evap_flux%index, &
                   deriv => this%model%ht%evap_flux%deriv)
          call dm%incr_face_diag(index, this%mesh%area(index)*deriv)
        end associate
      endif

      !! Internal HTC interface condition contribution.
      if (allocated(this%model%ht%ic_htc)) then
        call this%model%ht%ic_htc%compute(t, Tface)
        associate (index => this%model%ht%ic_htc%index, &
                   deriv => this%model%ht%ic_htc%deriv)
          if (associated(this%model%void_face)) then
            do j = 1, size(index,2) !FIXME? Bad form to modify deriv?
              if (any(this%model%void_face(index(:,j)))) deriv(:,j) = 0.0_r8
            end do
          end if
          call dm%incr_interface_flux3(index, deriv) !TODO: rename these methods
        end associate
      end if

      !! Internal gap radiation condition contribution.
      if (allocated(this%model%ht%ic_rad)) then
        call this%model%ht%ic_rad%compute(t, Tface)
        associate (index => this%model%ht%ic_rad%index, &
                   deriv => this%model%ht%ic_rad%deriv)
          if (associated(this%model%void_face)) then
            do j = 1, size(index,2) !FIXME? Bad form to modify deriv?
              if (any(this%model%void_face(index(:,j)))) deriv(:,j) = 0.0_r8
            end do
          end if
          call dm%incr_interface_flux3(index, deriv) !TODO: rename these methods
        end associate
      end if

      !! Dirichlet fix-ups for void faces.
      call dm%set_dir_faces (more_dir_faces)

      !! Enclosure radiation contributions to the preconditioner.
      !! TODO: what about factorization coupling?  Is this still correct?
      if (associated(this%model%ht%vf_rad_prob)) then
        do index = 1, size(this%model%ht%vf_rad_prob)
          faces => this%model%ht%vf_rad_prob(index)%faces
          allocate(values(size(faces)))
          call this%model%ht%vf_rad_prob(index)%rhs_deriv (t, Tface(faces), values)
          where (.not.this%model%ht%vf_rad_prob(index)%fmask) values = 0
          call dm%incr_face_diag (faces, this%mesh%area(faces) * values)
          deallocate(values)
        end do
      end if

      !! The matrix is now complete; re-compute the preconditioner.
      call this%hcprecon%compute

    end subroutine HT_precon_compute

    subroutine SD_precon_compute (index)

      integer, intent(in) :: index

      real(r8) :: values(this%mesh%ncell)
      type(mfd_diff_matrix), pointer :: matrix

      matrix => this%sdprecon(index)%matrix_ref()
      !! Jacobian of the diffusion operator that ignores nonlinearities.
      call this%model%sd(index)%diffusivity%compute_value(state, values)
      if (associated(this%model%void_cell)) where (this%model%void_cell) values = 0.0_r8
      call matrix%compute (values)
      !! Time derivative contribution to the diffusion equation Jacobian.
      values = this%mesh%volume / dt
      if (associated(this%model%void_cell)) where (this%model%void_cell) values = 1.0_r8
      call matrix%incr_cell_diag (values)
      !! Dirichlet BC fixups.
      if (allocated(this%model%sd(index)%bc_dir)) then
        call this%model%sd(index)%bc_dir%compute(t)
        call matrix%set_dir_faces(this%model%sd(index)%bc_dir%index)
      end if
      call matrix%set_dir_faces (more_dir_faces)
      !! The matrix is now complete; re-compute the preconditioner.
      call this%sdprecon(index)%compute

    end subroutine SD_precon_compute

  end subroutine HTSD_precon_compute

  subroutine HTSD_precon_apply (this, t, u, f)

    use mfd_disc_type

#ifdef G95_COMPILER_WORKAROUND
    type(HTSD_precon), intent(inout) :: this
#else
    type(HTSD_precon), intent(in) :: this
#endif
    real(r8), intent(in) :: t, u(:)
    real(r8), intent(inout) :: f(:)
    target :: u, f

    integer :: index
    real(r8), allocatable :: FTcell(:), FTface(:)
    real(r8), pointer :: state(:,:) => null()

    call start_timer ('HTSD precon apply')

    !! Precondition the HT component.
    if (associated(this%model%ht)) then
      call start_timer ('HT precon apply')
      call HT_precon_apply
      call stop_timer ('HT precon apply')
    end if

    !! Precondition the SD components.
    if (associated(this%model%sd)) then
      call start_timer ('SD precon apply')
      !! Initialize extra data needed to handle Soret coupling.
      if (this%have_soret_coupling) then
        state => HTSD_model_new_state_array(this%model, u)
        !! Off-process-extended copies of the preconditioned HT components.
        allocate(FTcell(this%mesh%ncell), FTface(this%mesh%nface))
        call HTSD_model_get_cell_temp_copy (this%model, f, FTcell)
        call gather_boundary (this%mesh%cell_ip, FTcell)
        call HTSD_model_get_face_temp_copy (this%model, f, FTface)
        call gather_boundary (this%mesh%face_ip, FTface)
        if (allocated(this%model%ht%bc_dir)) then
          FTface(this%model%ht%bc_dir%index) = 0.0_r8 ! temperature Dirichlet projection
        end if
        !TODO! void face dirichlet projection?
      end if
      do index = 1, this%model%num_comp
        call SD_precon_apply (index)
      end do
      if (associated(state)) deallocate(state)
      call stop_timer ('SD precon apply')
    end if

    call stop_timer ('HTSD precon apply')

  contains

    subroutine HT_precon_apply ()

      integer :: index, j, n
      real(r8), pointer :: f0(:), f1(:), f2(:), fq(:), Tface(:)
      real(r8) :: f1x(this%mesh%ncell), f2x(this%mesh%nface)
      real(r8), allocatable :: z(:)

      !! Glossary:
      !! Pointers into segments of the F array argument:
      !!  f0  -- cell enthalpy segment
      !!  f1  -- cell temperature segment
      !!  f2  -- face temperature segment
      !!  fq  -- a radiosity segment
      !! Local arrays:
      !!  f1x -- off-process extended copy of f1
      !!  f2x -- off-process extended copy of f2

      !! Precondition the radiosity components:
      !! Factorization and forward Gauss-Seidel coupling methods.
      if (associated(this%model%ht%vf_rad_prob)) then
        call start_timer ('VF rad precon')
        call HTSD_model_get_face_temp_view (this%model, f, f2)
        do index = 1, size(this%model%ht%vf_rad_prob)
          if (this%vfr_precon_coupling(index) == VFR_FGS .or. &
              this%vfr_precon_coupling(index) == VFR_FAC) then
            call HTSD_model_get_radiosity_view (this%model, index, f, fq)
            allocate(z(size(fq)))
            z = fq
            call this%model%ht%vf_rad_prob(index)%precon (t, z)
            if (this%vfr_precon_coupling(index) == VFR_FGS) fq = z
            !! Update the heat equation face residual.
            call this%model%ht%vf_rad_prob(index)%precon_matvec1 (t, z)
            do j = 1, size(z)
              if (this%model%ht%vf_rad_prob(index)%fmask(j)) then
                n = this%model%ht%vf_rad_prob(index)%faces(j)
                f2(n) = f2(n) + this%mesh%area(n) * z(j)
              end if
            end do
            deallocate(z)
          end if
        end do
        call stop_timer ('VF rad precon')
      end if

      !! Heat equation cell residual with the H/T relation residual eliminated.
      call HTSD_model_get_cell_heat_view (this%model, f, f0)
      call HTSD_model_get_cell_temp_view (this%model, f, f1)
      f1x(:this%mesh%ncell_onP) = f1 - (this%mesh%volume(:this%mesh%ncell_onP)/this%dt)*f0
      if (associated(this%model%void_cell)) then
        !! Restore heat equation cell residual on void cells.
        do j = 1, this%mesh%ncell_onP
          if (this%model%void_cell(j)) f1x(j) = f1(j)
        end do
      end if
      call gather_boundary (this%mesh%cell_ip, f1x)

      !! Heat equation face residual (with radiosity residuals optionally eliminated).
      call HTSD_model_get_face_temp_copy (this%model, f, f2x)
      call gather_boundary (this%mesh%face_ip, f2x)

      !! Precondition the heat equation.
      call this%hcprecon%apply(f1x, f2x)
      call HTSD_model_set_cell_temp (this%model, f1x, f)
      call HTSD_model_set_face_temp (this%model, f2x, f)

      !! Backsubstitute to obtain the preconditioned H/T-relation residual.
      f0 = f0 + this%dHdT(:this%mesh%ncell_onP)*f1

      !! Precondition the radiosity components:
      !! Factorization, backward Gauss-Seidel and Jacobi coupling methods.
      if (associated(this%model%ht%vf_rad_prob)) then
        call start_timer ('VF rad precon')
        call HTSD_model_get_face_temp_view (this%model, f, f2)
        call HTSD_model_get_face_temp_view (this%model, u, Tface)
        do index = 1, size(this%model%ht%vf_rad_prob)
          if (this%vfr_precon_coupling(index) == VFR_JAC .or. &
              this%vfr_precon_coupling(index) == VFR_BGS .or. &
              this%vfr_precon_coupling(index) == VFR_FAC) then
            !! Update the radiosity residual components.
            call HTSD_model_get_radiosity_view (this%model, index, f, fq)
            if (this%vfr_precon_coupling(index) /= VFR_JAC) then
              allocate(z(size(fq)))
              call this%model%ht%vf_rad_prob(index)%rhs_deriv (t, Tface(this%model%ht%vf_rad_prob(index)%faces), z)
              fq = fq + z * f2(this%model%ht%vf_rad_prob(index)%faces)
              deallocate(z)
            end if
            call this%model%ht%vf_rad_prob(index)%precon (t, fq)
          end if
        end do
        call stop_timer ('VF rad precon')
      end if

    end subroutine HT_precon_apply

    subroutine SD_precon_apply (index)

      integer, intent(in) :: index

      real(r8) :: Fcell(this%mesh%ncell), Fface(this%mesh%nface)
      real(r8), allocatable :: value(:)
      real(r8), pointer :: fview(:)

      !! Lower triangle coupling of HT and SD; forward elimination.
      if (associated(this%model%sd(index)%soret)) then
        !! Compute the update.
        allocate(value(this%mesh%ncell))
        call this%model%sd(index)%soret%compute_value(state, value)
        call this%model%sd(index)%diffusivity%compute_value(state, Fcell) ! Fcell used as temporary
        value = value * Fcell
        call this%model%disc%apply_diff (value, FTcell, FTface, Fcell, Fface)
        if (allocated(this%model%sd(index)%bc_dir)) then
          Fface(this%model%sd(index)%bc_dir%index) = 0.0_r8 ! concentration Dirichlet projection
        end if
        !TODO! void face dirichlet projection?
        !! Apply the update.
        call HTSD_model_get_cell_conc_view (this%model, index, f, fview)
        fview = fview - Fcell(:this%mesh%ncell_onP)
        call HTSD_model_get_face_conc_view (this%model, index, f, fview)
        fview = fview - Fface(:this%mesh%nface_onP)
      end if

      !! Off-process extended cell concentration components of F.
      call HTSD_model_get_cell_conc_copy (this%model, index, f, Fcell)
      call gather_boundary (this%mesh%cell_ip, Fcell)
      !! Off-process extended face concentration components of F.
      call HTSD_model_get_face_conc_copy (this%model, index, f, Fface)
      call gather_boundary (this%mesh%face_ip, Fface)
      !! Precondition the diffusion equation for this component.
      call this%sdprecon(index)%apply(Fcell, Fface)
      !! Return the on-process components of the result.
      call HTSD_model_set_cell_conc (this%model, index, Fcell, f)
      call HTSD_model_set_face_conc (this%model, index, Fface, f)

    end subroutine SD_precon_apply

  end subroutine HTSD_precon_apply

end module HTSD_precon_type
