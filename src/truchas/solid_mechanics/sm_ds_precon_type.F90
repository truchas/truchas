!!
!! Diagonal scaling preconditioner for solid mechanics.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! August 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_ds_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use sm_model_type
  use sm_bc_manager_type
  use truchas_timers
  implicit none
  private

  !! TODO extend a generic sm precon class
  type, public :: sm_ds_precon
    private
    type(sm_model), pointer, public :: model => null() ! unowned reference
    type(sm_bc_manager), pointer, public :: bc => null() ! unowned reference

    real(r8), allocatable :: diag(:,:), F(:,:,:), d1(:,:,:), d2(:,:,:)
    real(r8) :: omega, gamma
    integer :: niter
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type sm_ds_precon

contains

  subroutine init(this, model, params)

    use parameter_list_type

    class(sm_ds_precon), intent(inout) :: this
    type(sm_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params

    this%model => model
    this%bc => model%bc
    call params%get('num-iter', this%niter, default=1)
    call params%get('relaxation-parameter', this%omega, default=1.0_r8)
    call params%get('stress-relaxation-parameter', this%gamma, default=1.0_r8) ! legacy used 16/9

    allocate(this%diag(3,model%mesh%nnode_onP))
    this%diag = 0
    call compute_diagonals

  contains

    !! This computes the diagonal of the jacobian matrix by effectively
    !! multiplying that matrix by a series of matrices with only one nonzero
    !! element. This is done by alternating out local displacements and
    !! computing local stress contributions, so it should be fairly fast.
    subroutine compute_diagonals()

      use sm_bc_utilities, only: compute_stress

      integer :: n, d, p, k, xp, xnc, s
      real(r8) :: displ(3), lhs(3), strain(6), stress(6)

      associate (mesh => this%model%mesh, ig => this%model%ig)

        allocate(this%d1(3,3,mesh%nnode_onP), this%d2(3,3,mesh%nnode_onP), &
            this%F(3,3,mesh%nnode_onP))

        do n = 1, mesh%nnode_onP
          this%d1(:,:,n) = 0
          this%d2(:,:,n) = 0

          do d = 1, 3
            displ = 0
            displ(d) = 1

            do xp = ig%xnpoint(n), ig%xnpoint(n+1)-1
              k = xp - ig%xnpoint(n) + 1
              p = ig%npoint(xp)
              s = merge(-1, 1, ig%nppar(k,n))

              ! This is the local ID for node n according to the cell associated with p
              xnc = ig%xpxn(xp)
              call compute_strain_contribution(this, p, xnc, displ, strain)

              call compute_stress(1.0_r8, 0.0_r8, strain, stress)
              lhs = this%model%tensor_dot(stress, ig%n(:,p))
              this%d1(:,d,n) = this%d1(:,d,n) + s * lhs

              call compute_stress(0.0_r8, 1.0_r8, strain, stress)
              lhs = this%model%tensor_dot(stress, ig%n(:,p))
              this%d2(:,d,n) = this%d2(:,d,n) + s * lhs
            end do
          end do
        end do

      end associate

    end subroutine compute_diagonals

  end subroutine init


  !! Compute the strain contribution to the integration point p from the
  !! displacement at the node with ID xnc local to the cell containing p. This
  !! is the strain contribution from only the given displacement, and no others.
  !! This is a simplified version of sm_model_type::compute_total_strain.
  subroutine compute_strain_contribution(this, p, xnc, displ, strain)

    class(sm_ds_precon), intent(in) :: this
    integer, intent(in) :: p, xnc
    real(r8), intent(in) :: displ(:)
    real(r8), intent(out) :: strain(:)

    integer :: d
    real(r8) :: grad_displ(3,3)

    ASSERT(p > 0 .and. p <= this%model%ig%npt)

    do d = 1, 3
      grad_displ(:,d) = this%model%ig%grad_shape(p)%p(:,xnc) * displ(d)
    end do
    strain = this%model%strain_tensor(grad_displ)

  end subroutine compute_strain_contribution


  subroutine compute(this, t, dt, displ)

    class(sm_ds_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(inout) :: displ(:,:) ! need to update halo

    integer :: n
    real(r8) :: force(3,this%model%mesh%nnode)

    ASSERT(size(displ,dim=1) == 3 .and. size(displ,dim=2) >= this%model%mesh%nnode)

    call start_timer("precon-compute")

    ! For contact, we need the force.
    force = 0
    call this%model%mesh%node_imap%gather_offp(displ)
    if (this%bc%contact_active) then !.or. this%model%matl_model%viscoplasticity_enabled) then
      call this%model%compute_forces(t, displ, force)
      call this%model%mesh%node_imap%gather_offp(force)
    end if

    do n = 1, this%model%mesh%nnode_onP
      if (this%model%lame1_n(n) < 1e-6_r8 .and. this%model%lame2_n(n) < 1e-6_r8) then
        ! Set displacements to zero for empty or fluid filled cells
        this%diag(:,n) = 1
        this%F(:,:,n) = 0
        this%F(1,1,n) = 1
        this%F(2,2,n) = 1
        this%F(3,3,n) = 1
      else
        this%F(:,:,n) = this%model%lame1_n(n) * this%d1(:,:,n) + this%model%lame2_n(n) * this%d2(:,:,n)

        ! NB: Viscoplasticity is neglected
        ! if (this%model%matl_model%viscoplasticity_enabled) &
        !     call this%compute_viscoplasticity_precon_contribution(dt, n, this%F(:,:,n))

        this%diag(1,n) = this%F(1,1,n)
        this%diag(2,n) = this%F(2,2,n)
        this%diag(3,n) = this%F(3,3,n)

        ! under-relaxation
        this%F(:,:,n) = this%gamma * this%F(:,:,n)
        this%diag(:,n) = this%gamma * this%diag(:,n)
      end if
      ASSERT(all(this%diag(:,n) /= 0))
    end do

    call this%model%bc%apply_deriv_diagonal(t, this%model%scaling_factor, displ, force, this%diag, this%F)

    do n = 1, this%model%mesh%nnode_onP
      this%diag(:,n) = this%diag(:,n) / this%model%scaling_factor(n)
    end do

    call stop_timer("precon-compute")

  end subroutine compute


  subroutine apply(this, u, f)

    class(sm_ds_precon), intent(in), target :: this
    real(r8), intent(in), contiguous :: u(:,:) ! current displacement guess
    real(r8), intent(inout), contiguous :: f(:,:) ! in residual, out next displacement guess

    integer :: i, j
    real(r8) :: x(3)

    call start_timer("precon-apply")

    do j = 1, this%model%mesh%nnode_onP
      !f(:,j) = f(:,j) / this%diag(:,j)
      x = f(:,j) / this%diag(:,j)
      do i = 1, this%niter
        f(:,j) = f(:,j) + this%omega * (x - f(:,j))
        !f(:,j) = (f(:,j) - this%omega*f(:,j)) + this%omega*x
      end do
    end do

    call stop_timer("precon-apply")

  end subroutine apply

  ! TODO: Add the viscoplasticity contribution to the preconditioner.
  ! subroutine compute_viscoplasticity_precon_contribution(this, dt, vof, temperature, displ, F)

  !   class(sm_ds_precon), intent(inout) :: this
  !   real(r8), intent(in) :: dt
  !   real(r8), intent(inout) :: displ(:,:) ! need to update halo

  !   real(r8) :: precon(6,6,ig%npt)

  !   ! the total preconditioner is effectively:
  !   !   P_ij = dstress_i / ddispl_j
  !   !        = (dstress_i / dstrain_elastic_k) * (dstrain_elastic_k / ddispl_j
  !   !        = (dstress_i / dstrain_elastic_k) * [ (dstrain_total_k / ddispl_j)
  !   !                                            - (dstrain_plastic_k / ddispl_j)]
  !   ! since strain_elastic ~ strain_total - strain_plastic
  !   !   strain_plastic = strain_plastic_old + delta_strain_plastic
  !   !   delta_strain_plastic ~ dt * F(stressM)
  !   ! where stressM is a midpoint estimate of stress:
  !   !   stressM = stressM{(strain_total_old + strain_total) / 2 - strain_plastic_old - ...}
  !   ! Thus:
  !   !   dstrain_plastic_i / ddispl_j ~ dt * (dF_i / dstressM_k) * (dstressM_k / ddispl_j)
  !   ! But (dstressM_k / ddispl_j) = P_kj / 2. We estimate the viscoplastic contribution
  !   ! to the preconditioner as:
  !   !   P_ij <- P_ij - dt/2 * (dF_i / dstressM_k) * P_kj

  !   call this%model%compute_viscoplastic_precon(precon)
  !   call this%model%viscoplastic_model%compute_precon_stress(, viscoplastic_precon)
  !   viscoplastic_precon = matmul(viscoplastic_precon, F)
  !   do i = 1, 6
  !     do j = 1, 6
  !       F(i,j) = F(i,j) - dt/2 * viscoplastic_precon
  !     end do
  !   end do

  ! end subroutine compute_viscoplasticity_precon_contribution

end module sm_ds_precon_type
