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
  use sm_bc_type
  use truchas_timers
  implicit none
  private

  !! TODO extend a generic sm precon class
  type, public :: sm_ds_precon
    private
    type(sm_model), pointer, public :: model => null() ! unowned reference
    type(sm_bc), pointer, public :: bc => null() ! unowned reference

    real(r8), allocatable :: diag(:,:), d1(:,:), d2(:,:)
    real(r8) :: omega
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
    allocate(this%diag(3,model%mesh%nnode))
    call params%get('num-iter', this%niter, default=1)
    call params%get('precon-relaxation-parameter', this%omega, default=1.0_r8)

    this%diag = 0

    call compute_diagonals

  contains

    !! This computes the diagonal of the jacobian matrix by effectively
    !! multiplying that matrix by a series of matrices with only one nonzero
    !! element. This is done by alternating out local displacements and
    !! computing local stress contributions, so it should be fairly fast.
    subroutine compute_diagonals()

      integer :: n, d, p, k, xp, xnc, s
      real(r8) :: displ(3), lhs(3), strain(6), stress(6)

      associate (mesh => this%model%mesh, ig => this%model%ig)

        allocate(this%d1(3,mesh%nnode_onP), this%d2(3,mesh%nnode_onP))

        do n = 1, mesh%nnode_onP
          this%d1(:,n) = 0
          this%d2(:,n) = 0

          do d = 1, 3
            displ = 0
            displ(d) = 1

            do xp = ig%xnpoint(n), ig%xnpoint(n+1)-1
              k = xp - ig%xnpoint(n) + 1
              p = ig%npoint(xp)
              s = merge(-1, 1, btest(ig%nppar(n),k))

              ! This is the local ID for node n according to the cell associated with p
              xnc = ig%xpxn(xp)
              call compute_strain_contribution(this, p, xnc, displ, strain)

              call this%model%compute_stress(1.0_r8, 0.0_r8, strain, stress)
              lhs = this%model%tensor_dot(stress, ig%n(:,p))
              this%d1(d,n) = this%d1(d,n) + s * lhs(d)

              call this%model%compute_stress(0.0_r8, 1.0_r8, strain, stress)
              lhs = this%model%tensor_dot(stress, ig%n(:,p))
              this%d2(d,n) = this%d2(d,n) + s * lhs(d)
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

    strain(1) = grad_displ(1,1) ! exx
    strain(2) = grad_displ(2,2) ! eyy
    strain(3) = grad_displ(3,3) ! ezz
    strain(4) = (grad_displ(1,2) + grad_displ(2,1)) / 2 ! exy
    strain(5) = (grad_displ(1,3) + grad_displ(3,1)) / 2 ! exz
    strain(6) = (grad_displ(2,3) + grad_displ(3,2)) / 2 ! eyz

  end subroutine compute_strain_contribution


  subroutine compute(this, t, dt, displ)

    use index_partitioning, only: gather_boundary

    class(sm_ds_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt, displ(:,:)

    integer :: n, d, i, f, xn
    real(r8) :: force(3,this%model%mesh%nnode), displ_(3,this%model%mesh%nnode)

    call start_timer("precon-compute")

    do n = 1, this%model%mesh%nnode_onP
      if (this%model%lame1_n(n) < 1e-6_r8 .and. this%model%lame2_n(n) < 1e-6_r8) then
        ! Set displacements to zero for empty or fluid filled cells
        this%diag(:,n) = 1
      else
        this%diag(:,n) = this%model%lame1_n(n) * this%d1(:,n) + this%model%lame2_n(n) * this%d2(:,n)
      end if
      ASSERT(all(this%diag(:,n) /= 0))
    end do

    ! get off-rank halo
    displ_(:,:this%model%mesh%nnode_onP) = displ
    call gather_boundary(this%model%mesh%node_ip, displ_)
    call gather_boundary(this%model%mesh%node_ip, this%diag)
    if (this%bc%contact_active) then
      call this%model%compute_forces(t, displ_, force)
      call gather_boundary(this%model%mesh%node_ip, force)
    end if
    call this%model%bc%apply_deriv_diagonal(t, this%model%scaling_factor, displ_, force, this%diag)

    do n = 1, this%model%mesh%nnode_onP
      this%diag(:,n) = this%diag(:,n) / this%model%scaling_factor(n)
    end do

    call stop_timer("precon-compute")

  end subroutine compute


  subroutine apply(this, u, f)

    class(sm_ds_precon), intent(in), target :: this
    real(r8), intent(in), contiguous :: u(:,:) ! current displacement guess
    real(r8), intent(inout), contiguous :: f(:) ! in residual, out next displacement guess

    integer :: i, j
    real(r8) :: x
    real(r8), pointer :: diag(:) => null()

    call start_timer("precon-apply")

    !! TODO: reshape f outside this routine instead.
    diag(1:size(this%diag)) => this%diag

    do j = 1, size(f)
      !f(j) = f(j) / this%diag(j)
      x = f(j) / diag(j)
      do i = 1, this%niter
        f(j) = f(j) + this%omega*(x-f(j))
        !f(j) = (f(j) - this%omega*f(j)) + this%omega*x
      end do
    end do

    call stop_timer("precon-apply")

  end subroutine apply

end module sm_ds_precon_type
