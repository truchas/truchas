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

    real(r8), allocatable :: diag(:), d1(:), d2(:)
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
    allocate(this%diag(model%size()))
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

      integer :: n, d, j, p, k, xp, xnc, s
      real(r8) :: displ(3), lhs(3), strain(6), stress(6)

      associate (mesh => this%model%mesh, ig => this%model%ig)

        allocate(this%d1(3*mesh%nnode_onP), this%d2(3*mesh%nnode_onP))

        do n = 1, mesh%nnode_onP
          do d = 1, 3
            j = 3*(n-1) + d
            displ = 0
            displ(d) = 1

            this%d1(j) = 0
            this%d2(j) = 0
            do xp = ig%xnpoint(n), ig%xnpoint(n+1)-1
              k = xp - ig%xnpoint(n) + 1
              p = ig%npoint(xp)
              s = merge(-1, 1, btest(ig%nppar(n),k))

              ! This is the local ID for node n according to the cell associated with p
              xnc = ig%xpxn(xp)
              call compute_strain_contribution(this, p, xnc, displ, strain)

              call this%model%compute_stress(1.0_r8, 0.0_r8, strain, stress)
              lhs = this%model%tensor_dot(stress, ig%n(:,p))
              this%d1(j) = this%d1(j) + s * lhs(d)

              call this%model%compute_stress(0.0_r8, 1.0_r8, strain, stress)
              lhs = this%model%tensor_dot(stress, ig%n(:,p))
              this%d2(j) = this%d2(j) + s * lhs(d)
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

    class(sm_ds_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt, displ(:,:)

    integer :: n, d, j, i, f, xn
    real(r8) :: force(size(displ,dim=1),size(displ,dim=2))

    call start_timer("precon-compute")

    do n = 1, this%model%mesh%nnode_onP
      if (this%model%lame1_n(n) < 1e-6_r8 .and. this%model%lame2_n(n) < 1e-6_r8) then
        ! Set displacements to zero for empty or fluid filled cells
        j = 3*(n-1) + 1
        this%diag(j:j+2) = 1
      else
        do d = 1, 3
          j = 3*(n-1) + d
          this%diag(j) = (this%model%lame1_n(n) * this%d1(j) + this%model%lame2_n(n) * this%d2(j))
        end do
      end if
      ASSERT(this%diag(j) /= 0)
    end do

    do n = 1, this%model%mesh%nnode_onP
      j = 3*(n-1) + 1
      this%diag(j:j+2) = this%diag(j:j+2) / this%model%scaling_factor(n)
    end do

    call gap_conditions

    ! Dirichlet BCs
    do d = 1, 3
      call this%bc%displacement(d)%compute(t)
      associate (nodes => this%bc%displacement(d)%index)
        do i = 1, size(nodes)
          n = nodes(i)
          j = 3*(n-1) + d
          this%diag(j) = this%model%contact_penalty
        end do
      end associate
    end do

    ! The residual passed in will have, in the surface-aligned reference frame,
    ! the difference between the test normal displacement and the BC normal
    ! displacement. So we want the application of the preconditioner (inverse of
    ! diagonal) to result in this difference in that rotated reverence frame.
    call this%bc%displacementn%compute(t)
    associate (nodes => this%bc%displacementn%index, &
        rot => this%bc%displacementn%rotation_matrix, &
        normal => this%bc%displacementn%rotation_matrix(3,:,:))
      do i = 1, size(nodes)
        n = nodes(i)
        associate (rn => this%diag(3*(n-1)+1:3*(n-1)+3))
          rn = rn - rn * normal(:,i)**2
          rn = rn - this%model%contact_penalty * normal(:,i)**2
        end associate
      end do
    end associate

    call stop_timer("precon-compute")

  contains

    subroutine gap_conditions()

      associate (link => this%bc%gap_contact%index, dvalues => this%bc%gap_contact%dvalue)

        if (this%bc%gap_contact%enabled) then
          ! TODO-WARN: Need halo node displacements and stresses/residuals? For now just get
          !            everything working in serial.
          call this%model%compute_forces(t, displ, force)
          call this%bc%gap_contact%compute_deriv(t, displ, force, this%model%scaling_factor, this%diag)
#ifndef NDEBUG
          do i = 1, size(link, dim=2)
            block
              integer :: n1, n2
              n1 = link(1,i)
              n2 = link(2,i)
              if (n1 <= this%model%mesh%nnode_onP .or. n2 <= this%model%mesh%nnode_onP) then
                ASSERT(n1 <= this%model%mesh%nnode_onP .and. n2 <= this%model%mesh%nnode_onP)
              end if
            end block
          end do
#endif
        end if

        do i = 1, size(link, dim=2)
          block
            integer :: n1, n2
            real(r8) :: s, tn, l, dl(2)

            n1 = link(1,i)
            n2 = link(2,i)

            associate (r1 => this%diag(3*(n1-1)+1:3*(n1-1)+3), &
                r2 => this%diag(3*(n2-1)+1:3*(n2-1)+3))

              if (n1 <= this%model%mesh%nnode_onP) r1 = r1 + dvalues(:,1,i)
              if (n2 <= this%model%mesh%nnode_onP) r2 = r2 + dvalues(:,2,i)

              ! stress1 = matmul(rot(:,:,i), r(:,n1)) !+ this%rhs(:,n1)
              ! stress2 = matmul(rot(:,:,i), r(:,n2)) !+ this%rhs(:,n2)
              ! x1 = matmul(rot(:,:,i), displ(:,n1))
              ! x2 = matmul(rot(:,:,i), displ(:,n2))

              ! s = x1(3) - x2(3)
              ! tn = stress1(3)
              ! l = this%bc%gap_contact%contact_factor(s, tn)
              ! dl = this%bc%gap_contact%derivative_contact_factor(s, tn)
              ! dldu1 =  dl(1) + dl(2)*stress1(3)
              ! dldu2 = -dl(1) + dl(2)*stress2(3)

              ! ! In the first node we put the equal & opposite normal contact force constraint
              ! if (n1 <= this%model%mesh%nnode_onP) then
              !   !r1 = 1 / r1
              !   r1 = matmul(rot(:,:,i), r1)
              !   !r(1:2,n1) = stress1(1:2) ! If there is a sliding constraint... TODO: is this right?
              !   !r1(3) = stress1(3) + stress2(3)
              !   !r1(3) = r1(3) - l*1d3*1
              !   r1(3) = r1(3) - l*1d3*1 + dldu1*(stress2(3) + 1d3*(x2(3) - x1(3)))
              !   !r1(3) = 1
              !   r1 = matmul(transpose(rot(:,:,i)), r1)
              !   !r1 = 1 / r1
              ! end if

              ! ! In the second node we put the zero-displacement constraint
              ! if (n2 <= this%model%mesh%nnode_onP) then
              !   !r2 = 1 / r2
              !   r2 = matmul(rot(:,:,i), r2)
              !   !r(1:2,n2) = stress2(1:2) ! If there is a sliding constraint... TODO: is this right?
              !   !r2(3) = -1
              !   r2(3) = r2(3) - 1d3*1
              !   !r2(3) = -1
              !   r2 = matmul(transpose(rot(:,:,i)), r2)
              !   !r2 = 1 / r2
              ! end if
            end associate
          end block
        end do

      end associate

    end subroutine gap_conditions

  end subroutine compute


  subroutine apply(this, u, f)

    class(sm_ds_precon), intent(in) :: this
    real(r8), intent(in), contiguous :: u(:,:) ! current displacement guess
    real(r8), intent(inout), contiguous :: f(:) ! in residual, out next displacement step guess

    integer :: i, j
    real(r8) :: d

    call start_timer("precon-apply")

    do j = 1, size(this%diag)
      !f(j) = f(j) / this%diag(j)
      d = f(j) / this%diag(j)
      do i = 1, this%niter
        f(j) = f(j) + this%omega*(d-f(j))
        !f(j) = (f(j) - this%omega*f(j)) + this%omega*d
      end do
    end do

    call stop_timer("precon-apply")

  end subroutine apply

end module sm_ds_precon_type
