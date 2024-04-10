!!
!! Hypre preconditioner for solid mechanics.
!! Similar to mfd_diff_precon_type.F90.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! April 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_hypre_precon_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use sm_precon_class
  use sm_model_type
  use sm_bc_manager_type
  use pcsr_matrix_type
  use pcsr_precon_class
  use index_map_type
  use truchas_timers
  use truchas_logging_services
  implicit none
  private

  type matrix_box
    real(r8), allocatable :: p(:,:)
  end type matrix_box

  type, extends(sm_precon), public :: sm_hypre_precon
    private
    type(sm_model), pointer, public :: model => null() ! unowned reference
    type(sm_bc_manager), pointer, public :: bc => null() ! unowned reference

    type(pcsr_matrix), pointer :: A => null()
    type(index_map), pointer :: imap => null()
    class(pcsr_precon), allocatable :: precon

    type(matrix_box), allocatable :: displ_to_force(:)
    integer :: niter
    logical :: first = .true.
    real(r8), allocatable :: l1(:), l2(:)
  contains
    final :: sm_hypre_precon_delete
    procedure :: init
    procedure :: compute
    procedure :: apply
    procedure, private :: recompute_needed
  end type sm_hypre_precon

contains

  !! Final subroutine for emfd_nlsol_solver type objects.
  subroutine sm_hypre_precon_delete(this)
    type(sm_hypre_precon), intent(inout) :: this
    if (associated(this%imap)) deallocate(this%imap)
    if (associated(this%A)) deallocate(this%A)
  end subroutine sm_hypre_precon_delete


  subroutine init(this, model, params)

    use parameter_list_type
    use pcsr_precon_factory
    use parallel_communication, only: is_IOP

    class(sm_hypre_precon), intent(out) :: this
    type(sm_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params

    type(parameter_list), pointer :: plist, plist_params
    type(pcsr_graph), pointer :: g
    type(index_map), pointer :: row_imap
    integer, allocatable :: nvars(:)
    integer :: stat, c, clique(24), nnode
    integer :: n1, n1x, n2, n2x, ii, jj
    character(:), allocatable :: errmsg

    this%model => model
    this%bc => model%bc
    allocate(this%displ_to_force(this%model%ig%npt))

    !! Create a CSR matrix graph for the node-coupled system
    associate(mesh => model%mesh)
      allocate(g, this%imap, this%A)
      row_imap => mesh%node_imap
      allocate(nvars(merge(mesh%node_imap%global_size, 0, is_IOP)))
      nvars = 3
      call this%imap%init(row_imap, nvars)

      call g%init(this%imap)
      do c = 1, mesh%ncell
        associate (cn => mesh%cnode(mesh%xcnode(c):mesh%xcnode(c+1)-1))
          nnode = size(cn)
          clique(1:nnode) = 3*(cn - 1) + 1
          clique(nnode+1:2*nnode) = 3*(cn - 1) + 2
          clique(2*nnode+1:3*nnode) = 3*(cn - 1) + 3
          call g%add_clique(clique(:3*nnode))
          ! call g%add_clique(3*(cn - 1) + 1)
          ! call g%add_clique(3*(cn - 1) + 2)
          ! call g%add_clique(3*(cn - 1) + 3)

          ! do n1x = 1, size(cn)
          !   n1 = cn(n1x)
          !   do ii = 3*(n1 - 1) + 1, 3*(n1 - 1) + 3
          !     call g%add_edge(ii, ii)
          !   end do
          ! end do

          ! do n1x = 1, size(cn)
          !   n1 = cn(n1x)
          !   do n2x = 1, n1x
          !     n2 = cn(n2x)
          !     do ii = 3*(n1-1) + 1, 3*(n1-1) + 3
          !       do jj = 3*(n2-1) + 1, 3*(n2-1) + 3
          !         call g%add_edge(ii, jj)
          !         call g%add_edge(jj, ii)
          !       end do
          !     end do
          !   end do
          ! end do
        end associate
      end do
      call g%add_complete
      call this%A%init(g, take_graph=.true.)
    end associate

    !call params%get('num-iter', this%niter, default=1)
    plist => params%sublist("precon")
    plist_params => plist%sublist("params")
    call plist%set("method", "boomeramg")
    !call plist%set("method", "ssor")
    call plist_params%set("num-cycles", 2)
    ! call plist_params%set("print-level", 3)
    ! call plist_params%set("debug-level", 1)
    call alloc_pcsr_precon(this%precon, this%A, plist, stat, errmsg)
    if (stat /= 0) call tls_fatal("SOLID MECHANICS PRECON INIT: " // errmsg)

  end subroutine init


  subroutine compute(this, t, dt, displ)

    class(sm_hypre_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(inout) :: displ(:,:) ! need to update halo

    integer :: c, ii, jj, p, n, n1, n2, n3, xp, k, s, nn, nn1, nn2, nn3, i, j, xcn
    !real(r8) :: force(3,this%model%mesh%nnode)
    real(r8) :: strain_matrix(6,9), grad_shape(9,24), tensor_dot(3,6), stress_matrix(6,6)

    ASSERT(size(displ,dim=1) == 3 .and. size(displ,dim=2) >= this%model%mesh%nnode)

    if (.not.this%recompute_needed()) return

    call start_timer("precon-compute")

    call this%A%set_all(0.0_r8)

    ! For contact, we need the force.
    ! force = 0
    ! call this%model%mesh%node_imap%gather_offp(displ)
    ! if (this%bc%contact_active) then !.or. this%model%matl_model%viscoplasticity_enabled) then
    !   call this%model%compute_forces(t, displ, force)
    !   call this%model%mesh%node_imap%gather_offp(force)
    ! end if

    associate (mesh => this%model%mesh, ig => this%model%ig)
      strain_matrix = 0
      strain_matrix(1,1) = 1
      strain_matrix(2,5) = 1
      strain_matrix(3,9) = 1
      strain_matrix(4,[2,4]) = 0.5_r8
      strain_matrix(5,[3,7]) = 0.5_r8
      strain_matrix(6,[6,8]) = 0.5_r8

      do c = 1, mesh%ncell
        stress_matrix = 0
        do ii = 1, 6
          stress_matrix(ii,ii) = 2 * this%model%lame2(c)
        end do
        stress_matrix(:3,:3) = stress_matrix(:3,:3) + this%model%lame1(c)

        associate (cn => mesh%cnode(mesh%xcnode(c):mesh%xcnode(c+1)-1))
          do p = ig%xcpoint(c), ig%xcpoint(c+1)-1
            grad_shape = 0
            grad_shape(1:3,1:3*size(cn):3) = ig%grad_shape(p)%p
            grad_shape(4:6,2:3*size(cn):3) = ig%grad_shape(p)%p
            grad_shape(7:9,3:3*size(cn):3) = ig%grad_shape(p)%p

            tensor_dot = 0
            tensor_dot(1,[1,4,5]) = ig%n(:,p)
            tensor_dot(2,[4,2,6]) = ig%n(:,p)
            tensor_dot(3,[5,6,3]) = ig%n(:,p)

            ! displ_ = reshape(displ(:,cn), [3*size(cn)])
            ! grad_displ = matmul(grad_shape(:,:3*size(cn)), displ_)
            ! strain = matmul(strain_matrix, grad_displ)
            ! stress = matmul(stress_matrix, strain)

            this%displ_to_force(p)%p = matmul(tensor_dot, matmul(stress_matrix, &
                matmul(strain_matrix, grad_shape(:,:3*size(cn)))))

            ! NB: Viscoplasticity is neglected
            ! if (this%model%matl_model%viscoplasticity_enabled) &
            !     call this%compute_viscoplasticity_precon_contribution(dt, n, this%F(:,:,n))
          end do
        end associate
      end do

      do n = 1, mesh%nnode_onP
        n1 = 3*(n-1) + 1
        n2 = 3*(n-1) + 2
        n3 = 3*(n-1) + 3

        ! call this%A%set(n1, n1, 1.0_r8)
        ! call this%A%set(n2, n2, 1.0_r8)
        ! call this%A%set(n3, n3, 1.0_r8)
        ! cycle

        if (this%model%lame1_n(n) < 1e-6_r8 .and. this%model%lame2_n(n) < 1e-6_r8) then
          ! Set displacements to zero for empty or fluid filled cells
          call this%A%set(n1, n1, 1.0_r8 / this%model%scaling_factor(n))
          call this%A%set(n2, n2, 1.0_r8 / this%model%scaling_factor(n))
          call this%A%set(n3, n3, 1.0_r8 / this%model%scaling_factor(n))
        else
          do xp = ig%xnpoint(n), ig%xnpoint(n+1)-1
            k = xp - ig%xnpoint(n) + 1
            p = ig%npoint(xp)
            c = ig%pcell(p)
            s = merge(-1, 1, ig%nppar(k,n))

            associate (cn => mesh%cnode(mesh%xcnode(c):mesh%xcnode(c+1)-1))
              ! displ_ = reshape(displ(:,cn), [3*size(cn)])
              ! lhs(:,n) = lhs(:,n) + s * matmul(displ_to_force(p)%p, displ_)

              do xcn = 1, size(cn)
                nn = cn(xcn) ! displacement node
                nn1 = 3*(nn-1) + 1
                nn2 = 3*(nn-1) + 2
                nn3 = 3*(nn-1) + 3
                if (this%model%lame1_n(nn) < 1e-6_r8 .and. this%model%lame2_n(nn) < 1e-6_r8) cycle

                do i = n1, n3
                  ii = i-n1+1 ! length-3 residual
                  do j = nn1, nn3
                    jj = 3*(xcn-1) + j - nn1 + 1 ! length-3*size(cn) displacement
                    call this%A%add_to(i, j, s * this%displ_to_force(p)%p(ii, jj) &
                        / this%model%scaling_factor(n))
                  end do
                end do
              end do
            end associate
          end do
        end if
      end do
    end associate

    call this%model%bc%compute_deriv_full(t, this%model%scaling_factor, this%A)
    call this%precon%compute

    call stop_timer("precon-compute")

  end subroutine compute


  subroutine apply(this, u, f)
    class(sm_hypre_precon), intent(in), target :: this
    real(r8), intent(in), contiguous, target :: u(:,:) ! current displacement guess
    real(r8), intent(inout), contiguous, target :: f(:,:) ! in residual, out next displacement guess
    real(r8), pointer :: f_(:) => null()
    call start_timer("precon-apply")
    f_(1:3*size(f,dim=2)) => f
    call this%precon%apply(f_)
    call this%model%mesh%node_imap%gather_offp(f)
    call stop_timer("precon-apply")
  end subroutine apply

  ! TODO: Add the viscoplasticity contribution to the preconditioner.
  ! subroutine compute_viscoplasticity_precon_contribution(this, dt, vof, temperature, displ, F)

  !   class(sm_hypre_precon), intent(inout) :: this
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


  logical function recompute_needed(this)

    use parallel_communication, only: global_any

    class(sm_hypre_precon), intent(inout) :: this

    ! TODO: Figure out good parameters & expose to user input
    real(r8), parameter :: atol = 1d-3
    real(r8), parameter :: rtol = 1d-3
    real(r8) :: err1, err2
    integer :: c

    if (this%first) then
      recompute_needed = .true.
      this%first = .false.
    else
      recompute_needed = .false.
      do c = 1, this%model%mesh%ncell_onP
        err1 = abs(this%l1(c) - this%model%lame1(c)) / (atol + rtol*this%model%lame1(c))
        err2 = abs(this%l2(c) - this%model%lame2(c)) / (atol + rtol*this%model%lame2(c))
        if (err1 > 1 .or. err2 > 1) then
          recompute_needed = .true.
          exit
        end if
      end do
      recompute_needed = global_any(recompute_needed)
    end if

    if (recompute_needed) then
      this%l1 = this%model%lame1
      this%l2 = this%model%lame2
    end if

  end function recompute_needed

end module sm_hypre_precon_type
