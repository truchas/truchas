!! HT_2D_IC_SOLVER_TYPE
!!
!! Constructs a consistent initial state for the 2D heat-transfer DAE.

#include "f90_assert.fpp"

module ht_2d_ic_solver_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use ht_2d_model_type
  use ht_2d_vector_type
  use parameter_list_type
  use parallel_communication, only: global_dot_product
  use truchas_logging_services
  implicit none
  private

  type, public :: ht_2d_ic_solver
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(ht_2d_model), pointer :: model => null()
  contains
    procedure :: init
    procedure :: compute => compute_vector
  end type ht_2d_ic_solver

contains

  subroutine init(this, model)
    class(ht_2d_ic_solver), intent(out) :: this
    type(ht_2d_model), intent(in), target :: model
    this%model => model
    this%mesh => model%mesh
  end subroutine init


  !! Construct the initial state directly in the solver vector representation.
  subroutine compute_vector(this, t, dt, temp, u, udot, rel_tol, max_itr, stat, errmsg)
    class(ht_2d_ic_solver), intent(inout) :: this
    real(r8), intent(in) :: t, dt, temp(:), rel_tol
    type(ht_2d_vector), intent(inout) :: u, udot
    integer, intent(in) :: max_itr
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: state(:,:), Hcell(:)

    stat = 0
    errmsg = ''

    ASSERT(size(temp) == this%mesh%ncell_onP)

    call u%setval(0.0_r8)
    u%tc(:this%mesh%ncell_onP) = temp
    call this%mesh%cell_imap%gather_offp(u%tc)

    allocate(state(this%mesh%ncell_onP,0:0), Hcell(this%mesh%ncell_onP))
    state(:,0) = u%tc(:this%mesh%ncell_onP)
    call this%model%H_of_T%compute_value(state, Hcell)
    u%hc(:this%mesh%ncell_onP) = Hcell

    !! The zero face field is only an initial guess for the algebraic solve.
    call compute_face_temp(this, t, u, rel_tol, max_itr, stat, errmsg)
    if (stat /= 0) return

    call compute_udot(this, t, dt, u, udot, rel_tol, max_itr, stat, errmsg)
  end subroutine compute_vector


  !! Compute the initial time derivative by advancing enthalpy one small step,
  !! solving the algebraic variables at that advanced state, and differencing.
  subroutine compute_udot(this, t, dt, u, udot, rel_tol, max_itr, stat, errmsg)
    class(ht_2d_ic_solver), intent(inout) :: this
    real(r8), intent(in) :: t, dt, rel_tol
    type(ht_2d_vector), intent(inout) :: u, udot
    integer, intent(in) :: max_itr
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(ht_2d_vector) :: f, advanced
    real(r8) :: Tmin, Tmax
    integer :: j

    stat = 0
    errmsg = ''

    call udot%setval(0.0_r8)
    call f%init(u)
    call f%setval(0.0_r8)
    call this%model%residual(t, u, udot, f)
    udot%hc(:this%mesh%ncell_onP) = -f%tc(:this%mesh%ncell_onP) / &
        this%mesh%volume(:this%mesh%ncell_onP)

    call advanced%init(u)
    call advanced%copy(u)
    advanced%hc(:this%mesh%ncell_onP) = u%hc(:this%mesh%ncell_onP) + &
        dt * udot%hc(:this%mesh%ncell_onP)

    do j = 1, this%mesh%ncell_onP
      Tmin = u%tc(j) - 1.0_r8
      Tmax = u%tc(j) + 1.0_r8
      call this%model%T_of_H%compute(j, advanced%hc(j), Tmin, Tmax, advanced%tc(j))
    end do
    call this%mesh%cell_imap%gather_offp(advanced%tc)

    call compute_face_temp(this, t+dt, advanced, rel_tol, max_itr, stat, errmsg)
    if (stat /= 0) return

    udot%tc(:this%mesh%ncell_onP) = (advanced%tc(:this%mesh%ncell_onP) - &
        u%tc(:this%mesh%ncell_onP)) / dt
    udot%tf(:this%mesh%nface_onP) = (advanced%tf(:this%mesh%nface_onP) - &
        u%tf(:this%mesh%nface_onP)) / dt
    call udot%gather_offp()
  end subroutine compute_udot


  !! Solve the face-temperature algebraic equations using the A22 diffusion
  !! block.  This belongs to the initial-condition solver rather than the
  !! time-integration model because it constructs a consistent initial state.
  subroutine compute_face_temp(this, t, u, rel_tol, max_itr, stat, errmsg)
    use hypre_hybrid_type
    use mfd_2d_diff_matrix_type

    class(ht_2d_ic_solver), intent(inout) :: this
    real(r8), intent(in) :: t, rel_tol
    type(ht_2d_vector), intent(inout) :: u
    integer, intent(in) :: max_itr
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(mfd_2d_diff_matrix), target :: dm
    type(parameter_list), target :: params
    type(hypre_hybrid) :: solver
    type(ht_2d_vector) :: udot, f
    real(r8), allocatable :: coef(:), state(:,:), z(:)
    real(r8) :: norm
    integer :: num_itr, num_dscg_itr, num_pcg_itr
    character(80) :: msg

    stat = 0
    errmsg = ''

    call udot%init(u)
    call udot%setval(0.0_r8)
    call f%init(u)
    call f%setval(0.0_r8)
    call this%model%residual(t, u, udot, f)
    norm = sqrt(global_dot_product(f%tf(:this%mesh%nface_onP), f%tf(:this%mesh%nface_onP)))

    if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
      write (msg,'(a,es10.3)') 'ht_2d_ic_solver%compute_face_temp: initial ||rface||_2 = ', norm
      call TLS_info(trim(msg))
    end if
    if (norm == 0.0_r8) return

    allocate(coef(this%mesh%ncell), state(this%mesh%ncell_onP,0:0))
    state(:,0) = u%tc(:this%mesh%ncell_onP)
    call this%model%conductivity%compute_value(state, coef(:this%mesh%ncell_onP))
    call this%mesh%cell_imap%gather_offp(coef)
    call dm%init(this%model%disc)
    call dm%compute(coef)
    if (allocated(this%model%bc_dir)) then
      call this%model%bc_dir%compute(t)
      call dm%set_dir_faces(this%model%bc_dir%index)
    end if

    call params%set('krylov-method', 'cg')
    call params%set('max-ds-iter', max_itr)
    call params%set('max-amg-iter', max_itr)
    call params%set('rel-tol', rel_tol)
    if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
      call params%set('print-level', 1)
      call params%set('logging-level', 1)
    else
      call params%set('print-level', 0)
      call params%set('logging-level', 0)
    end if
    call solver%init(dm%a22, params)
    call solver%setup()

    allocate(z(this%mesh%nface_onP))
    z = 0.0_r8
    call solver%solve(f%tf(:this%mesh%nface_onP), z, stat)
    u%tf(:this%mesh%nface_onP) = u%tf(:this%mesh%nface_onP) - z

    if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
      call solver%get_metrics(num_itr, num_dscg_itr, num_pcg_itr, norm)
      write(msg,'(3(a,i0),a,es9.2)') 'solve: num_itr = ', num_itr, &
          ' (', num_dscg_itr, ', ', num_pcg_itr, '), ||r||/||b|| = ', norm
      call TLS_info(trim(msg))
    end if

    if (stat /= 0) then
      call solver%get_metrics(num_itr, num_dscg_itr, num_pcg_itr, norm)
      write(msg,'(3(a,i0),a,es9.2)') 'failed to converge: num_itr = ', num_itr, &
          ' (', num_dscg_itr, ', ', num_pcg_itr, '), ||r||/||b|| = ', norm
      errmsg = trim(msg)
    end if
  end subroutine compute_face_temp

end module ht_2d_ic_solver_type
