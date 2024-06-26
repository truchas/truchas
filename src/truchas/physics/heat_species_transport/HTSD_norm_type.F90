!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module HTSD_norm_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use HTSD_model_type
  use rad_problem_type
  use parallel_communication
  implicit none
  private

  type, public :: HTSD_norm
    private
    type(HTSD_model), pointer :: model
    logical :: verbose
    integer :: unit
    !! Heat transfer tolerances
    real(r8) :: abs_T_tol
    real(r8) :: rel_T_tol
    real(r8) :: abs_H_tol
    real(r8) :: rel_H_tol
    !! Species diffusion tolerances
    real(r8), allocatable :: abs_C_tol(:)
    real(r8), allocatable :: rel_C_tol(:)
  end type HTSD_norm

  public :: HTSD_norm_init
  public :: HTSD_norm_compute

contains

  subroutine HTSD_norm_init (this, model, params)

    use parameter_list_type
    use parallel_communication, only: is_IOP

    type(HTSD_norm), intent(out) :: this
    type(HTSD_model), intent(in), target :: model
    type(parameter_list), intent(inout) :: params
    integer :: stat
    character(:), allocatable :: errmsg
    this%model => model
    if (associated(model%ht)) then
      call params%get('abs-t-tol', this%abs_T_tol, stat, errmsg)
      if (stat /= 0) print *, params%path() // ': ' // errmsg
      call params%get('rel-t-tol', this%rel_T_tol)
      INSIST(valid_tol(this%abs_T_tol, this%rel_T_tol))
      call params%get('abs-h-tol', this%abs_h_tol)
      call params%get('rel-h-tol', this%rel_h_tol)
      INSIST(valid_tol(this%abs_H_tol, this%rel_H_tol))
    end if
    if (associated(model%sd)) then
      call params%get('abs-c-tol', this%abs_c_tol)
      call params%get('rel-c-tol', this%rel_c_tol)
      INSIST(size(this%abs_C_tol)==model%num_comp)
      INSIST(size(this%rel_C_tol)==model%num_comp)
      INSIST(all(valid_tol(this%abs_C_tol, this%rel_C_tol)))
    end if
    call params%get('verbose', this%verbose)
    this%verbose = (is_IOP .and. this%verbose)
    call params%get('unit', this%unit, default=-1)
  contains
    elemental logical function valid_tol (abs_tol, rel_tol)
      real(r8), intent(in) :: abs_tol, rel_tol
      valid_tol = (abs_tol >= 0.0_r8) .and. (rel_tol >= 0.0_r8) .and. &
                  ((abs_tol > 0.0_r8 .or. rel_tol > 0.0_r8))
    end function valid_tol
  end subroutine HTSD_norm_init

  subroutine HTSD_norm_compute (this, t, u, du, du_norm)

#ifdef G95_COMPILER_WORKAROUND
    type(HTSD_norm), intent(inout) :: this
#else
    type(HTSD_norm), intent(in) :: this
#endif
    real(r8), intent(in) :: t
    real(r8), intent(in), target :: u(:), du(:)
    real(r8), intent(out) :: du_norm

    integer :: n
    real(r8) :: ht_du_norm, sd_du_norm, qerror
    real(r8), pointer :: useg(:), duseg(:), qrad(:), temp(:)
    integer, pointer :: faces(:)
    real(r8), allocatable :: res(:), rhs(:)

    ASSERT(size(u) == size(du))
    ASSERT(size(u) == HTSD_model_size(this%model))

    if (associated(this%model%ht)) then
      ht_du_norm = 0.0_r8
      !! Cell temperature delta norm.
      call HTSD_model_get_cell_temp_view (this%model,  u,  useg)
      call HTSD_model_get_cell_temp_view (this%model, du, duseg)
      ht_du_norm = max(ht_du_norm, maxerr(useg, duseg, this%abs_T_tol, this%rel_T_tol, this%model%void_cell))
      !! Face temperature delta norm.
      call HTSD_model_get_face_temp_view (this%model,  u,  useg)
      call HTSD_model_get_face_temp_view (this%model, du, duseg)
      ht_du_norm = max(ht_du_norm, maxerr(useg, duseg, this%abs_T_tol, this%rel_T_tol, this%model%void_face))
      !! Cell enthalpy delta norm.
      call HTSD_model_get_cell_heat_view (this%model,  u,  useg)
      call HTSD_model_get_cell_heat_view (this%model, du, duseg)
      ht_du_norm = max(ht_du_norm, maxerr(useg, duseg, this%abs_H_tol, this%rel_H_tol, this%model%void_cell))
      ht_du_norm = global_maxval(ht_du_norm)
      !! Enclosure radiation system error norms
      !TODO! This is a quick hack that needs to be fixed.  The tolerance
      !TODO! is hardwired, and  we are (re)computing the residual.  The BDF2
      !TODO! integrator thinks it's computing the norm of the correction du but
      !TODO! in this case it is the actual residual norm.  We need more general
      !TODO! norms in the integrator.
      if (associated(this%model%ht%vf_rad_prob)) then
        !if (is_IOP) write(*,'(a)',advance='no')'ER error: ||res||/||rhs||='
        do n = 1, size(this%model%ht%vf_rad_prob)
          faces => this%model%ht%vf_rad_prob(n)%faces
          allocate(res(size(faces)), rhs(size(faces)))
          call HTSD_model_get_radiosity_view (this%model, n, u, qrad)
          call HTSD_model_get_face_temp_view (this%model, u, temp)
          call this%model%ht%vf_rad_prob(n)%residual (t, qrad, temp(faces), res)
          call this%model%ht%vf_rad_prob(n)%rhs (t, temp(faces), rhs)
          qerror = sqrt(global_sum(res**2)) / sqrt(global_sum(rhs**2))
          !if (is_IOP) write(*,'(e10.3)',advance='no') qerror
          ht_du_norm = max(ht_du_norm, qerror/1.0d-3)
          deallocate(res, rhs)
        end do
        !if (is_IOP) write(*,*)
      end if
    else
      ht_du_norm = 0.0_r8
    end if

    if (associated(this%model%sd)) then
      sd_du_norm = 0.0_r8
      do n = 1, this%model%num_comp
        !! Cell concentration delta norms.
        call HTSD_model_get_cell_conc_view (this%model, n,  u,  useg)
        call HTSD_model_get_cell_conc_view (this%model, n, du, duseg)
        sd_du_norm = max(sd_du_norm, maxerr(useg, duseg, this%abs_C_tol(n), this%rel_C_tol(n), this%model%void_cell))
        !! Face concentration delta norms.
        call HTSD_model_get_face_conc_view (this%model, n,  u,  useg)
        call HTSD_model_get_face_conc_view (this%model, n, du, duseg)
        sd_du_norm = max(sd_du_norm, maxerr(useg, duseg, this%abs_C_tol(n), this%rel_C_tol(n), this%model%void_face))
        sd_du_norm = global_maxval(sd_du_norm)
      end do
    else
      sd_du_norm = 0.0_r8
    end if

    du_norm = max(ht_du_norm, sd_du_norm)

  contains

    real(r8) function maxerr (u, du, atol, rtol, void)
      real(r8), intent(in) :: u(:), du(:), atol, rtol
      logical, pointer :: void(:)
      real(r8) :: array(size(du))
      if (associated(void)) then
        where (void(:size(du)))
          array = 0.0_r8
        elsewhere
          array = abs(du) / (atol + rtol*abs(u))
        end where
      else
        array = abs(du) / (atol + rtol*abs(u))
      end if
      maxerr = maxval(array)
    end function maxerr

  end subroutine HTSD_norm_compute

end module HTSD_norm_type
