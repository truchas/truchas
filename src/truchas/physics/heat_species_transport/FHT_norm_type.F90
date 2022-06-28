!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module FHT_norm_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use FHT_model_type
  use rad_problem_type
  use parallel_communication
  implicit none
  private

  type, public :: FHT_norm
    private
    type(FHT_model), pointer :: model
    real(r8) :: abs_tol
    real(r8) :: rel_tol
    real(r8), allocatable :: rad_tol(:)
    real(r8) :: err1_0, err2_0, err0
    logical  :: verbose
    integer  :: unit
  end type FHT_norm

  public :: FHT_norm_init, FHT_norm_fnorm

  interface FHT_norm_fnorm
    module procedure FHT_norm_f2norm
    !module procedure FHT_norm_fmaxnorm
    !module procedure FHT_norm_norm_old
  end interface

contains

  subroutine FHT_norm_init (this, model, params)

    use parameter_list_type

    type(FHT_norm), intent(out) :: this
    type(FHT_model), intent(in), target :: model
    type(parameter_list) :: params
    integer :: n
    this%model => model

    !TODO: add error checking to parameter list calls
    call params%get('abs-tol', this%abs_tol)
    call params%get('rel-tol', this%rel_tol)
    INSIST(valid_tol(this%abs_tol, this%rel_tol))
    call params%get('verbose', this%verbose)
    this%verbose = (is_IOP .and. this%verbose)
    call params%get('unit', this%unit, default=-1)

    !! Initialize the heat equation/view factor radiation
    if (associated(model%vf_rad_prob)) then
      n = size(model%vf_rad_prob)
      call params%get('rad-tol', this%rad_tol)
      INSIST(size(this%rad_tol) == n)
      INSIST(all(this%rad_tol >= 0.0_r8))
    end if
  contains
    logical function valid_tol (abs_tol, rel_tol)
      real(r8), intent(in) :: abs_tol, rel_tol
      valid_tol = (abs_tol >= 0.0_r8) .and. (rel_tol >= 0.0_r8) .and. &
                  ((abs_tol > 0.0_r8 .or. rel_tol > 0.0_r8))
    end function
  end subroutine FHT_norm_init

  subroutine FHT_norm_f2norm (this, t, u, hdot, f, error)

    type(FHT_norm), intent(inout) :: this
    real(r8), intent(in) :: t, u(:), hdot(:), f(:)
    real(r8), intent(out), optional :: error
    target :: u, f

    integer :: i, num1, num2
    real(r8), pointer :: fseg(:), qres(:), temp(:)
    integer, pointer :: faces(:)
    real(r8), allocatable :: rhs(:)
    real(r8) :: sum1, sum2, err

    ASSERT(size(u) == FHT_model_size(this%model))
    ASSERT(size(f) == FHT_model_size(this%model))

    call FHT_model_get_cell_temp_view (this%model, f, fseg)
    sum1 = global_sum(fseg**2, mask=.not.this%model%void_cell(:this%model%mesh%ncell_onP))
    num1 = global_count(.not.this%model%void_cell(:this%model%mesh%ncell_onP))
    call FHT_model_get_face_temp_view (this%model, f, fseg)
    sum2 = global_sum(fseg**2, mask=.not.this%model%void_face(:this%model%mesh%nface_onP))
    num2 = global_count(.not.this%model%void_face(:this%model%mesh%nface_onP))

    err = sqrt((sum1+sum2)/(num1+num2))

    if (present(error)) then

      error = err / (this%abs_tol + this%rel_tol * this%err0)
      if (this%verbose) then
        write(this%unit,'(2(a,es10.3))') '  HC error: ||F||_2 =', err, ', scaled =', error
      end if

      !! Enclosure radiation system error norms
      if (associated(this%model%vf_rad_prob)) then
        if (this%verbose) then
          write(this%unit,'(a)',advance='no')'  ER error: scaled ||res||_2/||rhs||_2 ='
        end if
        do i = 1, size(this%model%vf_rad_prob)
          faces => this%model%vf_rad_prob(i)%faces
          allocate(rhs(size(faces)))
          call FHT_model_get_radiosity_view (this%model, i, f, qres)
          call FHT_model_get_face_temp_view (this%model, u, temp)
          call this%model%vf_rad_prob(i)%rhs (t, temp(faces), rhs)
          err = sqrt(global_sum(qres**2)) / (this%rad_tol(i) * sqrt(global_sum(rhs**2)))
          if (this%verbose) then
            write(this%unit,'(es10.3)',advance='no') err
          end if
          error = max(error, err)
          deallocate(rhs)
        end do
        if (this%verbose) write(this%unit,*)
      end if

    else

      this%err0 = err
      if (this%verbose) then
      	write(this%unit,'(a)') 'HEAT TRANSFER FUNCTION NORMS'
        write(this%unit,'(a,es10.3,a)') '  HC error: ||F||_2 =', err, ' (initial)'
      end if

    end if

  end subroutine FHT_norm_f2norm


  subroutine FHT_norm_fmaxnorm (this, t, u, hdot, f, error)

    type(FHT_norm), intent(inout) :: this
    real(r8), intent(in) :: t, u(:), hdot(:), f(:)
    real(r8), intent(out), optional :: error
    target :: u, f

    integer :: i
    real(r8), pointer :: fseg(:), qres(:), temp(:)
    integer, pointer :: faces(:)
    real(r8), allocatable :: rhs(:)
    real(r8) :: err1, err2, err

    ASSERT(size(u) == FHT_model_size(this%model))
    ASSERT(size(f) == FHT_model_size(this%model))

    call FHT_model_get_cell_temp_view (this%model, f, fseg)
    err1 = global_maxval(abs(fseg), mask=.not.this%model%void_cell(:this%model%mesh%ncell_onP))
    call FHT_model_get_face_temp_view (this%model, f, fseg)
    err2 = global_maxval(abs(fseg), mask=.not.this%model%void_face(:this%model%mesh%nface_onP))
    err = max(err1,err2)

    if (present(error)) then

      error = err / (this%abs_tol + this%rel_tol * this%err0)
      if (this%verbose) then
        write(this%unit,'(2(a,es10.3))') '  HC error: ||F||_max =', err, ', scaled =', error
      end if

      !! Enclosure radiation system error norms
      !TODO! Fix the hardwired tolerance.
      if (associated(this%model%vf_rad_prob)) then
        if (this%verbose) then
          write(this%unit,'(a)',advance='no')'  ER error: scaled ||res||_2/||rhs||_2 ='
        end if
        do i = 1, size(this%model%vf_rad_prob)
          faces => this%model%vf_rad_prob(i)%faces
          allocate(rhs(size(faces)))
          call FHT_model_get_radiosity_view (this%model, i, f, qres)
          call FHT_model_get_face_temp_view (this%model, u, temp)
          call this%model%vf_rad_prob(i)%rhs (t, temp(faces), rhs)
          err = sqrt(global_sum(qres**2)) / (this%rad_tol(i) * sqrt(global_sum(rhs**2)))
          if (this%verbose) then
            write(this%unit,'(es10.3)',advance='no') err
          end if
          error = max(error, err)
          deallocate(rhs)
        end do
        if (this%verbose) write(this%unit,*)
      end if

    else

      this%err0 = err
      if (this%verbose) then
      	write(this%unit,'(a)') 'HEAT TRANSFER FUNCTION NORMS'
        write(this%unit,'(a,es10.3,a)') '  HC error: ||F||_max =', err, ' (initial)'
      end if

    end if

  end subroutine FHT_norm_fmaxnorm


  subroutine FHT_norm_fnorm_old (this, t, u, hdot, f, error)

    type(FHT_norm), intent(inout) :: this
    real(r8), intent(in) :: t, u(:), hdot(:), f(:)
    real(r8), intent(out), optional :: error
    target :: u, f

    integer :: i, loc1, loc2, pid1, pid2
    real(r8), pointer :: fseg(:), qres(:), temp(:)
    integer, pointer :: faces(:)
    real(r8), allocatable :: rhs(:)
    real(r8) :: qerror, err1, err2

    ASSERT(size(u) == FHT_model_size(this%model))
    ASSERT(size(f) == FHT_model_size(this%model))

    !! Heat equation residual norm.
    call FHT_model_get_cell_temp_view (this%model, f, fseg)
    call global_maxloc_sub (abs(fseg), pid1, loc1, mask=.not.this%model%void_cell(:this%model%mesh%ncell_onP))
    err1 = global_maxval(abs(fseg), mask=.not.this%model%void_cell(:this%model%mesh%ncell_onP))

    !! Flux matching residual norm.
    call FHT_model_get_face_temp_view (this%model, f, fseg)
    call global_maxloc_sub (abs(fseg), pid2, loc2, mask=.not.this%model%void_face(:this%model%mesh%nface_onP))
    err2 = global_maxval(abs(fseg), mask=.not.this%model%void_face(:this%model%mesh%nface_onP))

    if (present(error)) then

      err1 = err1 / (this%abs_tol + this%rel_tol * this%err1_0)
      err2 = err2 / (this%abs_tol + this%rel_tol * this%err2_0)
      error = max(err1, err2)

      if (this%verbose) write(this%unit,'(2es25.3)') err1, err2

!      !! Some extra diagnostics when F1 determines the error.
!      if (err1 == error .and. global_any(this%verbose)) then
!        call get_value (this%model%vfrac, pid1, loc1, bad_vfrac)
!        call FHT_model_get_cell_temp_view (this%model, u, useg)
!        call get_value (useg, pid1, loc1, bad_T)
!        call FHT_model_get_cell_temp_view (this%model, f, fseg)
!        call get_value (fseg, pid1, loc1, bad_F)
!        if (this%verbose) then
!          write(this%unit,'(4x,a,es13.5)') 'F1 cell data: vfrac=', bad_vfrac
!          write(this%unit,'(4x,a,es13.5)') 'F1 cell data: T=', bad_T
!          write(this%unit,'(4x,a,es13.5)') 'F1 cell data: F=', bad_F
!        end if
!      end if

      !! Enclosure radiation system error norms
      if (associated(this%model%vf_rad_prob)) then
        if (is_IOP) write(*,'(a)',advance='no')'ER error: ||res||/||rhs||='
        do i = 1, size(this%model%vf_rad_prob)
          faces => this%model%vf_rad_prob(i)%faces
          allocate(rhs(size(faces)))
          call FHT_model_get_radiosity_view (this%model, i, f, qres)
          call FHT_model_get_face_temp_view (this%model, u, temp)
          call this%model%vf_rad_prob(i)%rhs (t, temp(faces), rhs)
          qerror = sqrt(global_sum(qres**2)) / (this%rad_tol(i) * sqrt(global_sum(rhs**2)))
          if (is_IOP) write(*,'(e10.3)',advance='no') qerror
          error = max(error, qerror)
          deallocate(rhs)
        end do
        if (is_IOP) write(*,*)
      end if

    else

      this%err1_0 = err1
      this%err2_0 = err2
      if (this%verbose) then
      	write(this%unit,'(a)') 'HEAT TRANSFER FUNCTION NORMS'
	write(this%unit,'(a)') '  F1: heat conservation equation (cells)'
	write(this%unit,'(a)') '  F2: flux continuity equation (faces)'
	write(this%unit,'(2(a,es10.3))') '  ||F1_0||_max =', err1, ', ||F2_0||_max =', err2
	write(this%unit,'(a)') '  ||F1||_max/||F1_0||_max  ||F2||_max/||F2_0||_max'
      end if

    end if

  end subroutine FHT_norm_fnorm_old

end module FHT_norm_type
