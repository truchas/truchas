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
  use truchas_timers
  implicit none
  private

  !! TODO extend a generic sm precon class
  type, public :: sm_ds_precon
    private
    type(sm_model), pointer, public :: model => null() ! unowned reference

    real(r8), allocatable :: diag(:)
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

    character(:), allocatable :: context, errmsg
    integer :: stat

    this%model => model
    allocate(this%diag(model%size()))

    context = 'processing ' // params%name() // ': '

    call params%get('precon-num-iter', this%niter, default=1, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if

    call params%get('precon-relaxation-parameter', this%omega, default=1.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) then
      errmsg = context//errmsg
      return
    end if

  end subroutine init


  subroutine compute(this, t, dt, displ)
    class(sm_ds_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt, displ(:,:)
    call start_timer("precon-compute")
    INSIST(.false.) ! TODO -- update diag(:) from the model
    call stop_timer("precon-compute")
  end subroutine compute


  subroutine apply(this, u, f)

    class(sm_ds_precon), intent(in) :: this
    real(r8), intent(in), contiguous :: u(:,:) ! current displacement guess
    real(r8), intent(inout), contiguous :: f(:) ! in residual, out next displacement guess

    integer :: i, j
    real(r8) :: d

    call start_timer("precon-apply")

    do j = 1, size(this%diag)
      d = f(j) / this%diag(j)
      do i = 1, this%niter
        f(j) = f(j) + this%omega*(d-f(j))
      end do
    end do

    call stop_timer("precon-apply")

  end subroutine apply

end module sm_ds_precon_type
