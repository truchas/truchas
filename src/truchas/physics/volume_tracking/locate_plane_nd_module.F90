!!
!! LOCATE_PLANE_ND_MODULE
!!
!! This module provides a plane reconstruction
!! subroutine for the nested dissection method.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! March 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module locate_plane_nd_module

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use plane_type
  use polyhedron_type
  use brent_root_class
  use truchas_logging_services
  implicit none
  private

  public :: locate_plane_nd

  ! define type for error function
  type, extends(brent_root) :: vof_error_func
    real(r8)         :: tvol,norm(3),parvol
    type(polyhedron) :: poly
  contains
    procedure :: init => func_init
    procedure :: eval => func_eval
    procedure :: signed_eval => func_signed_eval
    procedure :: f => func_signed_eval
    !procedure :: f => func_eval
  end type vof_error_func

contains

  ! given a polyhedron, a normal, and a volume, calculate a plane constant
  ! and return a plane type
  ! Alternatively, should we return the two polyhedra? If they are already
  ! calculated, it would save resplitting the input polyhedron.

  type(plane) function locate_plane_nd (poly, norm, vol, cell_volume, cutvof, max_iterations)

    use polyhedron_type
    use plane_type

    type(polyhedron), intent(inout) :: poly
    real(r8), intent(in) :: norm(:), vol, cell_volume, cutvof
    integer, intent(in) :: max_iterations

    real(r8)             :: rho_min,rho_max
    integer              :: ierr
    type(vof_error_func) :: vof_error

    ASSERT(size(norm)==3)

    ! initialize error function
    call vof_error%init (norm, poly, vol, cell_volume)

    ! get bounds for Brent's method
    call rho_bracket (rho_min, rho_max, norm, poly, vof_error)

    ! start Brent's method
    locate_plane_nd%normal = norm
    vof_error%feps = 0.5_r8*cutvof; vof_error%maxitr = max_iterations
    call vof_error%find_root (rho_min, rho_max, locate_plane_nd%rho, ierr)
    ! note ~30 iterations seem to be necessary to pass current unit tests

  end function locate_plane_nd

  ! loop through all the polyhedron vertices, checking if a plane
  ! intersecting each one is above or below the target volume,
  ! thereby bracketing the allowed range of plane constants

  subroutine rho_bracket(rho_min, rho_max, norm, poly, volume_error)

    use near_zero_function

    real(r8), intent(out) :: rho_min, rho_max
    real(r8), intent(in) :: norm(:)
    type(polyhedron), intent(inout) :: poly
    type(vof_error_func), intent(inout) :: volume_error

    integer, parameter :: iter_max = 10
    real(r8) :: err_min,err_max,err, rho
    integer :: i, ierr, c
    type(plane) :: P

    err_min = -huge(1.0_r8); rho_min = -huge(1.0_r8)
    err_max =  huge(1.0_r8); rho_max =  huge(1.0_r8)

    ! find the outer bounds
    if (poly%nchildren == 0) then
      do i = 1,poly%parent%nVerts
        rho = dot_product(poly%parent%x(:,i),norm)

        if (rho <= rho_min .or. rho >= rho_max) cycle

        err = volume_error%signed_eval (rho)

        if (near_zero(err)) then
          err_min = err
          err_max = err

          rho_min = rho
          rho_max = rho
          return
        else if (0 < err .and. err < err_max) then
          err_max = err
          rho_max = rho
        else if (err_min < err .and. err < 0) then
          err_min = err
          rho_min = rho
        end if
      end do
    else
      do c = 1,poly%nchildren
        do i = 1,poly%convex_polyhedron(c)%nVerts
          rho = dot_product(poly%convex_polyhedron(c)%x(:,i),norm)

          if (rho <= rho_min .or. rho >= rho_max) cycle

          err = volume_error%signed_eval (rho)

          if (near_zero(err)) then
            err_min = err
            err_max = err

            rho_min = rho
            rho_max = rho
            return
          else if (0 < err .and. err < err_max) then
            err_max = err
            rho_max = rho
          else if (err_min < err .and. err < 0) then
            err_min = err
            rho_min = rho
          end if
        end do
      end do
    end if

    ! make sure the bounds were set
    if (rho_min==-huge(1.0_r8) .or. rho_max==huge(1.0_r8)) then
      P%normal = norm
      print '(a,3f10.3)', 'normal: ', norm
      call poly%print_data (normalized=.true.)
      !write(*,*) 'volume ',poly%volume ()
      do i = 1,poly%parent%nVerts
        rho = dot_product(poly%parent%x(:,i),norm)
        err = volume_error%signed_eval(rho)
        P%rho = rho
        print '(a,i3,3es14.4)', 'vertex, rho, error, vol: ', i, rho, err, &
            poly%volume_behind_plane (P,ierr)
      end do
      print '(a,2es14.4)', 'rho min,max: ', rho_min, rho_max
      call tls_fatal ('rho bounds not set')
    end if

  end subroutine rho_bracket

  ! error function procedures

  subroutine func_init (this, norm,poly,tvol,parvol)
    use polyhedron_type

    class(vof_error_func), intent(out) :: this
    real(r8),                 intent(in)  :: norm(:), tvol, parvol
    type(polyhedron),         intent(in)  :: poly

    this%norm   = norm
    this%poly   = poly
    this%tvol   = tvol
    this%parvol = parvol

  end subroutine func_init

  real(r8) function func_eval (this, x)
    class(vof_error_func), intent(inout) :: this
    real(r8), intent(in) :: x
    func_eval = abs(this%signed_eval (x))
  end function func_eval

  real(r8) function func_signed_eval (this, x)

    use plane_type

    class(vof_error_func), intent(inout) :: this
    real(r8),                 intent(in) :: x

    type(plane) :: P
    integer :: ierr

    P%rho = x; P%normal = this%norm
    func_signed_eval = (this%poly%volume_behind_plane (P,ierr) - this%tvol) / this%tvol

    if (ierr /= 0) call TLS_fatal ("func_signed_eval failed")

  end function func_signed_eval

end module locate_plane_nd_module
