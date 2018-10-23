!!
!! Provides the routine to locate a material interface
!! in a cell using the onion-skin infrastructure
!!
!! note: lots of duplication here from the nested dissection locate plane module
!!

#include "f90_assert.fpp"

module locate_plane_os_function

  use kinds, only: r8
  use plane_type
  use truncation_volume_type
  use brent_root_class
  use truchas_logging_services
  implicit none
  private

  public :: locate_plane_os

  type, extends(brent_root) :: vof_error_func
    real(r8) :: target_volume, cell_volume
    type(truncation_volume) :: trunc_vol
  contains
    procedure :: init
    procedure :: f => signed_eval
  end type vof_error_func

contains

  type(plane) function locate_plane_os(norm, vof, volume, node, cutoff, maxiter)

    real(r8), intent(in) :: norm(:), vof, volume, node(:,:), cutoff
    integer, intent(in) :: maxiter

    type(vof_error_func) :: vof_error
    real(r8) :: rho_min, rho_max
    integer :: ierr

    call vof_error%init(node, norm, vof, volume)
    call rho_bracket(rho_min, rho_max, norm, node, vof_error)

    locate_plane_os%normal = norm
    vof_error%feps = 0.5_r8 * cutoff; vof_error%maxitr = maxiter
    call vof_error%find_root(rho_min, rho_max, locate_plane_os%rho, ierr)

    ASSERT(ierr == 0)

  end function locate_plane_os

  subroutine init(this, nodex, normal, vof, cell_volume)
    class(vof_error_func), intent(out) :: this
    real(r8), intent(in) :: nodex(:,:), normal(:), vof, cell_volume
    call this%trunc_vol%init(nodex, normal)
    this%target_volume = vof * cell_volume
    this%cell_volume = cell_volume
  end subroutine init

  real(r8) function signed_eval(this, x)
    class(vof_error_func), intent(inout) :: this
    real(r8), intent(in) :: x
    signed_eval = (this%trunc_vol%volume(x) - this%target_volume) / this%cell_volume
  end function signed_eval

  ! loop through all the cell nodes, checking if a plane
  ! intersecting each one is above or below the target volume,
  ! bracketing the allowed range of plane constants
  subroutine rho_bracket(rho_min, rho_max, norm, node, vof_error)

    real(r8), intent(out) :: rho_min, rho_max
    real(r8), intent(in) :: norm(:), node(:,:)
    type(vof_error_func), intent(inout) :: vof_error

    real(r8) :: rho, err, err_min, err_max
    integer :: i

    err_min = -huge(1.0_r8); rho_min = -huge(1.0_r8)
    err_max =  huge(1.0_r8); rho_max =  huge(1.0_r8)

    ! find the outer bounds
    do i = 1,size(node, dim=2)
      rho = dot_product(node(:,i), norm)
      if (rho <= rho_min .or. rho >= rho_max) cycle

      err = vof_error%f(rho)

      if (0 < err .and. err < err_max) then
        err_max = err
        rho_max = rho
      else if (err_min < err .and. err < 0) then
        err_min = err
        rho_min = rho
      end if
    end do

    ! ensure the bounds were set
    INSIST(rho_min /= -huge(1.0_r8) .and. rho_max /= huge(1.0_r8))

  end subroutine rho_bracket

end module locate_plane_os_function
