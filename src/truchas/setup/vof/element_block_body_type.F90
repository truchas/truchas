!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module element_block_body_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use body_class
  use unstr_mesh_type
  implicit none
  private

  type, extends(body), public :: element_block_body
    private
    type(unstr_mesh), pointer :: mesh => null() ! do not own
    integer :: cblockid
  contains
    procedure :: eval
    procedure :: signed_distance
  end type element_block_body

  interface element_block_body
    procedure element_block_body_value
  end interface element_block_body

contains

  !! constructor for ELEMENT_BLOCK_BODY objects
  function element_block_body_value(mesh, cblockid) result(r)
    type(unstr_mesh), target, intent(in) :: mesh
    integer, intent(in) :: cblockid
    type(element_block_body) :: r
    r%mesh => mesh
    r%cblockid = cblockid
  end function element_block_body_value

  logical function eval(this, x, cellid)
    use bitfield_type, only: popcnt, trailz
    class(element_block_body), intent(in) :: this
    real(r8), intent(in) :: x(:)
    integer, intent(in) :: cellid
    INSIST(popcnt(this%mesh%cell_set_mask(cellid)) == 1)
    eval = this%cblockid == this%mesh%cell_set_id(trailz(this%mesh%cell_set_mask(cellid)))
  end function eval

  real(r8) function signed_distance(this, x)
    class(element_block_body), intent(in) :: this
    real(r8),             intent(in) :: x(:)
    signed_distance = 0
  end function signed_distance

end module element_block_body_type
