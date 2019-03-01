!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE REGION_DATA
  !=======================================================================
  ! Purpose(s):
  !
  !  Define advection module variables.
  !
  ! Contains: None
  !
  ! Author(s):  Ed Dendy     (dendy@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  implicit none
  private

  type, public :: REGION_type
    real(r8) :: x1(3), x2(3)
    logical :: flow_off
  end type REGION_type

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    real(r8), public :: x1,y1,z1,x2,y2,z2
    logical,  public :: flow_off
    integer,  public :: nregion

    type(REGION_type), public, dimension(50) :: Regions

END MODULE REGION_DATA
