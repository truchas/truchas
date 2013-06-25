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
  use kind_module,        only: int_kind, log_kind, real_kind

  implicit none

  public :: Region_type

  type REGION_type
    real(real_kind), dimension(3) :: x1
    real(real_kind), dimension(3) :: x2
    logical(log_kind)             :: flow_off
  end type REGION_type

  ! Private Module
  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    real(real_kind),   public   :: x1,y1,z1,x2,y2,z2
    logical(log_kind), public   :: flow_off
    integer(int_kind), public   :: nregion

    type(REGION_type), public, dimension(50) :: Regions

    

END MODULE REGION_DATA
