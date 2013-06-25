MODULE ADVECTION_DATA
  !=======================================================================
  ! Purpose(s):
  !
  !    Define advection module variables.
  !
  ! Contains: None
  !
  ! Author(s):  Ed Dendy     (dendy@lanl.gov)
  !
  !=======================================================================
  use kind_module, only: int_kind, real_kind

  implicit none

  private

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  character(LEN=80),  public, save                                :: limiter_type

  integer (int_kind), public, save                                :: advection_order_vol
  integer (int_kind), public, save                                :: advection_order_energy
  integer (int_kind), public, save                                :: advection_order_momentum
  integer (int_kind), public, save                                :: advection_order_species

  real(real_kind),    public, save, dimension(:,:),   allocatable :: Momentum_Delta
  real(real_kind),    public, save, dimension(:,:,:), allocatable :: Volume_Flux

END MODULE ADVECTION_DATA
