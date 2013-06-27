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
  use kinds, only: r8
  implicit none
  private

  character(80),  public, save :: limiter_type
  integer, public, save :: advection_order_vol
  integer, public, save :: advection_order_energy
  integer, public, save :: advection_order_momentum
  integer, public, save :: advection_order_species

  real(r8), public, save, dimension(:,:),   allocatable :: Momentum_Delta
  real(r8), public, save, dimension(:,:,:), allocatable :: Volume_Flux

END MODULE ADVECTION_DATA
