! Switch the comment to use the system random number generator
!#define SYSTEM_RANDOM_GEN
#define LOCAL_RANDOM_GEN
MODULE RANDOM_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Create a machine independent random number generation option
  !
  ! Public Interface(s):
  !
  !   * call INITIALIZE_RANDOM()
  !   * call GENERATE_RANDOM(Rlow, RHigh, N, Rand)
  !
  ! Contains: INITIALIZE_RANDOM, GENERATE_RANDOM
  !
  ! Author(s): Jim Sicilian, CCS-2, sicilian@lanl.gov
  !
  !=======================================================================
  use kind_module,    only: int_kind

  implicit none
  ! Private Module
  private
  ! Private Data
  integer(KIND = int_kind), save :: local_random_seed

  ! Public Procedures and Data
  public :: GENERATE_RANDOM, INITIALIZE_RANDOM
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
CONTAINS
    ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>
    SUBROUTINE INITIALIZE_RANDOM()
#ifdef SYSTEM_RANDOM_GEN
    use kind_module,    only: int_kind

    ! Local Variables
    integer(KIND = int_kind) :: kr, seedsize

    ! initialize the random number sequence using the FORTRAN 90 syntax
         call RANDOM_SEED (SIZE = seedsize)
         call RANDOM_SEED (PUT = (/(kr, kr = 1,seedsize)/))
#endif
#ifdef LOCAL_RANDOM_GEN
         local_random_seed = 61
#endif
    return
    END SUBROUTINE INITIALIZE_RANDOM

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    SUBROUTINE GENERATE_RANDOM(Rlow, RHigh, N, Rand)
    use kind_module,    only: int_kind, real_kind

    implicit none

    ! Local Variables
    integer(KIND= int_kind), parameter  :: k1 = 714025, &
                                           k2 = 150889,  &
                                           MaxSize = 1366
    ! rMaxSize is the reciprocal of MaxSize
    real(KIND= real_kind), parameter    :: rMaxSize = 7.32064421e-4
    integer(KIND= int_kind)             :: i

    ! Arguments
    integer(Kind= int_kind)             :: N
    real(Kind= real_kind)               :: Rlow, Rhigh
    real(Kind= real_kind), dimension(N) :: Rand
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

#ifdef SYSTEM_RANDOM_GEN
    ! generate the random number sequence using the FORTRAN 90 syntax
    ! these numbers will vary between compilers
         call RANDOM_NUMBER (Rand)
    ! scale the result to the desired interval
         Rand = (Rhigh-Rlow)*Rand + Rlow
#endif
#ifdef LOCAL_RANDOM_GEN
         do i = 1,N
             local_random_seed = MOD(local_random_seed*k1 + k2, MaxSize)
             Rand(i)=local_random_seed*(Rhigh-Rlow)*rMaxSize + Rlow
         enddo
#endif
    return
  END SUBROUTINE GENERATE_RANDOM
END MODULE RANDOM_MODULE
