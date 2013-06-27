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
  use kinds, only: r8
  implicit none
  private

  integer, save :: local_random_seed

  ! Public Procedures and Data
  public :: GENERATE_RANDOM, INITIALIZE_RANDOM
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
CONTAINS
    ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>
    SUBROUTINE INITIALIZE_RANDOM()
#ifdef SYSTEM_RANDOM_GEN

    ! Local Variables
    integer :: kr, seedsize

    ! initialize the random number sequence using the FORTRAN 90 syntax
         call RANDOM_SEED (SIZE = seedsize)
         call RANDOM_SEED (PUT = (/(kr, kr = 1,seedsize)/))
#endif
#ifdef LOCAL_RANDOM_GEN
         local_random_seed = 61
#endif
    END SUBROUTINE INITIALIZE_RANDOM

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    SUBROUTINE GENERATE_RANDOM(Rlow, RHigh, N, Rand)

    ! Local Variables
    integer, parameter :: k1 = 714025, k2 = 150889, MaxSize = 1366
    ! rMaxSize is the reciprocal of MaxSize
    real(r8), parameter :: rMaxSize = 7.32064421e-4
    integer :: i

    ! Arguments
    integer :: N
    real(r8) :: Rlow, Rhigh
    real(r8), dimension(N) :: Rand
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
         end do
#endif
  END SUBROUTINE GENERATE_RANDOM
END MODULE RANDOM_MODULE
