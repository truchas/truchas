MODULE LIMITER_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Encapsulate all scalars, arrays, and procedures
  !   related to the slope-limiting algorithms
  !
  ! Public Interface(s):
  !
  !   * call LIMITER (Extrapolated_Value, Cell_Value, Neighbor_Values, &
  !                   Slope_Limiter, element, location, method)
  !
  ! Contains: LIMITER
  !
  ! Author(s): Douglas B. Kothe (LANL Group T-3, dbk@lanl.gov)
  !
  !=======================================================================
  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: LIMITER

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Limiter namelist variables
  character(LEN = 80), public, save :: limiter_type

CONTAINS

  SUBROUTINE LIMITER (Extrapolated_Value, Cell_Value, Neighbor_Values, &
                      Slope_Limiter, element, location, method)
    !=======================================================================
    ! PURPOSE -
    !      Initialize parallel parameters, such as the number of PEs,
    !      the PE number of this process, and the I/O PE
    !=======================================================================
    use constants_module, only: one, one_third, zero, two
    use cutoffs_module,   only: alittle
    use kind_module,      only: int_kind, log_kind, real_kind
    use parameter_module, only: ncells, nfc
    use mesh_module,      only: Cell

    implicit none

    ! Arguments
    real(KIND = real_kind), dimension(ncells),     intent(IN)  :: Extrapolated_Value
    real(KIND = real_kind), dimension(ncells),     intent(IN)  :: Cell_Value
    real(KIND = real_kind), dimension(nfc,ncells), intent(IN)  :: Neighbor_Values
    real(KIND = real_kind), dimension(ncells),     intent(OUT) :: Slope_Limiter
    integer(KIND = int_kind),                      intent(IN)  :: element
    character(LEN = *), optional,                  intent(IN)  :: location
    character(LEN = *), optional,                  intent(IN)  :: method

    ! Local Variables
    logical(KIND = log_kind) :: specified_method, specified_location
    character(LEN = 80) :: limiter_method = 'Venkat', &
                           limiter_location = 'face'
    real(KIND = real_kind)                    :: Venkat_constant = one_third
    real(KIND = real_kind), dimension(ncells) :: Tmp1, Tmp2, Tmp3, Tmp4

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

    ! Process the incoming arguments
    specified_method   = PRESENT(method)
    specified_location = PRESENT(location)
    if (specified_method) limiter_method = method
    if (specified_location) limiter_location = location

    ! Compute the limiter according to the specified location and method
    select case(limiter_method)

       ! No limiting; set slope limiter to unity
       case default

          Slope_Limiter = one

       ! Venkatakrishnan's Limiter (AIAA Paper #AIAA-93-0880)
       case ('Venkat')

          Tmp1 = MAX(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp2 = MIN(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp3 = Extrapolated_Value - Cell_Value
          Tmp3 = SIGN(one,Tmp3)*(ABS(Tmp3) + alittle)
          Tmp4 = (Cell%Volume)**one_third  ! Approximation of cell length
          Tmp4 = (Venkat_constant*Tmp4)**3
          where (Tmp3 > alittle)
             Slope_Limiter = (Tmp3*(Tmp1**2 + Tmp4) + two*Tmp1*Tmp3**2) / &
                         (Tmp1**2 + two*Tmp3**2 + Tmp1*Tmp3 + Tmp4 + alittle)
          elsewhere
             Slope_Limiter = (Tmp3*(Tmp2**2 + Tmp4) + two*Tmp2*Tmp3**2) / &
                         (Tmp2**2 + two*Tmp3**2 + Tmp2*Tmp3 + Tmp4 + alittle)
          end where
          Slope_Limiter = MERGE(one, Slope_Limiter/Tmp3, ABS(Tmp3) <= alittle)

       ! Barth's Limiter (AIAA Paper #AIAA-89-0366)
       case ('Barth')

          Tmp1 = MAX(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp2 = MIN(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp3 = Extrapolated_Value - Cell_Value
          Tmp3 = MERGE(alittle, Tmp3, ABS(Tmp3) <= alittle)
          where (Tmp3 >= alittle)
             Slope_Limiter = MIN(one, Tmp1/Tmp3)
          elsewhere
             Slope_Limiter = MIN(one, Tmp2/Tmp3)
          end where
          Slope_Limiter = MERGE(one, Slope_Limiter, ABS(Tmp3) <= alittle)

    end select

    ! Make sure the limiter is bounded by zero and one
    Slope_Limiter = MIN(one, MAX(zero, Slope_Limiter))

    return

  END SUBROUTINE LIMITER

END MODULE LIMITER_MODULE
