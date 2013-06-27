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
  use kinds, only: r8
  implicit none
  private

  public :: LIMITER

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Limiter namelist variables
  character(80), public, save :: limiter_type

CONTAINS

  SUBROUTINE LIMITER (Extrapolated_Value, Cell_Value, Neighbor_Values, &
                      Slope_Limiter, element, location, method)
    !=======================================================================
    ! PURPOSE -
    !      Initialize parallel parameters, such as the number of PEs,
    !      the PE number of this process, and the I/O PE
    !=======================================================================
    use cutoffs_module,   only: alittle
    use parameter_module, only: ncells, nfc
    use mesh_module,      only: Cell

    ! Arguments
    real(r8), dimension(ncells),     intent(IN)  :: Extrapolated_Value
    real(r8), dimension(ncells),     intent(IN)  :: Cell_Value
    real(r8), dimension(nfc,ncells), intent(IN)  :: Neighbor_Values
    real(r8), dimension(ncells),     intent(OUT) :: Slope_Limiter
    integer, intent(IN)  :: element
    character(*), optional, intent(IN)  :: location
    character(*), optional, intent(IN)  :: method

    ! Local Variables
    logical :: specified_method, specified_location
    character(80) :: limiter_method = 'Venkat', limiter_location = 'face'
    real(r8) :: Venkat_constant = 1.0_r8 / 3.0_r8
    real(r8), dimension(ncells) :: Tmp1, Tmp2, Tmp3, Tmp4

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

          Slope_Limiter = 1.0_r8

       ! Venkatakrishnan's Limiter (AIAA Paper #AIAA-93-0880)
       case ('Venkat')

          Tmp1 = MAX(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp2 = MIN(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp3 = Extrapolated_Value - Cell_Value
          Tmp3 = SIGN(1.0_r8,Tmp3)*(ABS(Tmp3) + alittle)
          Tmp4 = (Cell%Volume)**(1.0d0/3.0d0)  ! Approximation of cell length
          Tmp4 = (Venkat_constant*Tmp4)**3
          where (Tmp3 > alittle)
             Slope_Limiter = (Tmp3*(Tmp1**2 + Tmp4) + 2.0_r8*Tmp1*Tmp3**2) / &
                         (Tmp1**2 + 2.0_r8*Tmp3**2 + Tmp1*Tmp3 + Tmp4 + alittle)
          elsewhere
             Slope_Limiter = (Tmp3*(Tmp2**2 + Tmp4) + 2.0_r8*Tmp2*Tmp3**2) / &
                         (Tmp2**2 + 2.0_r8*Tmp3**2 + Tmp2*Tmp3 + Tmp4 + alittle)
          end where
          Slope_Limiter = MERGE(1.0_r8, Slope_Limiter/Tmp3, ABS(Tmp3) <= alittle)

       ! Barth's Limiter (AIAA Paper #AIAA-89-0366)
       case ('Barth')

          Tmp1 = MAX(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp2 = MIN(Cell_Value, Neighbor_Values(element,:)) - Cell_Value
          Tmp3 = Extrapolated_Value - Cell_Value
          Tmp3 = MERGE(alittle, Tmp3, ABS(Tmp3) <= alittle)
          where (Tmp3 >= alittle)
             Slope_Limiter = MIN(1.0_r8, Tmp1/Tmp3)
          elsewhere
             Slope_Limiter = MIN(1.0_r8, Tmp2/Tmp3)
          end where
          Slope_Limiter = MERGE(1.0_r8, Slope_Limiter, ABS(Tmp3) <= alittle)

    end select

    ! Make sure the limiter is bounded by zero and one
    Slope_Limiter = MIN(1.0_r8, MAX(0.0_r8, Slope_Limiter))

  END SUBROUTINE LIMITER

END MODULE LIMITER_MODULE
