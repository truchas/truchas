MODULE TEMPGRAD_MODULE
  !======================================================================
  ! Purpose:
  !  Provide a Body_Temp_Gradient Type, array and function
  !
  ! Description:
  !  A cylindrical temperature gradient option is provided of the form
  !    T(dr,dz) = T0*(P(dr)+P(dz))
  !  where P(dr) = 1 + A*dr**L + B*dr**M + C*dr**N
  !  with (real) coefficients (A,B,C) = %R_Constants and
  !  (integer) powers (L,M,N) = %R_Exponents.
  !  P(dz) is similarly defined with coefficients %Z_Constants and
  !  powers %Z_Exponents.
  !
  !  The constant-term coefficients are defaulted to one. The temperature at
  !  the origin (T0) is taken from the current Temperature body namelist item,
  !  which is stored in Body_Temp()
  !
  !  dZ and dR in this gradient are measured from the point TG_Origin with
  !  the Z axis oriented in the TG_AXIS direciton.
  !
  !  A optional cutoff point for the gradient in Z is provided in TG_Z_LBound.
  !  The constant temperature is used for z < TG_Z_LBound (no variation
  !  in r or z).
  !
  ! Contains:
  !   Body_Temperature(ibody, zone, cell)
  !
  ! Author: Larry J. Cox (ljcox@lanl.gov
  ! Modified by Dave Korzekwa to add the two polynomials instead of 
  !  multiplying them.
  !======================================================================
  use kinds, only: r8
  USE parameter_module, ONLY: ndim, mbody
  IMPLICIT NONE
  PRIVATE

  ! Public Functions
  PUBLIC:: Body_Temperature

  ! Public TYPE
  PUBLIC:: TEMP_GRADIENT

  ! Internals:
  integer, PARAMETER            :: maxTerms = 3

  TYPE TEMP_GRADIENT
    real(r8), DIMENSION(ndim)      :: Origin = [0.0_r8, 0.0_r8, 0.0_r8]
    real(r8), DIMENSION(ndim)      :: Axis   = [0.0_r8, 0.0_r8, 1.0_r8]
    real(r8), DIMENSION(maxTerms)  :: Z_Constants = [0.0_r8, 0.0_r8, 0.0_r8]
    real(r8), DIMENSION(maxTerms)  :: R_Constants = [0.0_r8, 0.0_r8, 0.0_r8]
    real(r8), DIMENSION(2)         :: Z_Bound = [-huge(1.0_r8), huge(1.0_r8)]
    real(r8), DIMENSION(2)         :: R_Bound = [-huge(1.0_r8), huge(1.0_r8)]
    integer, DIMENSION(maxTerms):: Z_Exponents = [1, 2, 3]
    integer, DIMENSION(maxTerms):: R_Exponents = [1, 2, 3]
    logical                     :: On = .false.
  END TYPE TEMP_GRADIENT

  ! Public Variable
  TYPE(TEMP_GRADIENT), DIMENSION(mbody), PUBLIC  :: Body_Temp_Grad

CONTAINS

  PURE FUNCTION BODY_TEMPERATURE(ibody, TempOrigin, Zone, Cell) RESULT(NewTemp)

    ! Takes a body index, a single zone and its matching cell
    ! Returns the temperature value at the centroid of the zone/cell
    ! Note: This function assumes the input Zone and Cell
    !       refer to the same location

    USE cutoffs_module,    ONLY: alittle
    USE zone_module,       ONLY: CELL_AVG
    USE mesh_module,       ONLY: CELL_GEOMETRY

    ! Arguments
    integer, INTENT(in) :: ibody
    real(r8),   INTENT(in) :: TempOrigin
    TYPE(CELL_AVG), INTENT(in) :: Zone
    TYPE(CELL_GEOMETRY), INTENT(in) :: Cell

    ! Result
    real(r8) :: NewTemp

    ! Locals
    real(r8) :: cnt(ndim), tga(ndim)
    real(r8) :: dR, dZ, pR, pZ
    integer :: ip
    TYPE(TEMP_GRADIENT) :: BTG
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    BTG = Body_Temp_Grad(ibody)

    cnt = Cell%Centroid - BTG%Origin ! Centroid w.r.t. BTG Origin
    tga = BTG%Axis

    ! dZ: offset from gradient origin along gradient axis
    dZ = DOT_PRODUCT(tga,cnt)
    IF (ABS(dZ) < alittle) dZ = 0
    IF (dZ < BTG%Z_Bound(1)) dZ = BTG%Z_Bound(1)
    IF (dZ > BTG%Z_Bound(2)) dZ = BTG%Z_Bound(2)

    ! dR: perpendicular offset from gradient axis
    !     magnitude of cross_product (tga)X(cnt)
    dR =      ( tga(2)*cnt(3) - tga(3)*cnt(2) ) ** 2
    dR = dR + ( tga(1)*cnt(3) - tga(3)*cnt(1) ) ** 2
    dR = dR + ( tga(1)*cnt(2) - tga(2)*cnt(1) ) ** 2
    dR = SQRT(dR)
    IF (dR < alittle) dR = 0
    IF (dR < BTG%R_Bound(1)) dR = BTG%R_Bound(1)
    IF (dR > BTG%R_Bound(2)) dR = BTG%R_Bound(2)

    ! Polynomial values
    pZ = 1
    pR = 0
    DO ip = 1, maxTerms
      pZ = pZ + BTG%Z_Constants(ip) * dZ**BTG%Z_Exponents(ip)
      pR = pR + BTG%R_Constants(ip) * dR**BTG%R_Exponents(ip)
    END DO

    NewTemp = TempOrigin * (pZ + pR)
    IF (NewTemp < alittle) NewTemp = 0

  END FUNCTION BODY_TEMPERATURE

END MODULE TEMPGRAD_MODULE
