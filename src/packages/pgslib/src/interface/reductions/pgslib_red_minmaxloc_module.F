MODULE PGSLib_Red_MINMAXLOC_MODULE
  use pgslib_utility_module
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC:: PGSLib_Global_MINLOC
  PUBLIC:: PGSLib_Global_MAXLOC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines provided in this module
!       PGSLib_Global_MINLOC

      INTERFACE PGSLib_Global_MINLOC
         MODULE PROCEDURE PGS_Glbl_MINLOC_Int_1D_F
         MODULE PROCEDURE PGS_Glbl_MINLOC_Real_1D_F
         MODULE PROCEDURE PGS_Glbl_MINLOC_Double_1D_F
      END INTERFACE

      INTERFACE PGSLib_Global_MAXLOC
         MODULE PROCEDURE PGS_Glbl_MAXLOC_Int_1D_F
         MODULE PROCEDURE PGS_Glbl_MAXLOC_Real_1D_F
         MODULE PROCEDURE PGS_Glbl_MAXLOC_Double_1D_F
      END INTERFACE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Interfaces to the C routines

! MINLOC
      INTERFACE PGSLib_Global_MINLOC_C
         subroutine PGSLib_Global_MINLOC_INT_C(MinV, GlobalIndex)
           USE PGSLib_Type_MODULE
           integer (PGSLib_Int_Type), intent(IN OUT) :: GlobalIndex
           integer (PGSLib_Int_Type), intent(IN OUT) :: MinV
         END subroutine PGSLib_Global_MINLOC_INT_C

         subroutine PGSLib_Global_MINLOC_float_C(MinV, GlobalIndex)
           USE PGSLib_Type_MODULE
           integer (PGSLib_Int_Type), intent(IN OUT) :: GlobalIndex
           real    (PGSLib_Real_Type), intent(IN OUT) :: MinV
         END subroutine PGSLib_Global_MINLOC_float_C

         subroutine PGSLib_Global_MINLOC_Double_C(MinV, GlobalIndex)
           USE PGSLib_Type_MODULE
           integer (PGSLib_Int_Type), intent(IN OUT) :: GlobalIndex
           real    (PGSLib_Double_Type), intent(IN OUT) :: MinV
         END subroutine PGSLib_Global_MINLOC_Double_C

      END INTERFACE

! MAXLOC
      INTERFACE PGSLib_Global_MAXLOC_C
         subroutine PGSLib_Global_MAXLOC_INT_C(MinV, GlobalIndex)
           USE PGSLib_Type_MODULE
           integer (PGSLib_Int_Type), intent(IN OUT) :: GlobalIndex
           integer (PGSLib_Int_Type), intent(IN OUT) :: MinV
         END subroutine PGSLib_Global_MAXLOC_INT_C

         subroutine PGSLib_Global_MAXLOC_float_C(MinV, GlobalIndex)
           USE PGSLib_Type_MODULE
           integer (PGSLib_Int_Type), intent(IN OUT) :: GlobalIndex
           real    (PGSLib_Real_Type), intent(IN OUT) :: MinV
         END subroutine PGSLib_Global_MAXLOC_float_C

         subroutine PGSLib_Global_MAXLOC_Double_C(MinV, GlobalIndex)
           USE PGSLib_Type_MODULE
           integer (PGSLib_Int_Type), intent(IN OUT) :: GlobalIndex
           real    (PGSLib_Double_Type), intent(IN OUT) :: MinV
         END subroutine PGSLib_Global_MAXLOC_Double_C

      END INTERFACE

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! MINLOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function PGS_Glbl_MINLOC_Int_1D_F(Vector, MASK)
  USE PGSLib_Type_MODULE
  USE PGSLib_Index_GID_MODULE
  implicit none
  integer (PGSLib_Int_Type), dimension(1)             :: PGS_Glbl_MINLOC_Int_1D_F
  integer (PGSLib_Int_Type), intent(IN), dimension(:) :: Vector
  logical (PGSLib_Log_Type), intent(IN), OPTIONAL, dimension(:) :: MASK

  ! Local variables
  integer (PGSLib_Int_Type) :: MinV, LocalIndex(1), GlobalIndex
  type (PGSLib_GID), POINTER  :: GID
  
  ! First find the local MinLoc, and the associated MINVAL
  if (PRESENT(MASK)) then
     LocalIndex    = MINLOC(Vector, MASK=MASK)
     if (ANY(MASK)) then
        MinV       = Vector(LocalIndex(1))
     else
        MinV       = HUGE(Vector)
        LocalIndex = 1
     end if
  else
     LocalIndex    = MINLOC(Vector)
     MinV          = Vector(LocalIndex(1))
  end if

  ! Convert to global index
  ! We need a GID for that.  We will build one, although
  !    might eventually want to have a geometry manager and cache the GIDs.
  GID => PGSLib_Setup_GID(SIZE(Vector, 1))
  GlobalIndex = PGSLib_Index_Global_From_Local(LocalIndex(1), GID)
  call PGSLib_Deallocate_GID(GID)

  ! Global minloc returns with: global MinV and the corresponding GlobalIndex
  call pgslib_global_minloc_c(MinV, GlobalIndex)

  PGS_Glbl_MINLOC_Int_1D_F(1) = GlobalIndex

  RETURN
END function PGS_Glbl_MINLOC_Int_1D_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function PGS_Glbl_MINLOC_Real_1D_F(Vector, MASK)
  USE PGSLib_Type_MODULE
  USE PGSLib_Index_GID_MODULE
  implicit none
  real (PGSLib_Real_Type), dimension(1)           :: PGS_Glbl_MINLOC_Real_1D_F
  real (PGSLib_Real_Type), intent(IN), dimension(:) :: Vector
  logical (PGSLib_Log_Type), intent(IN), OPTIONAL, dimension(:) :: MASK

  ! Local variables
  integer (PGSLib_Int_Type) :: LocalIndex(1), GlobalIndex
  real (PGSLib_Real_Type) :: MinV
  type (PGSLib_GID), POINTER  :: GID
  
  ! First find the local MinLoc, and the associated MINVAL

  if (PRESENT(MASK)) then
     LocalIndex    = MINLOC(Vector, MASK)
     if (ANY(MASK)) then
        MinV       = Vector(LocalIndex(1))
     else
        MinV       = HUGE(Vector)
        LocalIndex = 1
     end if
  else
     LocalIndex    = MINLOC(Vector)
     MinV          = Vector(LocalIndex(1))
  end if
      

  ! Convert to global index
  ! We need a GID for that.  We will build one, although
  !    might eventually want to have a geometry manager and cache the GIDs.
  GID => PGSLib_Setup_GID(SIZE(Vector, 1))
  GlobalIndex = PGSLib_Index_Global_From_Local(LocalIndex(1), GID)
  call PGSLib_Deallocate_GID(GID)

  ! Global minloc returns with: global MinV and the corresponding GlobalIndex
  call pgslib_global_minloc_c(MinV, GlobalIndex)

  PGS_Glbl_MINLOC_Real_1D_F(1) = GlobalIndex

  RETURN
END function PGS_Glbl_MINLOC_Real_1D_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function PGS_Glbl_MINLOC_Double_1D_F(Vector, MASK)
  USE PGSLib_Type_MODULE
  USE PGSLib_Index_GID_MODULE
  implicit none
  real (PGSLib_Double_Type), dimension(1)             :: PGS_Glbl_MINLOC_Double_1D_F
  real (PGSLib_Double_Type), intent(IN), dimension(:) :: Vector
  logical (PGSLib_Log_Type), intent(IN), OPTIONAL, dimension(:) :: MASK

  ! Local variables
  integer (PGSLib_Int_Type) :: LocalIndex(1), GlobalIndex
  real (PGSLib_Double_Type) :: MinV
  type (PGSLib_GID), POINTER  :: GID
  
  ! First find the local MinLoc, and the associated MINVAL
  if (PRESENT(MASK)) then
     LocalIndex    = MINLOC(Vector, MASK)
     if (ANY(MASK)) then
        MinV       = Vector(LocalIndex(1))
     else
        MinV       = HUGE(Vector)
        LocalIndex = 1
     end if
  else
     LocalIndex    = MINLOC(Vector)
     MinV          = Vector(LocalIndex(1))
  end if
      


  ! Convert to global index
  ! We need a GID for that.  We will build one, although
  !    might eventually want to have a geometry manager and cache the GIDs.
  GID => PGSLib_Setup_GID(SIZE(Vector, 1))
  GlobalIndex = PGSLib_Index_Global_From_Local(LocalIndex(1), GID)
  call PGSLib_Deallocate_GID(GID)

  ! Global minloc returns with: global MinV and the corresponding GlobalIndex
  call pgslib_global_minloc_c(MinV, GlobalIndex)

  PGS_Glbl_MINLOC_Double_1D_F(1) = GlobalIndex

  RETURN
END function PGS_Glbl_MINLOC_Double_1D_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! MAXLOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function PGS_Glbl_MAXLOC_Int_1D_F(Vector, MASK)
  USE PGSLib_Type_MODULE
  USE PGSLib_Index_GID_MODULE
  implicit none
  integer (PGSLib_Int_Type), dimension(1)             :: PGS_Glbl_MAXLOC_Int_1D_F
  integer (PGSLib_Int_Type), intent(IN), dimension(:) :: Vector
  logical (PGSLib_Log_Type), intent(IN), OPTIONAL, dimension(:) :: MASK

  ! Local variables
  integer (PGSLib_Int_Type) :: MaxV, LocalIndex(1), GlobalIndex
  type (PGSLib_GID), POINTER  :: GID
  
  ! First find the local Maxloc, and the associated MAXVAL
  if (PRESENT(MASK)) then
     LocalIndex    = MAXLOC(Vector, MASK)
     if (ANY(MASK)) then
        MaxV       = Vector(LocalIndex(1))
     else
        MaxV       = HUGE(Vector)
        LocalIndex = 1
     end if
  else
     LocalIndex    = MAXLOC(Vector)
     MaxV          = Vector(LocalIndex(1))
  end if

  ! Convert to global index
  ! We need a GID for that.  We will build one, although
  !    might eventually want to have a geometry manager and cache the GIDs.
  GID => PGSLib_Setup_GID(SIZE(Vector, 1))
  GlobalIndex = PGSLib_Index_Global_From_Local(LocalIndex(1), GID)
  call PGSLib_Deallocate_GID(GID)

  ! Global maxloc returns with: global MaxV and the corresponding GlobalIndex
  call pgslib_global_maxloc_c(MaxV, GlobalIndex)

  PGS_Glbl_MAXLOC_Int_1D_F(1) = GlobalIndex

  RETURN
END function PGS_Glbl_MAXLOC_Int_1D_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function PGS_Glbl_MAXLOC_Real_1D_F(Vector, MASK)
  USE PGSLib_Type_MODULE
  USE PGSLib_Index_GID_MODULE
  implicit none
  real (PGSLib_Real_Type), DIMENSION(1)           :: PGS_Glbl_MAXLOC_Real_1D_F
  real (PGSLib_Real_Type), intent(IN), dimension(:) :: Vector
  logical (PGSLib_Log_Type), intent(IN), OPTIONAL, dimension(:) :: MASK

  ! Local variables
  integer (PGSLib_Int_Type) :: LocalIndex(1), GlobalIndex
  real (PGSLib_Real_Type) :: MaxV
  type (PGSLib_GID), POINTER  :: GID
  
  ! First find the local Maxloc, and the associated MAXVAL
  if (PRESENT(MASK)) then
     LocalIndex    = MAXLOC(Vector, MASK)
     if (ANY(MASK)) then
        MaxV       = Vector(LocalIndex(1))
     else
        MaxV       = HUGE(Vector)
        LocalIndex = 1
     end if
  else
     LocalIndex    = MAXLOC(Vector)
     MaxV          = Vector(LocalIndex(1))
  end if

  ! Convert to global index
  ! We need a GID for that.  We will build one, although
  !    might eventually want to have a geometry manager and cache the GIDs.
  GID => PGSLib_Setup_GID(SIZE(Vector, 1))
  GlobalIndex = PGSLib_Index_Global_From_Local(LocalIndex(1), GID)
  call PGSLib_Deallocate_GID(GID)

  ! Global maxloc returns with: global MaxV and the corresponding GlobalIndex
  call pgslib_global_maxloc_c(MaxV, GlobalIndex)

  PGS_Glbl_MAXLOC_Real_1D_F(1) = GlobalIndex

  RETURN
END function PGS_Glbl_MAXLOC_Real_1D_F

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function PGS_Glbl_MAXLOC_Double_1D_F(Vector, MASK)
  USE PGSLib_Type_MODULE
  USE PGSLib_Index_GID_MODULE
  implicit none
  real (PGSLib_Double_Type), dimension(1)             :: PGS_Glbl_MAXLOC_Double_1D_F
  real (PGSLib_Double_Type), intent(IN), dimension(:) :: Vector
  logical (PGSLib_Log_Type), intent(IN), OPTIONAL, dimension(:) :: MASK

  ! Local variables
  integer (PGSLib_Int_Type) :: LocalIndex(1), GlobalIndex
  real (PGSLib_Double_Type) :: MaxV
  type (PGSLib_GID), POINTER  :: GID
  
  ! First find the local Maxloc, and the associated MAXVAL
  if (PRESENT(MASK)) then
     LocalIndex    = MAXLOC(Vector, MASK)
     if (ANY(MASK)) then
        MaxV       = Vector(LocalIndex(1))
     else
        MaxV       = HUGE(Vector)
        LocalIndex = 1
     end if
  else
     LocalIndex    = MAXLOC(Vector)
     MaxV          = Vector(LocalIndex(1))
  end if

  ! Convert to global index
  ! We need a GID for that.  We will build one, although
  !    might eventually want to have a geometry manager and cache the GIDs.
  GID => PGSLib_Setup_GID(SIZE(Vector, 1))
  GlobalIndex = PGSLib_Index_Global_From_Local(LocalIndex(1), GID)
  call PGSLib_Deallocate_GID(GID)

  ! Global maxloc returns with: global MaxV and the corresponding GlobalIndex
  call pgslib_global_maxloc_c(MaxV, GlobalIndex)

  PGS_Glbl_MAXLOC_Double_1D_F(1) = GlobalIndex

  RETURN
END function PGS_Glbl_MAXLOC_Double_1D_F

END MODULE PGSLib_Red_MINMAXLOC_MODULE
