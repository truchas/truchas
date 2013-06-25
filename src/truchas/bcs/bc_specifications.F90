Module BC_Specifications
  !-----------------------------------------------------------------------------
  ! Purpose:
  !   Provide the data structures to support boundary condition specifiers.
  !   A SPECIFIER contains all the information to apply boundary
  !   conditions for a given variable.
  !
  ! Provides:
  !
  ! Documentation mostly in the documentation directory.
  !
  ! Author: Robert Ferrell (ferrell@cpca.com)
  !-----------------------------------------------------------------------------
  use bc_enum_types
  use bc_operators
  Implicit None
  Private

  PUBLIC :: BC_Specifier
  PUBLIC :: INITIALIZE
  PUBLIC :: BC_Spec_Get_Operator
  PUBLIC :: BC_Spec_Update_Operator
  PUBLIC :: BC_Spec_Retrieve_Value
  PUBLIC :: UpdateOperatorValue
  PUBLIC :: RetrieveOperatorValue
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Type definitions
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Type to store the data for the BC Specifier
  integer, parameter :: BC_SPEC_NAME_LEN = 256

  ! Type to identify an operator
  type BC_Specifier
     PRIVATE
     character (LEN=BC_SPEC_NAME_LEN) :: SPEC_NAME
     integer                          :: SPEC_ID
     type (BC_Operator), dimension(BC_MAX_OPERATORS) :: Operators
  end type BC_Specifier
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Generic Procedure Interfaces
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTERFACE INITIALIZE
     MODULE PROCEDURE InitSpecifier
  END INTERFACE

  INTERFACE BC_Spec_Get_Operator
     MODULE PROCEDURE GetOperator
  END INTERFACE

  INTERFACE BC_Spec_Update_Operator
     MODULE PROCEDURE UpdateOperatorValue
  END INTERFACE


  INTERFACE BC_Spec_Retrieve_Value
     MODULE PROCEDURE RetrieveOperatorValue
  END INTERFACE

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Data items provided by this module
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

  subroutine InitSpecifier(BC_SPec, NAME, ID)
    implicit none
    type(BC_Specifier), intent(INOUT), target :: BC_Spec
    character (LEN=*),  intent(IN   )         :: NAME
    integer,            intent(IN   )         :: ID

    type (BC_Operator),  POINTER     :: This_Operator
    integer :: Operator
    
    BC_Spec%SPEC_NAME = ''
    BC_Spec%SPEC_NAME = TRIM(NAME)
    BC_Spec%SPEC_ID   = ID
    ! Initialize all the operators
    do Operator = 1, BC_MAX_OPERATORS
       This_Operator => BC_Spec_Get_Operator(BC_Spec, Operator)
       call INITIALIZE(This_Operator, Operator)
    end do
    
  end subroutine InitSpecifier

  subroutine InvalidSpecifier(BC_Spec)
    ! Set the Specifier to invalid, so it will not be used
    ! Does not allocate or free any memory
    implicit none
    type(BC_Specifier), intent(INOUT) :: BC_Spec
    BC_Spec%SPEC_NAME = 'Invalid'
    BC_Spec%SPEC_ID   = BC_INVALID_ID
    return
  end subroutine InvalidSpecifier

  subroutine FreeSpecifier(BC_Spec)
    ! Free all the storage used by a specifier
    implicit none
    type(BC_Specifier), intent(INOUT), target :: BC_Spec
    
    ! Local variables
    type (BC_Operator),  POINTER     :: This_Operator
    integer :: Operator
    ! First free up each of the operators
    do Operator = 1, SIZE(BC_Spec%Operators)
       This_Operator => BC_Spec_Get_Operator(BC_Spec, Operator)
       call FREE(This_Operator)
    end do

    ! Now make this an invalid specifier
    call InvalidSpecifier(BC_Spec)
    return
  end subroutine FreeSpecifier


  function ValidSpecifier(BC_Spec)
    ! Return .TRUE. if this is a valid specifier, .FALSE. otherwise
    implicit none
    type(BC_Specifier), intent(IN   ) :: BC_Spec
    logical                           :: ValidSpecifier

    ValidSpecifier = BC_Spec%SPEC_ID /= BC_INVALID_ID
    return
  end function ValidSpecifier
    

    

  function GetOperator(BC_Spec, OP_ID) RESULT(Operator)
    implicit none
    type (BC_Specifier), intent(IN),  &
                         TARGET      :: BC_Spec
    integer,             intent(IN)  :: OP_ID
    type (BC_Operator),  POINTER     :: Operator
    
    Operator => BC_Spec%Operators(OP_ID)
    RETURN
  end function GetOperator

  subroutine RetrieveOperatorValue(BC_Spec, OPID, mask, Value, ivalue)

    !  Written by Sriram Swaminarayan
    !  This subroutine will take a BC specifier, 
    !  an input bc ID, and a new value and update the
    !  atlas to use the new value wherever that 
    !  particular BCID was used to specify a boundary condition
    use kind_module, only: real_kind, int_kind
    use parameter_module, only: ncells,   &
                                nfc
    use bc_atlases_data_types, only: bc_atlas,   &
                                     data_size,  &
                                     bc_get_face,&
                                     bc_get_cell,&
                                     bc_get_offset,&
                                     bc_get_values

    implicit none
    type (BC_Specifier), intent(IN),  &
                         TARGET      :: BC_Spec
    integer,             intent(IN)  :: OPID
    real(real_kind),     intent(OUT) :: Value
    integer(int_kind),   intent(IN)  :: iValue
    logical, dimension(nfc, ncells), intent(IN) :: mask
    

    integer                    :: i, n, j
    type(BC_OPERATOR), pointer :: Operator
    type(BC_ATLAS),    pointer :: Atlas
    integer, dimension(:), pointer :: bcells
    integer, dimension(:), pointer :: bfaces
    real(real_kind), dimension(:,:), pointer :: AtlasValues
    integer, dimension(:), pointer :: bdyOffsetList
    

    if (OPID == BC_NO_OP) return

    Operator => BC_Spec%Operators(OPID)
    Atlas    => BC_OP_GET_ATLAS(Operator)

    n = DATA_SIZE(Atlas)
    bfaces => BC_Get_Face(Atlas)
    bcells => BC_Get_Cell(Atlas)
    AtlasValues => BC_Get_Values(Atlas)
    bdyOffsetList => BC_Get_Offset(Atlas)
    do i = 1, n
       j = bdyOffsetList(i)
       if ( mask(bfaces(i), bcells(i))) then
          value = AtlasValues(ivalue,j)
          return
       end if
    end do
    
    RETURN
  end subroutine RetrieveOperatorValue


  subroutine UpdateOperatorValue(BC_Spec, OPID, mask, Value, ivalue)

    !  Written by Sriram Swaminarayan
    !  This subroutine will take a BC specifier, 
    !  an input bc ID, and a new value and update the
    !  atlas to use the new value wherever that 
    !  particular BCID was used to specify a boundary condition
    use kind_module, only: real_kind, int_kind
    use parameter_module, only: ncells,   &
                                nfc
    use bc_atlases_data_types, only: bc_atlas,   &
                                     data_size,  &
                                     bc_get_face,&
                                     bc_get_cell,&
                                     bc_get_offset,&
                                     bc_get_values

    implicit none
    type (BC_Specifier), intent(IN),  &
                         TARGET      :: BC_Spec
    integer,             intent(IN)  :: OPID
    real(real_kind),     intent(IN)  :: Value
    integer(int_kind),   intent(IN)  :: iValue
    logical, dimension(nfc, ncells), intent(IN) :: mask
    

    integer                    :: i, n, j
    type(BC_OPERATOR), pointer :: Operator
    type(BC_ATLAS),    pointer :: Atlas
    integer, dimension(:), pointer :: bcells
    integer, dimension(:), pointer :: bfaces
    real(real_kind), dimension(:,:), pointer :: AtlasValues
    integer, dimension(:), pointer :: bdyOffsetList
    
    if (OPID == BC_NO_OP) return

    Operator => BC_Spec%Operators(OPID)
    Atlas    => BC_OP_GET_ATLAS(Operator)

    n = DATA_SIZE(Atlas)
    bfaces => BC_Get_Face(Atlas)
    bcells => BC_Get_Cell(Atlas)
    AtlasValues => BC_Get_Values(Atlas)
    bdyOffsetList => BC_Get_Offset(Atlas)
    do i = 1, n
       j = bdyOffsetList(i)
       if ( mask(bfaces(i), bcells(i))) then
          AtlasValues(ivalue,j) = value
       end if
    end do
    
    RETURN
  end subroutine UpdateOperatorValue

END Module BC_SPECIFICATIONS
