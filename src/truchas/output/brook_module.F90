!
!  TRUCHAS DEVELOPERS:
!  
!  Warning: Do not use this module
#ifndef BYTES_PER_CHAR
#define BYTES_PER_CHAR     1
#endif

#ifndef BYTES_PER_SHORT
#define BYTES_PER_SHORT    2
#endif

#ifndef BYTES_PER_INT
#define BYTES_PER_INT      4
#endif

#ifndef BYTES_PER_LOGICAL
#define BYTES_PER_LOGICAL  4
#endif

#ifndef BYTES_PER_FLOAT
#define BYTES_PER_FLOAT    4
#endif

#ifndef BYTES_PER_DOUBLE
#define BYTES_PER_DOUBLE   8
#endif



! define the following if you want to run the local test program
! This will make it a local compile with main
! #define LOCALPROGRAM 1
! #define BROOK_DEBUG  1
! #define BROOK_DEBUG_UNIT 1
#ifdef BROOK_DEBUG
#ifndef BROOK_DEBUG_UNIT
#define BROOK_DEBUG_UNIT 1
#endif
#endif

module brook_module


!======================================================================
! Purpose(s):
! Provide support for output
!
! Public Interface:
! Parameters:   B_Strlen             ! Length of all internal strings
!
! Type :        Brook                ! Brook Structure
!
! Variables:
!               B_Stdout             ! Standard output brook
!               B_Stderr             ! Standard error brook
! 
! Subroutines:  
!               Brook_Set            ! Set brook prameters
!               Brook_Close          ! Close a brook file
!               Brook_Destroy        ! Destroy data associated with brook
!               Brook_Diagnostics    ! Diagnostic routine, not complete
!               Brook_Write          ! Write Data To Brooks
!               Brook_Write_XMLData  ! Write Complete XML Tag with data to xml brooks
!               Brook_OpenXMLTag     ! Open XML Tag in xml brooks
!               Brook_CloseXMLTag    ! Close XML Tag in xml brooks
!               Brook_WriteXMLComment! Add a comment to an xml brook
!               
! Query Functions:    
!               Brook_File
!               Brook_Unit
!               Brook_Form
!               Brook_Position
!               Brook_Closeable
!               Brook_[IDLA]Format
!               Brook_XMLRoot
!             
! Operators:    (B1) + (B2) [ + (B3) + (B5) ... ]
!
! Assignment:   (B1) = (B3) [ + (B4) + (B5) ... ]
!
! This module was written as a replacement for the streams module. 
! Basically it is a small subset of the stream module, hence the
! name brook!  
!
! I am hoping that this will not be as inscrutable as the stream
! module.  The basic premise here is that you open a number of
! brooks and pyss data into them.  You can combine brooks and the
! resultant brook will write data to both brooks.  
! 
! All entities in this module are prefixed with a 'Brook_' or a 'B_' 
! in the case of parameters so that they are easily identifiable in 
! the rest of the code.  
!
! The fundamental data type is Brook.  You will define a
! Brook (with attribute target) for every file that you want to
! write to.  
!   
! You can add brooks together as follows:
!   B1 = B1 + B3 + B4
! The addition will put copies of B3 and B4 into B1, so even if you
! destroy B1, you still have to destroy B3 and B4. 
!
! One important side effect of putting copies during addition is
! that if you should close the file, and then write to another
! brook that has that file name internally, you will lose all the
! data that was written to that file earlier since the file was
! reopened and written to.  Further, since these additions are COPIES, 
! any changes you make to the format will not be reflected in the 
! original entries.  Consequently, if you make a brook non-Closeable,
! then the original brook will still be Closeable. 
!
! The attributes of a Brook are set with Brook_Set.  Brook_Set has the
! following call frame:   
!   Subroutine Brook_Set(B,               & ! Brook, target
!                        File,            & ! File name to be associated with B
!                        Unit,            & ! Unit number for B
!                        Form,            & ! 'ascii', 'binary', or 'xml'
!                        Position,        & ! 'asis', 'rewind', or 'append', default: 'rewind'
!                        StandardHeaders, & ! will standard headers be written?
!                        Closeable,       & ! Can this file be closed?
!                        IFormat,         & ! Default integer format
!                        DFormat,         & ! Default real format
!                        LFormat,         & ! Default Logical Format
!                        AFormat,         & ! Default character Format
!                        XMLRoot,         & ! The XMLRoot for xml files
!                        iStatus)           ! Did an error occur?            
!  except for B, all other arguments are optional.
!  Remember, if you change form midbrook to binary from ascii, do not
!  forget to switch it back to ascii (and vice versa)
!
!  Subroutine Brook_Close(B,iStatus)
!   Close the file associated with struct B
!   Note that in case of composite Bs (i.e. where several brooks
!   have been strung together, this subroutine will close all files.
!
!  Subroutine Brook_Destroy(B, iStatus) 
!   Destroy data associated with struct B
!   You must call this when you are done with struct B or else
!   memory will leak.
!
!  Subroutine Brook_Write(B,                 & ! Brook, target
!                         Variable,          & ! Variable to be written, 
!                                            & ! real/integer/char/logical
!                                            & ! Constant or array 
!                         Format,            & ! Format for writing, default format 
!                                            & ! used if not present, ignored for 
!                                            & ! XML and binary files
!                         C_Array_Order,     & ! If present and true, C order , otherwise fortran order
!                                              ! for arrays
!                         XMLName,           & ! Name label for XML output
!                         XMLAttributes,     & ! Any attributes you want in the XML name tag
!                         XMLDataFile,       & ! Look-aside data file for XML output
!                         XMLDataFormat,     & ! Format of the look-aside data file
!                         iStatus            & ! Did an error occur?
!                         )

! Writes data in Variable to the different Brooks represented by B.
! All arguments except B are optional.
! The data is written based on the type of brook.  If it is an
! ascii brook, the data is written using Format.  If Format is not
! present, the data will be written using the default format.
! Argument Format is ignored for XML and binary brooks.  The
! arguments starting with 'XML' are only for XML brooks.  XMLName
! is the tag used for the XML field, and XMLAttributes are
! additional attributes you may want to send for the XMLTag.  If
! XMLDataFile is present, the variable data will be put into a file
! with the name contained therein.  Otherwise, the data will be
! contained inline in a tag called <Data> as a subset of variable.
! If XMLDataFormat is absent, the data will be written to
! XMLDataFile in ascii format.  If XMLDataFormat is present, the
! data will be written to XMLDataFile in that data format.  Valid
! values are 'xml', 'ascii', and 'binary'. iStatus returns an error
! status.  A value of 0 means no error.  If iStatus is absent, any 
! error will result in instant death.
! 
!
! 
! Subroutine Brook_WriteXMLData(B,             & ! Brook, target
!                               XMLTag,        & ! Name of tag to be opened
!                               XMLAttributes, & ! Attributes for tag to be opened
!                               XMLStringData, & ! String with data to be written inline
!                               iStatus        & ! Did an error occur?
!                              )
! This subroutine will write string data to any xml or ascii brook 
! in B and close the tag after the data is written. The data is 
! written inline, which is different from the openXMLTag and 
! closeXMLTag routines.  
!
! Subroutine Brook_OpenXMLTag(B,             & ! Brook, target
!                             XMLTag,        & ! Name of tag to be opened
!                             XMLAttributes, & ! Attributes for tag to be opened
!                             iStatus        & ! Did an error occur?
!                            )
! Opens an XML tag in any xml and ascii brooks present in B
!  
! Subroutine Brook_CloseXMLTag(B,             & ! Brook, target
!                              XMLTag,        & ! Name of tag to be opened
!                              XMLAttributes, & ! Attributes for tag to be opened
!                              iStatus        & ! Did an error occur?
!                             )
! Closes an XML tag to any xml and ascii brooks present in B
!
! Subroutine Brook_WriteXMLComment(B,         & ! Brook, target
!                                  Comment,   & ! Comment to be written
!                                  iStatus,   & ! did an error occur?
! Writes an xml commetn to all xml brooks in B
! If comment is a string array, each element
!
! Subroutine Brook_Flush(B, iStatus)
! Flushes all output pending to all brooks in B
!  
! Query Functions:    
!               Brook_File(B, iStatus)               Returns file name
!               Brook_Unit(B, iStatus)               Returns unit number
!               Brook_Form(B, iStatus)               Returns unit format
!               Brook_Position(B, iStatus)           Returns unit position
!               Brook_Closeable(B, iStatus)          Can Unit be closed  with Brook_Close
!               Brook_[IDLA]Format(B, iStatus)       Returns I, D, L, or A Format
!               Brook_XMLRoot(B, iStatus)            Returns XML Root
!               Brook_Next(B, iStatus)               Returns the next Brook
!
! Author(s): Sriram Swaminarayan (sriram@lanl.gov)
!            Februray 12, 2004
!
!======================================================================
IMPLICIT NONE

! Private Module
PRIVATE

! The external interface will only see the following public entities

public :: B_OP
PUBLIC :: BROOK
PUBLIC :: B_STDOUT
PUBLIC :: B_STDERR
!PUBLIC :: ASSIGNMENT(=)
PUBLIC :: OPERATOR(+)
! PUBLIC :: Brook_WRITE ! I use the individual routines from tbrook_Write()
PUBLIC :: Brook_OPENXMLTAG
PUBLIC :: Brook_CLOSEXMLTAG
PUBLIC :: BROOK_WRITEXMLDATA
PUBLIC :: BROOK_WRITEXMLTAG
PUBLIC :: BROOK_WRITEXMLCOMMENT
PUBLIC :: B_DATA
PUBLIC :: Brook_NEXT
PUBLIC :: Brook_CLOSE
PUBLIC :: BROOK_POINTER
PUBLIC :: Brook_DESTROY
PUBLIC :: Brook_SET
PUBLIC :: Brook_IFORMAT
PUBLIC :: Brook_DFORMAT
PUBLIC :: Brook_LFORMAT
PUBLIC :: Brook_AFORMAT
PUBLIC :: Brook_FILE
PUBLIC :: Brook_Flush
PUBLIC :: Brook_OFFSET
PUBLIC :: Brook_UNIT
PUBLIC :: Brook_FORM
PUBLIC :: Brook_POSITION
PUBLIC :: Brook_XMLROOT
PUBLIC :: Brook_CLOSEABLE
PUBLIC :: Brook_ADVANCE
PUBLIC :: Brook_Get_Error
PUBLIC :: Brook_Clear_Error
PUBLIC :: Brook_AssignUnit
PUBLIC :: Brook_ReleaseUnit
PUBLIC :: Brook_Remove_Duplicates
PUBLIC :: B_ToLowerCase
PUBLIC :: Brook_IForm
PUBLIC :: B_IFORM_ASCII
PUBLIC :: B_IFORM_XML
PUBLIC :: B_IFORM_BINARY
PUBLIC :: Brook_Write_File_Header


! Parameters
INTEGER,                 PARAMETER :: B_STDOUT_UNIT      = 6
INTEGER,                 PARAMETER :: B_STDERR_UNIT      = 7
INTEGER,                 PARAMETER :: B_STRLEN           = 4096
CHARACTER(LEN=B_STRLEN), PARAMETER :: BROOK_XML_VERSION  = 'Version="1.0"'
INTEGER,                 PARAMETER :: B_STATUS_NO_INIT   = -8
INTEGER,                 PARAMETER :: B_AUTO_UNIT        = -7
INTEGER,                 PARAMETER :: B_INVALID_UNIT     = -6
INTEGER,                 PARAMETER :: B_IFORM_UNKNOWN    = -5
INTEGER,                 PARAMETER :: B_ID_UNKNOWN       = -4
INTEGER,                 PARAMETER :: B_STATUS_UNKNOWN   = -2
INTEGER,                 PARAMETER :: B_UNIT_UNKNOWN     = -3
INTEGER,                 PARAMETER :: B_UNKNOWN          = -1
INTEGER,                 PARAMETER :: B_ISOPEN           = 1
INTEGER,                 PARAMETER :: B_ISCLOSED         = 2
INTEGER,                 PARAMETER :: B_TOPUNIT          = 3
INTEGER,                 PARAMETER :: B_ALLUNITS         = 4  
INTEGER,                 PARAMETER :: B_IFORM_ASCII      = 5
INTEGER,                 PARAMETER :: B_IFORM_BINARY     = 6
INTEGER,                 PARAMETER :: B_IFORM_XML        = 7
INTEGER,                 PARAMETER :: B_STATUS_INIT      = 8
INTEGER,                 PARAMETER :: B_UNIT_LOW         = 128
INTEGER,                 PARAMETER :: B_UNIT_HIGH        = 512
INTEGER,                 PARAMETER :: B_UNIT_INCREMENT   = 1
CHARACTER(len=B_Strlen), PARAMETER :: B_Position_Default = 'rewind'
CHARACTER(len=B_Strlen), PARAMETER :: B_File_Unknown     = 'Invalid File'
CHARACTER(len=B_Strlen), PARAMETER :: B_Form_Unknown     = 'Invalid Form'
CHARACTER(len=B_Strlen), PARAMETER :: B_Iformat_default  = '(I16)'
CHARACTER(len=B_Strlen), PARAMETER :: B_Dformat_default  = '(1es30.20)'
CHARACTER(len=B_Strlen), PARAMETER :: B_Lformat_default  = '(L1,x)'
CHARACTER(len=B_Strlen), PARAMETER :: B_Aformat_default  = '(a)'
CHARACTER(len=B_Strlen), PARAMETER :: B_XMLRoot_default  = 'BrookXML'


!   Brook Error Codes
INTEGER, PARAMETER :: BE_NONE              = 0
INTEGER, PARAMETER :: BE_INIT              = 1
INTEGER, PARAMETER :: BE_EXIT              = 2
INTEGER, PARAMETER :: BE_WRITE             = 3
INTEGER, PARAMETER :: BE_ASSIGN            = 4
INTEGER, PARAMETER :: BE_ADD               = 5
INTEGER, PARAMETER :: BE_SET               = 6
INTEGER, PARAMETER :: BE_OPEN              = 7
INTEGER, PARAMETER :: BE_CLOSE             = 8
INTEGER, PARAMETER :: BE_DESTROY           = 9
INTEGER, PARAMETER :: BE_QUERY             = 10
INTEGER, PARAMETER :: BE_RESHAPE           = 11
INTEGER, PARAMETER :: BE_ALLOCATE          = 12


! Error information holders
INTEGER,                 save :: BE_ID      = BE_NONE
CHARACTER(LEN=B_STRLEN), save :: BE_STRING  = ' ' 


TYPE :: B_Data
   PRIVATE
   ! Data in this type can only be modified by this module
   CHARACTER(LEN=B_Strlen) :: File
   INTEGER                 :: ID
   INTEGER                 :: Offset
   INTEGER                 :: Status
   INTEGER                 :: Unit
   INTEGER                 :: IForm
   LOGICAL                 :: StandardHeaders
   LOGICAL                 :: Closeable
   LOGICAL                 :: ADVANCE
   CHARACTER(LEN=B_Strlen) :: Form
   CHARACTER(LEN=B_Strlen) :: Position
   CHARACTER(LEN=B_Strlen) :: Iformat
   CHARACTER(LEN=B_Strlen) :: Dformat
   CHARACTER(LEN=B_Strlen) :: Lformat
   CHARACTER(LEN=B_Strlen) :: Aformat
   CHARACTER(LEN=B_Strlen) :: XMLRoot
END TYPE B_Data

TYPE :: Brook
   TYPE(B_Data),   POINTER :: DATA => Null()
   TYPE(Brook), POINTER :: Next => Null()
END TYPE Brook

TYPE :: Brook_Pointer
   TYPE(Brook),         POINTER :: S => NULL()
   TYPE(Brook_POINTER), POINTER :: Next => NULL()
END TYPE Brook_Pointer

INTERFACE OPERATOR (+)
   ! This will return a struct pointer
   MODULE PROCEDURE B_ADD_S_S
   MODULE PROCEDURE B_ADD_P_S
   MODULE PROCEDURE B_ADD_S_P
   MODULE PROCEDURE B_ADD_P_P
END INTERFACE

INTERFACE ASSIGNMENT (=)
   MODULE PROCEDURE B_ASSIGN
   MODULE PROCEDURE B_ASSIGN_STRUCT_POINTER
END INTERFACE

INTERFACE BROOK_WRITEXMLComment
   Module Procedure B_WXMLC_String
   Module Procedure B_WXMLC_StringArray
END INTERFACE

INTERFACE Brook_WRITE
   MODULE PROCEDURE B_WRITE_INTEGER_R0
   MODULE PROCEDURE B_WRITE_LOGICAL_R0
   MODULE PROCEDURE B_WRITE_STRING_R0
   MODULE PROCEDURE B_WRITE_DOUBLE_R0
   MODULE PROCEDURE B_WRITE_FLOAT_R0

   MODULE PROCEDURE B_WRITE_INTEGER_R1
   MODULE PROCEDURE B_WRITE_LOGICAL_R1
   MODULE PROCEDURE B_WRITE_STRING_R1
   MODULE PROCEDURE B_WRITE_DOUBLE_R1
   MODULE PROCEDURE B_WRITE_FLOAT_R1

   MODULE PROCEDURE B_WRITE_INTEGER_R2
   MODULE PROCEDURE B_WRITE_LOGICAL_R2
   MODULE PROCEDURE B_WRITE_DOUBLE_R2
   MODULE PROCEDURE B_WRITE_FLOAT_R2

   MODULE PROCEDURE B_WRITE_INTEGER_R3
   MODULE PROCEDURE B_WRITE_LOGICAL_R3
   MODULE PROCEDURE B_WRITE_DOUBLE_R3
   MODULE PROCEDURE B_WRITE_FLOAT_R3

   MODULE PROCEDURE B_WRITE_INTEGER_R4
   MODULE PROCEDURE B_WRITE_LOGICAL_R4
   MODULE PROCEDURE B_WRITE_DOUBLE_R4
   MODULE PROCEDURE B_WRITE_FLOAT_R4

   MODULE PROCEDURE B_WRITE_INTEGER_R5
   MODULE PROCEDURE B_WRITE_LOGICAL_R5
   MODULE PROCEDURE B_WRITE_DOUBLE_R5
   MODULE PROCEDURE B_WRITE_FLOAT_R5
END INTERFACE

! Variables
LOGICAL                          :: B_Initialized = .FALSE.
TYPE(Brook),              TARGET :: B_Tmp
TYPE(Brook),              TARGET :: B_Stdout
TYPE(Brook),              TARGET :: B_Stderr
INTEGER                          :: B_Tmp_Unit
logical,  dimension(B_UNIT_HIGH) :: B_AssignedUnits = .false.
TYPE(B_Data), PRIVATE, PARAMETER :: B_data_default = B_Data(                     &
                                     TRIM(B_File_Unknown),     & !File
                                     B_ID_Unknown,             & !ID
                                     -1,                       & !Offset
                                     B_Status_no_init,         & !Status
                                     B_Unit_Unknown,           & !Unit
                                     B_IForm_Unknown,          & !IForm
                                     .TRUE.,                   & !StandardHeaders
                                     .TRUE.,                   & !Closeable
                                     .FALSE.,                  & !ADVANCE
                                     TRIM(B_Form_Unknown),     & !Form
                                     TRIM(B_POSITION_default), & !Form
                                     TRIM(B_Iformat_default),  & !Iformat
                                     TRIM(B_Dformat_default),  & !Dformat
                                     TRIM(B_Lformat_default),  & !Lformat
                                     TRIM(B_Aformat_default),  & !Aformat
                                     TRIM(B_XMLRoot_default)   & !XMLRoot
                                     )


SAVE :: B_Tmp
SAVE :: B_Tmp_Unit
SAVE :: B_STDOUT
SAVE :: B_STDERR
SAVE :: B_Initialized
SAVE :: B_AssignedUnits

! Entities I have to make public because they belong to an
! interface.  These are not at the top because the average joe
! doesn't need to know about them.  
PUBLIC :: Brook_DIAGNOSTICS

PUBLIC :: B_WXMLC_String
PUBLIC :: B_WXMLC_StringArray

PUBLIC :: B_ASSIGN_STRUCT_POINTER
PUBLIC :: B_ASSIGN

PUBLIC :: B_ADD_S_S
PUBLIC :: B_ADD_S_P
PUBLIC :: B_ADD_P_S
PUBLIC :: B_ADD_P_P

PUBLIC :: B_WRITE_INTEGER_R0
PUBLIC :: B_WRITE_LOGICAL_R0
PUBLIC :: B_WRITE_DOUBLE_R0
PUBLIC :: B_WRITE_FLOAT_R0
PUBLIC :: B_WRITE_STRING_R0

PUBLIC :: B_WRITE_INTEGER_R1
PUBLIC :: B_WRITE_LOGICAL_R1
PUBLIC :: B_WRITE_DOUBLE_R1
PUBLIC :: B_WRITE_STRING_R1
PUBLIC :: B_WRITE_FLOAT_R1

PUBLIC :: B_WRITE_INTEGER_R2
PUBLIC :: B_WRITE_LOGICAL_R2
PUBLIC :: B_WRITE_DOUBLE_R2
PUBLIC :: B_WRITE_FLOAT_R2

PUBLIC :: B_WRITE_INTEGER_R3
PUBLIC :: B_WRITE_LOGICAL_R3
PUBLIC :: B_WRITE_DOUBLE_R3
PUBLIC :: B_WRITE_FLOAT_R3

PUBLIC :: B_WRITE_INTEGER_R4
PUBLIC :: B_WRITE_LOGICAL_R4
PUBLIC :: B_WRITE_DOUBLE_R4
PUBLIC :: B_WRITE_FLOAT_R4

PUBLIC :: B_WRITE_INTEGER_R5
PUBLIC :: B_WRITE_LOGICAL_R5
PUBLIC :: B_WRITE_DOUBLE_R5
PUBLIC :: B_WRITE_FLOAT_R5


PUBLIC :: Brook_Initialize
PUBLIC :: B_Initialize

PUBLIC :: Brook_Get_Standard_Attributes

CONTAINS

  SUBROUTINE Brook_Get_Standard_Attributes(str, name, type_code, rank, shape, C_Array_Order, map, mesh, offset)
    !
    ! Return a string formatted with the standard attributes:
    ! Name, DataType, Rank, Shape, ArrayOrder, Map, and Mesh
    implicit none
    character(len=*)                :: name
    character(len=*)                :: str
    character(len=*)                :: type_code
    integer                         :: rank
    integer, optional, dimension(:) :: shape
    logical, optional               :: C_Array_Order
    character(len=*), optional      :: map
    character(len=*), optional      :: mesh
    integer, optional, intent(in)   :: offset

    ! Local Variables
    character(len=256) :: tstring
    character(len=256) :: fstring
    character(len=1024) :: mstring
    character(len=1024) :: pstring
    integer :: i    

    fstring = 'FORTRAN'
    if ( present(C_Array_Order)) then
       if(C_Array_Order) then
          fstring = 'C'
       end if
    end if

1   FORMAT(100(i10,1x))
    if ( present(shape)) then
       write(tstring,1) shape
    else
       tstring = ''
    end if

    if(present(mesh)) then
       mstring = 'Mesh="'//TRIM(ADJUSTL(mesh))//'"'
    else
       mstring = ' '
    end if
    
    if(present(map)) then
       pstring = 'Map="'//TRIM(ADJUSTL(map))//'"'
    else
       pstring = ' '
    end if
    tstring = TRIM(ADJUSTL(tstring))
    fstring = TRIM(ADJUSTL(fstring))
    mstring = TRIM(ADJUSTL(mstring))
    pstring = TRIM(ADJUSTL(pstring))

    if ( present(offset)) then
       i = offset
    else 
       i = 0
    end if

    write(str,10) TRIM(name), TRIM(type_code), rank, TRIM(tstring), TRIM(mstring), TRIM(pstring), TRIM(fstring), i

10  FORMAT('Name="',a,'" DataType="',a1,'" Rank="',i2,'" Shape="',a,'" ',a,' ',a,' ArrayOrder="',a,'" Offset="',i20,'"')

    return
    
  END SUBROUTINE Brook_Get_Standard_Attributes
  SUBROUTINE Brook_WRITE_FILE_HEADER(B, iStatus)
    !---------------------------------------------------------------------------
    ! Purpose:
    ! Write headers so that truchas postprocessor can correctly interpret 
    ! FORTRAN files.  Currenty we only work on binary files
    !---------------------------------------------------------------------------
    implicit none
    
    type(Brook), intent(in), target :: B
    integer                         :: iStatus

    type(Brook), pointer :: P
    integer          :: ikind
    double precision :: rkind
    logical          :: lkind

    !---------------------------------------------------------------------------
    iStatus = BE_NONE

    ikind = 1
    rkind = 1.0d0
    lkind = .true.
    P => B

    if(Brook_StandardHeaders(P)) then
       IF ( iStatus == BE_NONE .and. Brook_Unit(P) > 0 ) THEN
          select case (Brook_IForm(P))
             ! Currently only xml and binary brooks have stuff written to them
          case (B_IFORM_BINARY)
             if ( iStatus == BE_NONE ) call Brook_Write(P, Variable=ikind, iStatus=iStatus)
             if ( iStatus == BE_NONE ) call Brook_Write(P, Variable=rkind, iStatus=iStatus)
             if ( iStatus == BE_NONE ) call Brook_Write(P, Variable=lkind, iStatus=iStatus)
          case (B_IFORM_XML)
10           FORMAT(a)
             WRITE(Brook_Unit(B),10) '<?xml version="1.0" ?>'
             CALL Unit_OpenXMLTag(Brook_Unit(P),XMLTag=TRIM(Brook_XMLROOT(P)), &
                  XMLAttributes=BROOK_XML_VERSION,iStatus=iStatus)
          end select
       END IF
    end if
    return

  END SUBROUTINE Brook_WRITE_FILE_HEADER

  SUBROUTINE Brook_WRITE_FILE_TRAILER(B, iStatus)
    !---------------------------------------------------------------------------
    ! Purpose:
    ! Write trailers so that truchas postprocessor can correctly interpret 
    ! FORTRAN files.  Currenty we only work on xml files
    !---------------------------------------------------------------------------
    implicit none
    
    type(Brook), intent(in), target :: B
    integer                         :: iStatus

    integer          :: ikind
    double precision :: rkind
    logical          :: lkind

    !---------------------------------------------------------------------------
    iStatus = BE_NONE

    ikind = 1
    rkind = 1.0d0
    lkind = .true.
    IF ( iStatus == BE_NONE ) THEN
       if(Brook_Unit(B) > 0 .and. Brook_StandardHeaders(B)) then
          select case (Brook_IForm(B))
             ! Currently only xml brooks have stuff written to them
          case (B_IFORM_XML)
             CALL Unit_CloseXMLTag(Brook_Unit(B), XMLTag=TRIM(Brook_XMLROOT(B)),iStatus=iStatus)
          end select
       end if
    END IF
    return

  END SUBROUTINE Brook_WRITE_FILE_TRAILER
  
  LOGICAL FUNCTION Brook_StandardHeaders(B)
    implicit none
    type(Brook), intent(in) :: B

    Brook_StandardHeaders = .true.
    if (Associated(B%Data)) then
       Brook_StandardHeaders = B%Data%StandardHeaders
    end if
    return
  END FUNCTION Brook_StandardHeaders

  SUBROUTINE Brook_Clear_Error(iStatus)
    integer, intent(inout) :: iStatus
    iStatus = BE_NONE
    BE_ID   = BE_NONE
    BE_STRING = ' ' 
    return
  END SUBROUTINE Brook_Clear_Error

  SUBROUTINE B_Set_Error(iStatus, EString)
    ! internal subroutine for setting error information
    integer, intent(in) :: iStatus
    character(len=*), intent(in) :: EString
    if ( iStatus /= BE_NONE ) then
       BE_ID = iStatus
       BE_STRING = EString
    end if
    write(*,*) 'BROOK ERROR: ', iStatus, EString
    return
  END SUBROUTINE B_Set_Error


  CHARACTER (LEN=B_STRLEN) FUNCTION BROOK_GET_ERROR(ERROR_ID)
    ! Return the error string corresponding to the error id
    INTEGER, intent(IN) :: ERROR_ID
    character(len=12) :: tmp
    
    Brook_Get_Error = 'NONE'
    if ( BE_ID /= BE_NONE) then
1      FORMAT('Brook Error: ',i2,' Explanation: ')
       write(tmp,1) BE_ID
       Brook_Get_Error = TRIM(tmp) // BE_STRING
    end if
    return
  END FUNCTION BROOK_GET_ERROR


  !=========================================
  ! Replace brook_abort() with your own abort
  SUBROUTINE Brook_Abort(s)
    character(len=*) :: s
    write(*,*) s
    stop 'Brook Panic'
  end SUBROUTINE Brook_Abort
    

  !=========================================
  ! SUBROUTINES Brook_WRITE_LOGICAL_*
  ! 
  ! Purpose: Write logical variables to output brooks
#define _TYPE_                  Logical
#define _TYPE_CODE_             'l'
#define _FORMAT_NAME_           Brook_LFormat
#define _FUNCTION_NAME_         B_WRITE_LOGICAL_R0
#define _FUNCTION_NAME_RESHAPE_ B_LOGICAL_RESHAPE
#define _NO_HIGHER_ORDER_XML_
#define _DIMENSION_
#define _DEFINE_RESHAPE_
#define _TYPE_SIZE_             BYTES_PER_LOGICAL
#include "brook_include.fpp"

#define _TYPE_                  Logical
#define _TYPE_CODE_             'l'
#define _RV_                    r1
#define _RP_                    pr1
#define _FORMAT_NAME_           Brook_LFormat
#define _FUNCTION_NAME_         B_WRITE_LOGICAL_R1
#define _FUNCTION_NAME_RESHAPE_ B_LOGICAL_RESHAPE
#define _DIMENSION_             dimension(:),
#define _TYPE_SIZE_             BYTES_PER_LOGICAL
#include "brook_include.fpp"

#define _TYPE_                  Logical
#define _TYPE_CODE_             'l'
#define _RV_                    r2
#define _RP_                    pr2
#define _FORMAT_NAME_           Brook_LFormat
#define _FUNCTION_NAME_         B_WRITE_LOGICAL_R2
#define _FUNCTION_NAME_RESHAPE_ B_LOGICAL_RESHAPE
#define _DIMENSION_             dimension(:,:),
#define _TYPE_SIZE_             BYTES_PER_LOGICAL
#include "brook_include.fpp"

#define _TYPE_                  Logical
#define _TYPE_CODE_             'l'
#define _RV_                    r3
#define _RP_                    pr3
#define _FORMAT_NAME_           Brook_LFormat
#define _FUNCTION_NAME_         B_WRITE_LOGICAL_R3
#define _FUNCTION_NAME_RESHAPE_ B_LOGICAL_RESHAPE
#define _DIMENSION_             dimension(:,:,:),
#define _TYPE_SIZE_             BYTES_PER_LOGICAL
#include "brook_include.fpp"

#define _TYPE_                  Logical
#define _TYPE_CODE_             'l'
#define _RV_                    r4
#define _RP_                    pr4
#define _FORMAT_NAME_           Brook_LFormat
#define _FUNCTION_NAME_         B_WRITE_LOGICAL_R4
#define _FUNCTION_NAME_RESHAPE_ B_LOGICAL_RESHAPE
#define _DIMENSION_             dimension(:,:,:,:),
#define _TYPE_SIZE_             BYTES_PER_LOGICAL
#include "brook_include.fpp"

#define _TYPE_                  Logical
#define _TYPE_CODE_             'l'
#define _RV_                    r5
#define _RP_                    pr5
#define _FORMAT_NAME_           Brook_LFormat
#define _FUNCTION_NAME_         B_WRITE_LOGICAL_R5
#define _FUNCTION_NAME_RESHAPE_ B_LOGICAL_RESHAPE
#define _DIMENSION_             dimension(:,:,:,:,:),
#define _TYPE_SIZE_             BYTES_PER_LOGICAL
#include "brook_include.fpp"

  ! END SUBROUTINES B_WRITE_LOGICAL_*
  !=========================================

  !=========================================
  ! SUBROUTINES B_WRITE_INTEGER_*
  ! Purpose: Write integer variables to output brooks
#define _TYPE_                  Integer
#define _TYPE_CODE_             'i'
#define _FORMAT_NAME_           Brook_IFormat
#define _FUNCTION_NAME_RESHAPE_ B_INTEGER_RESHAPE
#define _FUNCTION_NAME_         B_WRITE_INTEGER_R0
#define _NO_HIGHER_ORDER_XML_
#define _DIMENSION_
#define _DEFINE_RESHAPE_
#define _TYPE_SIZE_             BYTES_PER_INT
#include "brook_include.fpp"

#define _TYPE_                  Integer
#define _TYPE_CODE_             'i'
#define _RV_                    r1
#define _RP_                    pr1
#define _FORMAT_NAME_           Brook_IFormat
#define _FUNCTION_NAME_RESHAPE_ B_INTEGER_RESHAPE
#define _FUNCTION_NAME_         B_WRITE_INTEGER_R1
#define _DIMENSION_             dimension(:),
#define _TYPE_SIZE_             BYTES_PER_INT
#include "brook_include.fpp"

#define _TYPE_                  Integer
#define _TYPE_CODE_             'i'
#define _RV_                    r2
#define _RP_                    pr2
#define _FORMAT_NAME_           Brook_IFormat
#define _FUNCTION_NAME_         B_WRITE_INTEGER_R2
#define _FUNCTION_NAME_RESHAPE_ B_INTEGER_RESHAPE
#define _DIMENSION_             dimension(:,:),
#define _TYPE_SIZE_             BYTES_PER_INT
#include "brook_include.fpp"

#define _TYPE_                  Integer
#define _TYPE_CODE_             'i'
#define _RV_                    r3
#define _RP_                    pr3
#define _FORMAT_NAME_           Brook_IFormat
#define _FUNCTION_NAME_         B_WRITE_INTEGER_R3
#define _FUNCTION_NAME_RESHAPE_ B_INTEGER_RESHAPE
#define _DIMENSION_             dimension(:,:,:),
#define _TYPE_SIZE_             BYTES_PER_INT
#include "brook_include.fpp"

#define _TYPE_                  Integer
#define _TYPE_CODE_             'i'
#define _RV_                    r4
#define _RP_                    pr4
#define _FORMAT_NAME_           Brook_IFormat
#define _FUNCTION_NAME_         B_WRITE_INTEGER_R4
#define _FUNCTION_NAME_RESHAPE_ B_INTEGER_RESHAPE
#define _DIMENSION_             dimension(:,:,:,:),
#define _TYPE_SIZE_             BYTES_PER_INT
#include "brook_include.fpp"

#define _TYPE_                  Integer
#define _TYPE_CODE_             'i'
#define _RV_                    r5
#define _RP_                    pr5
#define _FORMAT_NAME_           Brook_IFormat
#define _FUNCTION_NAME_         B_WRITE_INTEGER_R5
#define _FUNCTION_NAME_RESHAPE_ B_INTEGER_RESHAPE
#define _DIMENSION_             dimension(:,:,:,:,:),
#define _TYPE_SIZE_             BYTES_PER_INT
#include "brook_include.fpp"
  ! END SUBROUTINES B_WRITE_INTEGER_*
  !=========================================

  !=========================================
  ! SUBROUTINES B_WRITE_DOUBLE_*
  ! Purpose: Write double precision variables to output brooks
#define _TYPE_                  double precision
#define _TYPE_CODE_             'd'
#define _FORMAT_NAME_           Brook_DFormat
#define _FUNCTION_NAME_         B_WRITE_DOUBLE_R0
#define _FUNCTION_NAME_RESHAPE_ B_DOUBLE_RESHAPE
#define _NO_HIGHER_ORDER_XML_
#define _DIMENSION_
#define _DEFINE_RESHAPE_
#define _TYPE_SIZE_             BYTES_PER_DOUBLE
#include "brook_include.fpp"

#define _TYPE_                  double precision
#define _TYPE_CODE_             'd'
#define _RV_                    r1
#define _RP_                    pr1
#define _FORMAT_NAME_           Brook_DFormat
#define _FUNCTION_NAME_RESHAPE_ B_DOUBLE_RESHAPE
#define _FUNCTION_NAME_         B_WRITE_DOUBLE_R1
#define _DIMENSION_             dimension(:),
#define _TYPE_SIZE_             BYTES_PER_DOUBLE
#include "brook_include.fpp"

#define _TYPE_                  double precision
#define _TYPE_CODE_             'd'
#define _RV_                    r2
#define _RP_                    pr2
#define _FORMAT_NAME_           Brook_DFormat
#define _FUNCTION_NAME_         B_WRITE_DOUBLE_R2
#define _FUNCTION_NAME_RESHAPE_ B_DOUBLE_RESHAPE
#define _DIMENSION_             dimension(:,:),
#define _TYPE_SIZE_             BYTES_PER_DOUBLE
#include "brook_include.fpp"

#define _TYPE_                  double precision
#define _TYPE_CODE_             'd'
#define _RV_                    r3
#define _RP_                    pr3
#define _FORMAT_NAME_           Brook_DFormat
#define _FUNCTION_NAME_         B_WRITE_DOUBLE_R3
#define _FUNCTION_NAME_RESHAPE_ B_DOUBLE_RESHAPE
#define _DIMENSION_             dimension(:,:,:),
#define _TYPE_SIZE_             BYTES_PER_DOUBLE
#include "brook_include.fpp"

#define _TYPE_                  double precision
#define _TYPE_CODE_             'd'
#define _RV_                    r4
#define _RP_                    pr4
#define _FORMAT_NAME_           Brook_DFormat
#define _FUNCTION_NAME_         B_WRITE_DOUBLE_R4
#define _FUNCTION_NAME_RESHAPE_ B_DOUBLE_RESHAPE
#define _DIMENSION_             dimension(:,:,:,:),
#define _TYPE_SIZE_             BYTES_PER_DOUBLE
#include "brook_include.fpp"

#define _TYPE_                  double precision
#define _TYPE_CODE_             'd'
#define _RV_                    r5
#define _RP_                    pr5
#define _FORMAT_NAME_           Brook_DFormat
#define _FUNCTION_NAME_         B_WRITE_DOUBLE_R5
#define _FUNCTION_NAME_RESHAPE_ B_DOUBLE_RESHAPE
#define _DIMENSION_             dimension(:,:,:,:,:),
#define _TYPE_SIZE_             BYTES_PER_DOUBLE
#include "brook_include.fpp"
  ! END SUBROUTINES B_WRITE_DOUBLE_*
  !=========================================
  !=========================================
  ! SUBROUTINES B_WRITE_FLOAT_*
  ! Purpose: Write real variables to output brooks
#define _TYPE_                  real
#define _TYPE_CODE_             'f'
#define _FORMAT_NAME_           Brook_DFormat
#define _FUNCTION_NAME_         B_WRITE_FLOAT_R0
#define _FUNCTION_NAME_RESHAPE_ B_FLOAT_RESHAPE
#define _NO_HIGHER_ORDER_XML_
#define _DIMENSION_
#define _DEFINE_RESHAPE_
#define _TYPE_SIZE_             BYTES_PER_FLOAT
#include "brook_include.fpp"

#define _TYPE_                  real
#define _TYPE_CODE_             'f'
#define _RV_                    r1
#define _RP_                    pr1
#define _FORMAT_NAME_           Brook_DFormat
#define _FUNCTION_NAME_RESHAPE_ B_FLOAT_RESHAPE
#define _FUNCTION_NAME_         B_WRITE_FLOAT_R1
#define _DIMENSION_             dimension(:),
#define _TYPE_SIZE_             BYTES_PER_FLOAT
#include "brook_include.fpp"

#define _TYPE_                  real
#define _TYPE_CODE_             'f'
#define _RV_                    r2
#define _RP_                    pr2
#define _FORMAT_NAME_           Brook_DFormat
#define _FUNCTION_NAME_         B_WRITE_FLOAT_R2
#define _FUNCTION_NAME_RESHAPE_ B_FLOAT_RESHAPE
#define _DIMENSION_             dimension(:,:),
#define _TYPE_SIZE_             BYTES_PER_FLOAT
#include "brook_include.fpp"

#define _TYPE_                  real
#define _TYPE_CODE_             'f'
#define _RV_                    r3
#define _RP_                    pr3
#define _FORMAT_NAME_           Brook_DFormat
#define _FUNCTION_NAME_         B_WRITE_FLOAT_R3
#define _FUNCTION_NAME_RESHAPE_ B_FLOAT_RESHAPE
#define _DIMENSION_             dimension(:,:,:),
#define _TYPE_SIZE_             BYTES_PER_FLOAT
#include "brook_include.fpp"

#define _TYPE_                  real
#define _TYPE_CODE_             'f'
#define _RV_                    r4
#define _RP_                    pr4
#define _FORMAT_NAME_           Brook_DFormat
#define _FUNCTION_NAME_         B_WRITE_FLOAT_R4
#define _FUNCTION_NAME_RESHAPE_ B_FLOAT_RESHAPE
#define _DIMENSION_             dimension(:,:,:,:),
#define _TYPE_SIZE_             BYTES_PER_FLOAT
#include "brook_include.fpp"

#define _TYPE_                  real
#define _TYPE_CODE_             'f'
#define _RV_                    r5
#define _RP_                    pr5
#define _FORMAT_NAME_           Brook_DFormat
#define _FUNCTION_NAME_         B_WRITE_FLOAT_R5
#define _FUNCTION_NAME_RESHAPE_ B_FLOAT_RESHAPE
#define _DIMENSION_             dimension(:,:,:,:,:),
#define _TYPE_SIZE_             BYTES_PER_FLOAT
#include "brook_include.fpp"
  ! END SUBROUTINES B_WRITE_FLOAT_*
  !=========================================

  !=========================================
  ! SUBROUTINES B_WRITE_STRING_*
  ! Purpose: Write string variables to output brooks
  ! Note only single arrays are allowed
#define _TYPE_                  Character(len=*)
#define _TYPE_CODE_             'c'
#define _FORMAT_NAME_           Brook_AFormat
#define _FUNCTION_NAME_         B_WRITE_STRING_R0
#define _IS_CHARACTER_TYPE_
#define _CHARLEN_               len(Variable)
#define _NO_HIGHER_ORDER_XML_
#define _DIMENSION_
#define _DEFINE_RESHAPE_
#define _FUNCTION_NAME_RESHAPE_ B_STRING_RESHAPE
#define _TYPE_SIZE_             BYTES_PER_CHAR
#include "brook_include.fpp"

#define _TYPE_                  Character(len=*)
#define _TYPE_CODE_             'c'
#define _RV_                    r1
#define _RP_                    pr1
#define _FORMAT_NAME_           Brook_AFormat
#define _FUNCTION_NAME_         B_WRITE_STRING_R1
#define _IS_CHARACTER_TYPE_
#define _CHARLEN_               len(Variable(1))
#define _NO_HIGHER_ORDER_XML_
#define _FUNCTION_NAME_RESHAPE_ B_STRING_RESHAPE
#define _DIMENSION_             dimension(:), 
#define _TYPE_SIZE_             BYTES_PER_CHAR
#include "brook_include.fpp"

  ! END SUBROUTINES B_WRITE_CHARACTER_*
  !=========================================

  ! Access Subroutines

  FUNCTION  Brook_Next(B) RESULT(p)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    TYPE(Brook), POINTER    :: P
    integer                 :: iStatus 
    iStatus = BE_NONE
    p=>NULL()
    IF ( .NOT. B_Initialized) THEN
       iStatus = BE_NONE
       CALL B_Initialize(iStatus)
    END IF
    P=>B%Next
    RETURN
  END FUNCTION Brook_Next

  CHARACTER(Len=B_Strlen) FUNCTION  Brook_DFormat(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_DFormat = TRIM(B%Data%DFormat)
    ELSE 
       Brook_DFormat = TRIM(B_DFormat_default)
    END IF
    RETURN
  END FUNCTION Brook_DFormat

  CHARACTER(Len=B_Strlen) FUNCTION  Brook_LFormat(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_LFormat = TRIM(B%Data%LFormat)
    ELSE 
       Brook_Lformat = TRIM(B_LFormat_default)
    END IF
    RETURN
  END FUNCTION Brook_LFormat

  CHARACTER(Len=B_Strlen) FUNCTION  Brook_AFormat(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    Brook_AFormat = TRIM(B_AFormat_default)
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_AFormat = TRIM(B%Data%AFormat)
    END IF
    RETURN
  END FUNCTION Brook_AFormat

  CHARACTER(Len=B_Strlen) FUNCTION  Brook_XMLRoot(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_XMLRoot = TRIM(B%Data%XMLRoot)
    ELSE 
       Brook_XMLRoot = TRIM(B_XMLRoot_default)
    END IF
    RETURN
  END FUNCTION Brook_XMLRoot

  CHARACTER(Len=B_Strlen) FUNCTION  Brook_IFormat(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_IFormat = TRIM(B%Data%IFormat)
    ELSE 
       Brook_IFormat = TRIM(B_IFormat_default)
    END IF
    RETURN
  END FUNCTION Brook_IFormat


  LOGICAL FUNCTION  Brook_Closeable(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_Closeable = B%Data%Closeable
    ELSE 
       Brook_Closeable = .FALSE.
    END IF
    RETURN
  END FUNCTION Brook_Closeable
  LOGICAL FUNCTION  Brook_Advance(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_Advance = B%Data%Advance
    ELSE 
       Brook_Advance = .TRUE.
    END IF
    RETURN
  END FUNCTION Brook_Advance

  INTEGER FUNCTION  Brook_UNIT(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_Unit = B%Data%Unit
       if (B%Data%Unit <= 0 .and. B%Data%Status  .ne. B_Status_no_init) then
          ! Inquire by file name if invalid unit
          INQUIRE(FILE=trim(BROOK_FILE(B)), number=Brook_Unit)
       end if
    ELSE 
       Brook_Unit = B_UNIT_UNKNOWN
    END IF
    RETURN
  END FUNCTION Brook_UNIT

  INTEGER FUNCTION  Brook_OFFSET(B)
    ! Access function for the data
    TYPE(Brook), target, INTENT(IN) :: B
    integer                 :: iStatus = 0
    Type(Brook), pointer :: P
    P=>B
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    Brook_Offset = 0
    IF ( ASSOCIATED(B%Data) ) THEN
       IF(BROOK_IFORM(B) == B_IFORM_BINARY) THEN
          Brook_Offset = B%Data%Offset
          if(Brook_Offset == 0 ) then
             call B_Open_TopLevel(P, iStatus=iStatus)
             Brook_Offset = B%Data%Offset
          end if
       END IF
    END IF
    RETURN
  END FUNCTION Brook_OFFSET

  INTEGER FUNCTION  Brook_STATUS(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_Status = B%Data%Status
    ELSE 
       Brook_Status = B_STATUS_UNKNOWN
    END IF
    RETURN
  END FUNCTION Brook_STATUS

  INTEGER FUNCTION  Brook_IFORM(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_IForm = B%Data%IForm
    ELSE 
       Brook_IForm = B_IFORM_UNKNOWN
    END IF
    RETURN
  END FUNCTION Brook_IFORM

  INTEGER FUNCTION  Brook_ID(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_Id = B%Data%Id
    ELSE 
       Brook_Id = B_ID_UNKNOWN
    END IF
    RETURN
  END FUNCTION Brook_ID


  CHARACTER(Len=B_Strlen) FUNCTION  Brook_Position(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_Position = TRIM(B%Data%Position)
    ELSE 
       Brook_Position = TRIM(B_Position_Default)
    END IF
    RETURN
  END FUNCTION Brook_Position

  CHARACTER(Len=B_Strlen) FUNCTION  Brook_Form(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_Form = TRIM(B%Data%Form)
    ELSE 
       Brook_Form = TRIM(B_Form_UNKNOWN)
    END IF
    RETURN
  END FUNCTION Brook_Form

  CHARACTER(Len=B_Strlen) FUNCTION  Brook_File(B)
    ! Access function for the data
    TYPE(Brook), INTENT(IN) :: B
    integer                 :: iStatus = 0
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
    END IF
    IF ( ASSOCIATED(B%Data) ) THEN
       Brook_File = TRIM(B%Data%File)
    ELSE 
       Brook_File = TRIM(B_FILE_UNKNOWN)
    END IF
    RETURN
  END FUNCTION Brook_File

  SUBROUTINE B_ASSIGN(Olhs, Orhs)
    TYPE(Brook), TARGET, INTENT(INOUT) :: Olhs
    TYPE(Brook), TARGET, INTENT(IN)    :: Orhs

    ! Local Variables
    TYPE(Brook), POINTER :: pr
    TYPE(Brook), POINTER :: pl
    TYPE(BROOK), TARGET  :: btmp
    INTEGER                 :: myStatus

    myStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(myStatus)
    END IF
    

    NULLIFY(btmp%Next)
    NULLIFY(btmp%Data)
    if (myStatus /= BE_NONE) return
       
    pl=>btmp
    
    ! Now copy the top level data
    pr=>Orhs
    CALL Brook_Initialize(B=pl, DATA=pr%Data, iStatus = myStatus)
    if (myStatus /= BE_NONE ) return
    
    ! Finally loop through the right hand side tail
    prLoop: DO
       pr=>Brook_Next(pr)
       IF (  .NOT. ASSOCIATED(pr) ) THEN
          EXIT prLoop
       END IF
       ALLOCATE(pl%Next, STAT=myStatus)
       IF ( myStatus /= BE_NONE ) then
          call B_Set_Error(BE_ASSIGN, 'Allocation Error in assign Brook = Brook')
          EXIT prLoop
       END IF
       
       pl=>Brook_Next(pl)
       CALL Brook_Initialize(B=pl, DATA=pr%Data, iStatus=myStatus)
!       write(*,*) 'assign unit: ', Brook_Unit(pr), Brook_Unit(pl)
       IF ( myStatus /= BE_NONE ) EXIT prLoop
       
    END DO prLoop
       
    IF ( myStatus /= BE_NONE ) return
    
    ! Destroy the lhs struct
    pl=>Olhs
    CALL Brook_Destroy(pl, iStatus=myStatus)
    IF ( myStatus /= BE_NONE ) return
    pl%Data=>btmp%Data
    pl%Next=>Brook_Next(btmp)
    
    btmp%Data=>Null()
    btmp%Next=>Null()
    
    RETURN
    
  END SUBROUTINE B_ASSIGN


  SUBROUTINE B_ALLOCATE_STRUCT_POINTER(P, S, iStatus)
    TYPE(BROOK_POINTER) :: P
    TYPE(BROOK), TARGET, OPTIONAL :: S
    INTEGER :: iStatus
    
    iStatus = 0
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       if ( iStatus /= BE_NONE) return
    END IF

    IF ( PRESENT(S) ) THEN
       IF ( ASSOCIATED(S%Data)) THEN
          ALLOCATE(p%S, STAT=iStatus)
          IF ( iStatus /= BE_NONE ) then
             call B_Set_Error(BE_ALLOCATE, 'Allocation error 1 in Brook allocate struct pointer')
             return
          end IF
          ALLOCATE(p%S%Data, STAT=iStatus)
          IF ( iStatus /= BE_NONE ) then
             DEALLOCATE(p%S)
             p%S=>NULL()
             call B_Set_Error(BE_ALLOCATE, 'Allocation error 2 in Brook allocate struct pointer')
             return
          end IF
          p%S%Data=S%Data
       ELSE 
          p%S => NULL()
       END IF
    ELSE 
       p%S => NULL()
    END IF
    p%Next => NULL()
    RETURN
  END SUBROUTINE B_ALLOCATE_STRUCT_POINTER

  SUBROUTINE B_DESTROY_STRUCT_POINTER(P, iStatus)
    TYPE(BROOK_POINTER), TARGET  :: P
    TYPE(BROOK_POINTER), POINTER :: Q
    TYPE(BROOK_POINTER), POINTER :: R
    integer                      :: iStatus
    LOGICAL :: first

    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       if ( iStatus /= BE_NONE) return
    END IF

    first = .TRUE.
    Q=>P
    DO
       IF ( .NOT. ASSOCIATED(Q)) EXIT
       R=>Q%Next
       IF (ASSOCIATED(Q%S)) THEN
          IF (ASSOCIATED(Q%S%Data)) THEN
             DEALLOCATE(Q%S%Data)
          END IF
          DEALLOCATE(Q%S)
          NULLIFY(Q%S)
       END IF
       NULLIFY(Q%Next)
       IF ( first ) THEN
          first = .FALSE.
       ELSE 
          DEALLOCATE(Q)
       END IF
       Q=>R
    END DO
    RETURN
  END SUBROUTINE B_DESTROY_STRUCT_POINTER
    

  SUBROUTINE B_ASSIGN_STRUCT_POINTER (OLhs, PRhs)
    TYPE(Brook),         TARGET, INTENT(INOUT) :: OLhs
    TYPE(BROOK_POINTER), TARGET, INTENT(IN)    :: PRHS

    TYPE(BROOK_POINTER), POINTER :: P
    TYPE(BROOK),         POINTER :: q

    INTEGER :: myStatus

    myStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(myStatus)
       if (myStatus /= BE_NONE) return
    END IF

    P=>PRhs
    q=>OLhs

    ! Destroy any data initially in OLhs
    CALL Brook_Destroy(q, iStatus=myStatus)
    if ( myStatus /= BE_NONE) return
    ! First the top level struct
    IF ( ASSOCIATED(P%S)) THEN
       CALL BROOK_INITIALIZE(B=q,DATA=P%S%Data, iStatus=myStatus)
    ELSE 
       CALL BROOK_INITIALIZE(B=q, iStatus=myStatus)
    END IF
    if ( myStatus /= BE_NONE) return

    P=>P%Next
    DO
       IF ( .NOT. ASSOCIATED(P) ) EXIT
       ALLOCATE(q%Next,stat=myStatus)
       IF ( myStatus /= BE_NONE ) then
          call B_Set_Error(BE_ASSIGN, 'allocation error 3 in assign struct pointer')
          return
       end IF
       q=>Brook_Next(q)
       IF ( ASSOCIATED(P%S)) THEN
          CALL BROOK_INITIALIZE(B=q,DATA=P%S%Data, iStatus=myStatus)
       ELSE 
          CALL BROOK_INITIALIZE(B=q, iStatus=myStatus)
       END IF
       IF ( myStatus /= BE_NONE ) return
       P=>P%Next
    END DO
          
    CALL B_Destroy_Struct_pointer(PRhs, iStatus=myStatus)
    RETURN
  END SUBROUTINE B_ASSIGN_STRUCT_POINTER


  FUNCTION B_ADD_S_S(Olhs, Orhs) RESULT(r)
    TYPE(Brook), TARGET, INTENT(IN) :: Olhs
    TYPE(Brook), TARGET, INTENT(IN) :: Orhs
    TYPE(Brook), POINTER :: pl
    TYPE(BROOK_POINTER), POINTER :: ptr2
    TYPE(BROOK_POINTER), TARGET  :: r
    INTEGER                          :: myStatus

    myStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(myStatus)
       if ( myStatus /= BE_NONE) return
    END IF

    NULLIFY(r%S)
    NULLIFY(r%Next)
    ptr2=>r

    ! Copy Top Level LHS to ptr
    pl=>Olhs
    CALL b_Allocate_Struct_Pointer(ptr2,pl, iStatus=myStatus)
    if ( myStatus /= BE_NONE) return

    ! Copy LHS Tail to ptr2
    plLoop: DO
       pl=>Brook_Next(pl)
       NULLIFY(ptr2%Next)
       IF ( .NOT. ASSOCIATED(pl)) THEN
          EXIT plLoop
       END IF
       ALLOCATE(ptr2%Next, STAT=myStatus)
       IF ( myStatus /= BE_NONE ) then
          call B_Set_Error(BE_ADD, 'allocation error 5 in Brook_add')
          return
       end IF
       CALL B_ALLOCATE_STRUCT_POINTER(ptr2%Next, pl, iStatus=myStatus)
       IF ( myStatus /= 0 ) return
       ptr2=>ptr2%Next
    END DO plLoop

    ! Now copy RHS to ptr2
    pl=>Orhs
    prLoop: DO
       NULLIFY(ptr2%Next)
       IF ( .NOT. ASSOCIATED(pl)) THEN
          EXIT prLoop
       END IF
       ALLOCATE(ptr2%Next, STAT=myStatus)
       IF ( myStatus /= BE_NONE ) then
          call B_Set_Error(BE_ADD, 'allocation error 6 in Brook_add')
          return
       end IF
       CALL B_ALLOCATE_STRUCT_POINTER(ptr2%Next, pl, iStatus=myStatus)
       IF ( myStatus /= BE_NONE ) then
          call B_Set_Error(BE_ADD, 'allocation error 7 in Brook_add')
          return
       end IF
       ptr2=>ptr2%Next
       pl=>Brook_Next(pl)
    END DO prLoop

    RETURN
  END FUNCTION B_ADD_S_S

  FUNCTION B_ADD_P_S(Plhs, Orhs) RESULT(r)
    TYPE(Brook_POINTER), TARGET, INTENT(IN) :: Plhs
    TYPE(Brook), TARGET, INTENT(IN) :: Orhs
    TYPE(Brook), POINTER :: pl
    TYPE(BROOK_POINTER), POINTER :: ptr2
    TYPE(BROOK_POINTER), TARGET  :: r
    INTEGER                          :: myStatus

    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(myStatus)
       if ( myStatus /= BE_NONE) return
    END IF

    NULLIFY(r%S)
    NULLIFY(r%Next)
    myStatus = 0
    ptr2=>r

    ! Copy LHS to ptr
    ptr2%S=>PLhs%S
    ptr2%Next =>PLhs%Next

    ! Find bottom of PLHS
    plLoop: DO
       IF ( .NOT. ASSOCIATED(ptr2%Next)) THEN
          EXIT plLoop
       END IF
       ptr2=>ptr2%Next
    END DO plLoop

    ! Now copy RHS to ptr2
    pl=>Orhs
    prLoop: DO
       NULLIFY(ptr2%Next)
       IF ( .NOT. ASSOCIATED(pl)) THEN
          EXIT prLoop
       END IF
       ALLOCATE(ptr2%Next, STAT=myStatus)
       IF ( myStatus /= BE_NONE ) then
          call B_Set_Error(BE_ADD, 'allocation error 8 in Brook_add')
          return
       end IF
       CALL B_ALLOCATE_STRUCT_POINTER(ptr2%Next, pl,iStatus=myStatus)
       IF ( myStatus /= BE_NONE ) then
          call B_Set_Error(BE_ADD, 'allocation error 9 in Brook_add')
          return
       end IF
       ptr2=>ptr2%Next
       pl=>Brook_Next(pl)
    END DO prLoop

    RETURN
  END FUNCTION B_ADD_P_S

  FUNCTION B_ADD_S_P(OLhs, PRhs) RESULT(r)
    TYPE(Brook), TARGET, INTENT(IN) :: OLhs
    TYPE(Brook_POINTER), TARGET, INTENT(IN) :: PRhs
    TYPE(Brook), POINTER :: pl
    TYPE(BROOK_POINTER), POINTER :: ptr2
    TYPE(BROOK_POINTER), TARGET  :: r
    INTEGER                         :: myStatus

    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(myStatus)
       if ( myStatus /= BE_NONE) return
    END IF
    NULLIFY(r%S)
    NULLIFY(r%Next)
    myStatus = 0
    ptr2=>r

    ! Copy Top Level LHS to ptr
    pl=>OLhs
    CALL B_Allocate_Struct_Pointer(ptr2,pl,iStatus=myStatus)
    IF ( myStatus /= BE_NONE) return
    
    ! Now the tail of the LHS
    pl=>Brook_Next(pl)
    plLoop: DO
       IF ( .NOT. ASSOCIATED(pl)) EXIT plLoop
       ALLOCATE(ptr2%Next, STAT=myStatus)
       IF ( myStatus /= BE_NONE ) then
          call B_Set_Error(BE_ADD, 'allocation error 10 in Brook_add')
          return
       end IF
       CALL B_Allocate_Struct_Pointer(ptr2%Next, pl, iStatus=myStatus)
       IF ( myStatus /= BE_NONE ) then
          call B_Set_Error(BE_ADD, 'allocation error 11 in Brook_add')
          return
       end IF
       ptr2=>ptr2%Next
       pl=>Brook_Next(pl)
    END DO plLoop


    ! Copy RHS Tail to ptr2
    ALLOCATE(ptr2%Next, STAT=myStatus)
    IF ( myStatus /= BE_NONE ) then
       call B_Set_Error(BE_ADD, 'allocation error 12 in Brook_add')
       return
    end IF
    ptr2=>ptr2%Next
    ptr2%S=>PRhs%S
    ptr2%Next=>PRhs%Next

    RETURN
  END FUNCTION B_ADD_S_P

  
  FUNCTION B_ADD_P_P(PLhs, PRhs) RESULT(r)
    TYPE(Brook_POINTER), TARGET, INTENT(IN) :: Plhs
    TYPE(Brook_POINTER), TARGET, INTENT(IN) :: Prhs
    TYPE(BROOK_POINTER), POINTER :: ptr2
    TYPE(BROOK_POINTER), TARGET  :: r
    INTEGER                         :: myStatus

    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(myStatus)
       if ( myStatus /= BE_NONE ) return
    END IF
    r%S=>PLhs%S
    r%Next=>PLhs%Next

    ! Now the tail of the LHS
    ptr2=>r
    plLoop: DO
       IF ( .NOT. ASSOCIATED(ptr2%Next)) EXIT plLoop
       ptr2=>ptr2%Next
    END DO plLoop

    ! Copy RHS to ptr2
    ALLOCATE(ptr2%Next, STAT=myStatus)
    IF ( myStatus /= BE_NONE ) then
       call B_Set_Error(BE_ADD, 'allocation error 12 in Brook_add')
       return
    end IF
    ptr2=>ptr2%Next
    ptr2%S=>PRhs%S
    ptr2%Next=>PRhs%Next

    RETURN
  END FUNCTION B_ADD_P_P

  
  SUBROUTINE Brook_Destroy(B, iStatus)
    ! Destroy any data associated with B and its tail
    ! 
    ! B_stdout and B_stderr will not be touched
    IMPLICIT NONE
    TYPE(Brook), TARGET,   INTENT(INOUT) :: B
    INTEGER,               INTENT(INOUT) :: iStatus
    
    LOGICAL :: first

    ! Local Variables
    TYPE(Brook),  POINTER :: ptr
    TYPE(Brook),  POINTER :: ptrNext

    iStatus = BE_NONE

    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       if ( iStatus /= BE_NONE ) return
    END IF

    first = .TRUE.
    ptrNext=>B
    PTRLOOP: DO
       IF ( .NOT. ASSOCIATED(ptrNext) ) THEN
          EXIT PTRLOOP
       END IF
       ptr => ptrNext
       ptrNext=>Brook_Next(ptr)
       IF (ASSOCIATED(ptr%Data)) THEN
          DEALLOCATE(ptr%Data)
       END IF
       ptr%Data => Null()
       ptr%Next => NULL()
       IF ( .NOT. first ) THEN
          ! In this case we will deallocate the 
          ! data pointers
          DEALLOCATE(ptr)
       ELSE 
          ! Do not deallocate the pointer
          first = .FALSE.
       END IF
    END DO PTRLOOP
  END SUBROUTINE Brook_Destroy

  SUBROUTINE B_Close_TopLevel(B, iStatus)
    ! B_stdout and B_stderr will not be touched
    IMPLICIT NONE
    TYPE(Brook), TARGET,   INTENT(IN)    :: B
    INTEGER,               INTENT(INOUT) :: iStatus
    
    ! Local Variables
    LOGICAL                   :: OpenP
    INTEGER                   :: iUnit
    character(LEN=B_STRLEN)   :: file
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN

       CALL B_Initialize(iStatus)
       ! Clearly there is nothing to be done if 
       ! we are not even initialized
       RETURN

    END IF

    ! Close XML Tags
    ! This is done in a loop because the root may be different for different files
    IF (ASSOCIATED(B%Data)) THEN
       iUnit = Brook_Unit(B)
       IF ( (iUnit /= B_STDOUT_UNIT)  .AND. &
            (iUnit /= B_STDERR_UNIT)  .AND. &
            (Brook_Closeable(B)) ) THEN
          openP = .FALSE.
          file='invalid b_unit'
          if ( iUnit .gt. 0 ) then
             INQUIRE(UNIT = iUnit, NAME=file, OPENED = openP)
          else 
             INQUIRE(FILE=trim(BROOK_FILE(B)), number=iUnit, OPENED = openP)
          end if
#ifdef BROOK_DEBUG_UNIT          
          write(*,*) 'CLOSETOP: NAME:',TRIM(BROOK_FILE(B)),' uFIle=',TRIM(file), ' iUnit=',iUnit
#endif
          IF ( openP ) THEN
             IF (Brook_IFORM(B) == B_IFORM_XML) THEN
                call Brook_WRITE_FILE_TRAILER(B, iStatus) 
                if (iStatus /= BE_NONE ) return
             END IF
             CLOSE(iUnit)
          END IF
          call Brook_ReleaseUnit(iUnit, iStatus)
       END IF
    END IF

    RETURN
  END SUBROUTINE B_Close_TopLevel

  RECURSIVE SUBROUTINE Brook_Close(B, iStatus)
    ! All files will be closed
    ! B_stdout and B_stderr will not be touched
    IMPLICIT NONE
    TYPE(Brook), TARGET,   INTENT(INOUT) :: B
    INTEGER,               INTENT(INOUT) :: iStatus
    
    ! Local Variables
    TYPE(Brook),  POINTER :: ptr

    iStatus = BE_NONE

    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       ! Clearly there is nothing to be done if 
       ! we are not even initialized
       RETURN    
    END IF

    ptr=>B
    PTRLOOP: DO
       IF ( .NOT. ASSOCIATED(ptr) ) EXIT PTRLOOP
       CALL B_Close_TopLevel(ptr,iStatus)
       IF (iStatus /= BE_NONE ) THEN
          RETURN
       END IF
       ptr =>Brook_Next( ptr)
    END DO PTRLOOP
    RETURN
  END SUBROUTINE Brook_Close

  SUBROUTINE BROOK_WRITEXMLDATA(B,             & ! Brook, target
                                XMLTag,        & ! Name of tag to be opened
                                XMLAttributes, & ! Attributes for tag to be opened
                                XMLStringData, & ! String with data to be written inline
                                iStatus        & ! Did an error occur?
                                )
    TYPE(BROOK),        TARGET, INTENT(IN)  :: B
    CHARACTER(LEN=*),           INTENT(IN)  :: XMLTag
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: XMLAttributes
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: XMLStringData
    INTEGER,                    INTENT(OUT) :: iStatus
    call Brook_WriteXMLTag(B,XMLTag,XMLAttributes,XMLStringData,iStatus)
  END SUBROUTINE BROOK_WRITEXMLDATA

  SUBROUTINE BROOK_WRITEXMLTAG(B,             & ! Brook, target
                               XMLTag,        & ! Name of tag to be opened
                               XMLAttributes, & ! Attributes for tag to be opened
                               XMLStringData, & ! String with data to be written inline
                               iStatus        & ! Did an error occur?
                               )
    TYPE(BROOK),        TARGET, INTENT(IN)  :: B
    CHARACTER(LEN=*),           INTENT(IN)  :: XMLTag
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: XMLAttributes
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: XMLStringData
    INTEGER,                    INTENT(OUT) :: iStatus

    TYPE(BROOK),      POINTER :: ptr
    INTEGER                   :: iUnit
    INTEGER                   :: iForm
    CHARACTER(LEN=3)          :: ADVANCE_STR 

    iStatus = BE_NONE

    ADVANCE_STR = 'no'
    if (present(XMLStringData)) then
       if(LEN_TRIM(XMLStringData)>0) then
          ADVANCE_STR = 'yes'
       end if
    end if

    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       ! Clearly there is nothing to be done if 
       ! we are not even initialized
       if ( iStatus /= BE_NONE ) return
    END IF

    ! First open the tag in all brooks
    CALL BROOK_OpenXMLTag(B,XMLTag,XMLAttributes,ADVANCE_STR=ADVANCE_STR, iStatus=iStatus)
    if ( iStatus /= BE_NONE ) return
    ! Now write to all ascii and xml brooks
    if (present(XMLStringData)) then
       if(LEN_TRIM(XMLStringData)>0) then
          ptr => B
          WRITELOOP: DO
             IF ( .NOT. ASSOCIATED (ptr) ) THEN
                EXIT WRITELOOP
             END IF
             iUnit = Brook_UNIT(ptr)
             IF ( iUnit /=  B_INVALID_UNIT) THEN
                iForm = Brook_IForm(ptr)
                IF ( iForm == B_IFORM_XML .OR. iForm == B_IFORM_ASCII) THEN
                   WRITE(iUnit,'(a)',ADVANCE=ADVANCE_STR) XMLStringData
                END IF
             END IF
             ptr => Brook_Next(ptr)
          END DO WRITELOOP
       end if
    end if
    ! Finally Close the tag in all brooks and exit
    CALL BROOK_CloseXMLTag(B,XMLTag, iStatus)
    RETURN
    
  END SUBROUTINE BROOK_WRITEXMLTAG

  SUBROUTINE Brook_OpenXMLTag(B, XMLTag, XMLAttributes, ADVANCE_STR, iStatus)
    TYPE(BROOK),        TARGET, INTENT(IN) :: B
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: XMLTag
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: XMLAttributes
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: ADVANCE_STR
    INTEGER,                    INTENT(INOUT) :: iStatus

    ! Local Variables
    TYPE(BROOK), POINTER :: ptr
    INTEGER              :: iUnit
    INTEGER              :: iForm

    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       ! Clearly there is nothing to be done if 
       ! we are not even initialized
       if ( iStatus /= BE_NONE ) return
    END IF
    ptr => B
    B_Loop: DO
       IF ( .NOT. ASSOCIATED(ptr) .OR. (iStatus /= 0) ) EXIT B_Loop
       iForm = Brook_IForm(ptr)
       IF ( iForm == B_IFORM_XML .OR. iForm == B_IFORM_ASCII) THEN
          CALL Brook_Open(ptr,iStatus=iStatus)
          if ( iStatus /= BE_NONE) return
          iUnit = Brook_Unit(ptr)
          CALL Unit_OpenXMLTag(iUnit, XMLTag=XMLTag, XMLAttributes=XMLAttributes, ADVANCE_STR=ADVANCE_STR, iStatus=iStatus)
          if ( iStatus /= BE_NONE) return
       END IF
       ptr =>Brook_Next( ptr)
    END DO B_Loop
    RETURN
  END SUBROUTINE Brook_OpenXMLTag

  SUBROUTINE Brook_CloseXMLTag(B, XMLTag, iStatus)
    TYPE(BROOK),        TARGET, INTENT(IN) :: B
    CHARACTER(LEN=*),           INTENT(IN) :: XMLTag
    INTEGER                                :: iStatus

    ! Local Variables
    TYPE(BROOK), POINTER :: ptr
    INTEGER              :: iUnit
    INTEGER              :: iForm

    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       ! Clearly there is nothing to be done if 
       ! we are not even initialized
       if ( iStatus /= BE_NONE ) return
    END IF
    ptr => B
    B_Loop: DO
       IF ( .NOT. ASSOCIATED(ptr) .OR. iStatus/= 0 ) EXIT B_Loop
       iForm = Brook_IForm(ptr)
       IF ( iForm == B_IFORM_XML .OR. iForm == B_IFORM_ASCII) THEN
          CALL Brook_Open(ptr,iStatus=iStatus)
          if ( iStatus /= BE_NONE ) return
          iUnit = Brook_Unit(ptr)
          CALL Unit_CloseXMLTag(iUnit, XMLTag=XMLTag, iStatus=iStatus)
          if ( iStatus /= BE_NONE ) return
       END IF
       ptr =>Brook_Next( ptr)
    END DO B_Loop
    RETURN
  END SUBROUTINE Brook_CloseXMLTag

  SUBROUTINE B_WXMLC_STRING(B, Comment, iStatus)
    implicit none
    type(Brook), target, intent(in) :: B
    character(LEN=*),    intent(in) :: Comment
    integer,             intent(inout) :: iStatus

    ! Local Variables
    TYPE(BROOK),      POINTER :: ptr
    INTEGER                   :: iUnit
    INTEGER                   :: iForm

    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       ! Clearly there is nothing to be done if 
       ! we are not even initialized
       if ( iStatus /= BE_NONE ) return
    END IF

    ! Open all Brooks
    call Brook_Open(B, iStatus)
    if ( iStatus /= BE_NONE ) return
    ! Now write to all ascii and xml brooks

1   FORMAT('  <!-- ')
2   FORMAT(a)
3   FORMAT(' -->')

    ptr => B
    WRITELOOP: DO
       IF ( .NOT. ASSOCIATED (ptr) ) THEN
          EXIT WRITELOOP
       END IF
       iUnit = Brook_UNIT(ptr)
       IF ( iUnit /=  B_INVALID_UNIT) THEN
          iForm = Brook_IForm(ptr)
          IF ( iForm == B_IFORM_XML .OR. iForm == B_IFORM_ASCII) THEN
             WRITE(iUnit,1, ADVANCE='NO')
             WRITE(iUnit,2, ADVANCE='NO') Comment
             WRITE(iUnit,3, ADVANCE='YES')
          END IF
       END IF
       ptr => Brook_Next(ptr)
    END DO WRITELOOP
       
    RETURN
    
  END SUBROUTINE B_WXMLC_STRING
  SUBROUTINE B_WXMLC_STRINGArray(B, Comment, iStatus)
    implicit none
    type(Brook),      target,       intent(in)    :: B
    character(LEN=*), dimension(:), intent(in)    :: Comment
    integer,                        intent(inout) :: iStatus

    ! Local Variables
    TYPE(BROOK),      POINTER :: ptr
    INTEGER                   :: iUnit
    INTEGER                   :: iForm

    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       ! Clearly there is nothing to be done if 
       ! we are not even initialized
       if ( iStatus /= BE_NONE ) return
    END IF

    ! Open all Brooks
    call Brook_Open(B, iStatus)
    if ( iStatus /= BE_NONE ) return
    ! Now write to all ascii and xml brooks

1   FORMAT('  <!-- ')
2   FORMAT('       ',a)
3   FORMAT(' -->')

    ptr => B
    WRITELOOP: DO
       IF ( .NOT. ASSOCIATED (ptr) ) THEN
          EXIT WRITELOOP
       END IF
       iUnit = Brook_UNIT(ptr)
       IF ( iUnit /=  B_INVALID_UNIT) THEN
          iForm = Brook_IForm(ptr)
          IF ( iForm == B_IFORM_XML .OR. iForm == B_IFORM_ASCII) THEN
             WRITE(iUnit,1, ADVANCE='NO')
             WRITE(iUnit,2, ADVANCE='YES') Comment
             WRITE(iUnit,3, ADVANCE='YES')
          END IF
       END IF
       ptr => Brook_Next(ptr)
    END DO WRITELOOP
       
    RETURN
    
  END SUBROUTINE B_WXMLC_STRINGARRAY


  SUBROUTINE Unit_OpenXMLTag(iUnit, XMLTag, XMLAttributes, ADVANCE_STR, iStatus)
    ! Open an xmlTag in unit iUnit.
    ! No Error check is done
    INTEGER,                    INTENT(IN) :: iUnit
    INTEGER                                :: iStatus
    CHARACTER(LEN=*),           INTENT(IN) :: XMLTag
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: ADVANCE_STR
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: XMLAttributes

    character(LEN=3) :: str

    if (present(ADVANCE_STR) ) then
       str = TRIM(ADVANCE_STR)
    else 
       str = 'YES'
    end if

#ifdef BROOK_DEBUG    
    WRITE(*,*) 'Opening: ', iUnit, TRIM(XMLTag)
#endif
    ! No Error because no brook functions
    iStatus = BE_NONE
10  FORMAT(/,'<',a)
11  FORMAT(' ', a)
12  FORMAT('>')
    WRITE(iUnit, 10, ADVANCE='NO') TRIM(XMLTag)
    IF (PRESENT(XMLAttributes)) THEN
       WRITE(iUnit, 11, ADVANCE='NO') TRIM(XMLAttributes)
    END IF
    WRITE(iUnit, 12, ADVANCE=str)
    RETURN
  END SUBROUTINE Unit_OpenXMLTag
  SUBROUTINE Unit_CloseXMLTag(iUnit, XMLTag, iStatus)
    ! Open an xmlTag in unit iUnit.
    ! No Error check is done
    INTEGER,             INTENT(IN) :: iUnit
    CHARACTER(LEN=*),    INTENT(IN) :: XMLTag
    INTEGER                         :: iStatus
    

    iStatus = BE_NONE
    
#ifdef BROOK_DEBUG    
    WRITE(*,*) 'Closing: ', iUnit, TRIM(XMLTag)
#endif
10  FORMAT('</',a,'>')
    WRITE(iUnit, 10) TRIM(XMLTag)
    RETURN
  END SUBROUTINE Unit_CloseXMLTag

  SUBROUTINE B_Initialize(iStatus)
    ! Initialize the basic CO data structures
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: iStatus

    iStatus = BE_NONE
    IF ( B_initialized ) THEN
       RETURN
    END IF

    B_initialized = .TRUE.

    CALL Brook_Set(B_Stdout, UNIT=B_STDOUT_UNIT, FILE='stdout', FORM='ascii', Closeable=.FALSE., iStatus=iStatus)
!!$    write(*,*) 'in b_initialzied, stdout is: ' , B_Stdout%Data%Unit, iStatus
    if (iStatus /= BE_NONE) return

    CALL Brook_Initialize(B = B_Stderr, iStatus = iStatus)
    if (iStatus /= BE_NONE) return

    CALL Brook_Set(B_Stderr, UNIT=B_STDERR_UNIT, FILE='stderr', FORM='ascii', Closeable=.FALSE., iStatus=iStatus)
    if (iStatus /= BE_NONE) return


    RETURN
  END SUBROUTINE B_Initialize

  SUBROUTINE Brook_Initialize(B, DATA, iStatus)
    IMPLICIT NONE
    ! initialize the struct B
    TYPE(Brook),            INTENT(INOUT) :: B
    TYPE(B_Data), OPTIONAL, POINTER       :: DATA
    INTEGER,                INTENT(INOUT) :: iStatus

    ! Local variables
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       if ( iStatus /= BE_NONE ) return
    END IF
    ALLOCATE(B%Data, STAT=iStatus)
    if ( iStatus /= 0 ) then
       call B_Set_Error(BE_INIT,'Allocation Error 1 in Brook_Initialzie')
       return
    end if

    IF ( PRESENT(DATA) ) THEN
       IF ( ASSOCIATED(DATA) ) THEN
          B%Data = DATA
       ELSE
          B%Data = B_data_Default
       END IF
    ELSE
       B%Data = B_data_Default
    END IF
    B%Next => Null()
    RETURN
    
  END SUBROUTINE Brook_Initialize

  SUBROUTINE B_ToLowerCase(S_in, S_out, iStatus)
    ! No Error Checking is done
    CHARACTER(len=*), INTENT(IN)  :: S_in
    CHARACTER(len=*), INTENT(OUT) :: S_out
    INTEGER                       :: iStatus
    INTEGER :: ILEN
    INTEGER :: i
    INTEGER :: isout
    INTEGER, SAVE :: iUpperA = -666
    INTEGER, SAVE :: iUpperZ = -666
    INTEGER, SAVE :: iOffset  = -666
    
    iStatus = BE_NONE
    IF ( iOffset == -666 ) THEN
       iUpperA = IACHAR('A')
       iUpperZ = IACHAR('Z')
       iOffset = IACHAR('a') - IACHAR('A')
    END IF

    S_out = S_in
    ILEN = LEN(TRIM(S_out))
    DO i = 1, ILEN
       isout = IACHAR(S_out(i:i))
       IF ( iUpperA<=isout .AND. isout <= iUpperZ ) THEN
          S_out(i:i) = ACHAR(iOffset + isout)
       END IF
    END DO

  END SUBROUTINE B_ToLowerCase

  RECURSIVE SUBROUTINE Brook_Set(B,               &
                                 File,            &
                                 Unit,            &
                                 Form,            &
                                 Position,        & 
                                 StandardHeaders, &
                                 Closeable,       &
                                 IFormat,         &
                                 DFormat,         &
                                 LFormat,         &
                                 AFormat,         &
                                 XMLRoot,         &
                                 iStatus)
    IMPLICIT NONE
    
    ! Argument List
    TYPE(Brook),      TARGET,   INTENT(INOUT) :: B
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: File
    INTEGER,          OPTIONAL, INTENT(IN)    :: Unit
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Form
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Position
    LOGICAL,          OPTIONAL, INTENT(IN)    :: StandardHeaders
    LOGICAL,          OPTIONAL, INTENT(IN)    :: Closeable
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: XMLRoot
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Iformat
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Dformat
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Lformat
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: Aformat
    INTEGER,                    INTENT(INOUT) :: iStatus
    ! Local Variables
    TYPE(Brook),       POINTER :: ptr
    CHARACTER(LEN=B_STRLEN)    :: SLocal

    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       if ( iStatus /= BE_NONE ) return
    END IF

    IF ( .NOT. ASSOCIATED(B%Data) ) THEN
       CALL Brook_Initialize(B, iStatus=iStatus)
       IF (iStatus /= BE_NONE ) THEN
          RETURN
       END IF
       B%Next=>Null()
    END IF

    ! 
    ! Only File, and unit apply to the top level struct. 
    ! All others cascade down all structs
    
    IF (PRESENT(File)) THEN
       B%Data%File = TRIM(File)
       B%Data%Status = B_STATUS_INIT  ! declare that this is a valid file
    END IF

    IF (PRESENT(UNIT)) THEN
       B%Data%Unit = UNIT
    END IF

    IF (PRESENT(Position)) THEN
       CALL B_ToLowerCase(Position, B%Data%Position, iStatus=iStatus)
       if (iStatus /= 0 ) return
    END IF


    IF (PRESENT(StandardHeaders)) THEN
       B%Data%StandardHeaders = StandardHeaders
    END IF

    ! Now cycle through the rest of the brooks
    ptr => B
    DO
       IF ( .NOT. ASSOCIATED (ptr) ) THEN
          EXIT
       END IF

       IF (PRESENT(Form)) THEN
          CALL B_ToLowerCase(Form,SLocal,iStatus)
          if (iStatus /= 0 ) return
          SELECT CASE (TRIM(SLocal))

          CASE ("ascii")
             PTR%Data%Form = 'ascii'
             PTR%Data%IForm = B_IFORM_ASCII
          CASE ("xml")
             PTR%Data%Form = 'xml'
             PTR%Data%IForm = B_IFORM_XML
          CASE ("binary")
             PTR%Data%Form = TRIM(Form)
             PTR%Data%IForm = B_IFORM_BINARY
             PTR%Data%Offset = 0 ! Brook data offset
          CASE default
             iStatus = BE_SET
             call B_Set_Error(istatus, 'Invalid Form sent to BE_Set')
             RETURN
          END SELECT
       END IF
       IF (PRESENT(Closeable)) THEN
          ptr%Data%Closeable = Closeable
       END IF
       
       IF (PRESENT(XMLRoot)) THEN
          ptr%Data%XMLRoot = TRIM(XMLRoot)
       END IF

       IF (PRESENT(Iformat)) THEN
          ptr%Data%Iformat = TRIM(Iformat)
       END IF
       
       IF (PRESENT(Dformat)) THEN
          ptr%Data%Dformat = TRIM(Dformat)
       END IF
       
       IF (PRESENT(Lformat)) THEN
          ptr%Data%Lformat = TRIM(Lformat)
       END IF
       
       IF (PRESENT(Aformat)) THEN
          ptr%Data%Aformat = TRIM(Aformat)
       END IF

       ptr =>Brook_Next(ptr)
    END DO
    
    RETURN
  END SUBROUTINE Brook_Set

  SUBROUTINE Brook_Diagnostics(B, iStatus)
    TYPE(Brook), TARGET, INTENT(IN) :: B
    TYPE(Brook), POINTER :: P
    integer :: iStatus
    INTEGER :: i
    INTEGER :: j
    INTEGER :: u

    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       if (iStatus /= BE_NONE ) return
    END IF

    i = 0
    j = 0

    P=>B
    DO
       IF ( .NOT. ASSOCIATED(P)) EXIT
       i = i + 1
10     FORMAT('Error: ',i4,' File :',i4,' ',a)
11     FORMAT('<File Name="',a,'">',/,'  <Error> ',i4,'  </Error>',/,'  <File_id>',i4,'  </File_id>',/,'</File>')
       IF ( ASSOCIATED(P%Data)) THEN
          CALL Brook_OPEN(p, iStatus=iStatus)
          if (iStatus /= 0 ) return
          u = Brook_UNIT(p)
          IF ( u /= B_INVALID_UNIT) THEN
             IF ( Brook_IFORM(P) == B_IFORM_ASCII) THEN
                WRITE(u, *) j,i,TRIM(Brook_File(P))
             ELSE IF ( Brook_IFORM(P) == B_IFORM_XML) THEN
                WRITE(u, 11) TRIM(Brook_File(P)),j,i
             ELSE
                WRITE(u) j,i,TRIM(Brook_File(P))
             END IF
          END IF
          WRITE(*,10) j, i, TRIM(Brook_File(P))
       ELSE 
          WRITE(*,10) 0, i, 'unassociated'
       END IF
       P=>Brook_Next(P)
    END DO
    WRITE(*,*) 'Totally ',i, ' Files'
    WRITE(*,*) ' '
  END SUBROUTINE Brook_Diagnostics

  LOGICAL FUNCTION B_OP(B, iStatus)
    TYPE(Brook), TARGET, INTENT(IN) :: B
    integer :: iStatus
    ! No Error Checking
    iStatus = BE_NONE
    B_OP = ASSOCIATED(B%Data)
    if ( B_OP ) B_OP = (B%Data%Status /= B_STATUS_NO_INIT)
  END FUNCTION B_OP

  INTEGER FUNCTION BROOK_AssignUnit()
    ! Find a free Unit and return it
    integer :: iUnit
    logical :: openP
#ifdef BROOK_DEBUG_UNIT
    integer :: iStatus
    CHARACTER(LEN=B_strlen) :: tmp
#endif

    FINDUNIT: DO iUnit = B_UNIT_LOW, B_UNIT_HIGH, B_UNIT_INCREMENT
       openP = .FALSE.
       if (B_AssignedUnits(iUnit)) CYCLE FINDUNIT
       INQUIRE(UNIT = iUnit, OPENED = openP)
       IF (.NOT. openP ) EXIT FINDUNIT
    END DO FINDUNIT
#ifdef BROOK_DEBUG_UNIT
    write(Tmp,*) 'Brook_AssignUnit:',iUnit
    call Brook_Write(B_stdout, Trim(Tmp), advance=.true., iStatus=iStatus)
#endif
    IF ( openP ) THEN
       ! Couldn't find a unit number
       BROOK_AssignUnit = B_INVALID_UNIT
    ELSE
       Brook_AssignUnit = iUnit
       B_AssignedUnits(iUnit) = .true.
    END IF
    return
  END FUNCTION BROOK_AssignUnit
  Subroutine BROOK_CleanUp()
    ! Close all units we assigned if they are open
    integer :: iUnit
    logical :: openP
#ifdef BROOK_DEBUG_UNIT
    integer :: iStatus
    CHARACTER(LEN=B_strlen) :: tmp
#endif

    FINDUNIT: DO iUnit = B_UNIT_LOW, B_UNIT_HIGH, B_UNIT_INCREMENT
       openP = .FALSE.
       if (B_AssignedUnits(iUnit)) CYCLE FINDUNIT
       INQUIRE(UNIT = iUnit, OPENED = openP)
       IF (.NOT. openP ) EXIT FINDUNIT
    END DO FINDUNIT
#ifdef BROOK_DEBUG_UNIT
    write(Tmp,*) 'Brook_AssignUnit:',iUnit
    call Brook_Write(B_stdout, Trim(Tmp), advance=.true., iStatus=iStatus)
#endif
    IF ( openP ) THEN
       ! Couldn't find a unit number
    ELSE
       B_AssignedUnits(iUnit) = .true.
    END IF
    return
  END Subroutine BROOK_CleanUp

  SUBROUTINE Brook_ReleaseUnit(iUnit, iStatus)
    integer, intent(in)  :: iUnit
    integer, intent(out) :: iStatus
#ifdef BROOK_DEBUG_UNIT
    integer :: iStatus2
    CHARACTER(LEN=B_strlen) :: tmp
#endif
    iStatus = 0
    if ( iUnit .le. B_UNIT_HIGH .and. iUnit .ge. B_UNIT_LOW ) then
       B_AssignedUnits(iUnit) = .false.
    end if
#ifdef BROOK_DEBUG_UNIT
    iStatus2 = 0
    write(Tmp,*) 'Brook_ReleaseUnit:',iUnit 
   call Brook_Write(B_stdout, Trim(Tmp), advance=.true., iStatus=iStatus2)
#endif
    return
  END SUBROUTINE Brook_ReleaseUnit

! SS iStatus
  RECURSIVE SUBROUTINE B_OPEN_TopLevel (B, inForm, inPosition, iStatus) 
    !==================================================================
    ! Purpose(s):
    !   Open file associated with B, if not opened already
    !   If no file, then error.
    !   If no error, then return 0.  
    !==================================================================
    IMPLICIT NONE

    ! Arguments
    TYPE(Brook),      TARGET,   INTENT(INOUT) :: B
    INTEGER,                    INTENT(INOUT) :: iStatus
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: inForm
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: inPosition

    ! Local variables
    LOGICAL                 :: openp
    INTEGER                 :: iUnit
    CHARACTER(LEN=B_strlen) :: File_Name
    CHARACTER(LEN=B_strlen) :: File_Form
    CHARACTER(LEN=B_strlen) :: File_Position
    CHARACTER(LEN=B_strlen) :: tmp
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       if ( iStatus /= BE_NONE ) return
    END IF

    IF (.NOT. ASSOCIATED(B%Data)) RETURN
    if ( B%Data%Status == B_STATUS_NO_INIT) RETURN

    iUnit = Brook_Unit(B)
    ! Don't try to open stdin or stdout
    IF (iUnit == B_stdout_unit  .OR. &
        iUnit == B_stderr_unit) THEN
       RETURN
    END IF
    IF (iUnit == B_AUTO_UNIT .or. iUnit == B_INVALID_UNIT .or. iUnit == B_UNIT_UNKNOWN .or. iUnit .le. 0) THEN
       ! Find a valid unit number between 128 and 512
       INQUIRE(FILE=TRIM(BROOK_FILE(B)), number=iUnit)
#ifdef BROOK_DEBUG_UNIT          
       write(*,*) 'OPENTOP: ', TRIM(BROOK_FILE(B)), ' iUnit=', iUnit
#endif
       if ( iUnit .le. 0 ) iUnit = Brook_AssignUnit()
       if ( iUnit == B_INVALID_UNIT) then
          ! Couldn't find a unit number
          iStatus = BE_OPEN
          call B_Set_Error(iStatus, 'Unable to get valid unit number')
          return
       else 
          call BROOK_SET(B, UNIT=iUNIT, iStatus=iStatus)
       end if
       call BROOK_SET(B=B, UNIT=iUnit, iStatus=iStatus)
       if ( iStatus /= BE_NONE) return
    END IF
          
    IF (iUnit == B_Unit_Unknown) THEN
       ! Not a valid unit in the unit number
       ! Should never be here
       iStatus = BE_OPEN
       call B_Set_Error(iStatus, 'Invalid unit number')
       return
    END IF

    IF (iStatus == BE_NONE ) THEN
       ! If valid unit number check, otherwise error
       INQUIRE(UNIT = iUnit, OPENED = openP, form = Tmp)
       ! If not opened, open it
       OPEN_FILE: IF (.NOT. openP) THEN
          ! Need a valid file 
          file_name = Brook_FILE (B)
          IF (file_name /= B_file_unknown) THEN
             IF ( PRESENT (inForm)) THEN
                File_Form = TRIM(inForm)
             ELSE 
                FILE_FORMAT: SELECT CASE (Brook_IForm (B))
                CASE (B_IForm_binary)
                   File_Form = 'unformatted'
                CASE (B_IForm_ascii)
                   File_Form = 'formatted'
                CASE (B_IForm_XML)
                   File_Form = 'formatted'
                CASE default FILE_FORMAT
                   ! Unknown file format
                   File_Form = 'formatted'
                END SELECT FILE_FORMAT
             END IF

             IF ( PRESENT(inPosition) ) THEN
                File_Position = TRIM(inPosition)
             ELSE 
                File_Position = TRIM(Brook_Position(B))
             END IF

             IF (iStatus == BE_NONE) THEN
                OPEN(UNIT     = iUnit,               &
                     FILE     = TRIM(file_name),     &
                     STATUS   = 'unknown',           &
                     POSITION = TRIM(File_Position), &
                     FORM     = TRIM(File_Form),     & 
#ifdef DARWIN_NAG_COMPILER_WORKAROUND
	             RECL    = 16777216,              &
#endif
                     IOSTAT=iStatus)
                if (iStatus == BE_NONE ) then
                   call Brook_Write_File_Header(B, iStatus)
                end if
             END IF
          ELSE
             ! Not a valid unit number
             iStatus = BE_OPEN
             call B_Set_Error(iStatus, 'Invalid Unit')
             return
          END IF
       ELSE 
          ! Make sure we have the correc tform, or else reopen and append
          OPEN_FORMAT: SELECT CASE(Tmp)
          CASE ('UNFORMATTED')
             IF (Brook_IFORM(B) /= B_IForm_Binary ) THEN
                CALL B_Close_TopLevel(B, iStatus)
                if ( iStatus /= BE_NONE) return
                CALL B_Open_TopLevel(B, inPosition = 'APPEND', iStatus=iStatus)
                if ( iStatus /= BE_NONE) return
                RETURN
             END IF
          CASE ('FORMATTED')
             IF (Brook_IFORM(B) == B_IForm_Binary ) THEN
                CALL B_Close_TopLevel(B, iStatus)
                if ( iStatus /= BE_NONE) return
                CALL B_Open_TopLevel(B, inPosition = 'APPEND', iStatus=iStatus)
                if ( iStatus /= BE_NONE) return
                RETURN
             END IF
          END SELECT OPEN_FORMAT
       END IF OPEN_FILE
    END IF
    
    RETURN

  END SUBROUTINE B_OPEN_TopLevel


  RECURSIVE SUBROUTINE Brook_OPEN (B, iStatus) 
    !==================================================================
    ! Purpose(s):
    !   Open file associated with B, if not opened already
    !   If no file, then error.
    !   If no error, then return 0.  
    !==================================================================
    IMPLICIT NONE

    ! Arguments
    TYPE(Brook), TARGET, INTENT(IN)    :: B
    INTEGER,             INTENT(INOUT) :: iStatus

    ! Local variables
    TYPE(Brook), POINTER  :: P

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    iStatus = BE_NONE
    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       if (iStatus /= BE_NONE ) return
    END IF

    ! Check all Brooks
    P => B
    BROOKS: DO
       IF ( iStatus /= BE_NONE ) EXIT BROOKS
       IF ( .NOT. ASSOCIATED(P) ) EXIT BROOKS
       IF (ASSOCIATED(P%Data)) then
          CALL B_Open_TopLevel(P,iStatus=iStatus)
          if (iStatus /= BE_NONE) return
       END IF
       P =>Brook_Next( P)
    END DO BROOKS

    RETURN

  END SUBROUTINE Brook_OPEN

  SUBROUTINE Brook_Flush (B, iStatus) 
    !==================================================================
    ! Purpose(s):
    !   Flush all brooks associated with B, if opened already
    !   Always return no error
    !==================================================================
#ifdef NAG_COMPILER
    use f90_unix, only: flush
#endif
    IMPLICIT NONE

    ! Arguments
    TYPE(Brook), TARGET, INTENT(IN)    :: B
    INTEGER,             INTENT(INOUT) :: iStatus

    ! Local variables
    TYPE(Brook), POINTER  :: P
    integer               :: iUnit
    logical               :: openP
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    iStatus = BE_NONE

    IF ( .NOT. B_Initialized) THEN
       CALL B_Initialize(iStatus)
       return ! Clearly we couldn't have an open brook
    END IF

    ! Flush all Brooks
    P => B
    BROOKS: DO
       IF ( iStatus /= BE_NONE ) EXIT BROOKS
       IF ( .NOT. ASSOCIATED(P) ) EXIT BROOKS

       openP = .false.
       iUnit = Brook_Unit(P)
       IF (iUnit > 0 ) then
          INQUIRE(UNIT = iUnit, OPENED = openP)
          if (openP) then
             CALL FLUSH(iUnit)
          END IF
       END IF
       P =>Brook_Next( P)
    END DO BROOKS
    RETURN

  END SUBROUTINE Brook_Flush

  
  SUBROUTINE Brook_Remove_Duplicates (B, iStatus)
    ! Remove duplicate brooks from B
    Type(Brook), target :: B
    integer, intent(inout) :: iStatus
    type(Brook), pointer :: P
    type(Brook), pointer :: Q
    type(Brook), pointer :: N
    integer              :: iUnitP

    logical :: removeN

    iStatus = 0

    P=>B
    OuterLoop: do
       if ( .not. associated(P) ) then
          ! We are done
          exit OuterLoop
       end if
       iUnitP = Brook_Unit(P)
       Q=>P
       N=>Brook_Next(P)
       removeN = .false.
       InnerLoop: do
          if ( .not. associated(N) ) then
             ! We are done
             exit InnerLoop
          end if
          ! 
          ! Compare unit number and File name
          if (Brook_Unit(N) > 0 .and. iUnitP > 0) then
             if (Brook_Unit(N) == Brook_Unit(P)) then
                removeN=.true.
                exit InnerLoop
             end if
          else if (TRIM(Brook_File(N)) == TRIM(Brook_File(P))) then
             removeN=.true.
             exit InnerLoop
          end if
          Q=>N
          N => Brook_Next(N)
       end do InnerLoop
       if ( removeN .and. associated(Q) .and. associated(N)) then
          ! remove N from our list
          ! point Q at the brook after next
          Q%Next => N%Next
          Nullify(N%Next)
          Deallocate(N%Data)
          Deallocate(N)
       end if
       P => Brook_Next(P)
    end do OuterLoop
  END SUBROUTINE Brook_Remove_Duplicates

end module brook_module

#ifdef LOCALPROGRAM

#define UTILIZE use

PROGRAM test
UTILIZE brook_module, only: BROOK,          &
                   B_STDOUT,       &
                   B_STDERR,       &
                   Brook_DESTROY,  &
                   Brook_WRITE,    &
                   Brook_CLOSE,    &
                   OPERATOR(+),    &
                   ASSIGNMENT(=),  &
                   Brook_SET,      &
                   Brook_AFormat,  &
                   Brook_Unit
  IMPLICIT NONE

  integer             :: iStatus
  TYPE(BROOK)         :: a
  TYPE(BROOK)         :: b
  TYPE(BROOK)         :: c
  REAL*8              :: r
  DOUBLE PRECISION    :: q(20)
  INTEGER             :: i2(5,2)
  INTEGER             :: i
  INTEGER             :: j
  LOGICAL             :: openP


  CALL Brook_Write(B_Stdout, Variable='Checking Brook_Set() on new brooks',iStatus=iStatus)
  CALL Brook_Set(a, form='Xml',    Closeable=.TRUE.,   file='a.xml',iStatus=iStatus)
  CALL Brook_Set(b, form='bINary', Closeable=.FALSE.,  UNIT=11, file='b.bin',iStatus=iStatus)
  CALL Brook_Set(c, form='ascII',  Closeable=.FALSE.,  UNIT=12, file='c.txt',iStatus=iStatus)



  c = c + a + b 

  q = 20.0d0
  DO i = 1, 10
     q(i) = i+1.123d0
  END DO
 
  DO i = 1, 5
     DO j = 1,2
        i2(i,j) = 10*i + j
     END DO
  END DO
  CALL Brook_Write(B_Stdout, Variable='Checking Brook_Set() on oldbrooks',iStatus=iStatus)
  CALL Brook_Set(c, Closeable=.TRUE.,iStatus=iStatus)

  CALL Brook_Write(B_Stdout, Variable='Checking Brook_Write()',iStatus=iStatus)
  ! Write a string variable to all files
  CALL Brook_Write(B=c,Variable='World',XMLName='Hello',iStatus=iStatus)

  ! Write a Logical variable to all files
  INQUIRE(UNIT = Brook_Unit(c), OPENED = openP)
  CALL Brook_Write(B=c,Variable=openP,XMLName='Is_C_Open',iStatus=iStatus)

  ! Write a constant integer, and a real variable to all files
  r = 15.123987981364761532467523847
  CALL Brook_Write(B=c,Variable=10,XMLName='Constant10',XMLAttributes='MYATTR="This is a Constant" ',iStatus=iStatus)
  CALL Brook_Write(B=c,Variable=r,XMLName='r',iStatus=iStatus)

  ! Write a 2D Real array to al lfiles with format specified
  CALL Brook_Write(B=c,Variable=q,XMLName='q',FORMAT='(10(E10.4,x))',iStatus=iStatus)

  ! Write a 2D integer array to all files with a look-aside file for some of them
  
  CALL Brook_Write(B=c,Variable=i2,XMLName='i2',XMLAttributes='Mesh="DefaultMesh"',iStatus=iStatus)
  CALL Brook_Write(B=c,Variable=i2,XMLName='i2',XMLAttributes='Mesh="EMMesh"', &
               XMLDataFile='corder.txt',XMLDataFormat='ascii',C_Array_Order=.FALSE.,iStatus=iStatus)
  CALL Brook_Write(B=c,Variable=i2,XMLName='i2',XMLAttributes='Mesh="None"',   &
               XMLDataFile='corder.xml',XMLDataFormat='xml',C_Array_Order=.TRUE.,iStatus=iStatus)

  CALL Brook_Write(B=c,Variable=i2,XMLName='i2',XMLAttributes='Mesh="None"',iStatus=iStatus)
  CALL Brook_Write(B=c,Variable=i2,XMLName='i2',XMLAttributes='Mesh="None"',C_ARRAY_ORDER=.true., iStatus=iStatus)

  ! Close all files
  CALL Brook_Close(c,iStatus=iStatus)

  ! Destroy all brooks.
  ! Note that we have to destroy a and b separately even though we
  ! had copied them into c.  The reason for this is that the addition
  ! places copies of the source brooks into the destination, instead
  ! of pointers.
  CALL Brook_Destroy(a,iStatus=iStatus)
  CALL Brook_Destroy(b,iStatus=iStatus)
  CALL Brook_Destroy(c,iStatus=iStatus)


  ! Create Brook a as an ascii brook
  CALL Brook_Set(a, form='ascii', UNIT=10, file='aa',iStatus=iStatus)

  ! Write i2 as an ASCII string
  CALL Brook_Set(a, form='ascii',iStatus=iStatus)
  CALL Brook_Write(a,Variable='i2 in ascii',iStatus=iStatus)
  CALL Brook_Write(a,Variable=i2,iStatus=iStatus)

  ! Change form to binary and write i2 in binary
  CALL Brook_Set(a, form='binary',iStatus=iStatus)
  CALL Brook_Write(a,Variable='i2 in binary',iStatus=iStatus)
  CALL Brook_Write(a,Variable=i2, XMLName='i2',iStatus=iStatus)

  ! Change form back to ASCII and write i2
  CALL Brook_Set(a, form='ascii',iStatus=iStatus)
  CALL Brook_Write(a,Variable='i2 in ascii',iStatus=iStatus)
  CALL Brook_Write(a,Variable=i2, C_Array_Order=.TRUE.,iStatus=iStatus)

  ! Close the file and destroy the struct
  CALL Brook_Close(a,iStatus=iStatus)
  CALL Brook_Destroy(a,iStatus=iStatus)

  CALL Brook_Write(B_Stdout, Variable='test done',iStatus=iStatus)
END PROGRAM test

#endif

#ifdef LOCALPROGRAM
#undef LOCALPROGRAM
#endif

#ifdef BROOK_DEBUG
#undef BROOK_DEBUG
#endif

#ifdef BROOK_DEBUG_UNIT
#undef BROOK_DEBUG_UNIT
#endif
