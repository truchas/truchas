!!CPP!! Include file for brook_module.F90
!!CPP!! Provides routines for output of native types

!!CPP!! 
!!CPP!! Note that character arrays cannot be more 
!!CPP!! than one dimension
!!CPP!!

#ifndef _FUNCTION_NAME_ 
#error "_FUNCTION_NAME_ must be defined before including this file"
#endif

#ifdef _IS_CHARACTER_TYPE_
#define _NO_HIGHER_ORDER_XML_
#else
#ifndef _FUNCTION_NAME_RESHAPE_ 
#error "_FUNCTION_NAME_RESHAPE_ must be defined before including this file"
#endif
#endif

#ifndef _FORMAT_NAME_
#error "_FORMAT_NAME_ must be defined before including this file"
#endif

#ifndef _TYPE_
#error "_TYPE_  must be defined before including this file"
#endif

#ifndef _DIMENSION_
#error "_DIMENSION_ must be defined before including this file"
#endif

#ifndef _TYPE_SIZE_
#error "_TYPE_SIZE_ must be defined before including this file"
#endif

RECURSIVE SUBROUTINE _FUNCTION_NAME_ (B, Variable, FORMAT,  C_Array_Order, ADVANCE, &
                                      XMLName, XMLAttributes, XMLDataFile, XMLDataFormat, &
                                      iStatus)
  !==================================================================
  ! Purpose(s):
  !   Write Variable to the output brook.
  !   Status is returned in iStatus
  !==================================================================
  !  use kind_module

  IMPLICIT NONE

  ! Arguments
  LOGICAL,          OPTIONAL :: ADVANCE
  TYPE(Brook), TARGET,   INTENT(IN) :: B
  _TYPE_ , _DIMENSION_  INTENT(IN)  :: Variable
  CHARACTER(LEN=*), OPTIONAL :: FORMAT
  CHARACTER(LEN=*), OPTIONAL :: XMLName
  CHARACTER(LEN=*), OPTIONAL :: XMLAttributes
  CHARACTER(LEN=*), OPTIONAL :: XMLDataFile
  CHARACTER(LEN=*), OPTIONAL :: XMLDataFormat
  LOGICAL,          OPTIONAL :: C_Array_Order
  INTEGER,          OPTIONAL :: iStatus

  ! Local variables
  CHARACTER(LEN=20)         :: a
  CHARACTER(LEN=B_STRLEN)   :: myFmt
  CHARACTER(LEN=B_STRLEN)   :: tmp
  TYPE(Brook),   POINTER    :: ThisB
  INTEGER                   :: error_code
  INTEGER                   :: iUnit
  INTEGER                   :: iDim
  LOGICAL                   :: OrderSwap
  LOGICAL                   :: askance
  INTEGER                   :: i1
  INTEGER                   :: iBytes
  INTEGER, ALLOCATABLE, DIMENSION(:) :: iPtr
#ifdef _NO_HIGHER_ORDER_XML_
   integer, pointer :: ptr
#else 
  _TYPE_ , _DIMENSION_  pointer :: ptr
#endif
  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  if ( PRESENT (iStatus) ) iStatus = BE_ID
  if ( BE_ID /= BE_NONE ) return
  error_code = BE_NONE
  if (present(ADVANCE)) then
     if (.not. ADVANCE) then
        a = 'NO'
     else 
        a = 'YES'
     end if
  else if (BROOK_ADVANCE(B)) then
     a = 'YES'
  else 
     a = 'NO'
  end if


#ifdef  _NO_HIGHER_ORDER_XML_

  OrderSwap = .false.   
  ptr => NULL()

#else

#ifndef _RV_
#error "_RV_ must be define before including this file"
#endif

#ifndef _RP_
#error "_RP_ must be define before including this file"
#endif

  if (PRESENT(C_Array_Order)) then
     OrderSwap = (C_Array_Order)
  else 
     OrderSwap = .false.
  end if

  if (OrderSwap) then

     call _FUNCTION_NAME_RESHAPE_ (iStatus=error_code, _RV_ = Variable, _RP_ = ptr)
     
     IF (error_code /= 0 ) THEN
        IF (PRESENT(iStatus)) THEN
           iStatus = error_code
        else 
           if ( error_code /= 0 ) then
              call Brook_Abort('error in output brook_include2.fpp')
           end if
        end IF
        RETURN
     END IF
  ELSE 
     ptr=>NULL()
  END IF
#undef _RV_
#undef _RP_
#endif


  ! Write to all brooks
  ThisB => B
  
  ! Open all brooks
  call Brook_Open(B,error_code)

  if ( present(XMLDataFile)) then
     Askance = .true.
  else 
     Askance = .false.
  end if

  WRITE_SCALAR: DO
     if (error_code /= 0 ) EXIT WRITE_SCALAR
     IF (.NOT. ASSOCIATED(ThisB)) THEN
#ifdef BROOK_DEBUG_1 
        WRITE(*,*) 'Unassociated ThisB, exiting WRITE_SCALAR '
#endif
        
        EXIT WRITE_SCALAR
     END IF

#ifdef BROOK_DEBUG 
     WRITE(*,*) 'Writing to Unit',Brook_Unit(ThisB),     &
                ' with File: ',TRIM(Brook_File(ThisB)),  &
                ' and form: ', TRIM(Brook_Form(ThisB)),  &
                ' (iForm=',Brook_IForm(ThisB),')'
#endif

     IF ( B_OP(ThisB, error_code) ) THEN

        iUnit = Brook_Unit(ThisB)

        OUT_FORMAT: SELECT CASE(Brook_IForm(ThisB))

        CASE (B_IForm_Binary)
           IF (OrderSwap) THEN
              WRITE(IUNIT) ptr
           ELSE 
              WRITE(IUnit) Variable
           END IF
           iBytes = _TYPE_SIZE_
           i1 = SIZE(SHAPE(VARIABLE))
           IF(i1 > 0 ) THEN
              allocate(iPTR(i1), STAT=iStatus)
              iPTR = SHAPE(VARIABLE)
              DO i1 = 1, SIZE(SHAPE(VARIABLE))
                 iBytes = iBytes*iPTR(i1)
              END DO
              DEALLOCATE(iPTR)
           END IF
#ifdef _IS_CHARACTER_TYPE_              
           iBytes = iBytes * _CHARLEN_ 
#endif
           ThisB%Data%Offset = ThisB%Data%Offset + iBytes 

        CASE (B_IForm_Ascii)

           IF (PRESENT(FORMAT)) THEN
              myFmt = FORMAT
           ELSE 
              myFmt = _FORMAT_NAME_(ThisB)
           END IF
           IF (OrderSwap) THEN
              WRITE(IUnit, FMT = myFMT, ADVANCE=a) ptr
           ELSE 
#ifdef _IS_CHARACTER_TYPE_              
              WRITE(IUnit, FMT = myFMT, ADVANCE=a) Variable
#else
              WRITE(IUnit, FMT = myFMT, ADVANCE=a) Variable
#endif
           END IF

        CASE (B_IForm_XML)
100        FORMAT('<Variable ',a,' ', a,'>')
11         FORMAT('</Variable>')
13         FORMAT(a)

           ! Open the tag XMLName
           iDim = SIZE(SHAPE(Variable))
           IF (PRESENT(XMLName)) THEN
              call Brook_Get_Standard_Attributes(tmp, TRIM(ADJUSTL(XMLName)), _TYPE_CODE_, iDim, SHAPE(Variable), C_Array_Order)
           ELSE
              call Brook_Get_Standard_Attributes(tmp, 'J_Doe', _TYPE_CODE_, iDim, SHAPE(Variable), C_Array_Order)
           END IF
           IF (PRESENT(XMLAttributes)) THEN
              write(IUnit, 100) TRIM(tmp), TRIM(ADJUSTL(XMLAttributes))
           ELSE
              write(IUnit, 100) TRIM(tmp), ' '
           END IF


           if ( .not. askance ) then
              ! Write the variable to the file
              call Unit_OpenXMLTag(iUnit,XMLTag='DataValues', iStatus=iStatus)
              IF (PRESENT(FORMAT)) THEN
                 myFmt = FORMAT
              ELSE 
                 myFmt = _FORMAT_NAME_(ThisB)
              END IF
              IF (OrderSwap) THEN
                 WRITE(IUnit, FMT = myFMT) ptr
              ELSE 
                 WRITE(IUnit, FMT = myFMT) Variable
              END IF
              call Unit_CloseXMLTag(iUnit,XMLTag='DataValues', iStatus=iStatus)
           else 
              if ( present(XMLDataFormat) ) then
                 write(tmp,*) 'Format="',TRIM(XMLDataFormat),'"'
              else 
                 write(tmp,*) 'Format="ascii"'
              end if
              call Unit_OpenXMLTag(iUnit,XMLTag='DataFile',XMLAttributes=TRIM(tmp), iStatus=iStatus)
              write(iUnit,13,ADVANCE='NO') TRIM(XMLDataFile)
              call Unit_CloseXMLTag(iUnit,XMLTag='DataFile', iStatus=iStatus)
           end if
           

           ! Close the Variable tag
           WRITE(iUnit,11) 
        CASE default

           ! Do Nothing

        END SELECT OUT_FORMAT

     END IF

     ThisB => Brook_Next(ThisB)

  END DO WRITE_SCALAR
 
  ! Write the look aside file, if requested
  if ( askance ) then
#ifdef BROOK_DEBUG
     write(*,*) 'Writing look-aside-file ...'
     if ( Present(XMLDataFormat) ) then
        write(*,*) '  Format is: ',TRIM(XMLDataFormat)
     end if
#endif
     call Brook_Set(B=B_Tmp,File=TRIM(XMLDataFile),FORM='ascii', Unit=B_Tmp_Unit,iStatus=iStatus)
     if ( present(XMLDataFormat) ) then
        call Brook_Set(B=B_Tmp,FORM=TRIM(XMLDataFormat), iStatus=iStatus)
     end if
     call Brook_Write(B                = B_Tmp,        &
                      Variable         = Variable,     &
                      Format           = Format,       &
                      C_Array_Order    = C_Array_Order,&
                      XMLName          = XMLName,      &
                      XMLAttributes    = XMLAttributes,&
                      iStatus          = iStatus       &
                      )
     call Brook_Close(B_Tmp, iStatus=iStatus)
     call Brook_Destroy(B_Tmp,iStatus=iStatus)
#ifdef BROOK_DEBUG
     write(*,*) 'Done'
#endif
  end if

  IF (PRESENT(iStatus)) THEN
     iStatus = error_code
  else 
     if ( error_code /= 0 ) then
        call Brook_Abort( 'error in output brook_include2.fpp')
     end if
  end IF

  if (Associated(ptr)) then
     deallocate(ptr)
  end if
  
  RETURN

END SUBROUTINE _FUNCTION_NAME_


#ifdef _DEFINE_RESHAPE_
SUBROUTINE _FUNCTION_NAME_RESHAPE_ (iStatus, &
                                    R1, PR1, &
                                    R2, PR2, &
                                    R3, PR3, &
                                    R4, PR4, &
                                    R5, PR5  &
                                    )
  ! Return a pointer to the reshaped variable

  integer, intent(INOUT) :: iStatus

  _TYPE_ , intent(IN), optional, dimension(:)          :: R1
  _TYPE_ , pointer, optional, dimension(:)                     :: PR1

  _TYPE_ , intent(IN), optional, dimension(:,:)        :: R2
  _TYPE_ , pointer, optional, dimension(:,:)                   :: PR2

  _TYPE_ , intent(IN), optional, dimension(:,:,:)      :: R3
  _TYPE_ , pointer, optional, dimension(:,:,:)                 :: PR3

  _TYPE_ , intent(IN), optional, dimension(:,:,:,:)    :: R4
  _TYPE_ , pointer, optional, dimension(:,:,:,:)               :: PR4

  _TYPE_ , intent(IN), optional, dimension(:,:,:,:,:)  :: R5
  _TYPE_ , pointer, optional, dimension(:,:,:,:,:)             :: PR5


  ! Local Variables
  integer :: s1, s2, s3, s4, s5
  integer :: i, j, k, l, m
  

  iStatus = BE_NONE
  
  if ( PRESENT (R1) .and. PRESENT(PR1) ) then
     s1 = size(R1,1)
     Allocate(PR1(s1), STAT=iStatus)
     PR1 = R1
  else if ( PRESENT (R2) .and. PRESENT(PR2) ) then
     s1 = size(R2,1)
     s2 = size(R2,2)
     Allocate(PR2(s2,s1), STAT=iStatus)
     if ( iStatus /= 0 ) then
        call B_Set_Error(BE_RESHAPE, 'Unable to allocate in reshape for C like array output')
        return
     end if
     do i = 1, s1
        do j = 1, s2
           PR2(j,i) = R2(i,j)
        end do
     end do
  else if ( PRESENT (R3) .and. PRESENT(PR3) ) then
     s1 = size(R3,1)
     s2 = size(R3,2)
     s3 = size(R3,3)
     Allocate(PR3(s3,s2,s1), STAT=iStatus)
     if ( iStatus /= 0 ) then
        call B_Set_Error(BE_RESHAPE, 'Unable to allocate in reshape for C like array output')
        return
     end if
     do i = 1, s1
        do j = 1, s2
           do k = 1, s3
              PR3(k,j,i) = R3(i,j,k)
           end do
        end do
     end do
  else if ( PRESENT (R4) .and. PRESENT(PR4) ) then
     s1 = size(R4,1)
     s2 = size(R4,2)
     s3 = size(R4,3)
     s4 = size(R4,4)
     Allocate(PR4(s4,s3,s2,s1), STAT=iStatus)
     call B_Set_Error(BE_RESHAPE, 'Unable to allocate in reshape for C like array output')
     if ( iStatus /= 0 ) then
           return
     end if
     do i = 1, s1
        do j = 1, s2
           do k = 1, s3
              do l = 1, s4
                 PR4(l,k,j,i) = R4(i,j,k,l)
              end do
           end do
        end do
     end do
  else if ( PRESENT (R5) .and. PRESENT(PR5) ) then
     s1 = size(R5,1)
     s2 = size(R5,2)
     s3 = size(R5,3)
     s4 = size(R5,4)
     s5 = size(R5,5)
     Allocate(PR5(s5,s4,s3,s2,s1), STAT=iStatus)
     call B_Set_Error(BE_RESHAPE, 'Unable to allocate in reshape for C like array output')
     if ( iStatus /= 0 ) then
           return
     end if
     do i = 1, s1
        do j = 1, s2
           do k = 1, s3
              do l = 1, s4
                 do m = 1, s5
                    PR5(m,l,k,j,i) = R5(i,j,k,l,m)
                 end do
              end do
           end do
        end do
     end do
  else 
     ! Error
     iStatus = BE_RESHAPE
     call B_Set_Error(BE_RESHAPE, 'Invalid call frame to reshape for C like array output')
  end if

  return

end SUBROUTINE _FUNCTION_NAME_RESHAPE_
#undef _DEFINE_RESHAPE_
#endif

#ifdef _IS_CHARACTER_TYPE_
#undef _IS_CHARACTER_TYPE_
#endif

#ifdef _NO_HIGHER_ORDER_XML_           
#undef _NO_HIGHER_ORDER_XML_
#endif

#ifdef _TYPE_CODE_
#undef _TYPE_CODE_
#endif

#ifdef _NO_HIGHER_ORDER_XML_
#undef _NO_HIGHER_ORDER_XML_
#endif

#undef _DIMENSION_
#undef _FUNCTION_NAME_RESHAPE_
#undef _FUNCTION_NAME_
#undef _FORMAT_NAME_
#undef _TYPE_SIZE_
#undef _TYPE_
