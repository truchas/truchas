!!CPP!! Include file for tbrook_module.F90
!!CPP!! Provides routines for output of native types

#ifdef TBROOK_COMMENT
!!CPP!! 
!!CPP!! Note that character arrays cannot be more 
!!CPP!! than one dimension
!!CPP!!

! This is the tbrook module include file.  This is where the bulk of
! the work gets done.  Unlike the streams module file which uses
! multiple inlclude statements for different dimensions,
! tbrook_module.F90 uses just a single include file, this one!  Do
! not be daunted by the large number of defines in this file.  They
! are there just to make my life easier.  Here is what the defines
! are, and what they do for you:
!
! _FUNCTION_NAME_ is the name of the subroutine that is being generated
! Typically this will be of the form TB_Write_DATATYPE_R[0-4]
!
! _FUNCTION_NAME_B_ is the name of the function that takes an integer
! as argument instead of a brook. This is needed currently for the
! internally defined brooks and is a kludge.  This will go away in
! later versions of this module.
!
! _BROOK_FUNC_ is the name of the function that needs to be invoked
! from the brooks module.
!
! _TYPE_ is the data type and is usually one of character, integer,
! logical, real or double precision.
!
! _DIMENSION_ is the dimension of the data being passed to above
! function.  Based on the dimension of the data we will decide
! whether to do a parallel collate or not. 
!
! _COLLATE_FUNC_ is the name of the function needed to do a collate. 
! The determination on collation or not is based on the
! dimensionality of the data being passed into this function.  If the
! dimension is 0, the _COLLATE_FUNC_ is defined and you will see the
! code at the bottom of this file.  All parallel collates are done
! within this one function defined by the 0 dimension include.  Note
! that we use pgslib to do parallel collates, and that collates are
! specific to the truchas code.  If you want to use tbrook_module.F90
! with your code, you will have to change the collate function and
! replace it with your own interface.  The advantage of using pgslib
! is that parallel and serial code can use the same flow and there is
! no branch required to separate them.  This makes programming at my
! end easier.  
!
! THere, thats all there is to it.  
!
! SOme other defines of interest:
! _PGSLIB_COLLATE_STRING_SIZE_ is the size of string for parallel
! character string collates.  For some reason truchas crashes at
! runtime if this number is not 70.  I shall track it down at some
! point in time, but do not deem it necessary to fix at this point
!
! _IS_CHARACTER_TYPE_ tells me whether the input is a character type
! or not. I'm not sure I need it, but have left it in there for
! historical reasons.  
!
! _DEFINE_COLLATE_FUNC_ is set with _DIMENSION_ = 0 and undefined
! othberwise.  The collations are defined within this collate
! function.  See note above for _PGSLIB_COLLATE_STRING_SIZE_
!
! _DO_PARALLEL_COLLATE_ is set if dimensionality is greater than 1
! and collations are done before printing if that is the case
!
! _RP_V[0-5], -RV_V[0-5] are defined for arguments to the mother
! collate function.
!
! -_end_-


#endif

#define _PGSLIB_COLLATE_STRING_SIZE_ 70

#ifndef _FUNCTION_NAME_ 
#error "_FUNCTION_NAME_ must be defined before including this file"
#endif

#ifndef _BROOK_FUNC_ 
#error "_BROOK_FUNC_ must be defined before including this file"
#endif

#ifndef _TYPE_
#error "_TYPE_  must be defined before including this file"
#endif

#define _DO_PARALLEL_COLLATE_

#if ( _DIMENSION_ == 0 )
#define _DIMENSION_S_ 
#define _DEFINE_COLLATE_FUNC_
#undef _DO_PARALLEL_COLLATE_

#elif ( _DIMENSION_ == 1 )
#define _RV_ v1
#define _RP_ pv1
#define _DIMENSION_S_ dimension(:),

#elif ( _DIMENSION_ == 2 )
#define _RV_ v2
#define _RP_ pv2
#define _DIMENSION_S_ dimension(:,:),

#elif ( _DIMENSION_ == 3 )
#define _RV_ v3
#define _RP_ pv3
#define _DIMENSION_S_ dimension(:,:,:),

#elif ( _DIMENSION_ == 4 )
#define _RV_ v4
#define _RP_ pv4
#define _DIMENSION_S_ dimension(:,:,:,:),

#elif ( _DIMENSION_ == 5 )
#define _RV_ v5
#define _RP_ pv5
#define _DIMENSION_S_ dimension(:,:,:,:,:),

#else
#error '_DIMENSION_ has invalid value'
#endif

#ifndef _COLLATE_FUNC_    
#error "_COLLATE_FUNC_ must be defined before calling this function"
#endif



  RECURSIVE SUBROUTINE _FUNCTION_NAME_ (TB_ID, Variable, FORMAT,  C_Array_Order, ADVANCE, OLDFILE, &
                                        XMLName, XMLAttributes, XMLDataFile, XMLDataFormat, &
                                        Scope, iStatus)
    !==================================================================
    ! Purpose(s):
    !   Write Variable to the output brook of structure TB(id)
    !   Status is returned in iStatus
    !==================================================================
    
    IMPLICIT NONE

    ! Arguments
    INTEGER,        INTENT(IN) :: TB_ID
    _TYPE_ , _DIMENSION_S_  INTENT(IN) :: Variable
    LOGICAL,          OPTIONAL :: ADVANCE
    CHARACTER(LEN=*), OPTIONAL :: FORMAT
    CHARACTER(LEN=*), OPTIONAL :: XMLName
    CHARACTER(LEN=*), OPTIONAL :: XMLAttributes
    CHARACTER(LEN=*), OPTIONAL :: XMLDataFile
    CHARACTER(LEN=*), OPTIONAL :: XMLDataFormat
    LOGICAL,          OPTIONAL :: C_Array_Order
    INTEGER,          OPTIONAL :: iStatus
    LOGICAL,          OPTIONAL :: OLDFILE
    INTEGER,          OPTIONAL :: SCOPE
    
    integer :: error_code 

    error_code = 0        
    if ( present (iStatus) ) error_code = iStatus

    IF ( .not. initialized ) then
       call TBrook_Initialize(iStatus)
       if ( iStatus /= 0 ) return
    END IF
        
    ! Insert the tag in the ascii base brook
    if ( TRIM(TBrook_Get_LastTag()) /= TRIM(TBrook_Get_Tag(TB_ID))) then
       call TBrook_Endline(BaseBrookAscii, iStatus=error_code)
       call TBrook_Write (BaseBrookAscii, Variable=TRIM(TBrook_Get_Tag(TB_ID)), Advance=.false., iStatus=error_code)
       call TBrook_Endline(BaseBrookAscii, iStatus=error_code)
       call TBrook_Set_LastTag(TB_ID)
    end if
    call _FUNCTION_NAME_B_ (TBrook_Get_Brook(TB_ID), Variable, Format, C_Array_Order, ADVANCE, OLDFILE, &
                            XMLName, XMLAttributes, XMLDataFile, XMLDataFormat, &
                            TB_ID, SCOPE, iStatus=error_code)
                         
    if ( present (iStatus) ) iStatus = error_code

    return

  END SUBROUTINE _FUNCTION_NAME_

  RECURSIVE SUBROUTINE _FUNCTION_NAME_B_ (B, Variable, FORMAT,  C_Array_Order, ADVANCE, OLDFILE, &
                                          XMLName, XMLAttributes, XMLDataFile, XMLDataFormat, &
                                          TB_ID, SCOPE, iStatus)
    !==================================================================
    ! Purpose(s):
    !   Write Variable to the output brook of structure TB(id)
    !   Status is returned in iStatus
    !==================================================================
    use brook_module, only : _BROOK_FUNC_ ,        &
                             brook_openxmltag,   &
                             b_stdout,           &
                             brook_closexmltag
    use parallel_info_module, only: p_info
    IMPLICIT NONE

    ! Arguments
    TYPE(Brook), TARGET, INTENT(IN) :: B
    INTEGER, OPTIONAL :: TB_ID
#ifdef _IS_CHARACTER_TYPE_
    _TYPE_ , _DIMENSION_S_  INTENT(IN) :: Variable
#else
    _TYPE_ , _DIMENSION_S_  INTENT(IN) :: Variable
#endif
    LOGICAL,          OPTIONAL :: ADVANCE
    CHARACTER(LEN=*), OPTIONAL :: FORMAT
    CHARACTER(LEN=*), OPTIONAL :: XMLName
    CHARACTER(LEN=*), OPTIONAL :: XMLAttributes
    CHARACTER(LEN=*), OPTIONAL :: XMLDataFile
    CHARACTER(LEN=*), OPTIONAL :: XMLDataFormat
    LOGICAL,          OPTIONAL :: C_Array_Order
    INTEGER,          OPTIONAL :: iStatus
    LOGICAL,          OPTIONAL :: OLDFILE
    INTEGER,          OPTIONAL :: SCOPE
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    TYPE(brook), pointer :: pB

#ifdef _IS_CHARACTER_TYPE_
#undef _TYPE_
#define _TYPE_ character(len = _PGSLIB_COLLATE_STRING_SIZE_ )
    _TYPE_ , _DIMENSION_S_  pointer :: VString
#endif

#ifdef _DO_PARALLEL_COLLATE_
#ifdef _IS_CHARACTER_TYPE_
    _TYPE_ , _DIMENSION_S_  pointer :: Var_Collate => NULL()
#else
    _TYPE_ , _DIMENSION_S_  pointer :: Var_Collate => NULL()
#endif
#else

#define Var_Collate Variable

#endif
    LOGICAL :: OLDFORMAT
    LOGICAL :: l_ADVANCE
    integer :: error_code
    INTEGER :: Local

    error_code = 0

    if (Present(SCOPE)) then
       Local = SCOPE
    else 
       Local = TB_SCOPE_GLOBAL
    end if


#ifdef _IS_CHARACTER_TYPE_
#if ( _DIMENSION_ == 0 )     
    allocate(VString)
#else
    allocate(VString(SIZE(Variable)))
#endif
    VString = Variable

#ifdef _DO_PARALLEL_COLLATE_
    if ( Local == TB_SCOPE_LOCAL ) then
       Var_Collate => VString
    else 
       Var_Collate => NULL()
       call _COLLATE_FUNC_ ( _RV_ = Variable, _RP_ = Var_Collate, Scope=Local, iStatus = error_code)
    endif
#endif

#else

#ifdef _DO_PARALLEL_COLLATE_
    Var_Collate => NULL()
    call _COLLATE_FUNC_ ( _RV_ = Variable, _RP_ = Var_Collate, Scope=Local, iStatus = error_code)
#endif

#endif
    if ( (Local == TB_SCOPE_LOCAL) .or. p_info%IOP) then
       
       if ( present(OLDFILE) ) then
          OLDFORMAT=OLDFILE
       else 
          OLDFORMAT = .false.
       end if

       pB => B
        if ( present(ADVANCE)) then
          l_ADVANCE = ADVANCE
       else 
          l_ADVANCE = .false.
       end if

111    FORMAT(a)
       if ( .not. initialized) then
          write(*,*) 'INITIALIZING FOR ', Variable
          pB=>B_Stdout
          IF ( .not. initialized ) then
             call TBrook_Initialize(iStatus)
             if ( iStatus /= 0 ) return
          END IF
       end if


       if ( .not. OLDFORMAT) then
          if ((.not.PRESENT(TB_ID)) .or. PRESENT(XMLName)) then
             call _BROOK_FUNC_ (B = pB,                         &
#if defined(_DO_PARALLEL_COLLATE_) || defined(_IS_CHARACTER_TYPE_)
                              Variable = Var_Collate,        &
#else
                              Variable = Variable,        &
#endif
                              Format=Format,                 &
                              C_Array_Order=C_Array_Order,   &
                              ADVANCE=ADVANCE,               &
                              XMLName=XMLName,               &
                              XMLAttributes=XMLAttributes,   &
                              XMLDataFile=XMLDataFile,       &
                              XMLDataFormat=XMLDataFormat,   &
                              iStatus=error_code)
          else 
             call _BROOK_FUNC_ (B = pB,                         &
#if defined(_DO_PARALLEL_COLLATE_) || defined(_IS_CHARACTER_TYPE_)
                              Variable = Var_Collate,        &
#else
                              Variable = Variable,        &
#endif
                              Format=Format,                 &
                              C_Array_Order=C_Array_Order,   &
                              ADVANCE=ADVANCE,               &
                              XMLName=TBrook_Get_Tag(TB_ID), &
                              XMLAttributes=XMLAttributes,   &
                              XMLDataFile=XMLDataFile,       &
                              XMLDataFormat=XMLDataFormat,   &
                              iStatus=error_code)
          end if
       else if ( no_file_set ) then 
          ! If not base fiel has been set, we return.
          ! This is a truchas specific if block and 
          ! will go away when we go to brooks as our 
          ! only means of output.
          return
       else 
          if ( present(TB_ID)) then
             if ( .not. TB_OPENTAG(TB_ID)) then
                call Brook_OpenXMLTag(TBrook_Get_BROOK(TB_ID), &
                                      XMLTag = TRIM(TBrook_Get_Tag(TB_ID)), &
                                      iStatus=error_code)
                TB_OPENTAG(TB_ID) = .true.
             end if
             call _BROOK_FUNC_ (B = TBrook_GET_BROOK(TB_ID),   &
                              Variable = Var_Collate,        & 
                              Format=Format,                 &
                              C_Array_Order=C_Array_Order,   &
                              ADVANCE=ADVANCE,               &
                              iStatus=error_code)
             if (l_ADVANCE) then
                call Brook_CloseXMLTag(TBrook_Get_BROOK(TB_ID), &
                                       XMLTag = TRIM(TBrook_Get_Tag(TB_ID)), &
                                       iStatus=error_code)
                TB_OPENTAG(TB_ID) = .false.
             end if
          end if
       end if
    end if

#ifdef _DO_PARALLEL_COLLATE_
    if ( associated(var_collate)) deallocate(Var_Collate)
#endif

#ifdef _IS_CHARACTER_TYPE_
    if ( Local /= TB_SCOPE_LOCAL ) then
       deallocate(VString)
    end if
#endif

    if ( present (iStatus) ) iStatus = error_code
  END SUBROUTINE _FUNCTION_NAME_B_

#ifdef _DEFINE_COLLATE_FUNC_
#undef _DEFINE_COLLATE_FUNC_

  ! The parallel collate
  SUBROUTINE _COLLATE_FUNC_ (iStatus, Scope, v1, pv1, v2, pv2, v3, pv3, v4, pv4, v5, pv5)

    use parallel_info_module, only: p_info
    use pgslib_module,        only: PGSLib_COLLATE, PGSLib_GLOBAL_SUM
    implicit none

    integer, intent(INOUT) :: iStatus
    integer :: scope

    _TYPE_ , dimension(:), optional          :: v1
    _TYPE_ , dimension(:), optional, pointer :: pv1

    _TYPE_ , dimension(:,:), optional          :: v2
    _TYPE_ , dimension(:,:), optional, pointer :: pv2

    _TYPE_ , dimension(:,:,:), optional          :: v3
    _TYPE_ , dimension(:,:,:), optional, pointer :: pv3

    _TYPE_ , dimension(:,:,:,:), optional          :: v4
    _TYPE_ , dimension(:,:,:,:), optional, pointer :: pv4

    _TYPE_ , dimension(:,:,:,:,:), optional          :: v5
    _TYPE_ , dimension(:,:,:,:,:), optional, pointer :: pv5

    ! Local Variables
    integer :: s1, s2, s3, s4, s5
    integer :: s1_tot, s2_tot, s3_tot, s4_tot, s5_tot
    integer :: i, j, k, l

    if ( iStatus /= 0 ) return

    if ( present (v1) .and. present(pv1) ) then
       s1 = SIZE(v1, DIM=1)
       if ( Scope == TB_SCOPE_LOCAL ) then
          s1_tot = s1
       else
          s1_tot = PGSLib_GLOBAL_SUM(s1)
       end if
       ! Allocate and collate the total vector

       if (p_info%IOP .or. (Scope == TB_SCOPE_LOCAL)) then
          ALLOCATE (pv1(s1_tot), STAT=iStatus)
       else
          ALLOCATE (pv1(0))
       end if
       if ( Scope == TB_SCOPE_LOCAL) then
          pv1 = v1
       else
          call PGSLib_COLLATE (pv1, v1)
       end if
    else if ( present (v2) .and. present(pv2) ) then
       s1 = SIZE(v2, DIM=1)
       s2 = SIZE(v2, DIM=2)

       if(Scope == TB_SCOPE_LOCAL) then
          s2_tot = s2
       else
          s2_tot = PGSLib_GLOBAL_SUM(s2)
       end if

       ! Allocate and collate the total vector
       if (p_info%IOP .or. (Scope == TB_SCOPE_LOCAL)) then
          ALLOCATE (pv2(s1, s2_tot), STAT=iStatus)
       else
          ALLOCATE (pv2(s1,0), STAT=iStatus)
       end if

       if ( Scope == TB_SCOPE_LOCAL) then
          pv2 = v2
       else
          do i = 1,s1
             call PGSLib_COLLATE (pv2(i,:), v2(i,:))
          end do
       end if

    else if ( present (v3) .and. present(pv3) ) then
       s1 = SIZE(v3, DIM=1)
       s2 = SIZE(v3, DIM=2)
       s3 = SIZE(v3, DIM=3)

       if(Scope == TB_SCOPE_LOCAL) then
          s3_tot = s3
       else
          s3_tot = PGSLib_GLOBAL_SUM(s3)
       end if

       ! Allocate and collate the total vector
       if (p_info%IOP .or. (Scope == TB_SCOPE_LOCAL)) then
          ALLOCATE (pv3(s1, s2, s3_tot), STAT=iStatus)
       else
          ALLOCATE (pv3(s1, s2, 0), STAT=iStatus)
       end if

       if ( Scope == TB_SCOPE_LOCAL) then
          pv3 = v3
       else
          do i = 1,s1
             do j = 1, s2
                call PGSLib_COLLATE (pv3(i,j,:), v3(i,j,:))
             end do
          end do
       end if

    else if ( present (v4) .and. present(pv4) ) then
       s1 = SIZE(v4, DIM=1)
       s2 = SIZE(v4, DIM=2)
       s3 = SIZE(v4, DIM=3)
       s4 = SIZE(v4, DIM=4)

       if(Scope == TB_SCOPE_LOCAL) then
          s4_tot = s4
       else
          s4_tot = PGSLib_GLOBAL_SUM(s4)
       end if

       ! Allocate and collate the total vector
       if (p_info%IOP .or. (Scope == TB_SCOPE_LOCAL)) then
          ALLOCATE (pv4(s1, s2, s3, s4_tot), STAT=iStatus)
       else
          ALLOCATE (pv4(s1, s2, s3, 0), STAT=iStatus)
       end if

       if ( Scope == TB_SCOPE_LOCAL) then
          pv4 = v4
       else
          do i = 1,s1
             do j = 1, s2
                do k = 1, s3
                   call PGSLib_COLLATE (pv4(i,j,k,:), v4(i,j,k,:))
                end do
             end do
          end do
       end if

    else if ( present (v5) .and. present(pv5) ) then

       s1 = SIZE(v5, DIM=1)
       s2 = SIZE(v5, DIM=2)
       s3 = SIZE(v5, DIM=3)
       s4 = SIZE(v5, DIM=4)
       s5 = SIZE(v5, DIM=5)

       if(Scope == TB_SCOPE_LOCAL) then
          s5_tot = s5
       else
          s5_tot = PGSLib_GLOBAL_SUM(s5)
       end if





       ! Allocate and collate the total vector
       if (p_info%IOP .or. (Scope == TB_SCOPE_LOCAL)) then
          ALLOCATE (pv5(s1, s2, s3, s4, s5_tot), STAT=iStatus)
       else
          ALLOCATE (pv5(s1, s2, s3, s4, 0), STAT=iStatus)
       end if

       if ( Scope == TB_SCOPE_LOCAL) then
          pv5 = v5
       else
          do i = 1,s1
             do j = 1, s2
                do k = 1, s3
                   do l = 1, s4
                      call PGSLib_COLLATE (pv5(i,j,k,l,:), v5(i,j,k,l,:))
                   end do
                end do
             end do
          end do
       end if
    else
       iStatus=1
    end if

    RETURN
  END SUBROUTINE _COLLATE_FUNC_
#endif

#ifdef _DO_PARALLEL_COLLATE_
#undef _DO_PARALLEL_COLLATE_
#else
#undef Var_Collate
#endif

#ifdef _IS_CHARACTER_TYPE_
#undef _IS_CHARACTER_TYPE_
#endif

#undef _RV_
#undef _RP_
#undef _TYPE_
#undef _DIMENSION_
#undef _DIMENSION_S_
#undef _BROOK_FUNC_
#undef _COLLATE_FUNC_
#undef _FUNCTION_NAME_
#undef _FUNCTION_NAME_B_
#undef _PGSLIB_COLLATE_STRING_SIZE_
