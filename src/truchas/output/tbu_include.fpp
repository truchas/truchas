#ifndef _TBU_FILEENTRY_FUNC_
#error "_TBU_FILEENTRY_FUNC_ must be defined before including this file"
#endif

#ifndef _TYPE_
#error "_TYPE_  must be defined before including this file"
#endif

#ifndef _DIMENSION_
#error "_DIMENSION_  must be defined before including this file"
#endif

#ifndef _TYPE_CODE_
#error "_TYPE_CODE_  must be defined before including this file"
#endif

#define _HIGH_DIM_ .true.

#if ( _DIMENSION_ == 0 )
#define _TYPE_DIM_ _TYPE_
#undef _HIGH_DIM_
#define _HIGH_DIM_ .false.

#elif ( _DIMENSION_ == 1 )
#define _TYPE_DIM_ _TYPE_, dimension(:)

#elif ( _DIMENSION_ == 2 )
#define _TYPE_DIM_ _TYPE_, dimension(:,:)

#elif ( _DIMENSION_ == 3 )
#define _TYPE_DIM_ _TYPE_, dimension(:,:,:)

#elif ( _DIMENSION_ == 4 )
#define _TYPE_DIM_ _TYPE_, dimension(:,:,:,:)

#elif ( _DIMENSION_ == 5 )
#define _TYPE_DIM_ _TYPE_, dimension(:,:,:,:,:)

#else
#error '_DIMENSION_ has invalid value'
#endif
    SUBROUTINE _TBU_FILEENTRY_FUNC_ (B_Base, B_LookAside, XMLName, Variable, iStatus, &
                                     ADVANCE, MESHLABEL, XMLATTRIBUTES, SCOPE, MAP, C_Array_Order, &
                                     str )
      use parallel_info_module, only: p_info

      type(Brook), target :: B_Base
      type(Brook), target :: B_LookAside
      character(LEN=*), intent(in)    :: XMLName
      _TYPE_DIM_ , intent(in) :: Variable
      integer, intent(inout) :: iStatus

      Logical,          optional, intent(in) :: Advance
      character(len=*), optional, intent(in) :: MESHLABEL
      character(len=*), optional, intent(in) :: XMLATTRIBUTES
      integer,          optional, intent(in) :: Scope
      character(len=*), optional, intent(in) :: map
      logical,          optional, intent(in) :: C_Array_Order
      character(len=*), optional, intent(out) :: str

      character(LEN=256) :: tstring
      character(LEN=1024) :: AttrString
      character(LEN=1024) :: BrookAttrString
      integer, dimension(1) :: s1
      integer               :: s1sum
      integer               :: myScope
      integer, dimension(:), pointer :: myShape
      character(len=256) :: mstring

      ! Do nothing if there is an error
      if ( iStatus /= 0 ) return

      mstring = 'DefaultMesh'
      if (present(MESHLABEL)) then
         mstring = TRIM(ADJUSTL(MESHLABEL))
      end if

      myScope = TBrook_Get_Scope(scope=scope)
      if ( _HIGH_DIM_ .and. myScope == TB_SCOPE_GLOBAL) then
         ! Global scope array output on all processors
         ! Get collated size in case of parallel writes
         allocate(myShape(size(SHAPE(Variable))), stat=iStatus)
         myshape = SHAPE(Variable)
         s1 = myshape(size(myshape))
         s1sum = PGSLib_GLOBAL_SUM(s1)
         myshape(size(myshape)) = s1sum

         if (p_info%IOP) then
            call Brook_Get_Standard_Attributes(str           = BrookAttrString,       &
                                               name          = XMLName,               &
                                               type_code     = _TYPE_CODE_,           &
                                               rank          = SIZE(SHAPE(Variable)), &
                                               shape         = myshape,               &
                                               C_Array_ORder = C_Array_Order,         &
                                               Map           = Map,                   &
                                               offset        = BROOK_OFFSET(B_LookAside), &
                                               Mesh          = TRIM(mstring))
         else
            BrookAttrString = " " 
         end if
                                            
         deallocate(myshape)
      else if (_HIGH_DIM_) then
         ! Local scope Array output on all processors
         call Brook_Get_Standard_Attributes(str           = BrookAttrString,       &
                                            name          = XMLName,               &
                                            type_code     = _TYPE_CODE_,           &
                                            rank          = SIZE(SHAPE(Variable)), &
                                            shape         = Shape(Variable),       &
                                            C_Array_ORder = C_Array_Order,         &
                                            Map           = Map,                   &
                                            offset        = BROOK_OFFSET(B_LookAside), &
                                            Mesh          = TRIM(mstring))
      else if ( myScope == TB_SCOPE_LOCAL .or. p_info%IOP) then
         ! Global scope constant output on an IO processor
         ! OR a local scope constant output on any processor
         call Brook_Get_Standard_Attributes(str           = BrookAttrString,       &
                                            name          = XMLName,               &
                                            type_code     = _TYPE_CODE_,           &
                                            rank          = SIZE(SHAPE(Variable)), &
                                            C_Array_ORder = C_Array_Order,         &
                                            Map           = Map,                   &
                                            offset        = BROOK_OFFSET(B_LookAside), &
                                            Mesh          = TRIM(mstring))
      else 
         ! Global scope constant on a non-IO processor
         BrookAttrString = " "
      end if

1     FORMAT(a,' ', a)
      if (present(XMLATTRIBUTES)) then
         write(AttrString,1) TRIM(ADJUSTL(BrookAttrString)), TRIM(ADJUSTL(XMLATTRIBUTES))
      else 
         write(AttrString,1) TRIM(ADJUSTL(BrookAttrString)), ' '
      end if

      if(p_info%iop .or. (myScope == TB_SCOPE_LOCAL)) then
         if (iStatus == 0 ) call TBrook_WriteXMLTag(B=B_Base, XMLTag="FILEVAR", XMLAttributes=AttrString, &
                                                    SCOPE=SCOPE, iStatus=iStatus)
      end if

      tstring = ' '
      if (present(MESHLABEL)) then
         tstring = TRIM(tstring) // ' MESH="' // TRIM(ADJUSTL(MESHLABEL)) //'"'
      else
         tstring = TRIM(tstring) // ' MESH="DefaultMesh"'
      end if
      if (present(Map)) then
         tstring = TRIM(tstring) // ' MAP="' // TRIM(ADJUSTL(Map)) //'"'
      else 
         tstring = TRIM(tstring) // ' MAP="Cells"'
      end if
      if (iStatus == 0 ) call TBrook_Write(B=B_LookAside, Variable=Variable, XMLName=XMLName, XMLAttributes=TRIM(tstring), &
                                           Advance=Advance, SCOPE=myScope, iStatus=iStatus)
      if(present(str)) then
         str = TRIM(ADJUSTL(BrookAttrString)) 
      end if
      
      return
    end SUBROUTINE _TBU_FILEENTRY_FUNC_

#undef _TYPE_
#undef _TYPE_DIM_
#undef _HIGH_DIM_
#undef _TYPE_CODE_
#undef _DIMENSION_
#undef _TBU_FILEENTRY_FUNC_    

