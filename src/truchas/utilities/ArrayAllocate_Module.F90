MODULE ARRAYALLOCATE_MODULE
  !=======================================================================
  ! Purpose(s):
  !   A module to handle array allocation and deallocation.
  !
  !   Usage:
  !      Use ArrayALLOCATE
  !      call ArrayCreate(ptr, lbn, ubn, ..., [string])
  !      call ArrayDestroy(ptr, [string])
  !
  !   Variables:
  !      Type, pointer, dimension(:) :: ptr
  !         ptr is a pointer to the array to be allocated or deallocated
  !      integer, intent(IN) :: lbn
  !         lbn is the lower bound of the array - multiple 'n's are used for 
  !         multidimensional arrays
  !      integer, intent(IN) :: ubn
  !         ubn is the upper bound of the array - multiple 'n's are used for 
  !         multidimensional arrays
  !      character(*), intent(IN), optional :: string
  !         string is an optional error string to be included in an error message
  !
  !   Example:
  !      allocate an array of rank 2, dimension 3x7
  !
  !      Use ArrayALLOCATE
  !      real, Pointer, Dimension(:,:) :: p2
  !      call ArrayCreate(p2,1,3,1,7,'P2_3x7')
  !          do stuff here
  !      call ArrayDestroy(p2,'P2')
  !
  !   Method:
  !      ArrayCreate is passed a pointer to an array, and the requisite number
  !      of array dimensions, and allocates the array, checking for errors.
  !      ArrayDestroy deallocates the array, again checking for errors.
  !      Generic external procedures are used.
  !
  !   Making New Routines To Handle New Types:
  !      To make Create and Destroy routines for new data types, you need to:
  !
  !      1)  Code up new Create and Destroy routines.  The easiest way to do this
  !          is to copy the Create and Destroy routines for an array of the desired
  !          rank from those already created, and change the type declarations.
  !
  !      2) Add the routine names to the interface blocks.
  !
  !      Please make the Destroy routine immediately follow the corresponding
  !      Create routine in this source file.
  !
  !      The names for the routines are of the form ArrayCreate_type_Rn, where
  !      type is replaced by the data type and n is replaced by the rank of the
  !      array.  For example, ArrayCreate_real_R3 is an array of type real,
  !      rank 3.
  !
  !      if the type is defined in a module by "use" association, use the only:
  !      qualifier to include only the definition of the needed type.
  !
  !   Notes:
  !      For a sample of code that was used to test these routines, see
  !      Subroutine ArrayALLOCATE_Test, at the end of this file.
  !
  !      New Create and Destroy routines will have to be coded for each new
  !      type and shape of array.  The coding is exceptionally easy,
  !      mostly cut and paste.
  !
  !      if you pass a pointer to ArrayCreate, and it is already pointing at an
  !      array, you are not going to like the results.  ArrayCreate will
  !      allocate a new array, and point the pointer at it.  You will lose the
  !      contents of the data space that the pointer was pointing at, unless
  !      you really know what you are doing.  Not recommended.
  !
  ! Contains: ArrayCreate_real_R1
  !           ArrayCreate_real_R2
  !           ArrayCreate_real_R3
  !           ArrayCreate_Log_Type_R1
  !           ArrayCreate_Log_Type_R2
  !           ArrayCreate_Log_Type_R3
  !           ArrayCreate_Int_Type_R1
  !           ArrayCreate_Int_Type_R2
  !           ArrayCreate_Int_Type_R3
  !           ArrayCreate_real_Type_R1
  !           ArrayCreate_real_Type_R2
  !           ArrayCreate_real_Type_R3
  !
  !           ArrayDestroy_real_R1
  !           ArrayDestroy_real_R2
  !           ArrayDestroy_real_R3
  !           ArrayDestroy_Log_Type_R1
  !           ArrayDestroy_Log_Type_R2
  !           ArrayDestroy_Log_Type_R3
  !           ArrayDestroy_Int_Type_R1
  !           ArrayDestroy_Int_Type_R2
  !           ArrayDestroy_Int_Type_R3
  !           ArrayDestroy_real_Type_R1
  !           ArrayDestroy_real_Type_R2
  !           ArrayDestroy_real_Type_R3
  !
  ! Author(s): Bryan R. Lally, LANL ESA-EPE (lally@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only : string_len
  use truchas_logging_services
  implicit none
  private

  ! Public Subroutines
  public :: ArrayCreate, ArrayDestroy

  ! Interface Procedures
  INTERFACE ArrayCreate
     MODULE PROCEDURE ArrayCreate_real_R1
     MODULE PROCEDURE ArrayCreate_real_R2
     MODULE PROCEDURE ArrayCreate_real_R3
     MODULE PROCEDURE ArrayCreate_Log_Type_R1
     MODULE PROCEDURE ArrayCreate_Log_Type_R2
     MODULE PROCEDURE ArrayCreate_Log_Type_R3
     MODULE PROCEDURE ArrayCreate_Int_Type_R1
     MODULE PROCEDURE ArrayCreate_Int_Type_R2
     MODULE PROCEDURE ArrayCreate_Int_Type_R3
     MODULE PROCEDURE ArrayCreate_real_Type_R1
     MODULE PROCEDURE ArrayCreate_real_Type_R2
     MODULE PROCEDURE ArrayCreate_real_Type_R3
  END INTERFACE

  INTERFACE ArrayDestroy
     MODULE PROCEDURE ArrayDestroy_real_R1
     MODULE PROCEDURE ArrayDestroy_real_R2
     MODULE PROCEDURE ArrayDestroy_real_R3
     MODULE PROCEDURE ArrayDestroy_Log_Type_R1
     MODULE PROCEDURE ArrayDestroy_Log_Type_R2
     MODULE PROCEDURE ArrayDestroy_Log_Type_R3
     MODULE PROCEDURE ArrayDestroy_Int_Type_R1
     MODULE PROCEDURE ArrayDestroy_Int_Type_R2
     MODULE PROCEDURE ArrayDestroy_Int_Type_R3
     MODULE PROCEDURE ArrayDestroy_real_Type_R1
     MODULE PROCEDURE ArrayDestroy_real_Type_R2
     MODULE PROCEDURE ArrayDestroy_real_Type_R3
  END INTERFACE

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE ARRAYCREATE_REAL_R1 (ArrayPtr, lb1, ub1, string)
    !=======================================================================
    ! Purpose(s):
    !   ALLOCATE real array of rank 1.
    !=======================================================================
    ! Argument List
    real, pointer, dimension(:) :: ArrayPtr
    integer, intent(IN) :: lb1
    integer, intent(IN) :: ub1
    character(*), intent(IN), optional :: string

    ! Local Variables
    character(64) :: routine = 'ArrayCreate_real_R1'
    character(string_len) :: error_string
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to allocate space
    ALLOCATE (ArrayPtr(lb1:ub1), STAT = status)
    

    if (PRESENT(string)) then
       error_string = string
    else
       error_string = routine
    end if

    if (status /= 0) call TLS_panic (trim(error_string) // ': allocate failed')

  END SUBROUTINE ARRAYCREATE_REAL_R1

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYDESTROY_REAL_R1 (ArrayPtr, string)
    !=======================================================================
    ! Purpose(s):
    !   DEALLOCATE real array of rank 1.
    !=======================================================================
    ! Argument List
    real, pointer, dimension(:) :: ArrayPtr
    character(*), intent(IN), optional :: string

    ! Local Variables
    character(64) :: routine = 'ArrayDestroy_real_R1'
    character(string_len) :: error_string
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to deallocate space
    DEALLOCATE (ArrayPtr)

  END SUBROUTINE ARRAYDESTROY_REAL_R1

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYCREATE_REAL_R2 (ArrayPtr, lb1, ub1, lb2, ub2, string)
    !=======================================================================
    ! Purpose(s):
    !   ALLOCATE real array of rank 2.
    !=======================================================================
    ! Argument List
    real, pointer, dimension(:,:) :: ArrayPtr
    integer, intent(IN) :: lb1
    integer, intent(IN) :: ub1
    integer, intent(IN) :: lb2
    integer, intent(IN) :: ub2
    character(*), intent(IN), optional :: string

    ! Local Variables
    character(64) :: routine = 'ArrayCreate_real_R2'
    character(string_len) :: error_string
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to allocate space
    ALLOCATE (ArrayPtr(lb1:ub1,lb2:ub2), STAT = status)

    if (PRESENT(string)) then
       error_string = string
    else
       error_string = routine
    end if

    if (status /= 0) call TLS_panic (trim(error_string) // ': allocate failed')

  END SUBROUTINE ARRAYCREATE_REAL_R2

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYDESTROY_REAL_R2 (ArrayPtr, string)
    !=======================================================================
    ! Purpose(s):
    !   DEALLOCATE real array of rank 2.
    !=======================================================================
    ! Argument List
    real, pointer, dimension(:,:) :: ArrayPtr
    character(*), intent(IN), optional :: string

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to deallocate space
    DEALLOCATE (ArrayPtr)

  END SUBROUTINE ARRAYDESTROY_REAL_R2

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYCREATE_REAL_R3 (ArrayPtr, lb1, ub1, lb2, ub2, lb3, ub3, string)
    !=======================================================================
    ! Purpose(s):
    !   ALLOCATE real array of rank 3.
    !=======================================================================
    ! Argument List
    real, pointer, dimension(:,:,:) :: ArrayPtr
    integer, intent(IN) :: lb1
    integer, intent(IN) :: ub1
    integer, intent(IN) :: lb2
    integer, intent(IN) :: ub2
    integer, intent(IN) :: lb3
    integer, intent(IN) :: ub3
    character(*), intent(IN), optional :: string

    ! Local Variables
    character(64) :: routine = 'ArrayCreate_real_R3'
    character(string_len) :: error_string
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to allocate space
    ALLOCATE (ArrayPtr(lb1:ub1,lb2:ub2,lb3:ub3), STAT = status)

    if (PRESENT(string)) then
       error_string = string
    else
       error_string = routine
    end if

    if (status /= 0) call TLS_panic (trim(error_string) // ': allocate failed')

  END SUBROUTINE ARRAYCREATE_REAL_R3

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYDESTROY_REAL_R3 (ArrayPtr, string)
    !=======================================================================
    ! Purpose(s):
    !   DEALLOCATE real array of rank 3.
    !=======================================================================
    ! Argument List
    real, pointer, dimension(:,:,:) :: ArrayPtr
    character(*), intent(IN), optional :: string

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to deallocate space
    DEALLOCATE (ArrayPtr)

  END SUBROUTINE ARRAYDESTROY_REAL_R3

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYCREATE_LOG_TYPE_R1 (ArrayPtr, lb1, ub1, string)
    !=======================================================================
    ! Purpose(s):
    !   ALLOCATE log_type array of rank 1.
    !=======================================================================

    ! Argument List
    logical, pointer, dimension(:) :: ArrayPtr
    integer, intent(IN) :: lb1
    integer, intent(IN) :: ub1
    character(*), intent(IN), optional :: string

    ! Local Variables
    character(64) :: routine = 'ArrayCreate_Log_Type_R1'
    character(string_len) :: error_string
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to allocate space
    ALLOCATE (ArrayPtr(lb1:ub1), STAT = status)

    if (PRESENT(string)) then
       error_string = string
    else
       error_string = routine
    end if

    if (status /= 0) call TLS_panic (trim(error_string) // ': allocate failed')

  END SUBROUTINE ARRAYCREATE_LOG_TYPE_R1

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYDESTROY_LOG_TYPE_R1 (ArrayPtr, string)
    !=======================================================================
    ! Purpose(s):
    !   DEALLOCATE log_type array of rank 1.
    !=======================================================================

    ! Argument List
    logical, pointer, dimension(:) :: ArrayPtr
    character(*), intent(IN), optional :: string

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to deallocate space
    DEALLOCATE (ArrayPtr)

  END SUBROUTINE ARRAYDESTROY_LOG_TYPE_R1

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYCREATE_LOG_TYPE_R2 (ArrayPtr, lb1, ub1, lb2, ub2, string)
    !=======================================================================
    ! Purpose(s):
    !   ALLOCATE log_type array of rank 2.
    !=======================================================================

    ! Argument List
    logical, pointer, dimension(:,:) :: ArrayPtr
    integer, intent(IN) :: lb1
    integer, intent(IN) :: ub1
    integer, intent(IN) :: lb2
    integer, intent(IN) :: ub2
    character(*), intent(IN), optional :: string

    ! Local Variables
    character(64) :: routine = 'ArrayCreate_Log_Type_R2'
    character(string_len) :: error_string
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to allocate space
    ALLOCATE (ArrayPtr(lb1:ub1,lb2:ub2), STAT = status)

    if (PRESENT(string)) then
       error_string = string
    else
       error_string = routine
    end if

    if (status /= 0) call TLS_panic (trim(error_string) // ': allocate failed')

  END SUBROUTINE ARRAYCREATE_LOG_TYPE_R2

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYDESTROY_LOG_TYPE_R2 (ArrayPtr, string)
    !=======================================================================
    ! Purpose(s):
    !   DEALLOCATE log_type array of rank 2.
    !=======================================================================

    ! Argument List
    logical, pointer, dimension(:,:) :: ArrayPtr
    character(*), intent(IN), optional :: string

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to deallocate space
    DEALLOCATE (ArrayPtr)

  END SUBROUTINE ARRAYDESTROY_LOG_TYPE_R2

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYCREATE_LOG_TYPE_R3 (ArrayPtr, lb1, ub1, lb2, ub2, lb3, ub3, string)
    !=======================================================================
    ! Purpose(s):
    !   ALLOCATE log_type array of rank 3.
    !=======================================================================

    ! Argument List
    logical, pointer, dimension(:,:,:) :: ArrayPtr
    integer, intent(IN) :: lb1
    integer, intent(IN) :: ub1
    integer, intent(IN) :: lb2
    integer, intent(IN) :: ub2
    integer, intent(IN) :: lb3
    integer, intent(IN) :: ub3
    character(*), intent(IN), optional :: string

    ! Local Variables
    character(64) :: routine = 'ArrayCreate_Log_Type_R3'
    character(string_len) :: error_string
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to allocate space
    ALLOCATE (ArrayPtr(lb1:ub1,lb2:ub2,lb3:ub3), STAT = status)

    if (PRESENT(string)) then
       error_string = string
    else
       error_string = routine
    end if

    if (status /= 0) call TLS_panic (trim(error_string) // ': allocate failed')

  END SUBROUTINE ARRAYCREATE_LOG_TYPE_R3

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYDESTROY_LOG_TYPE_R3 (ArrayPtr, string)
    !=======================================================================
    ! Purpose(s):
    !   DEALLOCATE log_type array of rank 3.
    !=======================================================================

    ! Argument List
    logical, pointer, dimension(:,:,:) :: ArrayPtr
    character(*), intent(IN), optional :: string

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to deallocate space
    DEALLOCATE (ArrayPtr)

  END SUBROUTINE ARRAYDESTROY_LOG_TYPE_R3

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYCREATE_INT_TYPE_R1 (ArrayPtr, lb1, ub1, string)
    !=======================================================================
    ! Purpose(s):
    !   ALLOCATE int_type array of rank 1.
    !=======================================================================

    ! Argument List
    integer, pointer, dimension(:) :: ArrayPtr
    integer, intent(IN) :: lb1
    integer, intent(IN) :: ub1
    character(*), intent(IN), optional :: string

    ! Local Variables
    character(64) :: routine = 'ArrayCreate_Int_Type_R1'
    character(string_len) :: error_string
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


    ! attempt to allocate space
    ALLOCATE (ArrayPtr(lb1:ub1), STAT = status)

    if (PRESENT(string)) then
       error_string = string
    else
       error_string = routine
    end if

    if (status /= 0) call TLS_panic (trim(error_string) // ': allocate failed')

  END SUBROUTINE ARRAYCREATE_INT_TYPE_R1

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYDESTROY_INT_TYPE_R1 (ArrayPtr, string)
    !=======================================================================
    ! Purpose(s):
    !   DEALLOCATE int_type array of rank 1.
    !=======================================================================

    ! Argument List
    integer, pointer, dimension(:) :: ArrayPtr
    character(*), intent(IN), optional :: string

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to deallocate space
    DEALLOCATE (ArrayPtr)

  END SUBROUTINE ARRAYDESTROY_INT_TYPE_R1

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYCREATE_INT_TYPE_R2 (ArrayPtr, lb1, ub1, lb2, ub2, string)
    !=======================================================================
    ! Purpose(s):
    !   ALLOCATE int_type array of rank 2.
    !=======================================================================

    ! Argument List
    integer, pointer, dimension(:,:) :: ArrayPtr
    integer, intent(IN) :: lb1
    integer, intent(IN) :: ub1
    integer, intent(IN) :: lb2
    integer, intent(IN) :: ub2
    character(*), intent(IN), optional :: string

    ! Local Variables
    character(64) :: routine = 'ArrayCreate_Int_Type_R2'
    character(string_len) :: error_string
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to allocate space
    ALLOCATE (ArrayPtr(lb1:ub1,lb2:ub2), STAT = status)

    if (PRESENT(string)) then
       error_string = string
    else
       error_string = routine
    end if

    if (status /= 0) call TLS_panic (trim(error_string) // ': allocate failed')

  END SUBROUTINE ARRAYCREATE_INT_TYPE_R2

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYDESTROY_INT_TYPE_R2 (ArrayPtr, string)
    !=======================================================================
    ! Purpose(s):
    !   DEALLOCATE int_type array of rank 2.
    !=======================================================================

    ! Argument List
    integer, pointer, dimension(:,:) :: ArrayPtr
    character(*), intent(IN), optional :: string

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to deallocate space
    DEALLOCATE (ArrayPtr)

  END SUBROUTINE ARRAYDESTROY_INT_TYPE_R2

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYCREATE_INT_TYPE_R3 (ArrayPtr, lb1, ub1, lb2, ub2, lb3, ub3, string)
    !=======================================================================
    ! Purpose(s):
    !   ALLOCATE int_type array of rank 3.
    !=======================================================================

    ! Argument List
    integer, pointer, dimension(:,:,:) :: ArrayPtr
    integer, intent(IN) :: lb1
    integer, intent(IN) :: ub1
    integer, intent(IN) :: lb2
    integer, intent(IN) :: ub2
    integer, intent(IN) :: lb3
    integer, intent(IN) :: ub3
    character(*), intent(IN), optional :: string

    ! Local Variables
    character(64) :: routine = 'ArrayCreate_Int_Type_R3'
    character(string_len) :: error_string
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to allocate space
    ALLOCATE (ArrayPtr(lb1:ub1,lb2:ub2,lb3:ub3), STAT = status)

    if (PRESENT(string)) then
       error_string = string
    else
       error_string = routine
    end if

    if (status /= 0) call TLS_panic (trim(error_string) // ': allocate failed')

  END SUBROUTINE ARRAYCREATE_INT_TYPE_R3

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYDESTROY_INT_TYPE_R3 (ArrayPtr, string)
    !=======================================================================
    ! Purpose(s):
    !   DEALLOCATE int_type array of rank 3.
    !=======================================================================

    ! Argument List
    integer, pointer, dimension(:,:,:) :: ArrayPtr
    character(*), intent(IN), optional :: string

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to deallocate space
    DEALLOCATE (ArrayPtr)

  END SUBROUTINE ARRAYDESTROY_INT_TYPE_R3

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYCREATE_REAL_TYPE_R1 (ArrayPtr, lb1, ub1, string)
    !=======================================================================
    ! Purpose(s):
    !   ALLOCATE real_type array of rank 1.
    !=======================================================================

    ! Argument List
    real(r8), pointer, dimension(:) :: ArrayPtr
    integer, intent(IN) :: lb1
    integer, intent(IN) :: ub1
    character(*), intent(IN), optional :: string

    ! Local Variables
    character(64) :: routine = 'ArrayCreate_real_Type_R1'
    character(string_len) :: error_string
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to allocate space
    ALLOCATE (ArrayPtr(lb1:ub1), STAT = status)

    if (PRESENT(string)) then
       error_string = string
    else
       error_string = routine
    end if

    if (status /= 0) call TLS_panic (trim(error_string) // ': allocate failed')

  END SUBROUTINE ARRAYCREATE_REAL_TYPE_R1

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYDESTROY_REAL_TYPE_R1 (ArrayPtr, string)
    !=======================================================================
    ! Purpose(s):
    !   DEALLOCATE real_type array of rank 1.
    !=======================================================================

    ! Argument List
    real(r8), pointer, dimension(:) :: ArrayPtr
    character(*), intent(IN), optional :: string

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to deallocate space
    DEALLOCATE (ArrayPtr)

  END SUBROUTINE ARRAYDESTROY_REAL_TYPE_R1

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYCREATE_REAL_TYPE_R2 (ArrayPtr, lb1, ub1, lb2, ub2, string)
    !=======================================================================
    ! Purpose(s):
    !   ALLOCATE real_type array of rank 2.
    !=======================================================================

    ! Argument List
    real(r8), pointer, dimension(:,:) :: ArrayPtr
    integer, intent(IN) :: lb1
    integer, intent(IN) :: ub1
    integer, intent(IN) :: lb2
    integer, intent(IN) :: ub2
    character(*), intent(IN), optional :: string

    ! Local Variables
    character(64) :: routine = 'ArrayCreate_real_Type_R2'
    character(string_len) :: error_string
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to allocate space
    ALLOCATE (ArrayPtr(lb1:ub1,lb2:ub2), STAT = status)

    if (PRESENT(string)) then
       error_string = string
    else
       error_string = routine
    end if

    if (status /= 0) call TLS_panic (trim(error_string) // ': allocate failed')

  END SUBROUTINE ARRAYCREATE_REAL_TYPE_R2

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYDESTROY_REAL_TYPE_R2 (ArrayPtr, string)
    !=======================================================================
    ! Purpose(s):
    !   DEALLOCATE real_type array of rank 2.
    !=======================================================================

    ! Argument List
    real(r8), pointer, dimension(:,:) :: ArrayPtr
    character(*), intent(IN), optional :: string

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to deallocate space
    DEALLOCATE (ArrayPtr)

  END SUBROUTINE ARRAYDESTROY_REAL_TYPE_R2

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYCREATE_REAL_TYPE_R3 (ArrayPtr, lb1, ub1, lb2, ub2, lb3, ub3, string)
    !=======================================================================
    ! Purpose(s):
    !   ALLOCATE real_type array of rank 3.
    !=======================================================================

    ! Argument List
    real(r8), pointer, dimension(:,:,:) :: ArrayPtr
    integer, intent(IN) :: lb1
    integer, intent(IN) :: ub1
    integer, intent(IN) :: lb2
    integer, intent(IN) :: ub2
    integer, intent(IN) :: lb3
    integer, intent(IN) :: ub3
    character(*), intent(IN), optional :: string

    ! Local Variables
    character(64) :: routine = 'ArrayCreate_real_Type_R3'
    character(string_len) :: error_string
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to allocate space
    ALLOCATE (ArrayPtr(lb1:ub1,lb2:ub2,lb3:ub3), STAT = status)

    if (PRESENT(string)) then
       error_string = string
    else
       error_string = routine
    end if

    if (status /= 0) call TLS_panic (trim(error_string) // ': allocate failed')

  END SUBROUTINE ARRAYCREATE_REAL_TYPE_R3

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE ARRAYDESTROY_REAL_TYPE_R3 (ArrayPtr, string)
    !=======================================================================
    ! Purpose(s):
    !   DEALLOCATE real_type array of rank 3.
    !=======================================================================

    ! Argument List
    real(r8), pointer, dimension(:,:,:) :: ArrayPtr
    character(*), intent(IN), optional :: string

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! attempt to deallocate space
    DEALLOCATE (ArrayPtr)

  END SUBROUTINE ARRAYDESTROY_REAL_TYPE_R3

END MODULE ARRAYALLOCATE_MODULE
