MODULE Output_Module
  !=======================================================================
  ! Purpose(s):
  !   Punt2: Alternate interface to Punt, useful for simple strings.
  !      Prints out the name of the routine that it was called from
  !      and n user strings. Does not return.
  !
  !   Usage:
  !      Use Output_Module
  !      Call Punt2('RoutineName',string1,string2,...)
  !
  !   Variables:
  !      Character(LEN=*), intent(IN)           :: Routine
  !      Character(LEN=*), intent(IN), optional :: string_n
  !
  !   Notes:
  !      The null string '' can be used to insert a blank line in the output.
  !      Currently written for up to 10 strings
  !
  ! Contains: Punt2
  !
  ! Author(s): Bryan R. Lally, LANL ESA-EPE (Lally@lanl.gov)
  !
  !=======================================================================
  implicit none

  ! Private Module
  private

  ! Public Variables

  ! Public Subroutines
  public :: Punt2

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  
CONTAINS

  SUBROUTINE Punt2 (Routine, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10)
  !=======================================================================
  ! Purpose(s):
  !   Alternate interface to Punt, useful for simple strings.
  !   Prints out the name of the routine that it was called from
  !   and n user strings. Does not return.
  !=======================================================================
    implicit none

    ! Argument List
    Character(LEN=*), intent(IN)           :: Routine
    Character(LEN=*), intent(IN), optional :: s1
    Character(LEN=*), intent(IN), optional :: s2
    Character(LEN=*), intent(IN), optional :: s3
    Character(LEN=*), intent(IN), optional :: s4
    Character(LEN=*), intent(IN), optional :: s5
    Character(LEN=*), intent(IN), optional :: s6
    Character(LEN=*), intent(IN), optional :: s7
    Character(LEN=*), intent(IN), optional :: s8
    Character(LEN=*), intent(IN), optional :: s9
    Character(LEN=*), intent(IN), optional :: s10

    ! Local Variables

    ! <><><><><><<><><><<><><><><><><><><><><><><<><><><><><><><><><><><><>

    Write (*,*) 'Punt2: failed in ',Routine
    If (PRESENT(s1))  Write (*,*) s1
    If (PRESENT(s2))  Write (*,*) s2
    If (PRESENT(s3))  Write (*,*) s3
    If (PRESENT(s4))  Write (*,*) s4
    If (PRESENT(s5))  Write (*,*) s5
    If (PRESENT(s6))  Write (*,*) s6
    If (PRESENT(s7))  Write (*,*) s7
    If (PRESENT(s8))  Write (*,*) s8
    If (PRESENT(s9))  Write (*,*) s9
    If (PRESENT(s10)) Write (*,*) s10
    Stop
    
  END SUBROUTINE Punt2

END MODULE Output_Module
