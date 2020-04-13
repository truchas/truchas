Module pgslib_process_info
  !-----------------------------------------------------------------------------
  ! purpose
  !
  ! gather with process information
  !
  ! Uses calls into the C file get_process_info.C
  !-----------------------------------------------------------------------------

  ! $Id: pgslib_process_info.F,v 1.1.1.1 2000/10/11 22:44:31 ferrell Exp $

  Implicit None
  Private

  Public PGSLib_Memory_Size

  interface 
    subroutine pgslib_get_process_id_c (pid) bind(c)
      use,intrinsic :: iso_c_binding, only: c_int
      integer(c_int) :: pid
    end subroutine
    subroutine pgslib_get_vm_size_c (vmsize) bind(c)
      use,intrinsic :: iso_c_binding, only: c_int
      integer(c_int) :: vmsize
    end subroutine
  end interface

Contains

  !-----------------------------------------------------------------------------

  function PGSLib_Memory_Size()
    !---------------------------------------------------------------------------
    ! purpose
    !
    ! return virtual memory size
    !
    !---------------------------------------------------------------------------

    implicit none

    integer :: PGSLib_Memory_Size
    ! local variables
    Integer :: vmsize

    !---------------------------------------------------------------------------

    Call pgslib_get_vm_size_C (vmsize)
    PGSLib_Memory_Size = vmsize


  End function PGSLib_Memory_Size

  !-----------------------------------------------------------------------------

End Module pgslib_process_info
