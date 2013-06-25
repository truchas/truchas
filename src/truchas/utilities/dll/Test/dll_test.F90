program dll_error_test

#ifdef NAG
  use f90_unix, only: exit
#endif

  use dynamic_linking_loader
  implicit none
  
  integer :: n = 4, stat
  character(len=200) :: errmsg
  integer(c_ptr_kind) :: so1
  
  select case (n)
  case (1)
    call dll_open ('./foo.so', so1, RTLD_NOW)
  case (2)
    call dll_open ('./foo.so', so1, RTLD_NOW, stat)
    if (stat /= 0) call exit(1)
  case (3)
    call dll_open ('./foo.so', so1, RTLD_NOW, stat, errmsg)
    if (stat /= 0) then
      write(0,*) trim(errmsg)
      call exit (1)
    end if
  case (4)
    call dll_open ('libm.so', so1, 4)
    call dll_close (so1)
    call dll_close (so1)
  end select
  
end program dll_error_test
