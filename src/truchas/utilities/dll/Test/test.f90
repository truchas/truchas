program main

  use dynamic_linking
  
  integer(c_ptr_kind) :: so_handle
  
  call dl_open ('./foo.so', so_handle, RTLD_NOW)
  
end program main


