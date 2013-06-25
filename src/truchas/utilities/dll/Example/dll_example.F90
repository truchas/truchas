program dll_example

  use dynamic_linking_loader
  implicit none
  
  integer(c_ptr_kind) :: so, funptr
  
  interface
    real function call_f (funptr, x)
      use dynamic_linking_loader, only: c_ptr_kind
#ifdef HAVE_PASS_BY_VALUE
      integer(c_ptr_kind), value :: funptr ! this is a lie
#else
      integer(c_ptr_kind) :: funptr ! this is a lie
#endif
      real, intent(in) :: x
    end function
  end interface
  
  !! First we load the shared object library that contains the procedure we
  !! we wish to call.  SO is the returned handle to library; we'll need this
  !! for later calls.  There are different modes for how the library is
  !! loaded; I don't understand the differences and I don't think it really
  !! matters too much.  See the documentation for dlopen.
  
  call dll_open ('./myfun.so', RTLD_NOW, so)
  
  !! Next we need to get the address of the procedure we wish to call in the
  !! shared library we just loaded.  The procedure is an external Fortran
  !! function named 'myfun2', however the compiler typically mangles this name,
  !! so the name that actually appears in the shared library will usually be
  !! different.  For external procedures (i.e., those not in a module) just an
  !! underscore is usually appended to the name.  For module procedures the
  !! name usually includes the name of the module; the nm utility can be used
  !! on the shared library to see a list of symbols it contains.

  call dll_symbol (so, 'myfun2_', funptr)
  
  !! Finally we want to actually call the function, however all we have is
  !! the C pointer to the function and we have no direct way of invoking the
  !! function in Fortran 95.  (Fortran 2003 solves this problem.)  Here we
  !! employ some trickery to invoke the function, which is technically illegal
  !! but which generally works.  Fortran does have a very restricted form of
  !! procedure pointer that is used when a procedure is passed as an argument;
  !! the dummy procedure is really a procedure pointer aliased to the actual
  !! procedure.  To exploit this feature, we use an intermediate Fortran
  !! procedure, here named CALL_F, whose first argument is a dummy procedure
  !! followed by the arguments of the function to be invoked.  However when
  !! this intermediate procedure is called, we pass it the C function pointer
  !! instead of a procedure name.  This requires us to lie about the interface
  !! to CALL_F, and to pass the pointer by value.  Many compilers permit
  !! pass-by-value for selected arguments as an extension, and here we have
  !! used the VALUE attribute from F2003.  If it is not possible to pass the
  !! pointer by value, then the intermediate procedure must be implemented
  !! in C, where the pointer-to-a-pointer that Fortran actually passes can
  !! be unpacked (see call_f.c).  However, with a C intermediate procedure
  !! one will be limited to function arguments that are interoperable with
  !! C; for example, it would probably not be possible to pass a derived
  !! type or F95 array (a F77 array is possible of course).
  
  print *, call_f(funptr, 1.0)
  
  !! When we are finished using the procedure we can unload the shared library.

  call dll_close (so)

end program dll_example
