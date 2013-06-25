program input_utilities_test_1

#ifdef NAG
  use f90_unix, only: exit
#endif
  use input_utilities
  use string_utilities
  implicit none
  
  integer :: x, status = 0
  namelist /foo/ x
  
  integer :: lun = 10
  character(len=13) :: inputfile = 'testinput.tmp'
  character(len=22) :: prog = 'input_utilities_test_1'
  
  call empty_file_test()
  call file1_test()
  call file2_test()
    
  call exit (status)

contains

  subroutine empty_file_test ()
  
    logical :: pass
    integer :: values(0)

    !! empty file
    open(lun,file=inputfile,action='write',position='rewind',status='replace')
    close(lun)
    
    open(lun,file=inputfile,action='read',position='rewind')
    call read_namelist_foo (lun, values, pass)
    close(lun)
    
    if (.not.pass) then
      status = 1
      write(0,*) prog // ': empty file test failed'
    end if
    
  end subroutine empty_file_test

  subroutine file1_test ()
  
    logical :: pass
    integer :: values(3) = (/1,2,6/)

    !! empty file
    open(lun,file=inputfile,action='write',position='rewind',status='replace')
    write(lun,fmt='(a)') '&foo x=1 /'
    write(lun,fmt='(a)') ' '
    write(lun,fmt='(a)') ''
    write(lun,fmt='(a)') '   &foo'
    write(lun,fmt='(a)') 'x=2/'
    write(lun,fmt='(a)') ' & foo x=3/'
    write(lun,fmt='(a)') '&'
    write(lun,fmt='(a)') 'random crap'
    write(lun,fmt='(a)') 'foo &foo x=4/'
    write(lun,fmt='(a)') '&bar x=5'
    write(lun,fmt='(a)') '/'
    write(lun,fmt='(a)') ' Thbbt'
    write(lun,fmt='(a)') '&foo x=6/'  
    close(lun)
    
    open(lun,file=inputfile,action='read',position='rewind')
    call read_namelist_foo (lun, values, pass)
    close(lun)
    
    if (.not.pass) then
      status = 1
      write(0,*) prog // ': file1 test failed'
    end if
    
  end subroutine file1_test

  subroutine file2_test ()
  
    logical :: pass
    integer :: values(2) = (/3,4/)

    !! empty file
    open(lun,file=inputfile,action='write',position='rewind',status='replace')
    write(lun,fmt='(a)') '&foobar x=1 /'
    write(lun,fmt='(a)') '&&foo'
    write(lun,fmt='(a)') ' x=2'
    write(lun,fmt='(a)') ' /'
    write(lun,fmt='(a)') ''
    write(lun,fmt='(a)') '                                     &Foo x=3/'
    write(lun,fmt='(a)') 'random crap'
    write(lun,fmt='(a)') '&fOO x=4/'
    write(lun,fmt='(a)') 'x=4'
    write(lun,fmt='(a)') '/'
    write(lun,fmt='(a)') ' Thbbt'
    write(lun,fmt='(a)') '                    &fu x=6/'  
    close(lun)
    
    open(lun,file=inputfile,action='read',position='rewind')
    call read_namelist_foo (lun, values, pass)
    close(lun)
    
    if (.not.pass) then
      status = 1
      write(0,*) prog // ': file2 test failed'
    end if
    
  end subroutine file2_test

  subroutine read_namelist_foo (lun, values, pass)
  
    integer, intent(in)  :: lun
    integer, intent(in)  :: values(:)
    logical, intent(out) :: pass
    
    integer :: n
    logical :: found
    
    pass = .true.
    n = 0
    
    do  ! until no more namelists
    
      call seek_to_namelist (lun, 'foO', found) ! check case insensitivity
      if (.not.found) exit
      
      n = n + 1
      
      read(lun,foo)
      
      if (n <= size(values)) then
        if (x /= values(n)) then
          write(0,*) 'error: foo namelist ' // i_to_c(n) // ': expected ' // i_to_c(values(n)) // ', got ' // i_to_c(x)
          pass = .false.
          exit
        end if
      end if
      
    end do
    
    if (pass .and. n /= size(values)) then
      write(0,*) 'error: wrong number of foo namelists found: expected ' // i_to_c(size(values)) // ', got ' // i_to_c(n)
      pass = .false.
    end if

  end subroutine read_namelist_foo

end program input_utilities_test_1
