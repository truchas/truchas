module command_line

  implicit none
  private

  public  :: parse_command_line
  private :: usage_halt, usage_summary

  !! Public variables initialized from the command line.
  character(len=32), public :: prog

contains

  subroutine parse_command_line (infile, outfile, seed)

    character(len=*), intent(out) :: infile, outfile
    integer, intent(out) :: seed

    integer :: i, n, num_arg, ios
    character(len=128) :: arg
    
    !! Default seed value
    seed = -1193317467

    call get_command_argument (0, arg)
    i = scan(arg, '/', back=.true.)
    prog = arg(i+1:)  ! Remove the leading path component, if any.

   !!!
   !!! PARSE COMMAND-LINE OPTIONS

    n = 0
    num_arg = command_argument_count()
    do while (n < num_arg)

      n = n + 1
      call get_command_argument (n, arg)

      select case (arg)

      case ('-h','--help')

        call usage_summary ()

      !!
      !! Example option with following numeric value
      !!
      !! case ('-n')
      !! 
      !!   n = n + 1
      !!   if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
      !!   call get_command_argument (n, arg)
      !!   read(unit=arg,fmt=*,iostat=ios) n_value
      !!   if (ios /= 0) call usage_halt ('invalid value for -n: ' // trim(arg))
      !!
      !! Continue to do sanity checking of the value.
      !!

      case ('-s')
      
        n = n + 1
        if (n > num_arg) call usage_halt ('option requires an argument: ' // trim(arg))
        call get_command_argument (n, arg)
        read(unit=arg,fmt=*,iostat=ios) seed
        if (ios /= 0) call usage_halt ('invalid value for -s: ' // trim(arg))
        
      case default

        if (arg(1:1) == '-') call usage_halt ('invalid option: ' // trim(arg))
        n = n - 1 ! Back-up.
        exit

      end select

    end do

   !!!
   !!! INTERPRET REMAINING ARGUMENTS AS FILE NAMES

    if (num_arg - n /= 2) call usage_halt ('exactly two file arguments are required')
    call get_command_argument (n+1, infile)
    call get_command_argument (n+2, outfile)

  end subroutine parse_command_line

  subroutine usage_halt (message)

    character(len=*), intent(in) :: message

    write(unit=*, fmt='(a)') trim(prog) // ': ' // message
    write(unit=*, fmt='(a)') 'Try `' // trim(prog) // ' --help'' for more information.'
    stop

  end subroutine usage_halt

  subroutine usage_summary ()

    write(unit=*,fmt='(a)') 'Usage: ' // trim(prog) // ' [options] INMESH OUTMESH'
    write(unit=*,fmt='(a)') ''
    write(unit=*,fmt='(a)') 'This program generates a modified mesh where each element has been'
    write(unit=*,fmt='(a)') 'replaced by an equivalent randomly-selected relabeling of the element'
    write(unit=*,fmt='(a)') 'nodes.  Only 8-node hex, 6-node wedge, 5-node pyramid, and 4-node tet'
    write(unit=*,fmt='(a)') 'elements are relabeled currently.'
    write(unit=*,fmt='(a)') ''
    write(unit=*,fmt='(a)') 'The output mesh is geometrically equivalent to the input mesh and'
    write(unit=*,fmt='(a)') 'thus should be numerically equivalent as well.'
    write(unit=*,fmt='(a)') ''
    write(unit=*,fmt='(a)') 'Options:'
    write(unit=*,fmt='(a)') '  -h, --help      Display this help and exit.'
    write(unit=*,fmt='(a)') '  -s, --seed      Integer seed for the random number generator.'
    stop

  end subroutine usage_summary

end module command_line
