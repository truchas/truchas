!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module command_line

  implicit none
  private

  public  :: parse_command_line
  private :: usage_halt, usage_summary

  !! Public variables initialized from the command line.
  character(len=32), public :: prog

  logical, public :: strict_ss_transf = .false.

contains

  subroutine parse_command_line (infile, outfile, ssid)

    character(len=*), intent(out) :: infile, outfile
    integer, pointer :: ssid(:)

    integer :: i, n, num_arg, ios
    character(len=128) :: arg

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
      
      case ('-s')
      
        strict_ss_transf = .true.

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

      case default

        if (arg(1:1) == '-') call usage_halt ('invalid option: ' // trim(arg))
        n = n - 1 ! Back-up.
        exit

      end select

    end do

    if (num_arg - n < 3) call usage_halt ('not enough arguments')

    !! Next two arguments are interpreted as file paths.
    call get_command_argument (n+1, infile)
    call get_command_argument (n+2, outfile)
    n = n + 2

    !! The remaining arguments are interpreted as side set numbers.
    allocate(ssid(num_arg-n))
    do i = 1, size(ssid)
      call get_command_argument (n+i, arg)
      read(arg,fmt=*,iostat=ios) ssid(i)
      if (ios /= 0) call usage_halt ('invalid sideset number: ' // trim(arg))
    end do

  end subroutine parse_command_line

  subroutine usage_halt (message)

    character(len=*), intent(in) :: message

    write(unit=*, fmt='(a)') trim(prog) // ': ' // message
    write(unit=*, fmt='(a)') 'Try `' // trim(prog) // ' -h'' for more information.'
    stop

  end subroutine usage_halt

  subroutine usage_summary ()

    write(unit=*,fmt='(a)') 'Usage: ' // trim(prog) // ' [options] infile outfile side-set-ID [side-set-ID ...]'
    write(unit=*,fmt='(a)') ' '
    write(unit=*,fmt='(a)') 'Adds gap elements to the Exodus-format mesh read from INFILE for each specified'
    write(unit=*,fmt='(a)') 'sideset.  The resulting mesh is written to the Exodus-format mesh file OUTFILE.'
    write(unit=*,fmt='(a)') ' '
    write(unit=*,fmt='(a)') 'Options:'
    write(unit=*,fmt='(a)') '  -s           Do not modify the side sets on output.  Normally the side'
    write(unit=*,fmt='(a)') '               sets specified on the command line are augmented if necessary'
    write(unit=*,fmt='(a)') '               to make them two-sided.' 
    write(unit=*,fmt='(a)') '  -h, --help   Display this help and exit.'
    stop

  end subroutine usage_summary

end module command_line
