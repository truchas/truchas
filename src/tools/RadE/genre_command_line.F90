!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module genre_command_line

  implicit none
  private

  public  :: parse_command_line

  character(:), allocatable :: prog

contains

  subroutine parse_command_line(infile, outfile, overwrite)

    character(:), allocatable, intent(out) :: infile, outfile
    logical, intent(out) :: overwrite

    integer :: i, n, num_arg
    character(255) :: arg

    call get_command_argument(0, arg)
    i = scan(arg, '/', back=.true.)
    prog = trim(arg(i+1:))  ! Remove the leading path component, if any.

   !!!
   !!! PARSE COMMAND-LINE OPTIONS

    n = 0
    num_arg = command_argument_count()
    overwrite = .false.

    do while (n < num_arg)

      n = n + 1
      call get_command_argument(n, arg)

      select case (arg)

      case ('-f')

        overwrite = .true.

      case ('-h','--help')

        call usage_summary

      !!
      !! Example option with following numeric value
      !!
      !! case ('-n')
      !!
      !!   n = n + 1
      !!   if (n > num_arg) call usage_halt('option requires an argument: ' // trim(arg))
      !!   call get_command_argument (n, arg)
      !!   read(unit=arg,fmt=*,iostat=ios) n_value
      !!   if (ios /= 0) call usage_halt('invalid value for -n: ' // trim(arg))
      !!
      !! Continue to do sanity checking of the value.
      !!

      case default

        if (arg(1:1) == '-') call usage_halt('invalid option: ' // trim(arg))
        n = n - 1 ! Back-up.
        exit

      end select

    end do

    if (num_arg - n /= 2) call usage_halt('exactly two file arguments are required')

    !! Next two arguments are interpreted as file paths.
    call get_command_argument(n+1, arg)
    infile = trim(arg)
    call get_command_argument(n+2, arg)
    outfile = trim(arg)

  end subroutine parse_command_line

  subroutine usage_halt(message)
    use,intrinsic :: iso_fortran_env, only: output_unit
    character(*), intent(in) :: message
    write(output_unit,fmt='(a)') prog // ': ' // message
    write(output_unit,fmt='(a)') 'Try `' // prog // ' -h'' for more information.'
    error stop
  end subroutine usage_halt

  subroutine usage_summary
    use,intrinsic :: iso_fortran_env, only: output_unit
    write(output_unit,fmt='(a)') 'Usage: ' // prog // ' [options] INFILE OUTFILE'
    write(output_unit,fmt='(a)') ' '
    write(output_unit,fmt='(a)') 'Generates a radiation enclosure data set OUTFILE according'
    write(output_unit,fmt='(a)') 'to the parameters read from INFILE.'
    write(output_unit,fmt='(a)') ' '
    write(output_unit,fmt='(a)') 'Options:'
    write(output_unit,fmt='(a)') '  -h, --help   Display this help and exit.'
    stop
  end subroutine usage_summary

end module genre_command_line
