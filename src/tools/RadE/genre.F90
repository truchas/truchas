program genre

  use re_utilities
  use re_dist_vf_type
  use re_encl_type
  use re_exodus_encl
  use re_chaparral_vf
  use genre_command_line
  use scl
  implicit none

  integer :: ios
  logical :: is_IOP, found
  type(encl) :: e
  type(encl_spec) :: spec
  type(chap_param) :: cpar
  type(dist_vf) :: vf
  character(len=512) :: infile, outfile

  call scl_init ()
  is_IOP = (scl_rank()==1)

  if (is_IOP) then
    call parse_command_line (infile, outfile)
    open(unit=10,file=trim(infile),status='old',action='read',iostat=ios)
  end if
  call scl_bcast (ios)
  if (ios /= 0) call re_halt ('Unable to open input file: ' // trim(infile))

  call read_enclosure_namelist (10, spec)
  call generate_encl (spec, e)
  call destroy (spec)

  call read_chaparral_namelist (10, cpar, found)
  if (found) then
    call calculate_vf (e, cpar, vf)
    if (is_IOP) call write_encl (e, outfile)
    call write_dist_vf (vf, outfile)
    call re_info ('Enclosure and VF data written to ' // trim(outfile))
    call destroy (vf)
    call destroy (e)
  else
    call write_encl (e, outfile)
    call re_info ('No CHAPARRAL namelist found')
    call re_info ('Enclosure data only written to ' // trim(outfile))
    call destroy (e)
  end if

  call scl_finalize()

end program genre
