!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program genre

  use re_utilities
  use re_dist_vf_type
  use re_encl_type
  use re_exodus_encl
  use re_chaparral_vf
  use re_patch_type
  use genre_command_line
  use scl
  implicit none

  integer :: ios
  logical :: is_IOP, found
  type(encl) :: e
  type(encl_spec) :: spec
  type(patch_param) :: ppar
  type(chap_param) :: cpar
  type(re_patch) :: ep
  type(dist_vf) :: vf
  character(len=512) :: infile, outfile
  character(len=32) :: string

  call scl_init ()
  is_IOP = (scl_rank()==1)

  if (is_IOP) then
    call parse_command_line (infile, outfile)
    open(unit=10,file=trim(infile),status='old',action='read',iostat=ios)
  end if
  call scl_bcast (ios)
  if (ios /= 0) call re_halt ('Unable to open input file: ' // trim(infile))

  !! Construct enclosure
  call read_enclosure_namelist (10, spec)
  call generate_encl (spec, e)
  call destroy (spec)

  string = 'Enclosure'
  call write_encl (e, trim(outfile))

  !! Generate patches
  call read_patches_namelist (10, ppar, found)
  if (found) then
    string = trim(string) // ' and patch'
  else
    call re_info ('No PATCHES namelist found')
  end if

  if (is_IOP) then
    call ep%generate_patches (e, ppar)
    call ep%write_patch_data (trim(outfile))
  end if

  !! Compute view factors
  call read_chaparral_namelist (10, cpar, found)
  if (found) then
    call calculate_vf (e, cpar, ep, vf)
    call write_dist_vf (vf, trim(outfile))
    string = trim(string) // ' and VF data'
  else
    call re_info ('No CHAPARRAL namelist found')
    string = trim(string) // ' data only'
  end if

  call re_info (trim(string) // ' written to ' // trim(outfile))

  call scl_finalize()

end program genre
