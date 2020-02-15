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
  use parameter_list_type
  use scl
  implicit none

  integer :: lun, ios, stat
  logical :: is_IOP, found
  type(encl) :: e
  type(re_patch) :: ep
  type(parameter_list) :: encl_params, chap_params, patch_params
  type(dist_vf) :: vf
  character(:), allocatable :: infile, outfile, msg, errmsg
  character(255) :: iom

  call scl_init
  is_IOP = (scl_rank()==1)

  if (is_IOP) call parse_command_line(infile, outfile)
  call scl_bcast_alloc(infile)
  call scl_bcast_alloc(outfile)

  if (is_IOP) open(newunit=lun,file=infile,status='old',action='read',iostat=ios,iomsg=iom)
  call scl_bcast(ios)
  if (ios /= 0) call re_halt('unable to open input file: ' // infile // ': ' // trim(iom))

  !! Construct enclosure
  call read_enclosure_namelist(lun, encl_params)
  call init_encl(e, encl_params, stat, errmsg)
  if (stat /= 0) call re_halt('error creating the enclosure: ' // errmsg)

  msg = 'Enclosure'
  call e%write(outfile)

  !! Generate patches
  call read_patches_namelist(lun, patch_params, found)
  if (found) then
    msg = msg // ' and patch'
  else
    call re_info('No PATCHES namelist found')
  end if

  if (is_IOP) then
    call ep%generate_patches(e, patch_params)
    call ep%write_patch_data(outfile)
  end if

  !! Compute view factors
  call read_chaparral_namelist(lun, chap_params, found)
  if (found) then
    call calculate_vf(e, chap_params, ep, vf)
    call write_dist_vf(vf, trim(outfile))
    msg = msg // ' and VF data'
  else
    call re_info('No CHAPARRAL namelist found')
    msg = msg // ' data only'
  end if

  call re_info(msg // ' written to ' // outfile)

  call scl_finalize

end program genre
