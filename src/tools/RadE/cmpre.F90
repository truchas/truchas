!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cmpre

  use re_utilities
  use re_dist_vf_type
  use re_vf_diff_type
  use scl
  use string_utilities, only: i_to_c
  implicit none

  integer :: n, rloc, cloc
  logical :: is_IOP
  real    :: max_norm, one_norm, two_norm, max_diff
  type(vf_diff) :: diff
  character(:), allocatable  :: prog, dataset1, dataset2
  character(255) :: arg

  call scl_init()
  is_IOP = (1 == scl_rank())  ! use process rank 1 for all IO

  !! Get the name by which the program was invoked (why?)
  call get_command_argument (0, arg)
  prog = trim(arg)
  n = scan(prog, '/', back=.true.)
  prog = prog(n+1:)  ! Remove the leading path component, if any.

  !! Verify we've got the right number of arguments.
  n = command_argument_count()
  if (n /= 2) then
    if (is_IOP) write(0,'(a)') 'Usage: ' // prog // ' dataset1 dataset2'
    call scl_finalize ()
    stop
  end if

  call get_command_argument (1, arg)
  dataset1 = trim(arg)
  call get_command_argument (2, arg)
  dataset2 = trim(arg)

  !! Read the enclosure dataset from the input files.
  call diff%init(dataset1, dataset2)

  call diff%compute_norms(max_norm, one_norm, two_norm, max_diff, rloc, cloc)
  if (is_IOP) write(*,'(3(a,es10.3))') 'max norm =', max_norm, '; one norm =', &
    one_norm, '; two norm =', two_norm
  if (is_IOP) write(*,'(a,es10.3,a)') 'max difference =', max_diff, ' at row ' &
    // i_to_c(rloc) // ', col ' // i_to_c(cloc)

  call scl_finalize()

end program cmpre
