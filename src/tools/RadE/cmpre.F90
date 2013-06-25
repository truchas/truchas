program cmpre

  use re_utilities
  use re_dist_vf_type
  use scl
  use string_utilities, only: i_to_c
  implicit none
  
  integer :: n, rloc, cloc
  logical :: is_IOP
  real    :: error0, error1, error2
  type(dist_vf) :: vf1, vf2
  character(len=32)  :: prog
  character(len=512) :: dataset1, dataset2
  
  call scl_init()
  is_IOP = (1 == scl_rank())  ! use process rank 1 for all IO

  !! Get the name by which the program was invoked (why?)
  call get_command_argument (0, prog)
  n = scan(prog, '/', back=.true.)
  prog = prog(n+1:)  ! Remove the leading path component, if any.

  !! Verify we've got the right number of arguments.
  n = command_argument_count()
  if (n /= 2) then
    if (is_IOP) write(0,'(a)') 'Usage: ' // trim(prog) // ' dataset1 dataset2'
    call scl_finalize ()
    stop
  end if
  
  call get_command_argument (1, dataset1)
  call get_command_argument (2, dataset2)
  
  !! Read the enclosure dataset from the infile.
  call read_dist_vf (vf1, trim(dataset1))
  call read_dist_vf (vf2, trim(dataset2))
  call vf_diff (vf1, vf2)
  call vf_diff_max_norm (vf2, error0)
  call vf_diff_one_norm (vf2, error1)
  call vf_diff_two_norm (vf2, error2)
  if (is_IOP) write(*,'(3(a,es10.3))') 'max norm =', error0, '; one norm =', error1, &
                                                             '; two norm =', error2
  call vf_diff_max (vf2, error0, rloc, cloc)
  if (is_IOP) write(*,'(a,es10.3,a)') 'max difference =', error0, &
      ' at row ' // i_to_c(rloc) // ', col ' // i_to_c(cloc)
  call destroy (vf1)
  call destroy (vf2)
  
  call scl_finalize()

end program cmpre
