program vizre

  use re_encl_type
  use re_dist_vf_type
  use re_graphics_gmv
  use vizre_command_line
  use string_utilities, only: i_to_c
  use scl
  implicit none

  integer :: n, nvar
  logical :: sym, has_vf, is_IOP
  integer, pointer :: row(:) => null(), col(:) => null()
  character(len=512) :: enclosure_file, gmv_file
  type(encl) :: e
  type(dist_vf) :: dvf
  real, pointer :: var(:) => null()

  call scl_init ()
  is_IOP = (scl_rank()==1)

  call parse_command_line (enclosure_file, gmv_file, row, col, sym)

  call read_encl (e, trim(enclosure_file), has_vf)
  if (is_IOP) then
    call gmv_open (trim(gmv_file))
    call gmv_write_encl (e, sym)
  end if

  if (has_vf) then ! we have a view factor matrix

    call read_dist_vf (dvf, enclosure_file)

    if (is_IOP) call gmv_begin_variables ()
    var => get_ambient_vf(dvf)
    if (is_IOP) call gmv_write_face_var (e, var, 'ambient', sym)
    deallocate(var)
    var => dvf_row_sum(dvf)
    if (is_IOP) call gmv_write_face_var (e, var, 'row sum', sym)
    deallocate(var)
    nvar = 2
    if (associated(col)) then ! we want to write certain vf matrix columns as face variables
      do n = 1, size(col)
        if (col(n) >= 1 .and. col(n) <= e%nface) then
          if (exceeded_variable_limit(nvar)) goto 99
          var => unpack_dvf_col(dvf, col(n))
          if (is_IOP) call gmv_write_face_var (e, var, 'col'//i_to_c(col(n)), sym)
          deallocate(var)
        else
          write(*,fmt='(a)') 'No such view factor matrix col: ' // i_to_c(col(n))
        end if
      end do
      deallocate(col)
    end if
    if (associated(row)) then ! we want to write certain vf matrix rows as face variables
      do n = 1, size(row)
        if (row(n) >= 1 .and. row(n) <= e%nface) then
          if (exceeded_variable_limit(nvar)) goto 99
          var => unpack_dvf_row(dvf, row(n))
          if (is_IOP) call gmv_write_face_var (e, var, 'row'//i_to_c(row(n)), sym)
          deallocate(var)
        else
          write(*,fmt='(a)') 'No such view factor matrix row: ' // i_to_c(row(n))
        end if
      end do
      deallocate(row)
    end if
99  if (is_IOP) call gmv_end_variables ()
  else if (associated(row) .or. associated(col)) then  ! we can't write the vf matrix
    write(*,'(a)') 'Enclosure contains no view factor data'
  end if

  if (is_IOP) call gmv_close ()

  call destroy (e)
  call destroy (dvf)

  call scl_finalize()

contains

  logical function exceeded_variable_limit (nvar)
    integer, intent(inout) :: nvar
    integer, parameter :: MAX_GMV_VARS = 200  ! GMV says 250 but it segfaults before then...
    nvar = nvar + 1
    exceeded_variable_limit = (nvar > MAX_GMV_VARS)
    if (exceeded_variable_limit) then
      if (is_IOP) write(*,'(a)') 'Reached the limit of GMV variables (' // &
          i_to_c(MAX_GMV_VARS) // ').  Skipping all further variable output.'
    end if
  end function

end program vizre
