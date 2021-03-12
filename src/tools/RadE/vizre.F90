!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

program vizre

  use re_encl_type
  use re_patch_type
  use re_dist_vf_type
  use re_graphics_gmv
  use re_utilities
  use vizre_command_line
  use string_utilities, only: i_to_c
  use scl
  implicit none

  integer, parameter :: MAX_GMV_VARS = 200  ! GMV says 250 but it segfaults before then...

  integer :: n, nvar, stat
  logical :: sym, has_vf, is_IOP
  integer, allocatable :: col(:), pcolor(:)
  character(len=512) :: enclosure_file, gmv_file
  character(:), allocatable :: errmsg
  type(encl) :: e
  type(dist_vf) :: dvf
  type(re_patch) :: ep
  real, allocatable :: patch_var(:)
  real, allocatable :: face_var(:)

  call scl_init()
  is_IOP = (scl_rank()==1)
  nvar = 0

  call parse_command_line(enclosure_file, gmv_file, col, sym)

  call e%read(trim(enclosure_file), has_vf)
  call ep%read(trim(enclosure_file))

  if (is_IOP) then
    call gmv_open(trim(gmv_file))
    call gmv_write_encl(e, sym)

    !! Initialize GMV if other data will be written
    if (has_vf .or. ep%has_patches) then
      allocate(face_var(e%nface))
      call gmv_begin_variables()
    end if
  end if

  !! We have a face-to-patch map
  if (is_IOP .and. ep%has_patches) then
    !! Write patch IDs
    nvar = nvar + 1
    call gmv_write_face_var(e, REAL(ep%f2p_map), 'pid', sym)

    !! Write patch coloring
    nvar = nvar + 1
    call ep%patch_coloring(e, pcolor, stat, errmsg)
    if (stat==0) then
      patch_var = REAL(pcolor)
      call ep%patch_to_face_array(patch_var, face_var)
      call gmv_write_face_var(e, face_var, 'pcolor', sym)
      deallocate(patch_var, pcolor)
    end if
  end if

  call scl_bcast(stat)
  if (stat/=0) call re_halt(errmsg)

  !! We have a view factor matrix
  if (has_vf) then

    call dvf%read(trim(enclosure_file))

    !! Write ambient view factors
    if (dvf%has_ambient) then
      nvar = nvar + 1
      patch_var = dvf%get_ambient_vf()
      if (is_IOP) then
          call ep%patch_to_face_array(patch_var, face_var)
          call gmv_write_face_var(e, face_var, 'ambient', sym)
      end if
      deallocate(patch_var)
    end if

    !! Write row sums
    nvar = nvar + 1
    patch_var = dvf%row_sum()
    if (is_IOP) then
        call ep%patch_to_face_array(patch_var, face_var)
        call gmv_write_face_var(e, face_var, 'row sum', sym)
    end if
    deallocate(patch_var)

    !! We want to write certain vf matrix columns as face variables
    if (allocated(col)) then
      do n = 1, size(col)
        if (col(n) >= 1 .and. col(n) <= e%nface) then
          nvar = nvar + 1
          if (nvar > MAX_GMV_VARS) exit
          patch_var = dvf%unpack_col(col(n))
          if (is_IOP) then
              call ep%patch_to_face_array(patch_var, face_var)
              !! If patch-based, scale column by face weight
              if (ep%has_patches) face_var = face_var * dvf%w(col(n))
              call gmv_write_face_var(e, face_var, 'col'//i_to_c(col(n)), sym)
          end if
          deallocate(patch_var)
        else
          write(*,fmt='(a)') 'No such view factor matrix col: ' // i_to_c(col(n))
        end if
      end do
      deallocate(col)
    end if

    if (is_IOP .and. nvar > MAX_GMV_VARS) then
      write(*,'(a)') 'Reached the limit of GMV variables (' // &
        i_to_c(MAX_GMV_VARS) // ').  Skipping all further variable output.'
    end if
  else if (allocated(col)) then  ! we can't write the vf matrix
    write(*,'(a)') 'Enclosure contains no view factor data'
  end if

  !! Finalize GMV
  if (is_IOP) then
    if (has_vf .or. ep%has_patches) then
      deallocate(face_var)
      call gmv_end_variables()
    end if
    call gmv_close()
  end if

  call scl_finalize()

end program vizre
