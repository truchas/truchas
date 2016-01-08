!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module truchas_danu_output_data
  use,intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR, c_associated
  implicit none
  private
  type(c_ptr), save, public :: fid = C_NULL_PTR ! h5 file id
  type(c_ptr), save, public :: sid = C_NULL_PTR ! Danu simulation id
  public :: close_danu_output_file
contains
  subroutine close_danu_output_file
    use danu_module, only: output_file_close
    if (c_associated(fid)) then
      call output_file_close (fid)
      fid = C_NULL_PTR
    end if
  end subroutine close_danu_output_file
end module truchas_danu_output_data
