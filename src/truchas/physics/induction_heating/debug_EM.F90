!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module debug_EM

  use kinds, only: r8
  use parallel_communication
  use index_partitioning
  
  implicit none
  private
  
  public :: verify_vector, verify_parallel_data
  
  integer, save :: lun = -1, mode = 0
  
contains

  subroutine verify_vector (vec, gsd, label)
  
    real(kind=r8), intent(in) :: vec(:)
    type(ip_desc), intent(in) :: gsd
    character(len=*), intent(in) :: label
    
    integer :: n
    real(kind=r8) :: err1, err2, err3
    real(kind=r8), pointer :: ref(:), g_ref(:)
    
    if (lun == -1) then
      open(newunit=lun,file='data.dat',form='unformatted')
    end if
    
    if (nPE == 1) then
      write(unit=lun) size(vec)
      write(unit=lun) vec
    else
      if (is_IOP) read(unit=lun) n
      call allocate_collated_array (g_ref, n)
      if (is_IOP) read(unit=lun) g_ref
      allocate(ref(size(vec)))
      call distribute (ref, g_ref)
      call gather_boundary (gsd, ref)
      err1 = global_maxval(abs(ref))
      err2 = global_maxval(abs(vec-ref))
      err3 = global_maxval(abs(vec-ref)/abs(ref),(abs(ref)>0.0_r8))
      if (is_IOP) print *, label, err1, err2, err3
    end if
    
  end subroutine verify_vector
  
  subroutine verify_parallel_data (array, label)
  
    use string_utilities, only: i_to_c
  
    real(kind=r8), intent(in) :: array(:)
    character(len=*), intent(in) :: label
    
    integer :: n
    logical :: exists
    character(len=8) :: fname
    real(kind=r8) :: err1, err2, err3
    real(kind=r8), allocatable :: ref_array(:)
    
    if (lun == -1) then
      fname = 'data.' // i_to_c(this_PE)
      inquire(file=fname, exist=exists)
      exists = global_all(exists)
      if (exists) then
        mode = 2
        open(newunit=lun, file=fname, form='unformatted', action='read', position='rewind')
      else
        mode = 1
        open(newunit=lun, file=fname, form='unformatted', action='write', position='rewind')
      end if
    end if
    
    select case (mode)
    case(1) ! Writing data to disk
      write(unit=lun) size(array)
      write(unit=lun) array
    case(2) ! Comparing data against disk
      read(unit=lun) n
      allocate(ref_array(n))
      read(unit=lun) ref_array
      err1 = global_maxval(abs(ref_array))
      err2 = global_maxval(abs(array-ref_array))
      err3 = global_maxval(abs(array-ref_array)/abs(ref_array),(abs(ref_array)>0.0_r8))
      write(0,*) label//'(', this_PE, ')', err1, err2, err3
      deallocate(ref_array)
    end select
    
  end subroutine verify_parallel_data
      
    
end module debug_EM
