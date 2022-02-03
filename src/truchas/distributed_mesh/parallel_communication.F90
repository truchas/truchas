!!
!! PARALLEL_COMMUNICATION
!!
!! This module provides data that describe the configuration of the parallel
!! machine, and some basic parallel communication procedures.  For the most
!! part these are just renamed PGSLib procedures, but some have been extended
!! to handle multidimensional arrays.  Application code is expected to use
!! the facilities provided by this module, and not use PGSLib directly.  This
!! will facilitate replacing PGSLib with another parallel communication library.
!!
!! Neil N. Carlson <nnc@newmexico.com>
!! Last revised 24 Mar 2004
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module parallel_communication

  use,intrinsic :: iso_fortran_env, only: int8
  use pgslib_module, halt_parallel_communication => PGSLib_Finalize, &
                     broadcast  => PGSLib_BCast, &
                     distribute => PGSLIB_Dist,  &
                     collate    => PGSLib_Collate, &
                     global_any => PGSLib_Global_Any, &
                     global_all => PGSLib_Global_All, &
                     global_count => PGSLib_Global_Count, &
                     global_sum => PGSLib_Global_Sum, &
                     global_minval => PGSLib_Global_Minval, &
                     global_maxval => PGSLib_Global_Maxval, &
                     global_dot_product => PGSLib_Global_Dot_Product
  
  implicit none
  private
  
  public :: init_parallel_communication, halt_parallel_communication
  public :: broadcast, distribute, collate
  public :: global_any, global_all, global_count, global_sum, global_minval, global_maxval, global_dot_product
  public :: allocate_collated_array
  public :: global_maxloc, get_value
  
  !! Public 
  integer, public, save :: nPE = 1          ! total number of processing entities (PE)
  integer, public, save :: this_PE = 1      ! the number of this PE
  integer, public, save :: IO_PE = 1        ! the number of the PE allowed to do I/O
  logical, public, save :: is_IOP = .true.  ! True if this PE is the IO PE.
  integer, public, save :: delta_IOP = 1    ! 1 if this PE is the IO PE, 0 otherwise
  
  interface distribute
    module procedure distribute_L2, distribute_I2, distribute_S2, distribute_D2, distribute_int8_2
  end interface
  
  interface collate
    module procedure collate_L2, collate_I2, collate_S2, collate_D2, collate_int8_2
  end interface
  
  interface allocate_collated_array
    module procedure allocate_CA_L1, allocate_CA_L2
    module procedure allocate_CA_I1, allocate_CA_I2
    module procedure allocate_CA_S1, allocate_CA_S2
    module procedure allocate_CA_D1, allocate_CA_D2
  end interface
  
  interface get_value
    module procedure get_value_i, get_value_r
  end interface
  
contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! 
 !! INIT_PARALLEL_COMMUNICATION
 !!
 !! Initialize the module data that describe the configuration of the parallel
 !! machine.  For now we just piggy-back off of Truchas' initialization of
 !! PGSLib.
 !!
 
  subroutine init_parallel_communication (argv)
  
    use pgslib_module

    character(PGSLib_CL_MAX_TOKEN_LENGTH), dimension(:), pointer :: argv
    
    call PGSLib_INITIALIZE (0, argv=argv)
    nPE     = PGSLib_Inquire_NPE()
    this_PE = PGSLib_Inquire_ThisPE_Actual()
    IO_PE   = PGSLib_Inquire_IO_Root_PE()
    is_IOP  = PGSLIB_Inquire_IO_P()
    delta_IOP = 0
    if (is_IOP) delta_IOP = 1

  end subroutine init_parallel_communication
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Specific procedures for DISTRIBUTE: extend PGSLib_Dist to rank-2 arrays
 !!
 !!

#define _LOGICAL_DATA_
#include "distribute.fpp"
 
#define _INT8_DATA_
#include "distribute.fpp"
 
#define _INTEGER_DATA_
#include "distribute.fpp"
 
#define _SINGLE_DATA_
#include "distribute.fpp"
 
#define _DOUBLE_DATA_
#include "distribute.fpp"
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Specific procedures for COLLATE: extend PGSLib_Collate to rank-2 arrays
 !!

#define _LOGICAL_DATA_
#include "collate.fpp"
 
#define _INT8_DATA_
#include "collate.fpp"
 
#define _INTEGER_DATA_
#include "collate.fpp"
 
#define _SINGLE_DATA_
#include "collate.fpp"
 
#define _DOUBLE_DATA_
#include "collate.fpp"

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Specific procedures for ALLOCATE_COLLATED_ARRAY
 !!
 !! These routines allocate an array of the specified size on the IO processor,
 !! but on all other processors simply allocate a zero-sized array.  Because
 !! it is illegal to pass an unassociated pointer to a procedure when the
 !! corresponding actual argument is not itself a pointer, these routines are
 !! useful for allocating arrays to be passed to a global procedure which will
 !! use the array on the IO processor only; on all other processors the array
 !! is simply a place holder.
 !!
 !! NOTES
 !!
 !! (1) These routines have the same behavior as the intrinsic ALLOCATE; any
 !! storage the pointer might be associated with is _not_ deallocated before
 !! allocating the new storage.
 !!
 !! (2) These routine return an optional, _global_ status value.  We assume
 !! there is no way a zero-sized allocate can fail, so we simply broadcast
 !! the status value on the IO processor to the others instead of the more
 !! complex 'global_any'.
 !!
 !! (3) The unusual semantics of the STAT specifier for the ALLOCATE statement
 !! do not permit passing the procedure's optional STAT argument directly to
 !! the ALLOCATE.  Thus the dual allocate forms delineated by the if-else.
 !!

#define _LOGICAL_DATA_
#include "allocate_collated_array.fpp"
 
#define _INTEGER_DATA_
#include "allocate_collated_array.fpp"
 
#define _SINGLE_DATA_
#include "allocate_collated_array.fpp"
 
#define _DOUBLE_DATA_
#include "allocate_collated_array.fpp"

  subroutine global_maxloc (array, pid, lindex, mask)
  
    double precision, intent(in) :: array(:)
    integer, intent(out) :: pid, lindex
    logical, intent(in), optional :: mask(:)
    
    integer :: maxlocs(nPE)
    double precision :: maxvals(nPE)
    
    if (present(mask)) then
      ASSERT(size(mask) == size(array))
    end if
    
    call collate (maxvals, maxval(array, mask))
    call collate (maxlocs, maxloc(array, 1, mask))
    
    if (is_IOP) then
      pid = maxloc(maxvals,dim=1)
      lindex = maxlocs(pid)
      !TODO! There is a border case that this routine doesn't handle correctly.
      !TODO! If the first PE has a 0-sized array and the values of all other
      !TODO! array elements are equal to -huge(1.0d0) (the result maxval returns
      !TODO! for a 0-sized array), this routine returns lindex=0, pid=1 when it
      !TODO! really ought to return lindex=1 and the first pid that has an array
      !TODO! with a positive size.
    end if
    
    call broadcast (pid)
    call broadcast (lindex)
    
  end subroutine global_maxloc
  
  subroutine get_value_i (array, pid, lindex, value)
  
    integer, intent(in) :: array(:)
    integer, intent(in) :: pid, lindex
    integer, intent(out) :: value
    
    integer :: values(nPE)
    
    if (this_PE == pid) then
      call collate (values, array(lindex))
    else
      call collate (values, 0)  ! send a dummy value
    end if
    
    if (is_IOP) value = values(pid)
    call broadcast (value)
    
  end subroutine get_value_i
  
  subroutine get_value_r (array, pid, lindex, value)
  
    double precision, intent(in) :: array(:)
    integer, intent(in) :: pid, lindex
    double precision, intent(out) :: value
    
    double precision :: values(nPE)
    
    ASSERT(pid >= 1 .and. pid <= nPE)
    
    if (this_PE == pid) then
      ASSERT(lindex >= 1 .and. lindex <= size(array))
      call collate (values, array(lindex))
    else
      call collate (values, 0.0d0)  ! send a dummy value
    end if
    
    if (is_IOP) value = values(pid)
    call broadcast (value)
    
  end subroutine get_value_r

end module parallel_communication
