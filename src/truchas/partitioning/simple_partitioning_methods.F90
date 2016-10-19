!!
!! SIMPLE_PARTITIONING_METHODS
!!
!! This module provides several procedures for generating a partition of an
!! index set that employ simple methods.  These are alternatives to the more
!! sophisticated graph partitioning methods from GRAPH_PARTITIONER_FACTORY.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! October 2016
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module simple_partitioning_methods

  implicit none
  private

  public :: read_partition, get_block_partition

contains

  !! This subroutine reads the partitioning from a disk file. The file format
  !! is simple, consisting of a sequence of integer values (one or multiple
  !! per line -- it does not matter); the first value is the partition number
  !! of the first index, the second value the partition number of the second
  !! index, and so forth. The partition numbers are returned in the array PART.
  !! The number of values must match the size of the array.  1-based numbering
  !! of partitions is used internally, however the numbering used in the file
  !! may differ (e.g., 0-based); PFIRST specifies the number of the first
  !! partition (0 or 1).  The subroutine performs the necessary translation.
  !! The number of partitions is NPART and the read values must be consistent
  !! with that value and PFIRST.  The subroutine does not require that every
  !! partition number is assigned to at least one index, only that the numbers
  !! belong to the expected range defined by PFIRST and NPART.  STAT returns
  !! a nonzero value if an error occurs, and the allocatable character ERRMSG
  !! an explanatory error message.

  subroutine read_partition(file, pfirst, npart, part, stat, errmsg)

    use string_utilities, only: i_to_c

    character(*), intent(in) :: file
    integer, intent(in)  :: pfirst, npart
    integer, intent(out) :: part(:), stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: lun, pmin, pmax, plast

    open(newunit=lun,file=file,status='old',position='rewind',action='read',iostat=stat)
    if (stat /= 0) then
      errmsg = 'unable to open file ' // file // ': iostat=' // i_to_c(stat)
      stat = 1
      return
    end if

    read(lun,*,iostat=stat) part
    close(lun)
    if (stat /= 0) then
      errmsg = 'error reading file ' // file // ' (iostat=' // i_to_c(stat) &
          // '); attempted to read ' // i_to_c(size(part)) // ' integers'
      stat = 1
      return
    end if

    pmin = minval(part)
    pmax = maxval(part)

    plast = pfirst + npart - 1

    if (pmin < pfirst .or. pmax > plast) then
      errmsg = 'invalid partition data range [' // i_to_c(pmin) // ',' // i_to_c(pmax) &
          // ']; expected range [' // i_to_c(pfirst) // ',' // i_to_c(plast) // ']'
      stat = 1
      return
    end if

    if (pfirst /= 1) part = part - pfirst + 1 ! map to [1,npart]
    stat = 0

  end subroutine read_partition

  !! This subroutine computes a simple block partition.  The indices, which
  !! correspond to the elements of the PART array, are partitioned into NPART
  !! equal-sized blocks.  The first index block is assigned to partition 1,
  !! the next to partition 2, and so forth.  The partitioning is returned in
  !! the PART array.  The partition sizes differ by at most 1.

  subroutine get_block_partition(npart, part)

    integer, intent(in)  :: npart   ! number of partitions
    integer, intent(out) :: part(:) ! partition vector

    integer :: j, j1, n, n1, offset

    ASSERT(npart >= 1)

    offset = 0
    n1 = size(part)/npart
    j1 = modulo(size(part), npart)
    do j = 1, npart
      n = n1
      if (j <= j1) n = n + 1
      part(offset+1:offset+n) = j
      offset = offset + n
    end do
    ASSERT(offset == size(part))

  end subroutine get_block_partition

end module simple_partitioning_methods
