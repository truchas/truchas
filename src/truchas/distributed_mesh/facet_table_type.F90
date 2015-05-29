!!
!! FACET_TABLE_TYPE
!!
!! This module defines an efficient data structure for enumerating the facets
!! in a mesh (faces or edges) and establishing their orientation.  This is
!! intended to be a short-lived structure used while instantiating the final
!! mesh data structure.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for Fortran 2008, Jan 2014
!!
!! PROGRAMMING INTERFACE
!!
!!  The module defines the FACET_TABLE derived type.  It has the following type
!!  bound procedures.
!!
!!  INIT (MAX_FACET, NODE_MAX) initialized the table.  MAX_FACET should be
!!    somewhat larger that the number of facets that will be inserted into
!!    the table.  This determines the size of the internal hash table used
!!    in the implementation; too small and there is insufficient space or
!!    excessive hash collisions, too large and there is wasted memory use.
!!    The facets will be input as lists of node numbers, and NODE_MAX should
!!    be the maximum (or larger) of those numbers.  Since every node belongs
!!    to a face or edge, this is typically the number of nodes in the mesh.
!!
!!  NUMBER_OF_FACETS() returns the number of facets stored in the table.
!!
!!  FILL_FACTOR() returns the fraction of the table that is used.
!!
!!  GET_FACET_LABEL (FACET, LABEL, INSERT) returns the integer LABEL assigned
!!    to the facet specified by the rank-1 integer array FACET.  If the facet
!!    is not found in the table and INSERT is true, the facet is inserted and
!!    assigned the next numeric label in sequence starting with 1.  When the
!!    facet is found in the table, its orientation is compared with the stored
!!    facet.  If it is the same, the (positive) numeric label of the stored
!!    facet is returned; otherwise it has the opposite orientation and the
!!    negative of the numeric label is returned. If INSERT is false, a 0 label
!!    indicates the face was not found; otherwise a 0 label indicates a full
!!    table and is an error.
!!
!!  WRITE_HASH_PERFORMANCE (LUN) writes a summary of the hash table performance
!!    to logical unit LUN, which must be open for writing.  This can be helpful
!!    in development and debugging to identify issues in right-sizing the hash
!!    table.
!!
!!  DUMP_FACET_TABLE (UNIT) dumps the contents of the hash table to the
!!    specified logical unit in a readable format.  For debugging.
!!

#include "f90_assert.fpp"

module facet_table_type

  use facet_hash_type
  implicit none
  private

  type, private :: table_record
    integer, allocatable :: facet(:)
    integer :: label = 0
  end type table_record

  type, public :: facet_table
    private
    integer :: nfacet = 0   ! number of facets stored in the table
    integer :: tsize = 0    ! Table size
    type(facet_hash) :: fh
    type(table_record), allocatable :: record(:)
    !! Performance counters
    integer :: nss = 0  ! Number of successful searches
    integer :: nus = 0  ! Number of unsuccessful searches
    integer :: nsp = 0  ! Number of probes in successful searches
    integer :: nup = 0  ! Number of probes in unsuccessful searches
    integer :: msp = 0  ! Maximum probes in a successful probe
    integer :: mup = 0  ! Maximum probes in an unsuccessful probe
  contains
    procedure :: init
    procedure :: get_facet_label
    procedure :: number_of_facets
    procedure :: fill_factor
    procedure :: write_hash_performance
    procedure :: dump_facet_table
  end type facet_table

contains

  subroutine init (this, max_facet, node_max)
    class(facet_table), intent(out) :: this
    integer, intent(in) :: max_facet  ! maximum number of facets to be stored
    integer, intent(in) :: node_max   ! maximum node label referenced by a facet
    ASSERT(max_facet > 0)
    ASSERT(node_max > 0)
    this%tsize = max_facet
    call this%fh%init (this%tsize, node_max)
    allocate(this%record(0:this%tsize-1))
  end subroutine init

  pure integer function number_of_facets (this)
    class(facet_table), intent(in) :: this
    number_of_facets = this%nfacet
  end function number_of_facets

  pure real function fill_factor (this)
    class(facet_table), intent(in) :: this
    fill_factor = real(this%nfacet) / real(this%tsize)
  end function fill_factor

  subroutine get_facet_label (this, facet, label, insert)

    use cell_topology, only: normalize_facet, parity

    class(facet_table), intent(inout) :: this
    integer, intent(inout) :: facet(:)
    integer, intent(out) :: label
    logical, intent(in)  :: insert

    integer :: p, np, i, j

    ASSERT(size(facet) > 1)
    ASSERT(all(facet > 0))
    ASSERT(.not.degenerate(facet))

    call normalize_facet (facet)
    call this%fh%hash (facet, i, j)

    np = 0
    label = 0

    !! Search the hash table for the facet
    do while (allocated(this%record(i)%facet))
      np = np + 1 ! update the number of probes
      p = parity(facet, this%record(i)%facet)
      if (p /= 0) then  ! we located the facet
        label = p * this%record(i)%label
        !! Update performance counters
        this%nss = this%nss + 1
        this%nsp = this%nsp + np
        this%msp = max(np, this%msp)
        return
      end if
      !! Next probe address.
      i = i - j
      if (i < 0) i = i + this%tsize
    end do

    if (.not.insert) return

    !! Facet not found; insert it at this location.
    if (this%nfacet >= this%tsize - 1) return ! table is full
    this%nfacet = 1 + this%nfacet
    label = this%nfacet
    this%record(i)%label = label
    this%record(i)%facet = facet
    !! Update performance counters
    this%nus = this%nus + 1
    this%nup = this%nup + np
    this%mup = max(np, this%mup)

  end subroutine get_facet_label


  subroutine dump_facet_table (this, unit)
    class(facet_table), intent(in) :: this
    integer, intent(in) :: unit
    integer :: j
    write(unit,*) 'TSIZE=', this%tsize
    write(unit,*) 'NFACET=', this%nfacet
    if (.not.allocated(this%record)) then
      write(unit,*) 'RECORD NOT ASSOCIATED'
      return
    end if
    do j = lbound(this%record,dim=1), ubound(this%record,dim=1)
      if (allocated(this%record(j)%facet)) then
        write(unit,*) 'RECORD(', j, ')%LABEL=', this%record(j)%label, &
                      ', %FACET=', this%record(j)%facet
      else
        write(unit,*) 'RECORD(', j, ') IS UNUSED'
      end if
    end do
  end subroutine dump_facet_table

  subroutine write_hash_performance (this, lun)
    class(facet_table), intent(in) :: this
    integer, intent(in) :: lun
    write(lun,'(/,a)') 'HASH TABLE PERFORMANCE'
    write(lun,'(a)') '-----------------------------------------'
    write(lun,*) '                Fill factor:', real(this%nfacet) / real(this%tsize)
    write(lun,*) '          Number of entries:', this%nfacet
    write(lun,*) '    Avg probes, succ search:', real(this%nsp) / real(this%nss)
    write(lun,*) '  Avg probes, unsucc search:', real(this%nup) / real(this%nus)
    write(lun,*) '    Max probes, succ search:', this%msp
    write(lun,*) '  Max probes, unsucc search:', this%mup
    write(lun,*) 'Succ to unsucc search ratio:', real(this%nss) / real(this%nus)
  end subroutine write_hash_performance

  !! Auxillary function returns true if the facet is degnerate.
  pure logical function degenerate (facet)
    integer, intent(in) :: facet(:)
    integer :: j, k
    degenerate = .true.
    do j = 1, size(facet)-1
      do k = j+1, size(facet)
        if (facet(j) == facet(k)) return
      end do
    end do
    degenerate = .false.
  end function degenerate

end module facet_table_type
