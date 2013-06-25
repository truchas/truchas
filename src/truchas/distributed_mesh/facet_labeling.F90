#include "f90_assert.fpp"

module facet_labeling

  use hashing
  implicit none
  private
  
  public :: create_facet_table, destroy_facet_table, get_facet_label, get_facet_label2
  public :: number_of_facets, fill_factor, write_hash_performance
  public :: dump_facet_table
  
  type, private :: table_record
    integer, pointer :: facet(:) => null()
    integer :: label = 0
  end type table_record
  
  type, public :: facet_table
    private
    integer :: nfacet = 0   ! number of facets stored in the table
    integer :: tsize = 0    ! Table size
    type(hash_param) :: hpar
    type(table_record), pointer :: record(:) => null()
    !! Performance counters
    integer :: nss = 0  ! Number of successful searches
    integer :: nus = 0  ! Number of unsuccessful searches
    integer :: nsp = 0  ! Number of probes in successful searches
    integer :: nup = 0  ! Number of probes in unsuccessful searches
    integer :: msp = 0  ! Maximum probes in a successful probe
    integer :: mup = 0  ! Maximum probes in an unsuccessful probe
  end type facet_table

contains

  subroutine get_facet_label (table, facet, label)
  
    use cell_topology

    type(facet_table), intent(inout) :: table
    integer, pointer :: facet(:)
    integer, intent(out) :: label
    
    integer :: p, np, i, j

    ASSERT( associated(facet) )
    ASSERT( size(facet) > 1 )
    ASSERT( all(facet > 0) )
    ASSERT( .not.degenerate(facet) )
    
    call normalize_facet (facet)
    
    call hash (table%hpar, facet, i, j)
    
    np = 0
    label = 0
    
    !! Search the hash table for the facet
    do while (associated(table%record(i)%facet))
      np = np + 1 ! update the number of probes
      p = parity(facet, table%record(i)%facet)
      if (p /= 0) then  ! we located the facet
        label = p * table%record(i)%label
        !! Update performance counters
        table%nss = table%nss + 1
        table%nsp = table%nsp + np
        table%msp = max(np, table%msp)
        deallocate(facet)
        return
      end if
      !! Next probe address.
      i = i - j
      if (i < 0) i = i + table%tsize
    end do
  
    !! Facet not found; insert it at this location.
    if (table%nfacet >= table%tsize - 1) return ! table is full
    table%nfacet = 1 + table%nfacet
    label = table%nfacet
    table%record(i)%label = label
    table%record(i)%facet => facet
    facet => null()
    !! Update performance counters
    table%nus = table%nus + 1
    table%nup = table%nup + np
    table%mup = max(np, table%mup)

  end subroutine get_facet_label
  
  subroutine get_facet_label2 (table, facet, label)
  
    use cell_topology

    type(facet_table), intent(inout) :: table
    integer, intent(in) :: facet(:)
    integer, intent(out) :: label
    
    integer :: p, np, i, j

    ASSERT( size(facet) > 1 )
    ASSERT( all(facet > 0) )
    ASSERT( .not.degenerate(facet) )
    
    !!! the input facet must be normalized !!!
    !call normalize_facet (facet)
    
    call hash (table%hpar, facet, i, j)
    
    np = 0
    label = 0
    
    !! Search the hash table for the facet
    do while (associated(table%record(i)%facet))
      np = np + 1 ! update the number of probes
      p = parity(facet, table%record(i)%facet)
      if (p /= 0) then  ! we located the facet
        label = p * table%record(i)%label
        !! Update performance counters
        table%nss = table%nss + 1
        table%nsp = table%nsp + np
        table%msp = max(np, table%msp)
        return
      end if
      !! Next probe address.
      i = i - j
      if (i < 0) i = i + table%tsize
    end do
  
  end subroutine get_facet_label2
  
  subroutine create_facet_table (table, max_facet, node_max)
  
    type(facet_table), intent(out) :: table
    integer, intent(in) :: max_facet  ! maximum number of facets to be stored
    integer, intent(in) :: node_max   ! maximum node label referenced by a facet
    
    ASSERT( max_facet > 0 )
    ASSERT( node_max > 0 )
    
    table%tsize = max_facet
    call initialize_hash_param (table%hpar, table%tsize, node_max)
    
    allocate(table%record(0:table%tsize-1))
    
  end subroutine create_facet_table


  subroutine destroy_facet_table (table)
  
    type(facet_table), intent(inout) :: table
    
    integer :: j
    type(facet_table) :: default
    
    if (associated(table%record)) then
      do j = lbound(table%record,dim=1), ubound(table%record,dim=1)
        if (associated(table%record(j)%facet)) deallocate(table%record(j)%facet)
      end do
      deallocate(table%record)
    end if
    
    table = default ! set default values
    
  end subroutine destroy_facet_table

  subroutine dump_facet_table (table, unit)
    type(facet_table), intent(in) :: table
    integer, intent(in) :: unit
    integer :: j
    write(unit,*) 'TSIZE=', table%tsize
    write(unit,*) 'NFACET=', table%nfacet
    if (.not.associated(table%record)) then
      write(unit,*) 'RECORD NOT ASSOCIATED'
      return
    end if
    do j = lbound(table%record,dim=1), ubound(table%record,dim=1)
      if (associated(table%record(j)%facet)) then
        write(unit,*) 'RECORD(', j, ')%LABEL=', table%record(j)%label, ', %FACET=', table%record(j)%facet
      else
        write(unit,*) 'RECORD(', j, ') IS UNUSED'
      end if
    end do
  end subroutine dump_facet_table

  pure integer function number_of_facets (table)
    type(facet_table), intent(in) :: table
    number_of_facets = table%nfacet
  end function number_of_facets
  
  pure real function fill_factor (table)
    type(facet_table), intent(in) :: table
    fill_factor = real(table%nfacet) / real(table%tsize)
  end function fill_factor
  
  subroutine write_hash_performance (table, lun)
    type(facet_table), intent(in) :: table
    integer, intent(in) :: lun
    write(lun,'(/,a)') 'HASH TABLE PERFORMANCE'
    write(lun,'(a)') '-----------------------------------------'
    write(lun,*) '                Fill factor:', real(table%nfacet) / real(table%tsize)
    write(lun,*) '          Number of entries:', table%nfacet
    write(lun,*) '    Avg probes, succ search:', real(table%nsp) / real(table%nss)
    write(lun,*) '  Avg probes, unsucc search:', real(table%nup) / real(table%nus)
    write(lun,*) '    Max probes, succ search:', table%msp
    write(lun,*) '  Max probes, unsucc search:', table%mup
    write(lun,*) 'Succ to unsucc search ratio:', real(table%nss) / real(table%nus)
  end subroutine write_hash_performance
  
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

end module facet_labeling
