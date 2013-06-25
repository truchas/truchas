!! MATERIAL_MESH_FUNCTION
!!
!! The derived type MAT_MF is an opaque structure that describes how materials
!! are distributed across a mesh.  Conceptually, the mesh cells are partitioned
!! into regions, and a set of one or more materials associated with each region.
!! A cell may contain any mixture of the materials associated with the region
!! to which the cell belongs.  The cells in a single-material region obviously
!! contain only that material.  Multi-material regions include material volume
!! fraction data for each cell that describes the mixture of materials.
!!
!! ACCESSOR PROCEDURES
!!
!!  MMF_NUM_REG(THIS) returns the number of regions in the MAT_MF object THIS.
!!
!!  MMF_NUM_REG_MAT(THIS, N) returns the number of materials that may be
!!    present in region N of the MAT_MF object THIS.
!!
!!  MMF_REG_MATIDS(THIS, N) returns a pointer to a rank-1 integer array
!!    containing the IDs of the materials that may be present in region N of
!!    the MAT_MF object THIS.  The size of the target array is
!!    MMF_NUM_REG_MAT(THIS, N).  The ID list is ordered.
!!
!!  MMF_REG_CELLS(THIS, N) returns a pointer to a rank-1 integer array
!!    containing the indices of the cells that comprise region N of the
!!    MAT_MF object THIS.  The index list is ordered.
!!
!!  MMF_REG_VOL_FRAC(THIS, N) returns a pointer to a rank-2 real array
!!    containing the material volume fractions for region N of the MAT_MF
!!    object THIS.  For single-material regions a null pointer is returned.
!     If
!!
!!      MATID => MMF_REG_MATIDS(THIS, N),
!!      CELL  => MMF_REG_CELLS(THIS, N),
!!      VF    => MMV_REG_VOL_FRAC(THIS, N),
!!
!!    then VF(j,k) is the volume fraction of MATID(k) in CELL(j).  VF has
!!    shape [SIZE(CELL), SIZE(MATID)]. (Comment: should we swap dimensions?)
!!    Application code is expected to define the values of the target array;
!!    the MAT_MF object merely manages the storage.
!!
!!  MMF_MESH(THIS) returns a pointer to the distributed mesh object on which
!!    the MAT_MF object THIS is based.
!!
!! CREATION/DESTRUCTION PROCEDURES
!!
!!  CALL MMF_PREP (THIS, MESH) prepares the MAT_MF object THIS to begin
!!    receiving the data associated with one or more mesh material regions.
!!    MESH is the distributed mesh.  The structure THIS maintains a pointer
!!    to this mesh (read-only), and so if the actual argument is not a pointer
!!    it must be given the target attribute.
!!
!!  CALL MMF_DEFINE_REGION (THIS, SETID, MATID, STAT, ERRMSG) defines a
!!    material region in the MAT_MF object THIS.  The region consists of the
!!    group of cells belonging to the cell set IDs given in the rank-1 integer
!!    array SETID, and the IDs of the materials that may be present in the
!!    region are given in the rank-1 integer array MATID.  A nonzero value is
!!    assigned to the integer STAT if an error condition occurs, and an
!!    explanatory message is assigned to the character string ERRMSG.
!!    Possible errors are referring to an unknown cell set ID, or attempting
!!    to define overlapping material regions.  This routine can be called
!!    multiple times after the initial call to MMF_PREP.
!!
!!  CALL MMF_DONE (THIS, STAT, ERRMSG) performs the final configuration of
!!    material mesh function object THIS after all the desired calls to
!!    MMF_DEFINE_REGION have been made.  A nonzero value is assigned to
!!    the integer STAT if an error condition occurs, and an explanatory
!!    message is assigned to the character string ERRMSG.  It is an error
!!    if the defined regions do not cover the entire mesh.
!!
!!  CALL DESTROY (THIS) deallocates any storage associated with the MMF
!!    object THIS, returning it to its default initialization state.
!!

#include "f90_assert.fpp"

module material_mesh_function

  use kinds, only: r8
  use distributed_mesh, only: dist_mesh
  use string_utilities, only: i_to_c

  implicit none
  private
  
  public :: mmf_num_reg, mmf_num_reg_mat, mmf_reg_matID, mmf_reg_cell, mmf_reg_vol_frac, mmf_mesh
  public :: mmf_prep, mmf_define_region, mmf_done, mmf_destroy
  public :: mmf_dump
  public :: mmf_get_all_matID
  public :: mmf_reg_has_matID, mmf_has_matID
  public :: mmf_get_mat_vol_frac, mmf_get_void_vol_frac
  
  type, public :: mat_mf
    private
    type(dist_mesh), pointer :: mesh => null()
    integer :: nreg = -1
    type(region), pointer :: reg(:) => null()
    !! Temporary components used during setup.
    integer, pointer :: tag(:) => null()
    type(list), pointer :: matID_list => null()
  end type

  type :: region
    integer,  pointer :: matID(:)   => null() ! list of material IDs for the region
    integer,  pointer :: cell(:)    => null() ! list of cells comprising the region
    real(r8), pointer :: vfrac(:,:) => null() ! volume fraction array
  end type
  
  type :: list
    integer, pointer :: matID(:) => null()
    type(list), pointer :: next => null()
  end type list
  
contains

  integer function mmf_num_reg (this)
    type(mat_mf), intent(in) :: this
    mmf_num_reg = size(this%reg)
  end function mmf_num_reg
  
  integer function mmf_num_reg_mat (this, n)
    type(mat_mf), intent(in) :: this
    integer, intent(in) :: n
    mmf_num_reg_mat = size(this%reg(n)%matID)
  end function mmf_num_reg_mat
  
  function mmf_reg_matID (this, n) result (list)
    type(mat_mf), intent(in) :: this
    integer, intent(in) :: n
    integer, pointer :: list(:)
    list => this%reg(n)%matID
  end function mmf_reg_matID
  
  logical function mmf_reg_has_matID (this, n, matID)
    type(mat_mf), intent(in) :: this
    integer, intent(in) :: n, matID
    integer :: j
    mmf_reg_has_matID = .true.
    do j = 1, size(this%reg(n)%matID)
      if (matID == this%reg(n)%matID(j)) return
      if (matID < this%reg(n)%matID(j)) exit ! list is sorted
    end do
    mmf_reg_has_matID = .false.
  end function mmf_reg_has_matID
  
  logical function mmf_has_matID (this, matID)
    type(mat_mf), intent(in) :: this
    integer, intent(in) :: matID
    integer :: n
    mmf_has_matID = .true.
    do n = 1, this%nreg
      if (mmf_reg_has_matID(this, n, matID)) return
    end do
    mmf_has_matID = .false.
  end function mmf_has_matID
  
  subroutine mmf_get_all_matID (this, list, drop_void)
    type(mat_mf), intent(in) :: this
    integer, pointer :: list(:)
    logical, optional :: drop_void
    integer :: i, n, offset
    integer, allocatable :: temp(:)
    logical :: drop_zero
    drop_zero = .false.
    if (present(drop_void)) drop_zero = drop_void
    n = 0
    do i = 1, this%nreg
      n = n + size(this%reg(i)%matID)
    end do
    allocate(temp(n))
    offset = 0
    do i = 1, this%nreg
      n = size(this%reg(i)%matID)
      temp(offset+1:offset+n) = this%reg(i)%matID
      offset = offset + n
    end do
    call sort_list (temp, n)
    if (temp(1) == 0 .and. drop_zero) then
      allocate(list(n-1))
      list = temp(2:n)
    else
      allocate(list(n))
      list = temp(:n)
    end if
    deallocate(temp)
  end subroutine mmf_get_all_matID
  
  function mmf_reg_cell (this, n) result (list)
    type(mat_mf), intent(in) :: this
    integer, intent(in) :: n
    integer, pointer :: list(:)
    list => this%reg(n)%cell
  end function mmf_reg_cell
  
  function mmf_reg_vol_frac (this, n) result (vf)
    type(mat_mf), intent(in) :: this
    integer, intent(in) :: n
    real(r8), pointer :: vf(:,:)
    vf => this%reg(n)%vfrac
  end function mmf_reg_vol_frac
  
  function mmf_mesh (this) result (mesh)
    type(mat_mf), intent(in) :: this
    type(dist_mesh), pointer :: mesh
    mesh => this%mesh
  end function mmf_mesh
  
  subroutine mmf_prep (this, mesh)
  
    type(mat_mf), intent(out) :: this
    type(dist_mesh), intent(in), target :: mesh
    
    this%mesh => mesh
    this%nreg = 0
    allocate(this%tag(mesh%ncell))
    this%tag = 0
    
  end subroutine mmf_prep
  
  subroutine mmf_define_region (this, setID, matID, stat, errmsg)
  
    type(mat_mf), intent(inout) :: this
    integer,   intent(in) :: setID(:)
    integer,   intent(in) :: matID(:)
    integer,   intent(out) :: stat
    character(len=*), intent(out) :: errmsg
    
    integer :: n, sorted(size(matID))
    type(list), pointer :: first
    
    !! Verify that THIS is in the correct state.
    INSIST( this%nreg>=0 .and. associated(this%tag) )
    
    call set_tag_array (this, setID, stat, errmsg)
    if (stat /= 0) return
    
    INSIST( size(matID) > 0 )
    INSIST( all(matID >= 0) )
    
    !! Sort the supplied material ID array and remove any duplicates.
    sorted = matID
    call sort_list (sorted, n)
    
    !! Prepend the processed material ID array to the MATID_LIST list.
    allocate(first)
    allocate(first%matID(n))
    first%matID = sorted(:n)
    first%next => this%matID_list
    this%matID_list => first
    
  end subroutine mmf_define_region

  !!
  !! This auxillary subroutine tags the cells identified by the list of
  !! cell set IDs, checking for error conditions in the process.
  !!

  subroutine set_tag_array (this, setID, stat, errmsg)

    type(mat_mf), intent(inout) :: this
    integer, intent(in) :: setID(:)
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: i, j
    logical :: mask(this%mesh%ncell)
    integer(kind(this%mesh%cell_set_mask)) :: bitmask

    !! Create the bitmask corresponding to setID.
    bitmask = 0
    do i = 1, size(setID)
      do j = size(this%mesh%cell_set_ID), 1, -1
        if (setID(i) == this%mesh%cell_set_ID(j)) exit
      end do
      if (j == 0) then
        stat = 1
        errmsg = 'unknown cell set ID: ' // i_to_c(setID(i))
        return
      end if
      bitmask = ibset(bitmask, j)
    end do

    !! Identify the cells specified by setID.
    mask = (iand(bitmask, this%mesh%cell_set_mask) /= 0)

    !! Check that these cells don't overlap those from preceding calls.
    if (any(mask .and. this%tag /= 0)) then
      stat = 1
      errmsg = 'region overlaps a preceding region'
      return
    end if

    !! Set the tag array.
    this%nreg = 1 + this%nreg
    this%tag = merge(this%nreg, this%tag, mask)

    stat = 0
    errmsg = ''

  end subroutine set_tag_array
  
  !!
  !! Auxillary routine to (insertion) sort the short material ID list,
  !! with duplicate values dropped; N is the final list length.
  !!
  
  subroutine sort_list (list, n)
  
    integer, intent(inout) :: list(:)
    integer, intent(out) :: n
    
    integer :: j, k, next

    !! Insertion sort of the list.
    do j = 2, size(list)
      k = j
      next = list(j)
      do while (k > 1)  ! find location k of next value
        if (list(k-1) <= next) exit
        list(k) = list(k-1)
        k = k - 1
      end do
      list(k) = next
    end do
    
    !! Cull duplicates.
    n = min(1, size(list))
    do j = 2, size(list)
     if (list(j) == list(n)) cycle
     n = n + 1
     if (j /= n) list(n) = list(j)
    end do
    
  end subroutine sort_list
  
  subroutine mmf_done (this, stat, errmsg)
  
    type(mat_mf), intent(inout) :: this
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg
    
    integer :: n, ncell, nmat
    integer, allocatable :: index(:)
    type(list), pointer :: first
    
    !! Verify that THIS is in the correct state.
    INSIST( this%nreg>= 0 .and. associated(this%tag) )
    
    ASSERT( minval(this%tag)>=0 .and. maxval(this%tag)<=this%nreg )
    
    !! Verify that every cell has been associated with a region.
    n = count(this%tag <= 0)
    if (n /= 0) then
      stat = 1
      errmsg = i_to_c(n) // ' cells not associated with a region'
      return
    end if
    
    !! Initialize the REG structure-array component of THIS using the
    !! accumulated data in the temporary MATID_LIST and TAG components.
    !! Note that MATID_LIST is treated as a stack (LIFO) so in order to avoid
    !! unnecessary surprises to the user we fill the REG array from the end
    !! as we pop data off MATID_LIST which preserves the order of regions as
    !! defined by the user (not strictly necessary).

    allocate(this%reg(this%nreg))
    
    allocate(index(size(this%tag))) ! used for generating cell index arrays
    forall (n = 1:size(index)) index(n) = n
    
    do n = this%nreg, 1, -1
    
      !! Unlink the first element in the region list.
      first => this%matID_list
      ASSERT( associated(first) )
      this%matID_list => first%next
      
      !! Take the material ID list from FIRST.
      this%reg(n)%matID => first%matID
      deallocate(first) ! we're finished with it now
      
      !! Initialize the cell index array for this region.
      ncell = count(this%tag==n)
      allocate(this%reg(n)%cell(ncell))
      this%reg(n)%cell = pack(index, mask=(this%tag==n))
      
      !! Allocate the volume fraction array for this region if appropriate.
      nmat = size(this%reg(n)%matID)
      if (nmat > 1) allocate(this%reg(n)%vfrac(ncell,nmat))
      
    end do
    
    ASSERT( .not.associated(this%matID_list) )
    deallocate(index, this%tag)
    
    stat = 0
    errmsg = ''
    
  end subroutine mmf_done
  
  subroutine mmf_destroy (this)
    type(mat_mf), intent(inout) :: this
    type(mat_mf) :: default
    integer :: n
    type(list), pointer :: first
    if (associated(this%reg)) then
      do n = 1, size(this%reg)
        if (associated(this%reg(n)%matID)) deallocate(this%reg(n)%matID)
        if (associated(this%reg(n)%cell))  deallocate(this%reg(n)%cell)
        if (associated(this%reg(n)%vfrac)) deallocate(this%reg(n)%vfrac)
      end do
      deallocate(this%reg)
    end if
    if (associated(this%tag)) deallocate(this%tag)
    do while (associated(this%matID_list))
      first => this%matID_list
      this%matID_list => first%next
      if (associated(first%matID)) deallocate(first%matID)
      deallocate(first)
    end do
    this = default  ! assign default initialization values
  end subroutine mmf_destroy
  
  subroutine mmf_dump (this, unit)
  
    type(mat_mf), intent(in) :: this
    integer, intent(in) :: unit
    
    integer :: n, j
    type(list), pointer :: l
    
    write(unit,'(/,a,i0)') '%nreg=', this%nreg
    
    if (associated(this%reg)) then
      do n = 1, size(this%reg)
        write(unit,'(a,i0,a)') '%reg(', n, ')='
        if (associated(this%reg(n)%matID)) then
          write(unit,'(2x,a,(t11,9(i0,:,1x)))') '%matID=', this%reg(n)%matID
        else
          write(unit,'(2x,a)') '%matID => null'
        end if
        if (associated(this%reg(n)%cell)) then
          write(unit,'(2x,a,(t11,9(i0,:,1x)))') '%cell=', this%reg(n)%cell
        else
          write(unit,'(2x,a)') '%cell => null'
        end if
        if (associated(this%reg(n)%vfrac)) then
          do j = 1, size(this%reg(n)%vfrac,dim=1)
            write(unit,'(2x,a,i0,a,(t17,9f7.4))') '%vfrac(', j, ',:)=', this%reg(n)%vfrac(j,:)
          end do
        else
          write(unit,'(2x,a)') '%vfrac => null'
        end if
      end do
    else
      write(unit,'(a)') '%reg => null'
    end if
    
    if (associated(this%tag)) then
      write(unit,'(a,(t11,9(i0,:,1x)))') '%tag=', this%tag
    else
      write(unit,'(a)') '%tag => null'
    end if
    
    if (associated(this%matID_list)) then
      write(unit,'(a)') '%matID_list =>'
      l => this%matID_list
      do while (associated(l))
        write(unit,'(2x,a,(t11,9(i0,:,1x)))') '%matID=', l%matID
        l => l%next
      end do
    else
      write(unit,'(a)') '%matID_list => null'
    end if
    
  end subroutine mmf_dump
  
  subroutine mmf_get_void_vol_frac (mmf, void_vol_frac)
    type(mat_mf), intent(in) :: mmf
    real(r8), intent(out) :: void_vol_frac(:)
    call mmf_get_mat_vol_frac (mmf, 0, void_vol_frac)
  end subroutine mmf_get_void_vol_frac
  
  subroutine mmf_get_mat_vol_frac (mmf, matid, vol_frac)
  
    type(mat_mf), intent(in) :: mmf
    integer, intent(in) :: matid
    real(r8), intent(out) :: vol_frac(:)
    
    integer :: n, i
    integer, pointer :: matids(:), cell(:)
    real(r8), pointer :: vfrac(:,:)
    type(dist_mesh), pointer :: mesh
    
    mesh => mmf_mesh(mmf)
    ASSERT(size(vol_frac) == mesh%ncell)
    
    vol_frac = 0.0_r8
    do n = 1, mmf_num_reg(mmf)
      matids => mmf_reg_matID(mmf, n)
      do i = 1, size(matids)
        if (matids(i) == matid) then
          vfrac => mmf_reg_vol_frac(mmf, n)
          cell  => mmf_reg_cell(mmf, n)
          if (associated(vfrac)) then
            vol_frac(cell) = vfrac(:,i)
          else  ! single material region
            vol_frac(cell) = 1.0_r8
          end if
          exit
        end if
      end do
    end do
    
  end subroutine mmf_get_mat_vol_frac

end module material_mesh_function
