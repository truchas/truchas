!!
!! EXODUS_MESH_TYPE
!!
!! Neil N. Carlson <nnc@lanl.gov> 16 Sep 2004
!! Last revised 10 Nov 2004
!!
!! This module provides derived data types for encapsulating the mesh data
!! described by an Exodus II file, and several basic procedures that operate
!! on those types.
!!
!! This is a private module; application code should only use the top-level
!! module EXODUS.
!!
!! The derived types are designed to mirror the structure of the data as
!! stored in the Exodus II file (see Appendix A in [1]).  Moreover, the
!! derived types only provide for that data which is currently useful.
!! For example, they do not support
!!  o optional node and element number maps that specify a mapping from the
!!    internal numbering to a user-space numbering;
!!  o an optional element order map that specifies a 'good' order to process
!!    the elements;
!!  o distribution factors for node sets and side sets;
!!  o element attributes used for some esoteric element types.
!!
!! [1] L.A.Schoof and V.R.Yarberry, "Exodus II: A Finite Element Data Model",
!!     Sandia report SAND92-2137.  This can be obtained at
!!     http://endo.sandia.gov/SEACAS/Documentation/exodusII.pdf
!!

module exodus_mesh_type

  implicit none
  private

  public :: defined, destroy

  !! These are really only useful in testing situations.
  public :: operator(.eq.), operator(.ne.), dump_type

  integer, parameter :: r64 = selected_real_kind(10,50) ! 64-bit IEEE float (hopefully!)

  !! Values documented in [1], and consistent with Cubit 8/9 output.
  integer, parameter, public :: MAX_STR_LENGTH  = 32  ! Maximum character string length
  integer, parameter, public :: MAX_LINE_LENGTH = 80  ! Maximum character line length

  type, public :: exodus_mesh
    integer :: num_dim  = 0   ! spatial dimension of the mesh
    integer :: num_node = 0   ! number of nodes in the mesh
    integer :: num_elem = 0   ! number of elements in the mesh
    integer :: num_eblk = 0   ! number of element blocks
    integer :: num_nset = 0   ! number of node sets
    integer :: num_sset = 0   ! number of side sets
    type(elem_blk), pointer :: eblk(:) => null()    ! list of element blocks
    type(node_set), pointer :: nset(:) => null()    ! list of node sets
    type(side_set), pointer :: sset(:) => null()    ! list of side sets
    real(kind=r64), pointer :: coord(:,:) => null() ! node coordinates
    character(len=MAX_LINE_LENGTH) :: title
  end type exodus_mesh

  type, public :: elem_blk
    integer :: ID = 0                           ! element block ID
    integer :: num_elem = 0                     ! number of elements in the element block
    character(len=MAX_STR_LENGTH) :: elem_type  ! type of elements in the element block
    integer, pointer :: connect(:,:) => null()  ! list of nodes that define each element
  end type elem_blk

  type, public :: node_set
    integer :: ID = 0                      ! node set ID
    integer :: num_node = 0                ! number of nodes in the node set
    integer, pointer :: node(:) => null()  ! list of nodes
  end type node_set

  type, public :: side_set
    integer :: ID = 0                      ! side set ID
    integer :: num_side = 0                ! number of sides in the side set
    integer, pointer :: elem(:) => null()  ! list of sides described as
    integer, pointer :: face(:) => null()  !   (element, local face) pairs
  end type side_set

  interface defined
    module procedure defined_exodus_mesh, defined_elem_blk, defined_node_set, defined_side_set
  end interface

  interface destroy
    module procedure destroy_exodus_mesh, destroy_elem_blk, destroy_node_set, destroy_side_set
  end interface

  interface operator(.eq.)
    module procedure eq_exodus_mesh, eq_elem_blk, eq_node_set, eq_side_set
  end interface

  interface operator(.ne.)
    module procedure ne_exodus_mesh, ne_elem_blk, ne_node_set, ne_side_set
  end interface

  interface dump_type
    module procedure dump_exodus_mesh, dump_elem_blk, dump_node_set, dump_side_set
  end interface dump_type

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Specific procedures for the generic DESTROY
 !!
 !! These procedures deallocate any allocated component of the argument,
 !! and assign default values to the remaining components, returning the
 !! argument to the default initialization state.
 !!

  elemental subroutine destroy_exodus_mesh (mesh)
    type(exodus_mesh), intent(inout) :: mesh
    type(exodus_mesh) :: default
    if (associated(mesh%coord)) deallocate(mesh%coord)
    if (associated(mesh%eblk)) then
      call destroy (mesh%eblk)
      deallocate(mesh%eblk)
    end if
    if (associated(mesh%nset)) then
      call destroy (mesh%nset)
      deallocate(mesh%nset)
    end if
    if (associated(mesh%sset)) then
      call destroy (mesh%sset)
      deallocate(mesh%sset)
    end if
    mesh = default  ! assign default initialization value
  end subroutine destroy_exodus_mesh


  elemental subroutine destroy_elem_blk (eb)
    type(elem_blk), intent(inout) :: eb
    type(elem_blk) :: default
    if (associated(eb%connect)) deallocate(eb%connect)
    eb = default  ! assign default initialization value
  end subroutine destroy_elem_blk


  elemental subroutine destroy_node_set (ns)
    type(node_set), intent(inout) :: ns
    type(node_set) :: default
    if (associated(ns%node)) deallocate(ns%node)
    ns = default  ! assign default initialization value
  end subroutine destroy_node_set


  elemental subroutine destroy_side_set (ss)
    type(side_set), intent(inout) :: ss
    type(side_set) :: default
    if (associated(ss%elem)) deallocate(ss%elem)
    if (associated(ss%face)) deallocate(ss%face)
    ss = default  ! assign default initialization value
  end subroutine destroy_side_set

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DEFINED_EXODUS_MESH (specific procedure for generic DEFINED)
 !!
 !! This function returns the value true if the exodus mesh data structure MESH
 !! is defined; otherwise it returns the value false.  Defined means that the
 !! components of the variable are properly defined and at least superficially
 !! reasonable; however, it does not examine the topology or geometry of the
 !! mesh to ensure it is valid.
 !!
 !! For element blocks, in particular, the element type and number of nodes per
 !! element are not examined for consistency as Exodus II does not prescribe
 !! a set of allowed element types.  These values are left open for application
 !! interpretation, and so the application must ensure that those values are
 !! consistent.
 !!
 !! This function is expected to be used mostly in assertion checks.
 !!

  elemental logical function defined_exodus_mesh (mesh)

    type(exodus_mesh), intent(in) :: mesh

    integer :: n

    CHECKLIST: do

      defined_exodus_mesh = .false.

      !! Check the node coordinate info.
      if (.not.associated(mesh%coord)) exit CHECKLIST
      if (mesh%num_dim /= size(mesh%coord,dim=1)) exit CHECKLIST
      if (mesh%num_dim < 1 .or. mesh%num_dim > 3) exit CHECKLIST
      if (mesh%num_node /= size(mesh%coord,dim=2)) exit CHECKLIST

      !! Check the element block info.
      if (.not.associated(mesh%eblk)) exit CHECKLIST
      if (.not.all(defined(mesh%eblk))) exit CHECKLIST
      if (.not.unique(mesh%eblk%ID)) exit CHECKLIST
      if (mesh%num_elem /= sum(mesh%eblk%num_elem)) exit CHECKLIST
      do n = 1, size(mesh%eblk) ! Check that the CONNECT values are in-range.
        if (minval(mesh%eblk(n)%connect) < 1) exit CHECKLIST
        if (maxval(mesh%eblk(n)%connect) > mesh%num_node) exit CHECKLIST
      end do

      !! If present, check the node set info.
      if (associated(mesh%nset)) then
        if (mesh%num_nset /= size(mesh%nset)) exit CHECKLIST
        if (.not.all(defined(mesh%nset))) exit CHECKLIST
        if (.not.unique(mesh%nset%ID)) exit CHECKLIST
        do n = 1, size(mesh%nset) ! Check that the NODE values are in-range.
          if (minval(mesh%nset(n)%node) < 1) exit CHECKLIST
          if (maxval(mesh%nset(n)%node) > mesh%num_node) exit CHECKLIST
        end do
      else
        if (mesh%num_nset /= 0) exit CHECKLIST
      end if

      !! If present, check the side set info.
      if (associated(mesh%sset)) then
        if (mesh%num_sset /= size(mesh%sset)) exit CHECKLIST
        if (.not.all(defined(mesh%sset))) exit CHECKLIST
        if (.not.unique(mesh%sset%ID)) exit CHECKLIST
        do n = 1, size(mesh%sset) ! Check that the ELEM and FACE values are in-range.
          if (minval(mesh%sset(n)%elem) < 1) exit CHECKLIST
          if (maxval(mesh%sset(n)%elem) > mesh%num_elem) exit CHECKLIST
          if (minval(mesh%sset(n)%face) < 1) exit CHECKLIST
          !if (maxval(mesh%sset(n)%face) > ???) exit CHECKLIST
        end do
      else
        if (mesh%num_sset /= 0) exit CHECKLIST
      end if

      defined_exodus_mesh = .true.
      exit

    end do CHECKLIST

  end function defined_exodus_mesh

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Specific procedures for the generic DEFINED
 !!
 !! These functions return the value true if the derived type argument is
 !! defined; otherwise they return the value false.  Defined means that the
 !! variable has data and that the data is superficially consistent.  The
 !! values of the array components themselves are not examined.
 !!
 !! These functions may take an array as an argument (returning a conformable
 !! logical result), however each array element examined independently of the
 !! others.
 !!

  elemental logical function defined_elem_blk (eb)

    type(elem_blk), intent(in) :: eb

    CHECKLIST: do
      defined_elem_blk = .false.
      if (eb%ID <= 0) exit
      if (.not.associated(eb%connect)) exit
      if (eb%num_elem /= size(eb%connect,dim=2)) exit
      defined_elem_blk = .true.
      exit
    end do CHECKLIST

  end function defined_elem_blk

  elemental logical function defined_node_set (ns)

    type(node_set), intent(in) :: ns

    CHECKLIST: do
      defined_node_set = .false.
      if (ns%ID <= 0) exit
      if (.not.associated(ns%node)) exit
      if (ns%num_node /= size(ns%node)) exit
      defined_node_set = .true.
      exit
    end do CHECKLIST

  end function defined_node_set

  elemental logical function defined_side_set (ss)

    type(side_set), intent(in) :: ss

    CHECKLIST: do
      defined_side_set = .false.
      if (ss%ID <= 0) exit
      if (.not.associated(ss%elem)) exit
      if (.not.associated(ss%face)) exit
      if (size(ss%elem) /= size(ss%face)) exit
      if (ss%num_side /= size(ss%elem)) exit
      defined_side_set = .true.
      exit
    end do CHECKLIST

  end function defined_side_set

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! UNIQUE
 !!
 !! This auxillary function returns the value true if the integer values in
 !! the vector LIST are unique; otherwise it returns the value false.
 !!
 !! A naive O(n^2) algorithm is used; please don't use it on long lists!
 !!

  pure logical function unique (list)

    integer, intent(in) :: list(:)

    integer :: i, j

    unique = .false.
    do i = 1, size(list)-1
      do j = i+1, size(list)
        if (list(i) == list(j)) return
      end do
    end do
    unique = .true.

  end function unique

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! EQ_* and NE_* functions for defined operators .EQ. and .NE.
 !!
 !! These functions test the equality and inequality of the derived data types
 !! defined by this module.
 !!

  logical function eq_exodus_mesh (a, b)
    type(exodus_mesh), intent(in) :: a, b
    eq_exodus_mesh = .false.
    if (a%num_dim  /= b%num_dim)  return
    if (a%num_node /= b%num_node) return
    if (a%num_elem /= b%num_elem) return
    if (a%num_eblk /= b%num_eblk) return
    if (a%num_nset /= b%num_nset) return
    if (a%num_sset /= b%num_sset) return
    if (associated(a%eblk) .neqv. associated(b%eblk)) return
    if (associated(a%eblk)) then
      if (size(a%eblk) /= size(b%eblk)) return
      if (any(a%eblk /= b%eblk)) return
    end if
    if (associated(a%nset) .neqv. associated(b%nset)) return
    if (associated(a%nset)) then
      if (size(a%nset) /= size(b%nset)) return
      if (any(a%nset /= b%nset)) return
    end if
    if (associated(a%sset) .neqv. associated(b%sset)) return
    if (associated(a%sset)) then
      if (size(a%sset) /= size(b%sset)) return
      if (any(a%sset /= b%sset)) return
    end if
    if (associated(a%coord) .neqv. associated(b%coord)) return
    if (associated(a%coord)) then
      if (any(shape(a%coord) /= shape(b%coord))) return
      if (any(a%coord /= b%coord)) return
    end if
    eq_exodus_mesh = .true.
  end function eq_exodus_mesh

  logical function ne_exodus_mesh (a, b)
    type(exodus_mesh), intent(in) :: a, b
    ne_exodus_mesh = .not. (a == b)
  end function ne_exodus_mesh

  elemental logical function eq_elem_blk (a, b)
    type(elem_blk), intent(in) :: a, b
    eq_elem_blk = .false.
    if (a%ID /= b%ID) return
    if (a%num_elem /= b%num_elem) return
    if (a%elem_type /= b%elem_type) return
    if (associated(a%connect) .neqv. associated(b%connect)) return
    if (associated(a%connect)) then
      if (any(shape(a%connect) /= shape(b%connect))) return
      if (any(a%connect /= b%connect)) return
    end if
    eq_elem_blk = .true.
  end function eq_elem_blk

  elemental logical function ne_elem_blk (a, b)
    type(elem_blk), intent(in) :: a, b
    ne_elem_blk = .not. (a == b)
  end function ne_elem_blk

  elemental logical function eq_node_set (a, b)
    type(node_set), intent(in) :: a, b
    eq_node_set = .false.
    if (a%ID /= b%ID) return
    if (a%num_node /= b%num_node) return
    if (associated(a%node) .neqv. associated(b%node)) return
    if (associated(a%node)) then
      if (size(a%node) /= size(b%node)) return
      if (any(a%node /= b%node)) return
    end if
    eq_node_set = .true.
  end function eq_node_set

  elemental logical function ne_node_set (a, b)
    type(node_set), intent(in) :: a, b
    ne_node_set = .not. (a == b)
  end function ne_node_set

  elemental logical function eq_side_set (a, b)
    type(side_set), intent(in) :: a, b
    eq_side_set = .false.
    if (a%ID /= b%ID) return
    if (a%num_side /= b%num_side) return
    if (associated(a%elem) .neqv. associated(b%elem)) return
    if (associated(a%elem)) then
      if (size(a%elem) /= size(b%elem)) return
      if (any(a%elem /= b%elem)) return
    end if
    if (associated(a%face) .neqv. associated(b%face)) return
    if (associated(a%face)) then
      if (size(a%face) /= size(b%face)) return
      if (any(a%face /= b%face)) return
    end if
    eq_side_set = .true.
  end function eq_side_set

  elemental logical function ne_side_set (a, b)
    type(side_set), intent(in) :: a, b
    ne_side_set = .not. (a == b)
  end function ne_side_set

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! Specific procedures for the generic DUMP_TYPE
 !!
 !! Testing/debugging routines that print the contents of the derived data
 !! types.  These are not intended for general use.
 !!

  subroutine dump_exodus_mesh (mesh, lun)

    type(exodus_mesh), intent(in) :: mesh
    integer, intent(in) :: lun

    integer :: j

    write(unit=lun,fmt='(a)') 'EXODUS_MESH('
    write(unit=lun,fmt='(t4,a,i6)') 'NUM_DIM= ', mesh%num_dim
    write(unit=lun,fmt='(t4,a,i6)') 'NUM_NODE=', mesh%num_node
    write(unit=lun,fmt='(t4,a,i6)') 'NUM_ELEM=', mesh%num_elem
    write(unit=lun,fmt='(t4,a,i6)') 'NUM_EBLK=', mesh%num_eblk
    write(unit=lun,fmt='(t4,a,i6)') 'NUM_NSET=', mesh%num_nset
    write(unit=lun,fmt='(t4,a,i6)') 'NUM_SSET=', mesh%num_sset
    if (associated(mesh%eblk)) then
      write(unit=lun,fmt='(t4,a)') 'EBLK= ...'
      do j = 1, size(mesh%eblk)
        call dump_elem_blk (mesh%eblk(j), lun)
      end do
    else
      write(unit=lun,fmt='(t4,a)') 'EBLK => NULL()'
    end if
    if (associated(mesh%nset)) then
      write(unit=lun,fmt='(t4,a)') 'NSET= ...'
      do j = 1, size(mesh%nset)
        call dump_node_set (mesh%nset(j), lun)
      end do
    else
      write(unit=lun,fmt='(t4,a)') 'NSET => NULL()'
    end if
    if (associated(mesh%sset)) then
      write(unit=lun,fmt='(t4,a)') 'SSET= ...'
      do j = 1, size(mesh%sset)
        call dump_side_set (mesh%sset(j), lun)
      end do
    else
      write(unit=lun,fmt='(t4,a)') 'SSET => NULL()'
    end if
    if (associated(mesh%coord)) then
      write(unit=lun,fmt='(t4,a,(t12,3(1x,es22.14)))') 'COORD=', mesh%coord(:,1)
      do j = 2, size(mesh%coord,dim=2)
        write(unit=lun,fmt='((t12,3(1x,es22.14)))') mesh%coord(:,j)
      end do
    else
      write(unit=lun,fmt='(t4,a)') 'COORD => NULL()'
    end if

  end subroutine dump_exodus_mesh

  subroutine dump_node_set (ns, lun)

    type(node_set), intent(in) :: ns
    integer, intent(in) :: lun

    write(unit=lun,fmt='(a)') 'NODE_SET('
    write(unit=lun,fmt='(t4,a,i6)') 'ID=', ns%ID
    write(unit=lun,fmt='(t4,a,i6)') 'NUM_NODE=', ns%num_node
    if (associated(ns%node)) then
      write(unit=lun,fmt='(t4,a,(t10,10(1x,i6)))') 'NODE=', ns%node
    else
      write(unit=lun,fmt='(t4,a)') 'NODE => NULL()'
    end if
    write(unit=lun,fmt='(t4,a)') ')'

  end subroutine dump_node_set

  subroutine dump_side_set (ss, lun)

    type(side_set), intent(in) :: ss
    integer, intent(in) :: lun

    write(unit=lun,fmt='(a)') 'SIDE_SET('
    write(unit=lun,fmt='(t4,a,i6)') 'ID=', ss%ID
    write(unit=lun,fmt='(t4,a,i6)') 'NUM_SIDE=', ss%num_side
    if (associated(ss%elem)) then
      write(unit=lun,fmt='(t4,a,(t10,10(1x,i6)))') 'ELEM=', ss%elem
    else
      write(unit=lun,fmt='(t4,a)') 'ELEM => NULL()'
    end if
    if (associated(ss%face)) then
      write(unit=lun,fmt='(t4,a,(t10,10(1x,i6)))') 'FACE=', ss%face
    else
      write(unit=lun,fmt='(t4,a)') 'FACE => NULL()'
    end if
    write(unit=lun,fmt='(t4,a)') ')'

  end subroutine dump_side_set

  subroutine dump_elem_blk (eb, lun)

    type(elem_blk), intent(in) :: eb
    integer, intent(in) :: lun

    integer :: j

    write(unit=lun,fmt='(a)') 'ELEMENT_BLOCK('
    write(unit=lun,fmt='(t4,a,i6)') 'ID=', eb%ID
    write(unit=lun,fmt='(t4,a,i6)') 'NUM_ELEM=', eb%num_elem
    write(unit=lun,fmt='(t4,a)')  'ELEM_TYPE= "' // trim(eb%elem_type) // '"'
    if (associated(eb%connect)) then
      write(unit=lun,fmt='(t4,a,(t12,8(1x,i6)))') 'CONNECT=', eb%connect(:,1)
      do j = 2, size(eb%connect,dim=2)
          write(unit=lun,fmt='((t12,8(1x,i6)))') eb%connect(:,j)
      end do
    else
      write(unit=lun,fmt='(t4,a)') 'CONNECT => NULL()'
    end if
    write(unit=lun,fmt='(t4,a)') ')'

  end subroutine dump_elem_blk

end module exodus_mesh_type
