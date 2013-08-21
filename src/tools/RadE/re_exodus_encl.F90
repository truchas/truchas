!!
!! RE_EXODUS_ENCL
!!
!! This provides a method for generating an enclosure by extracting a surface
!! from an Exodus mesh file.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 4 April 2008
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL READ_ENCLOSURE_NAMELIST (LUN, SPEC) reads the first occurrence of an
!!    ENCLOSURE namelist from the file opened on unit LUN.  The values read
!!    (and defaults for any others) are returned in the enclosure specification
!!    structure SPEC.  It is an error if the namelist is not found; this and
!!    any IO errors are handled by the subroutine.  The namelist values are
!!    also checked for correctness.  Execution of the program is gracefully
!!    terminated if any errors are encountered.  This is a collective
!!    procedure, but only for the purposes of error handling; input occurs on
!!    process rank 1 and the returned SPEC is only valid on process rank 1.
!!
!!    The ENCLOSURE namelist contains the following variables:
!!      name ~ user-supplied name for the enclosure (req)
!!      mesh_file ~ path of the source Exodus mesh file (req)
!!      coord_scale_factor ~ scale factor to apply to the surface coordinates
!!      ignore_block_IDs ~ list of mesh element blocks to ignore (opt)
!!      side_set_IDs ~ list of mesh side sets defining the surface (req)
!!      symmetries ~ list of up to 3 symmetry operations (opt):
!!                   Mirror<a>, a = X, Y, Z; e.g. MirrorZ
!!                   Rot<a><n>, a = X, Y, Z, n integer; e.g., RotZ3, RotX16
!!      displace_side_set_IDs ~ list of surface side sets to displace;
!!      displacement ~ the constant (x,y,z) displacement amount.
!!
!!  CALL GENERATE_ENCL (SPEC, E) generates the enclosure E using the enclosure
!!    specification SPEC.  This is a collective procedure.  SPEC is only
!!    relevant on process rank 1, and the generation of the enclosure takes
!!    place there too and then replicated across all processes.
!!
!!  CALL DESTROY (SPEC) deallocates any allocated storage associated with
!!    the enclosure specification SPEC.
!!
!! IMPLEMENTATION NOTES
!!
!!  (1) The extraction of the enclosure surface from the Exodus mesh is done
!!  in serial, but is very quick.  Ultimately the enclosure is to be replicated
!!  on all processes, so one could imagine for simplicity just having each
!!  process do the extraction itself.  However this would require each process
!!  to read the source mesh, or to replicate it, having one process read it.
!!  Since the mesh is potential very large -- much larger than the output
!!  enclosure mesh -- this isn't a good idea.  Better to extract the enclosure
!!  on one process, leaving the others idle, and then replicate the much smaller
!!  results.
!!

#include "f90_assert.fpp"

module re_exodus_encl

  use kinds, only: r8
  use scl
  use re_utilities
  implicit none
  private

  public :: encl_spec, read_enclosure_namelist, generate_encl, destroy

  integer, parameter :: MAX_NAME_LEN = 32, MAX_FILE_LEN = 256, MAX_IDS = 128

  type :: encl_spec
    character(len=MAX_NAME_LEN) :: name
    character(len=MAX_FILE_LEN) :: mesh_file
    real(r8) :: scale
    integer, pointer :: ssid(:) => null()
    integer, pointer :: ebid(:) => null()
    logical :: mirror(3)
    integer :: rot_axis, num_rot
    integer, pointer :: disp_ssid(:) => null()
    real(r8) :: disp(3)
  end type encl_spec

  interface destroy
    module procedure destroy_encl_spec
  end interface

contains

  subroutine read_enclosure_namelist (lun, spec)

    use string_utilities
    use input_utilities

    integer, intent(in) :: lun
    type(encl_spec), intent(out) :: spec

    !! The ENCLOSURE namelist variables; user visible.
    character(len=MAX_NAME_LEN) :: name
    character(len=MAX_FILE_LEN) :: mesh_file
    character(len=7) :: symmetries(3)
    integer :: side_set_ids(MAX_IDS), ignore_block_ids(MAX_IDS), displace_side_set_ids(MAX_IDS)
    real(r8) :: coord_scale_factor, displacement(3)
    namelist /enclosure/ name, mesh_file, side_set_ids, ignore_block_ids, symmetries, &
        coord_scale_factor, displace_side_set_ids, displacement

    integer :: j, n, ios, stat, rot_axis, num_rot
    logical :: is_IOP, found, file_exists, mirror(3)
    character(len=len(symmetries)) :: sym

    is_IOP = (scl_rank()==1)  ! process rank 1 does the reading

    !! Seek to the first instance of the ENCLOSURE namelist.
    if (is_IOP) then
      rewind (lun)
      call seek_to_namelist (lun, 'ENCLOSURE', found, iostat=ios)
    end if
    call scl_bcast (ios)
    if (ios /= 0) call re_halt ('error reading file connected to unit ' // &
                                i_to_c(lun) // ': iostat=' // i_to_c(ios))

    !! This is a required namelist.
    call scl_bcast (found)
    if (.not.found) call re_halt ('ENCLOSURE namelist not found')

    !! Read the namelist, assigning default values first.
    call re_info ('Reading ENCLOSURE namelist ...')
    if (is_IOP) then
      name = NULL_C
      mesh_file = NULL_C
      side_set_ids = NULL_I
      ignore_block_ids = NULL_I
      symmetries = NULL_C
      coord_scale_factor = 1.0_r8
      displace_side_set_ids = NULL_I
      displacement = 0.0_r8
      read(lun,nml=enclosure,iostat=ios)
    end if
    call scl_bcast (ios)
    if (ios /= 0) call re_halt ('Error reading ENCLOSURE namelist: iostat=' // i_to_c(ios))

    if (is_IOP) then
      stat = 0

      !! Check the user-supplied NAME for the namelist.
      if (name == NULL_C) call data_err ('NAME must be assigned a value')
      name = raise_case(name)

      !! Check the MESH_FILE path.
      if (mesh_file == NULL_C) call data_err ('MESH_FILE must be assigned a value')
      inquire(file=mesh_file,exist=file_exists)
      if (.not.file_exists) call data_err ('no such MESH_FILE: ' // trim(mesh_file))

      !! Check for a non-empty SIDE_SET_IDS.
      if (count(side_set_ids /= NULL_I) == 0) &
          call data_err ('SIDE_SET_IDS must contain at least one value')

      !! Unpack SYMMETRIES and check the values.
      mirror = .false.
      rot_axis = 0
      num_rot  = 0
      do j = 1, size(symmetries)
        if (symmetries(j) == NULL_C) cycle
        sym = raise_case(symmetries(j))
        if (sym(1:6) == 'MIRROR' .and. len_trim(sym) == 7) then
          n = scan('XYZ',set=sym(7:7))
          if (n == 0) call data_err ('unknown symmetry value: '//trim(symmetries(j)))
          mirror(n) = .true.
        else if (sym(1:3) == 'ROT') then
          n = scan('XYZ',set=sym(4:4))
          if (n == 0) call data_err ('unknown symmetry value: '//trim(symmetries(j)))
          rot_axis = n
          read(unit=sym(5:),fmt=*,iostat=ios) num_rot
          if (ios /= 0) then
            call data_err ('unknown symmetry value: '//trim(symmetries(j)))
          else if (num_rot < 2) then
            call data_err ('number of rotations must be > 1: ' // trim(symmetries(j)))
          end if
        else
          call data_err ('unknown symmetry value: ' // trim(symmetries(j)))
        end if
      end do

      !! Check COORD_SCALE_FACTOR.
      if (coord_scale_factor <= 0.0_r8) call data_err ('COORD_SCALE_FACTOR must be > 0.0')
    end if

    call scl_bcast (stat)
    if (stat /= 0) call re_halt ('errors found in ENCLOSURE namelist variables')

    !! Everything checks out; stuff the values into the return data structure.
    if (is_IOP) then
      spec%name = name
      spec%mesh_file = mesh_file
      spec%scale = coord_scale_factor
      spec%ssid => ptr_to_packed(side_set_ids, mask=(side_set_ids /= NULL_I))
      spec%ebid => ptr_to_packed(ignore_block_ids, mask=(ignore_block_ids /= NULL_I))
      spec%mirror = mirror
      spec%rot_axis = rot_axis
      spec%num_rot  = num_rot
      spec%disp_ssid => ptr_to_packed(displace_side_set_ids, mask=(displace_side_set_ids /= NULL_I))
      spec%disp = displacement
    end if

  contains

    function ptr_to_packed (source, mask) result (ptr)
      integer, intent(in) :: source(:)
      logical, intent(in) :: mask(:)
      integer, pointer :: ptr(:)
      ASSERT( size(source) == size(mask) )
      allocate(ptr(count(mask)))
      ptr = pack(source, mask)
    end function ptr_to_packed

    subroutine data_err (errmsg)
      character(len=*), intent(in) :: errmsg
      stat = 1
      call re_info ('  ERROR: ' // errmsg)
    end subroutine data_err

  end subroutine read_enclosure_namelist

  subroutine destroy_encl_spec (spec)
    type(encl_spec), intent(inout) :: spec
    if (associated(spec%ssid)) deallocate(spec%ssid)
    if (associated(spec%ebid)) deallocate(spec%ebid)
  end subroutine destroy_encl_spec

  subroutine generate_encl (spec, e)
    use exodus
    use re_encl_type
    type(encl_spec), intent(in)  :: spec
    type(encl), intent(out) :: e
    type(exodus_mesh) :: mesh
    integer :: stat
    character(len=128) :: errmsg
    if (scl_rank() == 1) then
      call read_exodus_mesh (spec%mesh_file, mesh)
      call extract_surface_from_exodus (mesh, spec%ssid, spec%ebid, e, stat, errmsg)
      call destroy (mesh)
      if (stat == 0) then
        e%name     = spec%name
        e%mirror   = spec%mirror
        e%rot_axis = spec%rot_axis
        e%num_rot  = spec%num_rot
        if (spec%scale /= 1.0_r8) e%x = spec%scale * e%x
        call displace_surfaces (e, spec%disp, spec%disp_ssid, stat, errmsg)
      end if
    end if
    call scl_bcast (stat)
    if (stat /= 0) call re_halt ('GENERATE_ENCL: ' // trim(errmsg))
    call bcast_encl (e)
  end subroutine generate_encl

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! CALL EXTRACT_SURFACE_FROM_EXODUS
 !!
 !! This subroutine does all the serious work of extracting a surface mesh
 !! from an Exodus mesh.  The surface is specified as a list of side set IDs.
 !! There are a couple of restrictions: all surface faces must be on the
 !! boundary of the mesh, and each surface face belongs to exactly one of the
 !! specified side sets.  The code checks that these conditions are satisfied,
 !! which leads to some of the complexity of the code.  Certain element blocks
 !! of the mesh may be ignored by specifying their IDs.  This changes the
 !! effective mesh boundary, and is useful when void regions contained in the
 !! volume of a radiation enclosure have been meshed.
 !!
 !! The extracted surface mesh is a typical unstructured 2D mesh (with 3D
 !! coordinates) consisting of tri and quad faces, all consistently oriented
 !! so that the top of the surface points into the enclosure.  Exodus meshes
 !! consisting of tet, hex and wedge elements can be handled.  Each face is
 !! tagged with the (unique) side set ID that identified it.  In addition,
 !! each face is tagged with the (unique) corresponding cell/side index pair
 !! of the Exodus mesh.  This enables the surface face to be associated with
 !! the corresponding face of the volume mesh at a later time.
 !!
 !! Considerable care is taken to ensure that the extraction scales well with
 !! increasing mesh size, by restricting attention to the part of the mesh near
 !! the specified surface, and through the use of hash tables for searching.
 !!

  subroutine extract_surface_from_exodus (mesh, ssid, ebid, surf, stat, errmsg)

    use exodus
    use string_utilities, only: i_to_c
    use re_encl_type
    use hashing

    type(exodus_mesh), intent(in)  :: mesh
    integer, intent(in) :: ssid(:)  ! IDs of side sets describing the surface
    integer, intent(in) :: ebid(:)  ! IDs of element blocks to omit from the mesh
    type(encl), intent(out) :: surf  ! mesh of the extracted surface
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: i, j, k, n, n1, n2, b, l, bitmask, tsize
    integer, pointer :: list(:), side_sig(:)
    logical, allocatable :: surf_sset(:), omit_elem(:)
    integer, allocatable :: surf_side(:), nhbr_side(:)  ! bitmask arrays
    integer, allocatable :: nmap(:), gmap(:)

    !! Hash table structure
    type :: side_row
      integer :: j=0, k=0 ! key values
      integer :: ssn=0    ! data value
      integer, pointer :: list(:) => null() ! more data values
    end type
    type :: side_table
      integer :: size=0, nentry=0
      type(hash_param) :: hpar
      type(side_row), pointer :: row(:) => null()
    end type
    type(side_table) :: stab

    !! Another hash table structure
    type :: nhbr_row
      integer, pointer :: list(:) => null() ! key values, no data
    end type
    type :: nhbr_table
      integer :: size=0, nentry=0
      type(hash_param) :: hpar
      type(nhbr_row), pointer :: row(:) => null()
    end type
    type(nhbr_table) :: ntab

    integer, parameter :: MAX_SIDE = 6  ! works for tet, hex and wedge elements.

    integer, target, save :: TETRA4_SIDE_SIG(4)
    data TETRA4_SIDE_SIG/b'1011', b'1110', b'1101', b'0111'/

    integer, target, save :: WEDGE6_SIDE_SIG(5)
    data WEDGE6_SIDE_SIG/b'011011', b'110110', b'101101', b'000111', b'111000'/

    integer, target, save :: HEX8_SIDE_SIG(6)
    data HEX8_SIDE_SIG/b'00110011', b'01100110', b'11001100', b'10011001', b'00001111', b'11110000'/

    !! Generate the mask array SURF_SSET that tags the surface side sets.
    allocate(surf_sset(mesh%num_sset))
    surf_sset = .false.
    do i = 1, size(ssid)
      n = mesh%num_sset
      do while (n > 0)
        if (mesh%sset(n)%ID == ssid(i)) exit
        n = n - 1
      end do
      if (n == 0) then
        stat = -1
        errmsg = 'no side set with ID ' // i_to_c(ssid(i))
        return
      end if
      surf_sset(n) = .true.
    end do

    !! The side set IDs become the surface face group IDs.
    allocate(surf%group_id_list(count(surf_sset)), gmap(mesh%num_sset))
    surf%group_id_list = pack(mesh%sset%ID, mask=surf_sset)
    gmap = unpack((/(j, j=1,size(gmap))/), mask=surf_sset, field=0)

    !! Generate the mask arrays OMIT_ELEM.
    allocate(omit_elem(mesh%num_elem))
    omit_elem = .false.
    do i = 1, size(ebid)
      n = mesh%num_eblk
      do while (n > 0)
        if (mesh%eblk(n)%ID == ebid(i)) exit
        n = n - 1
      end do
      if (n == 0) then
        stat = -1
        errmsg = 'no element block with ID ' // i_to_c(ebid(i))
        return
      end if
      n2 = sum(mesh%eblk(1:n)%num_elem)
      n1 = n2 - mesh%eblk(n)%num_elem + 1
      omit_elem(n1:n2) = .true.
    end do

    !! Store the side information temporarily in a hash table.  This enables
    !! us to skip duplicate sides in a side set and detect conflicting sides
    !! from different side sets.  The element-based bitmask array SURF_SIDE
    !! keeps track of the sides stored in the table.  In the process we mark
    !! the surface nodes with a nonzero value in NMAP.

    tsize = 1.25 * sum(mesh%sset%num_side,mask=surf_sset)  ! ensures 20% free space
    call create_side_table (stab, tsize, mesh%num_node)

    allocate(surf_side(mesh%num_elem), nmap(mesh%num_node))
    nmap = 0
    surf_side = 0

    n1 = 0  ! surface face count
    n2 = 0  ! surface face node count
    do n = 1, mesh%num_sset
      if (.not.surf_sset(n)) cycle
      do i = 1, mesh%sset(n)%num_side
        j = mesh%sset(n)%elem(i)
        k = mesh%sset(n)%face(i)
        if (omit_elem(j)) cycle
        list => side_node_list(mesh, j, k)
        INSIST( associated(list) )
        call store_side (stab, j, k, n, list, stat)
        select case (stat)
        case (0)
          n1 = n1 + 1
          n2 = n2 + size(list)
          surf_side(j) = ibset(surf_side(j),k-1)
          nmap(list) = 1
        case (1) ! duplicate side; okay but might want to warn
          ! set some sort of warning flag
        case default  ! side exists from different side set; fatal
          stat = -1
          errmsg = 'face belongs to multiple side sets'
          return
        end select
      end do
    end do

    !! Surface nodes have been tagged with a nonzero value in NMAP.
    !! Now modify NMAP to give the node-to-surface node mapping.
    n = 0 ! surface node index
    do j = 1, size(nmap)
      if (nmap(j) == 0) cycle
      n = n + 1
      nmap(j) = n
    end do
    surf%nnode = n

    !! Create the surface node coordinate array.
    allocate(surf%x(size(mesh%coord,dim=1),surf%nnode))
    do j = 1, mesh%num_node
      if (nmap(j) > 0) surf%x(:,nmap(j)) = mesh%coord(:,j)
    end do

    !! Identify all mesh element sides that could possibly be a neighbor of
    !! one of the surface cell sides specified above and store the results
    !! in the element-based bitmask array NHBR_SIDE.  These would be sides
    !! whose nodes are all surface nodes.

    allocate(nhbr_side(mesh%num_elem))
    nhbr_side = 0
    j = 0 ! element index
    tsize = 0 ! neighbor count
    do b = 1, mesh%num_eblk

      select case (mesh%eblk(b)%elem_type)
      case ('TETRA', 'TETRA4')
        side_sig => TETRA4_SIDE_SIG
      case ('WEDGE', 'WEDGE6')
        side_sig => WEDGE6_SIDE_SIG
      case ('HEX', 'HEX8')
        side_sig => HEX8_SIDE_SIG
      case default
        stat = -1
        errmsg = 'unable to handle element type: ' // trim(mesh%eblk(b)%elem_type)
        return
      end select

      do l = 1, mesh%eblk(b)%num_elem
        j = j + 1
        if (omit_elem(j)) cycle

        !! Mark the element vertices that are surface nodes.
        list => mesh%eblk(b)%connect(:,l)
        bitmask = 0
        do k = 1, size(list)
          if (nmap(list(k))>0) bitmask = ibset(bitmask, k-1)
        end do

        !! Mark sides having all surface nodes.
        if (bitmask == 0) cycle ! nothing to see on this element
        do k = 1, size(side_sig)
          if (iand(bitmask, side_sig(k)) == bitmask) then
            nhbr_side(j) = ibset(nhbr_side(j), k-1)
            tsize = tsize + 1
          end if
        end do
      end do
    end do

    !! Put the hypothetical neighbor face to each of neighbor sides identified
    !! above into a hash table.  If the face of a surface side appears in this
    !! table it isn't a boundary face.  The reason for not simply including
    !! this in the above loop is because we can't reliably and reasonably
    !! estimate the size of the needed table.

    tsize = max(4.0,1.25 * tsize)  ! at least 20% free space in hash table
    call create_nhbr_table (ntab, tsize, mesh%num_node)
    do j = 1, mesh%num_elem
      if (nhbr_side(j) == 0) cycle
      do k = 1, MAX_SIDE
        if (btest(nhbr_side(j), k-1)) then
          list => side_node_list(mesh, j, k, reverse=.true.)
          ASSERT( associated(list) )
          call add_nhbr (ntab, list)
        end if
      end do
    end do

    !! Finally we can create the surface face connection arrays.
    !! We pull the data out of the side table

    surf%nface = n1
    allocate(surf%xface(n1+1), surf%fnode(n2), surf%gnum(n1))
    allocate(surf%src_elem(n1), surf%src_side(n1))

    n = 0 ! surface face index
    surf%xface(1) = 1
    do j = 1, mesh%num_elem
      if (surf_side(j) == 0) cycle
      do k = 1, MAX_SIDE
        if (btest(surf_side(j),k-1)) then
          n = n + 1
          surf%src_elem(n) = j
          surf%src_side(n) = k
          call retrieve_side (stab, j, k, i, list)
          surf%xface(n+1) = surf%xface(n) + size(list)
          surf%fnode(surf%xface(n):surf%xface(n+1)-1) = nmap(list)
          surf%gnum(n) = gmap(i)
          if (nhbr_exists(ntab, list)) then
            stat = -1
            errmsg = 'surface face not on the mesh boundary'
            return
          end if
        end if
      end do
    end do
    ASSERT( n == surf%nface )
    ASSERT( surf%xface(surf%nface+1) == size(surf%fnode)+1 )

    !! Clean up.
    call destroy_side_table (stab)
    call destroy_nhbr_table (ntab)
    deallocate(surf_sset, omit_elem, nhbr_side, surf_side, nmap, gmap)

    stat = 0

  contains

    !! These auxillary routines implement the methods of the surface side
    !! hash table.  Though each is called in only one location and could
    !! be inlined, they've been encapsulated for clarity.

    subroutine create_side_table (table, size, key_max)
      type(side_table), intent(out) :: table
      integer, intent(in) :: size, key_max
      table%size = size
      call initialize_hash_param (table%hpar, hsize=table%size, kmax=key_max)
      allocate(table%row(0:table%size-1))
    end subroutine

    subroutine destroy_side_table (table)
      type(side_table), intent(inout) :: table
      if (associated(table%row)) then
        do j = lbound(table%row,dim=1), ubound(table%row,dim=1)
          if (associated(table%row(j)%list)) deallocate(table%row(j)%list)
        end do
        deallocate(table%row)
      end if
    end subroutine

    subroutine store_side (table, j, k, ssn, list, status)
      type(side_table), intent(inout) :: table
      integer, intent(in) :: j, k, ssn
      integer, pointer :: list(:)
      integer, intent(out) :: status
      integer :: n, incr
      call hash (table%hpar, (/j, k/), n, incr)
      do while (table%row(n)%j /= 0)
        if (j == table%row(n)%j) then
          if (k == table%row(n)%k) exit
        end if
        n = n - incr
        if (n < 0) n = n + table%size
      end do
      if (table%row(n)%j == 0) then
        !! Store new side into empty row.
        INSIST( table%nentry < table%size - 1 )
        table%row(n)%j = j
        table%row(n)%k = k
        table%row(n)%ssn = ssn
        table%row(n)%list => list
        status = 0
      else if (table%row(n)%ssn == ssn) then
        !! Side already stored but with consistent side set ID -- okay, but skip it.
        status = 1
      else
        !! Side already stored and with conflicting side set ID -- not good.
        status = -1
      end if
    end subroutine

    subroutine retrieve_side (table, j, k, ssn, list)
      type(side_table), intent(in) :: table
      integer, intent(in) :: j, k
      integer, intent(out) :: ssn
      integer, pointer :: list(:)
      integer :: n, incr
      call hash (table%hpar, (/j, k/), n, incr)
      do while (table%row(n)%j /= 0)
        if (j == table%row(n)%j) then
          if (k == table%row(n)%k) exit
        end if
        n = n - incr
        if (n < 0) n = n + table%size
      end do
      INSIST( table%row(n)%j /= 0 )
      ssn = table%row(n)%ssn
      list => table%row(n)%list
    end subroutine

    !! These auxillary routines implement the methods of the neighbor side
    !! hash table.  Though each is called in only one location and could
    !! be inlined, they've been encapsulated for clarity.

    subroutine create_nhbr_table (table, size, key_max)
      type(nhbr_table), intent(out) :: table
      integer, intent(in) :: size, key_max
      table%size = size
      call initialize_hash_param (table%hpar, hsize=table%size, kmax=key_max)
      allocate(table%row(0:table%size-1))
    end subroutine

    subroutine destroy_nhbr_table (table)
      type(nhbr_table), intent(inout) :: table
      if (associated(table%row)) then
        do j = lbound(table%row,dim=1), ubound(table%row,dim=1)
          if (associated(table%row(j)%list)) deallocate(table%row(j)%list)
        end do
        deallocate(table%row)
      end if
    end subroutine

    subroutine add_nhbr (table, list)
      type(nhbr_table), intent(inout) :: table
      integer, pointer :: list(:)
      integer :: n, inc
      call hash (table%hpar, list, n, inc)
      do while (associated(table%row(n)%list))
        n = n - inc
        if (n < 0) n = n + table%size
      end do
      INSIST( table%nentry < table%size - 1 )
      table%row(n)%list => list
    end subroutine

    logical function nhbr_exists (table, list)
      type(nhbr_table), intent(in) :: table
      integer, intent(in) :: list(:)
      integer :: n, inc
      call hash (table%hpar, list, n, inc)
      nhbr_exists = .true.
      do while (associated(table%row(n)%list))
        if (size(table%row(n)%list) == size(list)) then
          if (all(table%row(n)%list == list)) return
        end if
        n = n - inc
        if (n < 0) n = n + table%size
      end do
      nhbr_exists = .false.
    end function

  end subroutine extract_surface_from_exodus

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! DISPLACE_SURFACES
 !!
 !! Displaces the position of the of the given enclosure surfaces by a constant
 !! amount.  Those surfaces must be disconnected from the remaining enclosure
 !! surfaces.
 !!

  subroutine displace_surfaces (surf, disp, ssid, stat, errmsg)

    use re_encl_type
    use string_utilities, only: i_to_c

    type(encl), intent(inout) :: surf
    real(r8), intent(in) :: disp(:)
    integer, intent(in) :: ssid(:)
    integer, intent(out) :: stat
    character(len=*), intent(out) :: errmsg

    integer :: j
    logical, allocatable :: gmask(:), mask(:)

    stat = 0

    if (size(ssid) == 0) return ! nothing to do

    !! Tag the surface face groups (== side sets) that are to be displaced.
    allocate(gmask(size(surf%group_id_list)))
    do j = 1, size(surf%group_id_list)
      gmask(j) = any(surf%group_id_list(j) == ssid)
    end do

    if (.not.any(gmask)) then ! nothing to do
      call re_info ('Warning: directed to displace side sets but none found!')
      return
    end if

    !! Tag all nodes belonging to the displaced surface faces.
    allocate(mask(surf%nnode))
    mask = .false.
    do j = 1, surf%nface
      if (gmask(surf%gnum(j))) mask(surf%fnode(surf%xface(j):surf%xface(j+1)-1)) = .true.
    end do

    !! Check that no nodes belonging to undisplaced surface faces are tagged.
    do j = 1, surf%nface
      if (gmask(surf%gnum(j))) cycle
      if (any(mask(surf%fnode(surf%xface(j):surf%xface(j+1)-1)))) then
        stat = -1
        errmsg = 'undisplaced side set ' // i_to_c(surf%gnum(j)) // &
                 ' shares nodes with a displaced side set'
        return
      end if
    end do

    !! Displace the tagged nodes.
    do j = 1, surf%nnode
      if (mask(j)) surf%x(:,j) = surf%x(:,j) + disp
    end do

  end subroutine displace_surfaces

end module re_exodus_encl
