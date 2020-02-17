!!
!! RE_EXODUS_ENCL
!!
!! This provides a method for generating an enclosure by extracting a surface
!! from an Exodus mesh file.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 4 April 2008; revised February 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL READ_ENCLOSURE_NAMELIST(LUN, PARAMS) reads the first occurrence of
!!    an ENCLOSURE namelist from the file opened on unit LUN. The values read
!!    (and defaults for any others) are returned in the parameter list PARAMS.
!!    It is an error if the namelist is not found; this and any IO errors are
!!    handled by the subroutine.  The namelist values are also checked for
!!    correctness.  Execution of the program is gracefully terminated if any
!!    errors are encountered.  This is a collective procedure; input occurs on
!!    rank 1 but the returned PARAMS is replicated on all ranks.
!!
!!    The ENCLOSURE namelist contains the following variables:
!!      name ~ user-supplied name for the enclosure (req)
!!      mesh_file ~ path of the source Exodus mesh file (req)
!!      coord_scale_factor ~ scale factor to apply to the surface coordinates (opt)
!!      exodus_block_modulus - remap block IDs to their value modulo this modulus (opt)
!!      ignore_block_IDs ~ list of mesh element blocks to ignore (opt)
!!      side_set_IDs ~ list of mesh side sets defining the surface (req)
!!      symmetries ~ list of up to 3 symmetry operations (opt):
!!                   Mirror<a>, a = X, Y, Z; e.g. MirrorZ
!!                   Rot<a><n>, a = X, Y, Z, n integer; e.g., RotZ3, RotX16
!!      displace_block_IDs ~ list of element blocks to displace;
!!      displacements ~ a list of (x,y,z) displacements to apply (opt).
!!
!!  CALL INIT_ENCL(THIS, PARAMS, STAT, ERRMSG) initializes the ENCL object
!!    THIS as specified by the parameter list PARAMS. The object is replicated
!!    on all ranks, however PARAMS is only significant on rank 1. The integer
!!    STAT returns a non-zero value if an error occurs, and the deferred length
!!    allocatable character ERRMSG returns an explanatory message. The
!!    following parameters, with the meanings above, are used:
!!
!!      'name'
!!      'mesh-file'
!!      'coord-scale-factor'
!!      'exodus-block-modulus'
!!      'side-set-ids'
!!      'ignore-block-ids'
!!      'symmetries'
!!
!!  CALL INIT_ENCL_LIST(THIS, PARAMS, STAT, ERRMSG) initializes the ENCL_LIST
!!    object THIS as specified by the parameter list PARAMS.  The object is
!!    replicated on all ranks, however PARAMS is only significant on rank 1.
!!    The integer STAT returns a non-zero value if an error occurs, and the
!!    deferred length allocatable character ERRMSG returns an explanatory
!!    message. The following parameters are recognized in addtion to those
!!    above:
!!
!!      'displace-block-ids'
!!      'displacements'
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

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use scl
  use re_utilities
  use parameter_list_type
  implicit none
  private

  public :: read_enclosure_namelist, init_encl, init_encl_list

  integer, parameter :: MAX_NAME_LEN = 32, MAX_FILE_LEN = 256, MAX_IDS = 128, MAX_DISPL = 1000

contains

  subroutine read_enclosure_namelist(lun, params)

    use string_utilities, only: raise_case, i_to_c
    use input_utilities
    use toolpath_table, only: known_toolpath

    integer, intent(in) :: lun
    type(parameter_list), intent(out) :: params

    !! The ENCLOSURE namelist variables; user visible.
    character(MAX_NAME_LEN) :: name, displacement_toolpath
    character(MAX_FILE_LEN) :: mesh_file
    character(7) :: symmetries(3)
    integer :: side_set_ids(MAX_IDS), ignore_block_ids(MAX_IDS)
    integer :: displace_block_ids(MAX_IDS), exodus_block_modulus
    real(r8) :: coord_scale_factor, displacements(3,MAX_DISPL)
    namelist /enclosure/ name, mesh_file, coord_scale_factor, exodus_block_modulus, &
        side_set_ids, ignore_block_ids, symmetries, &
        displace_block_ids, displacements, displacement_toolpath

    integer :: j, n, ios, stat, rot_axis, num_rot, ndispl
    logical :: is_IOP, found, mirror(3)
    character(len=len(symmetries)) :: sym
    character(255) :: iom

    is_IOP = (scl_rank()==1)  ! process rank 1 does the reading

    call re_info ('Reading ENCLOSURE namelist ...')

    !! Locate the ENCLOSURE namelist (required)
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'ENCLOSURE', found, iostat=ios)
    end if
    call scl_bcast(ios)
    if (ios /= 0) call re_halt('error reading input file: iostat=' // i_to_c(ios))
    call scl_bcast(found)
    if (.not.found) call re_halt('ENCLOSURE namelist not found')

    !! Default values
    name = NULL_C
    mesh_file = NULL_C
    coord_scale_factor = 1.0_r8
    exodus_block_modulus = NULL_I
    side_set_ids = NULL_I
    ignore_block_ids = NULL_I
    symmetries = NULL_C
    displace_block_ids = NULL_I
    displacements = NULL_R
    displacement_toolpath = NULL_C

    if (is_IOP) read(lun,nml=enclosure,iostat=ios,iomsg=iom)
    call scl_bcast(ios)
    if (ios /= 0) call re_halt('error reading ENCLOSURE namelist: ' // trim(iom))

    !! Broadcast the namelist variables
    call scl_bcast(name)
    call scl_bcast(mesh_file)
    call scl_bcast(coord_scale_factor)
    call scl_bcast(exodus_block_modulus)
    call scl_bcast(side_set_ids)
    call scl_bcast(ignore_block_ids)
    call scl_bcast(symmetries)
    call scl_bcast(displace_block_ids)
    call scl_bcast(displacements)
    call scl_bcast(displacement_toolpath)

    !! Check the user-supplied NAME for the namelist.
    if (name == NULL_C) call re_halt('NAME not specified')
    call params%set('name', raise_case(trim(name)))

    !! Check the MESH_FILE path.
    if (mesh_file == NULL_C) call re_halt('MESH_FILE not specified')
    inquire(file=mesh_file,exist=found)
    if (.not.found) call re_halt('no such MESH_FILE: ' // trim(mesh_file))
    call params%set('mesh-file', trim(mesh_file))

    !! Check COORD_SCALE_FACTOR.
    if (coord_scale_factor <= 0.0_r8) call re_halt('COORD_SCALE_FACTOR must be > 0.0')
    call params%set('coord-scale-factor', coord_scale_factor)

    !! Check EXODUS_BLOCK_MODULUS.
    if (exodus_block_modulus /= NULL_I) then
      if (exodus_block_modulus < 0) call re_halt('EXODUS_BLOCK_MODULUS must be >= 0')
      call params%set('exodus-block-modulus', exodus_block_modulus)
    end if

    !! Check for a non-empty SIDE_SET_IDS.
    if (count(side_set_ids /= NULL_I) == 0) call re_halt('SIDE_SET_IDS not specified')
    call params%set('side-set-ids', pack(side_set_ids, mask=(side_set_ids /= NULL_I)))

    if (any(ignore_block_ids /= NULL_I)) &
    call params%set('ignore-block-ids', pack(ignore_block_ids, mask=(ignore_block_ids /= NULL_I)))

    !! Unpack SYMMETRIES and check the values.
    mirror = .false.
    rot_axis = 0
    num_rot  = 0
    do j = 1, size(symmetries)
      if (symmetries(j) == NULL_C) cycle
      sym = raise_case(symmetries(j))
      if (sym(1:6) == 'MIRROR' .and. len_trim(sym) == 7) then
        n = scan('XYZ',set=sym(7:7))
        if (n == 0) call re_halt('unknown symmetry value: ' // trim(symmetries(j)))
        mirror(n) = .true.
      else if (sym(1:3) == 'ROT') then
        n = scan('XYZ',set=sym(4:4))
        if (n == 0) call re_halt('unknown symmetry value: ' // trim(symmetries(j)))
        rot_axis = n
        read(unit=sym(5:),fmt=*,iostat=ios) num_rot
        if (ios /= 0) then
          call re_halt('unknown symmetry value: ' // trim(symmetries(j)))
        else if (num_rot < 2) then
          call re_halt('number of rotations must be > 1: ' // trim(symmetries(j)))
        end if
      else
        call re_halt('unknown symmetry value: ' // trim(symmetries(j)))
      end if
    end do
    call params%set('mirror', mirror)
    call params%set('rot-axis', rot_axis)
    call params%set('num-rot', num_rot)

    if (count(displace_block_ids /= NULL_I) > 0) then
      call params%set('displace-block-ids', &
          pack(displace_block_ids, mask=(displace_block_ids /= NULL_I)))
      !! Check DISPLACEMENTS
      do ndispl = size(displacements,dim=2), 1, -1
        if (any(displacements(:,ndispl) /= NULL_R)) exit
      end do
      if (ndispl == 0 .and. displacement_toolpath == NULL_C) &
          call re_halt('neither DISPLACEMENTS nor DISPLACEMENT_TOOLPATH specified')
      if (ndispl > 0 .and. displacement_toolpath /= NULL_C) &
          call re_halt('both DISPLACEMENTS and DISPLACEMENT_TOOLPATH specified')
      if (ndispl > 0) then
        if (any(displacements(:,:ndispl) == NULL_R)) &
            call re_halt('DISPLACEMENTS not fully specified')
        call params%set('displacements', displacements(:,:ndispl))
      else
        if (.not.known_toolpath(displacement_toolpath)) &
            call re_halt('unknown toolpath: ' // trim(displacement_toolpath))
        call params%set('displacement-toolpath', trim(displacement_toolpath))
      end if
   end if

  end subroutine read_enclosure_namelist


  subroutine init_encl_list(this, params, stat, errmsg)

    use exodus_mesh_type
    use exodus_mesh_io, only: read_exodus_mesh
    use re_encl_type
    use toolpath_type
    use toolpath_table

    type(encl_list), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, n
    type(exodus_mesh) :: mesh
    character(:), allocatable :: mesh_file, name
    logical, allocatable :: mask(:)
    type(toolpath), pointer :: tp => null()

    !! Initialize the base enclosure held by the parent ENCL type.
    if (scl_rank() == 1) then
      !! Read the volume mesh.
      call params%get('mesh-file', mesh_file)
      call re_info('Reading mesh from ' // mesh_file)
      call read_exodus_mesh(mesh_file, mesh)
      call merge_block_ids(mesh, params)
      call re_info('Generating enclosure surface')
      call check_for_altered_surfaces(mesh, params)
      !! Initialize the parent ENCL type object.
      call encl_init_aux(this, mesh, params, stat, errmsg)
    end if
    call scl_bcast(stat)
    if (stat /= 0) then
      call scl_bcast_alloc(errmsg)
      return
    end if

    !! Generate a mask of all cells being displaced, if any.
    if (scl_rank() == 1) call get_displaced_elem_mask(mesh, params, mask, stat, errmsg)
    call scl_bcast(stat)
    if (stat /= 0) then
      call scl_bcast_alloc(errmsg)
      return
    end if

    !! Initialize the additional components of the ENCL_LIST type object.
    if (scl_rank() == 1) then
      if (allocated(mask)) then ! have displacements
        !! Mask of enclosure surface nodes being moved.
        allocate(this%mask(this%nnode), source=.false.)
        do j = 1, this%nface
          if (mask(this%src_elem(j))) then
            associate (face => this%fnode(this%xface(j):this%xface(j+1)-1))
              this%mask(face) = .true.
            end associate
          end if
        end do
        !! List of displacements and associated labels.
        if (params%is_parameter('displacements')) then
          !! Get displacements directly from the input.
          call params%get('displacements', this%dx)
          this%n = size(this%dx,dim=2)
          block
            character(8) :: fmt
            write(fmt,'(i0)') this%n
            n = max(3, len_trim(fmt))
            write(fmt,'("(i",i0,".",i0,")")') n, n
            allocate(character(n)::this%label(this%n))
            write(this%label,fmt) (j, j=1, this%n)
          end block
        else
          !! Extract displacements and labels from a partitioned toolpath.
          call params%get('displacement-toolpath', name)
          tp => toolpath_ptr(name)
          if (tp%has_partition()) then
            call tp%get_partition(coord=this%dx, hash=this%label)
            !! Cull any duplicates (possible)
            n = 1 ! top of list of uniques
            do j = 2, size(this%dx,dim=2)
              if (any(this%label(j) == this%label(1:n))) cycle
              n = n + 1
              if (j == n) cycle
              this%dx(:,n) = this%dx(:,j)
              this%label(n) = this%label(j)
            end do
            if (n < size(this%dx,dim=2)) then
              this%dx = this%dx(:,:n)
              this%label = this%label(:n)
            end if
            this%n = size(this%dx,dim=2)
          else
            stat = -1
            errmsg = 'toolpath is not partitioned'
          end if
        end if
      end if
    end if
    call scl_bcast(stat)
    if (stat /= 0) then
      call scl_bcast_alloc(errmsg)
      return
    end if

    if (scl_rank() == 1) then
      if (allocated(mask)) then ! have displacements
        !! Define the initial enclosure in the sequence.
        this%index = 1
        if (this%n > 1) this%x0 = this%x
        do j = 1, this%nnode
          if (this%mask(j)) this%x(:,j) = this%x(:,j) + this%dx(:,1)
        end do
      else  ! no displacements -- just the single base enclosure.
        this%n = 1
        this%index = 1
      end if
    end if

    call this%bcast

  end subroutine init_encl_list

  subroutine init_encl(this, params, stat, errmsg)

    use exodus_mesh_type
    use exodus_mesh_io, only: read_exodus_mesh
    use re_encl_type

    type(encl), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(exodus_mesh) :: mesh
    character(:), allocatable :: mesh_file
    logical,  allocatable :: mask(:)
    real(r8), allocatable :: disp(:)

    if (scl_rank() == 1) then
      call params%get('mesh-file', mesh_file)
      call read_exodus_mesh(mesh_file, mesh)
      call merge_block_ids(mesh, params)
      call check_for_altered_surfaces(mesh, params)
      call encl_init_aux(this, mesh, params, stat, errmsg)
    end if
    call scl_bcast(stat)
    if (stat /= 0) then
      call scl_bcast_alloc(errmsg)
      return
    end if

    call this%bcast

  end subroutine init_encl

  !! This auxiliary subroutine initializes an ENCL object, extracting the
  !! surface from the passed MESH object according to the passed parameter
  !! list. [SERIAL PROCEDURE]

  subroutine encl_init_aux(this, mesh, params, stat, errmsg)

    use re_encl_type
    use exodus_mesh_type

    class(encl), intent(out) :: this
    type(exodus_mesh), intent(in) :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    integer, allocatable :: ssid(:), ebid(:), disp_ssid(:)
    logical, allocatable :: mirror(:)
    real(r8) :: scale

    !! Extract the surface from the mesh
    call params%get('side-set-ids', ssid)
    call params%get('ignore-block-ids', ebid, default=[integer::])
    call extract_surface_from_exodus(mesh, ssid, ebid, this, stat, errmsg)
    if (stat /= 0) return
    call params%get('coord-scale-factor', scale, default=1.0_r8)
    if (scale /= 1.0_r8) this%x = scale * this%x

    !! Other parameters destined for Chaparral.
    call params%get('name', this%name)
    call params%get('mirror', mirror)
    this%mirror = mirror
    call params%get('rot-axis', this%rot_axis)
    call params%get('num-rot', this%num_rot)

  end subroutine encl_init_aux

  !! This auxiliary procedure handles the optional merging of exodus block IDs;
  !! IDs of secondary blocks replaced with congruent values. [SERIAL PROCEDURE]

  subroutine merge_block_ids(mesh, params)
    use exodus_mesh_type
    use string_utilities, only: i_to_c
    type(exodus_mesh), intent(inout) :: mesh
    type(parameter_list), intent(inout) :: params
    integer :: n, exodus_block_modulus, new_id
    character(:), allocatable :: msg
    !! Overwrite block IDs with their congruent values (10000 Cubit default for 3D meshes)
    call params%get('exodus-block-modulus', exodus_block_modulus, default=10000)
    if (exodus_block_modulus > 0) then
      do n = 1, mesh%num_eblk
        associate (id => mesh%eblk(n)%id)
          new_id = modulo(id, exodus_block_modulus)
          if (new_id /= id) then
            msg = 'Element block ' // i_to_c(id) // ' merged with block ' // i_to_c(new_id)
            call re_info(msg)
            id = new_id
          end if
        end associate
      end do
    end if
  end subroutine merge_block_ids

  !! Ignoring an element block may alter the surface defined by a side set that
  !! references it.  This auxiliary subroutine examines the side sets used to
  !! define the enclosure surface for any reference to an ignored element block
  !! and writes a warning message if any are found. [SERIAL PROCEDURE]

  subroutine check_for_altered_surfaces(mesh, params)

    use exodus_mesh_type
    use string_utilities, only: i_to_c

    type(exodus_mesh), intent(in) :: mesh
    type(parameter_list) :: params

    integer :: i, j, n, m
    integer, allocatable :: list(:), ssid(:)
    logical :: mask(mesh%num_eblk), warned
    integer :: refcount(mesh%num_eblk)
    character(:), allocatable :: string

    if (.not.params%is_parameter('ignore-block-ids')) return ! nothing to check

    !! Ignored element block mask.
    call params%get('ignore-block-ids', list)
    do j = 1, mesh%num_eblk
      mask(j) = any(mesh%eblk(j)%id == list)
    end do

    !! Look for surface side sets that use an ignored element block.
    warned = .false.
    call params%get('side-set-ids', ssid)
    do j = 1, mesh%num_sset
      if (all(mesh%sset(j)%id /= ssid)) cycle ! unused side set
      !! Count the references to ignored blocks by this side set.
      refcount = 0  ! block reference count
      do i = 1, mesh%sset(j)%num_side ! iterate over sides in the side set
        !! Compute the block index N referenced by the side
        n = 1
        m = mesh%sset(j)%elem(i)
        do while (m > mesh%eblk(n)%num_elem)
          m = m - mesh%eblk(n)%num_elem
          n = n + 1
        end do
        if (mask(n)) refcount(n) = refcount(n) + 1
      end do
      !! Warn if any ignored block has a non-zero reference count.
      list = pack(mesh%eblk%id, mask=(refcount>0))
      if (size(list) > 0) then
        string = i_to_c(list(1))
        do i = 2, size(list)
          if (any(list(i) == list(:i-1))) cycle ! skip duplicates
          string = string // ', ' // i_to_c(list(i))
        end do
        call re_info('WARNING: side set ' // i_to_c(mesh%sset(j)%id) // &
                     ' surface may be altered by ignoring blocks ' // string // ';')
        warned = .true.
      end if
    end do
    if (warned) call re_info('         please visually verify the computed enclosure surface.')

  end subroutine check_for_altered_surfaces

  !! This auxiliary subroutine returns an element-based mask that identifies
  !! the elements that are to be displaced.  An unallocated mask is returned
  !! if no displacement was specified. This also verifies that the displaced
  !! elements are not connected to the remaining elements. [SERIAL PROCEDURE]

  subroutine get_displaced_elem_mask(mesh, params, mask, stat, errmsg)

    use exodus_mesh_type
    use string_utilities, only: i_to_c

    type(exodus_mesh), intent(in) :: mesh
    type(parameter_list), intent(inout) :: params
    logical, allocatable, intent(out) :: mask(:)
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    integer :: j, n, offset
    integer, allocatable :: dbid(:), ibid(:)
    logical, allocatable :: node_mask(:)

    stat = 0

    call params%get('ignore-block-ids', ibid, default=[integer::])
    do n = 1, size(ibid)  ! verify ignored blocks are valid blocks
      if (.not.any(ibid(n) == mesh%eblk%id)) then
        stat = -1
        errmsg = 'invalid ignored block ID: ' // i_to_c(ibid(n))
        return
      end if
    end do

    call params%get('displace-block-ids', dbid, default=[integer::])
    do n = 1, size(dbid) ! verify displaced blocks are valid blocks
      if (.not.any(dbid(n) == mesh%eblk%id)) then
        stat = -1
        errmsg = 'invalid displaced block ID: ' // i_to_c(dbid(n))
        return
      end if
      if (any(dbid(n) == ibid)) then
        stat = -1
        errmsg = 'invalid displaced block ID: ' // i_to_c(dbid(n))
        return
      end if
    end do

    if (size(dbid) == 0) return ! mask not allocated
    allocate(mask(mesh%num_elem), source=.false.)

    !! Tag the displaced elements and nodes.
    allocate(node_mask(mesh%num_node), source=.false.)
    offset = 0
    do n = 1, mesh%num_eblk
      if (any(mesh%eblk(n)%id == dbid)) then
        mask(offset+1:offset+mesh%eblk(n)%num_elem) = .true.
        do j = 1, mesh%eblk(n)%num_elem
          node_mask(mesh%eblk(n)%connect(:,j)) = .true.
        end do
      end if
      offset = offset + mesh%eblk(n)%num_elem
    end do

    !! Verify the displaced blocks are disconnected from the others.
    do n = 1, mesh%num_eblk
      if (any(mesh%eblk(n)%id == dbid)) cycle
      if (any(mesh%eblk(n)%id == ibid)) cycle
      do j = 1, mesh%eblk(n)%num_elem
        if (any(node_mask(mesh%eblk(n)%connect(:,j)))) then
          stat = -1
          errmsg = 'displaced blocks connected to undisplaced blocks'
          return
        end if
      end do
    end do

  end subroutine get_displaced_elem_mask

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

    use exodus_mesh_type
    use string_utilities, only: i_to_c
    use re_encl_type
    use facet_hash_type

    type(exodus_mesh), intent(in)  :: mesh
    integer, intent(in) :: ssid(:)  ! IDs of side sets describing the surface
    integer, intent(in) :: ebid(:)  ! IDs of element blocks to omit from the mesh
    type(encl), intent(out) :: surf  ! mesh of the extracted surface
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

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
      type(facet_hash) :: hpar
      type(side_row), pointer :: row(:) => null()
    end type
    type(side_table) :: stab

    !! Another hash table structure
    type :: nhbr_row
      integer, pointer :: list(:) => null() ! key values, no data
    end type
    type :: nhbr_table
      integer :: size=0, nentry=0
      type(facet_hash) :: hpar
      type(nhbr_row), pointer :: row(:) => null()
    end type
    type(nhbr_table) :: ntab

    integer, parameter :: MAX_SIDE = 6  ! works for tet, hex and wedge elements.

    integer, target, save :: TETRA4_SIDE_SIG(4)
    data TETRA4_SIDE_SIG/b'1011', b'1110', b'1101', b'0111'/

    integer, target, save :: PYRAMID5_SIDE_SIG(5)
    data PYRAMID5_SIDE_SIG/b'10011', b'10110', b'11100', b'11001', b'01111'/

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
        list => mesh%side_node_list(j, k)
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

      select case (mesh%eblk(b)%elem_type(1:3))
      case ('TET')
        side_sig => TETRA4_SIDE_SIG
      case ('PYR')
        side_sig => PYRAMID5_SIDE_SIG
      case ('WED')
        side_sig => WEDGE6_SIDE_SIG
      case ('HEX')
        side_sig => HEX8_SIDE_SIG
      case default
        stat = -1
        errmsg = 'unable to handle element type: ' // mesh%eblk(b)%elem_type
        return
      end select

      do l = 1, mesh%eblk(b)%num_elem
        j = j + 1
        if (omit_elem(j)) cycle

        !! Mark the element vertices that are surface nodes.
        associate(list => mesh%eblk(b)%connect(:,l))
          bitmask = 0
          do k = 1, size(list)
            if (nmap(list(k))>0) bitmask = ibset(bitmask, k-1)
          end do
        end associate

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
          list => mesh%side_node_list(j, k, reverse=.true.)
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
      call table%hpar%init (hsize=table%size, kmax=key_max)
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
      call table%hpar%hash ([j, k], n, incr)
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
      call table%hpar%hash ([j, k], n, incr)
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
      call table%hpar%init (hsize=table%size, kmax=key_max)
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
      call table%hpar%hash (list, n, inc)
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
      call table%hpar%hash (list, n, inc)
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

end module re_exodus_encl
