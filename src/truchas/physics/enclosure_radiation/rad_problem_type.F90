!!
!! RAD_PROBLEM_TYPE
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

!#define GMV_SURFACE_DIAGNOSTICS

module rad_problem_type

  use kinds, only: r8
  use parallel_communication
  use parallel_permutations
  use rad_solver_type
  use encl_vf_class
  use static_vf_class
  use moving_vf_type
  use rad_encl_type
  use sim_event_queue_type, only: event_action
  implicit none
  private

  integer, parameter :: PC_JACOBI = 1
  integer, parameter :: PC_CHEBY  = 2

  type, public :: rad_problem
    integer, allocatable :: faces(:)
    logical, allocatable :: fmask(:)
    !! The rest are private
    type(rad_solver) :: sol
    integer, allocatable :: ge_faces(:)
    integer :: nface_hc
    integer :: nface_er
    type(par_perm) :: perm_er_to_hc, perm_hc_to_er
    !! Radiosity system preconditioner parameters
    integer :: pc_numitr = 1
    integer :: pc_method = PC_JACOBI
    type(rad_encl), allocatable :: encl ! radiation enclosure surface
    class(encl_vf), allocatable :: vf   ! enclosure view factors
  contains
    procedure :: init
    procedure :: residual
    procedure :: heat_flux
    procedure :: solve_radiosity
    procedure :: precon
    procedure :: precon_matvec1
    procedure :: rhs
    procedure :: rhs_deriv
    procedure :: update_moving_vf
    procedure :: add_moving_vf_events
  end type rad_problem

  integer, parameter, public :: DT_POLICY_NONE   = 0
  integer, parameter, public :: DT_POLICY_NEXT   = 1
  integer, parameter, public :: DT_POLICY_FACTOR = 2
  integer, parameter, public :: DT_POLICY_VALUE  = 3

  type, extends(event_action), public :: vf_event
    private
    integer  :: dt_policy = DT_POLICY_NONE
    real(r8) :: c = 0.0_r8
  contains
    procedure :: init_dt => vf_event_init_dt
  end type

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! INIT
 !!

  subroutine init (this, mesh, name, params, tinit)

    use unstr_mesh_type
    use physical_constants, only: stefan_boltzmann, absolute_zero
    use scalar_func_class
    use scalar_func_factories, only: alloc_scalar_func
    use rad_encl_file_type
    use rad_encl_func_type
    use parameter_list_type
    use truchas_logging_services

    class(rad_problem), target, intent(out) :: this
    type(unstr_mesh),   intent(in)  :: mesh
    character(len=*),   intent(in)  :: name
    type(parameter_list), intent(inout) :: params
    real(r8), intent(in) :: tinit

    integer :: stat
    logical :: flag
    character(:), allocatable :: filename, method, errmsg
    type(rad_encl_file) :: file
    character(len=127) :: string
    integer, allocatable :: setids(:)
    real(r8) :: csf, tol
    class(scalar_func), allocatable :: tamb
    type(rad_encl_func), allocatable :: eps
    type(parameter_list), pointer :: plist

    call TLS_info ('  initializing enclosure radiation problem "' // trim(name) // '" ...')

    !! Load the view factor data
    if (params%is_parameter('toolpath')) then
      call alloc_moving_vf(tinit, params, this%vf, filename)
    else  ! static view factors
      call alloc_static_vf(params, this%vf, filename)
    end if

    !! Load the distributed enclosure surface
    allocate(this%encl)
    call params%get('coord-scale-factor', csf, default=1.0d0)
    if (is_IOP) call file%open_ro(filename)
    call this%encl%init(file, this%vf%face_map, csf)
    if (is_IOP) call file%close
    this%nface_er = this%encl%nface_onP
    ASSERT(this%nface_er == this%vf%nface)

    !! Identify the HC faces that correspond to the ER surface faces.
    call connect_to_mesh(mesh, filename, this%faces, this%ge_faces, stat)
    if (stat /= 0) then
      write(string,'(a,i0)') 'no matching mesh face for enclosure face ', stat
      call TLS_fatal(string)
    end if
    this%nface_hc = size(this%faces)
    ASSERT(this%vf%nface_tot == global_sum(this%nface_hc))
    allocate(this%fmask(this%nface_hc), source=.true.)

    !! Verify that the identified faces are boundary faces.
    call boundary_face_check (mesh, this%faces, stat, setids)
    if (stat /= 0) then
      write(string,'(4x,a,i0,a)') 'Error: ', stat, ' enclosure faces are not mesh boundary faces;'
      call TLS_info (trim(string))
      if (size(setids) > 0) then
        write(string,'(11x,a,i0,99(:,1x,i0))') 'faces belong to face sets ', setids
      else
        write(string,'(11x,a)') 'faces belong to no mesh face sets.'
      end if
      call TLS_info (trim(string))
      write(string,'(11x,a)') 'Perhaps an internal interface should have been defined?'
      call TLS_info (trim(string))
      deallocate(setids)
      call TLS_fatal ('error initializing enclosure radiation problem "' // trim(name) // '"')
    end if

    !! Create the parallel permutations between the HC and ER partitions.
    call create_par_perm (this%ge_faces, this%vf%face_map, this%perm_hc_to_er, this%perm_er_to_hc)
    INSIST(defined(this%perm_er_to_hc))
    INSIST(defined(this%perm_hc_to_er))

#ifdef GMV_SURFACE_DIAGNOSTICS
    call write_mesh_surface (trim(name)//'-mesh-surf.gmv', mesh, this%faces, this%ge_faces)
    call write_encl_surface (trim(name)//'-encl-surf.gmv', this%sol)
#endif

    !! Check that HC radiation faces are geometrically equal to the enclosure surface.
    call params%get('check-geometry', flag, default=.false.)
    if (flag) then
      call TLS_warn ('Not able to check matching enclosure surface geometry for mixed cell meshes')
      !FIXME -- REWRITE CHECK_SURFACE TO OPERATE ON A UNSTR_MESH TYPE MESH
      !call check_surface (mesh, this%encl, this%faces, this%perm_er_to_hc, stat)
      !if (stat /= 0) then
      !  write(string,'(4x,a,i0,a)') 'enclosure surface doesn''t match mesh: ', stat, ' faces differ'
      !  call TLS_fatal (string)
      !end if
    end if

    !! Create the enclosure radiation solver.
    call this%sol%init(this%vf)

    !! Set ER system parameters.
    call alloc_scalar_func(params, 'ambient-temp', tamb, stat, errmsg)
    INSIST(stat == 0)
    call this%sol%set_ambient (tamb)

    allocate(eps)
    plist => params%sublist('surfaces')
    call define_emissivity(this%encl, plist, eps, stat, errmsg)
    if (stat /= 0) call TLS_fatal('error defining emissivity: ' // errmsg)
    call this%sol%set_emissivity (eps)

    call this%sol%set_stefan_boltzmann (stefan_boltzmann)
    call this%sol%set_absolute_zero (absolute_zero)

    !! Go ahead and set the Chebyshev iteration parameters.  This assumes
    !! that the emissivities are time-independent so that we can compute
    !! the parameters once.  Really these need to be computed whenever the
    !! emissivities change significantly, but this involves computing
    !! eigenvalues of the radiosity system and isn't cheap.
    call this%sol%set_cheby_param (time=tinit)

    call params%get('error-tol', tol, default=1e-3_r8)
    call this%sol%set_solver_controls (tol, maxitr=500) !TODO! input maxitr value

    call params%get('precon-method', method, default='JACOBI')
    call params%get('precon-iter', this%pc_numitr, default=1)
    select case (method)
    case ('JACOBI')
      this%pc_method = PC_JACOBI
    case ('CHEBYSHEV')
      this%pc_method = PC_CHEBY
    case default
      INSIST(.false.)
    end select

  end subroutine init

  subroutine define_emissivity(encl, plist, eps, stat, errmsg)

    use rad_encl_func_type
    use parameter_list_type
    use scalar_func_factories

    type(rad_encl), intent(in), target :: encl
    type(parameter_list), intent(inout) :: plist
    type(rad_encl_func), intent(out) :: eps
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: params
    class(scalar_func), allocatable :: f
    integer, allocatable :: block_ids(:)

    call eps%prep(encl)
    piter = parameter_list_iterator(plist, sublists_only=.true.)
    do while (.not.piter%at_end())
      params => piter%sublist()
      call alloc_scalar_func(params, 'emissivity', f, stat, errmsg)
      if (stat /= 0) exit
      call params%get('face-block-ids', block_ids)
      call eps%add_function(f, block_ids, stat, errmsg)
      if (stat /= 0) exit
      call piter%next
    end do
    if (stat /= 0) then
      errmsg = 'error defining emissivity for surface ' // piter%name() // ': ' // errmsg
      return
    end if
    call eps%done(stat, errmsg)

  end subroutine

  subroutine alloc_static_vf(params, vf, filename)

    use parameter_list_type
    use static_vf_class
    use facet_vf_type
    use patch_vf_type
    use rad_encl_file_type
    use truchas_logging_services

    type(parameter_list), intent(inout) :: params
    class(encl_vf), allocatable, intent(out) :: vf
    character(:), allocatable, intent(out) :: filename

    type(rad_encl_file) :: file
    class(static_vf), allocatable :: svf
    logical :: has_patches

    call params%get('enclosure-filename', filename)
    call TLS_info('    reading enclosure radiation view factors from ' // filename)
    if (is_IOP) call file%open_ro(filename)
    if (is_IOP) has_patches = file%has_patches()
    call broadcast(has_patches)

    if (has_patches) then
      allocate(patch_vf :: svf)
    else
      allocate(facet_vf :: svf)
    end if
    call svf%init(file)
    call move_alloc(svf, vf)
    if (is_IOP) call file%close

  end subroutine

  subroutine alloc_moving_vf(tinit, params, vf, filename)

    use parameter_list_type
    use moving_vf_type
    use toolpath_type
    use toolpath_table, only: toolpath_ptr

    real(r8), intent(in) :: tinit
    type(parameter_list), intent(inout) :: params
    class(encl_vf), allocatable, intent(out) :: vf
    character(:), allocatable, intent(out) :: filename

    integer :: n, j
    character(:), allocatable :: tpname, hash(:), ext, basename, filenames(:)
    real(r8), allocatable :: times(:)
    type(toolpath), pointer :: tp
    type(moving_vf), allocatable :: mvf

    call params%get('toolpath', tpname)
    tp => toolpath_ptr(tpname)
    INSIST(tp%has_partition())
    call tp%get_partition(hash=hash, time=times)

    call params%get('enclosure-filename', filename)
    ext = file_extension(filename)
    n = index(filename,ext,back=.true.) - 1 ! strip file extension
    basename = filename(:n) // '.'

    n = len(basename) + len(hash) + len(ext)
    allocate(character(n) :: filenames(size(hash)))
    do j = 1, size(hash)
      filenames(j) = basename // hash(j) // ext
    end do

    allocate(mvf)
    call mvf%init(tinit, times, filenames)
    filename = mvf%filename()
    call move_alloc(mvf, vf)

  contains

    !! Return the file extension, or empty string if none
    function file_extension(filename) result(ext)
      character(*), intent(in) :: filename
      character(:), allocatable :: ext
      integer :: n
      n = scan(filename, '/', back=.true.)
      ext = filename(n+1:)
      n = scan(ext, '.', back=.true.)
      if (n > 1) then
        ext = ext(n:)
      else
        ext = ''
      end if
    end function

  end subroutine alloc_moving_vf

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! CONNECT_TO_MESH
 !!
 !! This auxillary routine re-establishes the correspondence between faces of
 !! the radiation enclosure and faces of the distributed mesh.  FILE is the
 !! path to the radiation enclosure file.  LM_FACES and GE_FACES are integer
 !! rank-1 array pointers that are allocated by this routine.  LM_FACE(:) is
 !! a per-process list of local mesh face indices and GE_FACE(:) is the
 !! corresponding list of global enclosure face indices; that is, for each j,
 !! local face LM_FACE(j) corresponds to global enclosure face GE_FACE(j).
 !! The routine uses the source info stored in the enclosure file to establish
 !! the correspondence.  If an enclosure face cannot be matched to a mesh face,
 !! the (global) index of that enclosure face is returned in STAT; STAT = 0
 !! indicates success.
 !!
 !! NB: this routine doesn't check for correctly matching geometry.
 !!

  subroutine connect_to_mesh (mesh, path, lm_faces, ge_faces, stat)

    use unstr_mesh_type
    use permutations
    use rad_encl_file_type

    type(unstr_mesh), intent(in) :: mesh
    character(*), intent(in) :: path
    integer, allocatable :: lm_faces(:) ! returned list of local HC mesh faces
    integer, allocatable :: ge_faces(:) ! returned list of global ER faces
    integer, intent(out) :: stat

    integer :: j, n, offset, nface
    integer :: last(nPE), bsize(nPE)
    integer, allocatable :: fcell(:), fcell_l(:), fside(:), fside_l(:), perm(:), perm2(:)
    integer, allocatable :: map(:)
    type(rad_encl_file) :: file

    !! Read the source info from the enclosure file: FCELL and FSIDE.
    if (is_IOP) then
      call file%open_ro(path)
      call file%get_source_info(fcell, fside)
      call file%close
      nface = size(fcell)
    else
      allocate(fcell(0), fside(0))
    end if

    !! Mapping from external cell numbers to internal (global) cell numbers.
    allocate(map(merge(mesh%cell_imap%global_size, 0, is_IOP)))
    call collate (mesh%xcell(:mesh%ncell_onP), map)
    if (is_IOP) then
      ASSERT(is_perm(map))
      call invert_perm (map)
    end if

    !! Map the source cell numbers to internal numbers.
    if (is_IOP) then
      do j = nface, 1, -1
        if (fcell(j) < 1 .or. fcell(j) > size(map)) exit  ! no matching mesh face
        fcell(j) = map(fcell(j))  ! map to internal cell index
      end do
    end if
    deallocate(map)
    call broadcast (j)
    if (j /= 0) then
      stat = j
      return
    end if

    !! Sort the cell/side pairs so that they are ordered by cell process rank.
    call collate (mesh%cell_imap%last_gid, last) ! last global cell index on the processes
    if (is_IOP) then
      allocate(perm(nface))
      call partition_sort (last, fcell, bsize, perm)
      call reorder (fcell, perm)
      call reorder (fside, perm)
    end if

    !! Distribute the cell/side pairs to the processes owning the cell.
    call distribute (bsize, n)
    allocate(fcell_l(n), fside_l(n))
    call distribute (fcell, fcell_l)
    call distribute (fside, fside_l)
    deallocate(fcell, fside)
    offset = mesh%cell_imap%first_gid - 1
    fcell_l = fcell_l - offset ! local cell index
    ASSERT(all(fcell_l >= 1))
    ASSERT(all(fcell_l <= mesh%ncell_onP))
    ASSERT(all(fside_l >= 1))

    !! Find the global mesh face that corresponds to each of the enclosure
    !! faces.  FCELL_L holds the results temporarily.
    do j = n, 1, -1
      associate (faces => mesh%cface(mesh%xcface(fcell_l(j)):mesh%xcface(fcell_l(j)+1)-1))
        if (fside_l(j) > size(faces)) exit ! no matching mesh face
        fcell_l(j) = mesh%face_imap%global_index(faces(fside_l(j)))
      end associate
    end do
    stat = global_maxval(j) ! get one of the unmatched faces, if any
    if (stat /= 0) return

    !! Generate the global mapping MAP of enclosure faces to mesh faces.
    allocate(map(merge(nface, 0, is_IOP)))
    call collate (fcell_l, map)
    if (is_IOP) then  ! undo the partition sort
      call reorder (map, perm, forward=.true.)
      ASSERT(all(map >= 1))
      ASSERT(all(map <= mesh%face_imap%global_size))
    end if
    deallocate(fcell_l, fside_l)

    !! Sort the mesh face list MAP so that it is ordered by face process rank.
    call collate (mesh%face_imap%last_gid, last)
    if (is_IOP) then
      call partition_sort (last, map, bsize, perm)
      call reorder (map, perm)
    end if

    !! Distribute the face list MAP, apply local index offset, and sort; this
    !! gives a per-process list LM_FACES of mesh faces that are enclosure faces.
    call distribute (bsize, n)
    allocate(lm_faces(n), perm2(n))
    call distribute (map, lm_faces)
    offset = mesh%face_imap%first_gid - 1
    lm_faces = lm_faces - offset  ! local face index
    call heapsort (lm_faces, perm2)
    call reorder (lm_faces, perm2)
    ASSERT(all(lm_faces >= 1))
    ASSERT(all(lm_faces <= mesh%face_imap%onp_size))

    !! Generate the mapping from mesh enclosure faces to enclosure faces:
    !! it is the composition of the partition sort and the heap sort.
    call broadcast (bsize)
    offset = sum(bsize(:this_PE-1))
    call collate (perm2+offset, map)  ! the global heapsort permutation
    if (is_IOP) then
      ASSERT(is_perm(map))
      do j = 1, size(map)
        map(j) = perm(map(j))
      end do
      ASSERT(is_perm(map))
      deallocate(perm)
    end if
    deallocate(perm2)

    !! Distribute the mapping from mesh enclosure faces to enclosure faces; this
    !! gives a per-process list GE_FACES of enclosure faces corresponding to LM_FACES.
    allocate(ge_faces(n))
    call distribute (map, ge_faces)
    deallocate(map)

    stat = 0  ! success

  end subroutine connect_to_mesh

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! PARTITION_SORT
 !!
 !! This auxillary routine computes the permutation PERM(:) and block sizes
 !! BSIZE(:) such that the reordered list LIST(PERM(:)) will have the property
 !! that the first BSIZE(1) values lie in the interval (0, LAST(1)], the next
 !! BSIZE(2) values belong to the interval (LAST(1), LAST(2)], and so on
 !! through the number of processes.  The expectation is that LAST(j) is the
 !! last global index belonging to process rank j, which describes a block
 !! partition of an index set that the values of LIST belong to.  LIST is not
 !! modified.  PERM is a mapping from new numbering to original, and it
 !! preserves the relative order of elements whose LIST values belong to the
 !! same partition.
 !!

  subroutine partition_sort (last, list, bsize, perm)

    use permutations

    integer, intent(in)  :: last(:)
    integer, intent(in)  :: list(:)
    integer, intent(out) :: bsize(:)
    integer, intent(out) :: perm(:)

    integer :: i, i1, i2, j, n
    integer :: pnum(size(list)), next(size(last))

    ASSERT(size(last) == size(bsize))
    ASSERT(size(list) == size(perm))
    ASSERT(size(list) >= 1)
    ASSERT(all(last(2:)>=last(:size(last)-1)))
    ASSERT(all(list >= 1))
    ASSERT(all(list <= last(size(last))))

    !! Identify the partition to which each list value N belongs using a binary
    !! search on the ordered array LAST: I1, I2 satisfy LAST(I1)<N<=LAST(I2);
    !! search terminates when I2-I1==1, with I2 the partition number.
    bsize = 0
    do j = 1, size(list)
      i1 = 0; i2 = size(last)  ! initialization of search loop
      do while (i2-i1 > 1)
        i = (i1 + i2)/2
        if (list(j) > last(i)) then
          i1 = i
        else
          i2 = i
        end if
      end do
      pnum(j) = i2
      bsize(i2) = bsize(i2) + 1
    end do

    !! NEXT(j) is the next index for partition j.
    next(1) = 1
    do j = 2, size(bsize)
      next(j) = next(j-1) + bsize(j-1)
    end do

    !! Generate the permutation.
    do j = 1, size(list)
      n = next(pnum(j))
      perm(n) = j
      next(pnum(j)) = n + 1
    end do

    ASSERT(is_perm(perm))

  end subroutine partition_sort

  subroutine heapsort (array, perm)

    integer, intent(in)  :: array(:)
    integer, intent(out) :: perm(:)

    integer :: j, tmp

    ASSERT(size(array) == size(perm))

    do j = 1, size(perm)
      perm(j) = j
    end do

    !! Create the heap from the bottom up.
    do j = size(perm)/2, 1, -1 ! start with the last parent
      call sift_down (j, size(perm))
    end do

    do j = size(perm), 2, -1
      tmp = perm(j)
      perm(j) = perm(1)
      perm(1) = tmp
      if (j > 2) call sift_down (1, j-1)
    end do

  contains

    subroutine sift_down (first, last)
      integer, intent(in) :: first, last
      integer :: root, child, proot
      root = first
      proot = perm(root)
      child = 2*root
      do while (child <= last)
        if (child < last) then  ! there is a sibling, take the larger
          if (array(perm(child)) < array(perm(child+1))) child = child + 1
        end if
        if (array(proot) >= array(perm(child))) exit  ! tree at root is a heap
        perm(root) = perm(child)
        root = child
        child = 2*root
      end do
      perm(root) = proot
    end subroutine

  end subroutine heapsort

  subroutine boundary_face_check (mesh, faces, stat, setids)

    use bitfield_type
    use unstr_mesh_type

    type(unstr_mesh), intent(in) :: mesh
    integer, intent(in) :: faces(:)
    integer, intent(out) :: stat
    integer, allocatable, intent(out) :: setids(:)

    integer :: j, n
    type(bitfield) :: bitmask

    stat = 0  ! count of non-boundary faces
    bitmask = ZERO_BITFIELD
    do j = 1, size(faces)
      if (btest(mesh%face_set_mask(faces(j)),pos=0)) cycle  ! a boundary face
      stat = stat + 1
      bitmask = ior(bitmask, mesh%face_set_mask(faces(j)))
    end do

    stat = global_sum(stat)
    if (stat == 0) return

    bitmask = global_ior(bitmask)

    !! Create the list of involved side set IDS.
    n = 0 ! count first to allocate
    do j = 1, size(mesh%face_set_id)
      if (btest(bitmask,j)) n = n + 1
    end do
    allocate(setids(n))
    n = 0 ! now store the data
    do j = 1, size(mesh%face_set_id)
      if (btest(bitmask,j)) then
        n = n + 1
        setids(n) = mesh%face_set_id(j)
      end if
    end do

  end subroutine boundary_face_check

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! CHECK_SURFACE
 !!
 !! This routine verifies that the mesh faces identified as radiation faces
 !! are geometrically identical to the enclosure surface faces as required.
 !! STAT returns 0 if the two sets of faces are identical; otherwise it returns
 !! the number of faces that were found to differ.  MESH is the distributed
 !! heat conduction mesh and ENCL is the distributed radiation enclosure.
 !! FACES(:) is the per-process list of mesh faces that are radiation faces.
 !! PERM_ER_TO_HC is the parallel permutation mapping radiation enclosure faces
 !! to mesh radiation faces.
 !!
 !! IMPLEMENTATION NOTE
 !!
 !! This relies on face vertices being ordered around the face and the property
 !! that mesh boundary faces (a radiating face must be on the boundary) have an
 !! outward orientation -- the same orientation that the enclosure surface face
 !! should have.
 !!

!  subroutine check_surface (mesh, encl, faces, perm_er_to_hc, stat)
!
!    use unstr_mesh_type
!    use rad_encl_type
!
!    type(unstr_mesh), intent(in) :: mesh
!    type(rad_encl),  intent(in) :: encl
!    integer,         intent(in) :: faces(:)
!    type(par_perm),  intent(in) :: perm_er_to_hc
!    integer, intent(out) :: stat
!
!    integer :: i, j, k, i1, i2, dimen, nvert, badcnt
!    real(r8), allocatable :: mesh_vert(:,:,:)
!    real(r8) :: tolsq
!
!    !! Spatial tolerance (squared); I don't like this -- should be user-specified.
!    tolsq = 1.0d-8*sum((maxval(encl%coord,2) - minval(encl%coord,2))**2)
!
!    dimen = size(mesh%x,dim=1)
!    nvert = size(mesh%fnode,dim=1) ! number of vertices per face in mesh
!
!    !! Verify that the spatial dimensions are the same.
!    if (dimen /= size(encl%coord,dim=1)) then
!      stat = -1
!      return
!    end if
!
!    !! Gather the mesh face vertex coordinates.
!    allocate(mesh_vert(dimen,nvert,encl%nface))
!    do k = 1, nvert
!      do i = 1, dimen
!        call reorder (perm_er_to_hc, mesh_vert(i,k,:), mesh%x(i,mesh%fnode(k,faces)))
!      end do
!    end do
!
!    badcnt = 0  ! count of mismatched faces
!
!    do j = 1, encl%nface
!
!      associate (fnode => encl%fnode(encl%xface(j):encl%xface(j+1)-1))
!
!        if (size(fnode) /= nvert) then  ! different type faces
!          badcnt = badcnt + 1
!          cycle
!        end if
!
!        !! Locate the mesh face vertex that matches the first enclosure face vertex.
!        do i1 = nvert, 1, -1
!          if (vertices_match(mesh_vert(:,i1,j), encl%coord(:,fnode(1)))) exit
!        end do
!        if (i1 == 0) then ! no matching vertex
!          badcnt = badcnt + 1
!          cycle
!        end if
!
!        !! Verify that the remaining vertices match, in order.
!        do i2 = 2, nvert
!          i1 = modulo(i1,nvert) + 1
!          if (vertices_match(mesh_vert(:,i1,j), encl%coord(:,fnode(i2)))) cycle
!          badcnt = badcnt + 1
!          exit
!        end do
!
!      end associate
!
!    end do
!
!    deallocate(mesh_vert)
!    stat = global_sum(badcnt)
!
!  contains
!
!    logical function vertices_match (a, b)
!      real(r8), intent(in) :: a(:), b(:)
!      vertices_match = (sum((a-b)**2) < tolsq)
!    end function vertices_match
!
!  end subroutine check_surface

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SOLVE_RADIOSITY
 !!
 !! This subroutine solves the radiosity system for the radiosity from the
 !! enclosure surface given the temperature on the surface.  The surface
 !! emissivities and ambient temperature parameters of the system may be
 !! time-dependent.
 !!

  subroutine solve_radiosity (this, time, temp, qrad, stat, numitr, error)

    class(rad_problem), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: temp(:)
    real(r8), intent(inout) :: qrad(:)
    integer,  intent(out) :: stat
    integer,  intent(out) :: numitr
    real(r8), intent(out) :: error

    real(r8), dimension(this%nface_er) :: qrad_er, temp_er

    ASSERT(size(temp) == this%nface_hc)
    ASSERT(size(qrad) == this%nface_hc)

    !! Form the local temperature and radiosity vectors in the ER ordering.
    call reorder (this%perm_er_to_hc, temp_er, temp)
    call reorder (this%perm_er_to_hc, qrad_er, qrad)

    select type (vf => this%vf)
    type is (moving_vf)
      call vf%set_time(time)
    end select

    call this%sol%solve_radiosity (time, temp_er, qrad_er, stat, numitr, error)

    !! Form the local radiosity vector in the HC ordering.
    call reorder (this%perm_hc_to_er, qrad, qrad_er)

  end subroutine solve_radiosity

  subroutine precon (this, time, z)

    class(rad_problem), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(inout) :: z(:)

    real(r8) :: z_er(this%nface_er)

    ASSERT(size(z) == this%nface_hc)

    select type (vf => this%vf)
    type is (moving_vf)
      call vf%set_time(time)
    end select

    select case (this%pc_method)
    case (PC_JACOBI)
      if (this%pc_numitr == 1) return ! the effect of 1 iteration
      call reorder (this%perm_er_to_hc, z_er, z)  ! form the Z vector in the ER ordering
      call this%sol%precon_jacobi (time, this%pc_numitr, z_er)
      call reorder (this%perm_hc_to_er, z, z_er)  ! form the Z vector in the HC ordering
    case (PC_CHEBY)
      call reorder (this%perm_er_to_hc, z_er, z)  ! form the Z vector in the ER ordering
      call this%sol%precon_cheby (time, this%pc_numitr, z_er)
      call reorder (this%perm_hc_to_er, z, z_er)  ! form the Z vector in the HC ordering
    end select

  end subroutine precon

  subroutine precon_matvec1 (this, time, z)

    class(rad_problem), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(inout) :: z(:)

    real(r8) :: z_er(this%nface_er)

    ASSERT(size(z) == this%nface_hc)

    !! Form the Z vector in the ER ordering.
    call reorder (this%perm_er_to_hc, z_er, z)

    select type (vf => this%vf)
    type is (moving_vf)
      call vf%set_time(time)
    end select

    call this%sol%precon_matvec1 (time, z_er)

    !! Form the Z vector in the HC ordering.
    call reorder (this%perm_hc_to_er, z, z_er)

  end subroutine precon_matvec1

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! HEAT_FLUX
 !!
 !! Calculates the net heat flux through the enclosure surface given the
 !! radiosity from the enclosure surface.  The flux may depend on a possibly
 !! time-dependent ambient temperature.
 !!

  subroutine heat_flux (this, time, qrad, temp, flux)

    class(rad_problem), intent(inout)  :: this
    real(r8), intent(in)  :: time
    real(r8), intent(in)  :: qrad(:)  ! local radiosity vector
    real(r8), intent(in)  :: temp(:)
    real(r8), intent(out) :: flux(:)  ! local heat flux vector

    real(r8), dimension(this%nface_er) :: qrad_er, temp_er, flux_er

    ASSERT(size(qrad) == this%nface_hc)
    ASSERT(size(temp) == this%nface_hc)
    ASSERT(size(flux) == this%nface_hc)

    !! Form the local radiosity vector in the ER ordering.
    call reorder (this%perm_er_to_hc, qrad_er, qrad)
    call reorder (this%perm_er_to_hc, temp_er, temp)

    select type (vf => this%vf)
    type is (moving_vf)
      call vf%set_time(time)
    end select

    call this%sol%heat_flux (time, qrad_er, temp_er, flux_er)
!    call ERS_compute_heat_flux (this%sol, time, qrad_er, flux_er)

    !! Form the local heat flux vector in the HC ordering.
    call reorder (this%perm_hc_to_er, flux, flux_er)

  end subroutine heat_flux

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! RESIDUAL
 !!
 !! This subroutine calculates the residual of the radiosity system given the
 !! temperature and radiosity on the enclosure surface.  The emissivities and
 !! ambient temperature parameters of the system may be time-dependent.
 !!

  subroutine residual (this, time, qrad, temp, res)

    class(rad_problem), intent(inout) :: this
    real(r8), intent(in)  :: time
    real(r8), intent(in)  :: qrad(:)
    real(r8), intent(in)  :: temp(:)
    real(r8), intent(out) :: res(:)

    real(r8), dimension(this%nface_er) :: qrad_er, temp_er, res_er

    ASSERT(size(qrad) == this%nface_hc)
    ASSERT(size(temp) == this%nface_hc)
    ASSERT(size(res)  == this%nface_hc)

    !! Form the local radiosity and temperature vectors in the ER ordering.
    call reorder (this%perm_er_to_hc, temp_er, temp)
    call reorder (this%perm_er_to_hc, qrad_er, qrad)

    select type (vf => this%vf)
    type is (moving_vf)
      call vf%set_time(time)
    end select

    call this%sol%residual (time, qrad_er, temp_er, res_er)

    !! Form the local residual vector in the HC ordering.
    call reorder (this%perm_hc_to_er, res, res_er)

  end subroutine residual

  subroutine rhs_deriv (this, time, temp, drhs)

    class(rad_problem), intent(inout) :: this
    real(r8), intent(in)  :: time
    real(r8), intent(in)  :: temp(:)
    real(r8), intent(out) :: drhs(:)

    real(r8), dimension(this%nface_er) :: temp_er, drhs_er

    ASSERT(size(temp) == this%nface_hc)
    ASSERT(size(drhs) == this%nface_hc)

    !! Form the local temperature vector in the ER ordering.
    call reorder (this%perm_er_to_hc, temp_er, temp)

    select type (vf => this%vf)
    type is (moving_vf)
      call vf%set_time(time)
    end select

    call this%sol%rhs_deriv (time, temp_er, drhs_er)

    !! Form the local derivative vector in the HC ordering.
    call reorder (this%perm_hc_to_er, drhs, drhs_er)

  end subroutine rhs_deriv

  subroutine rhs (this, time, temp, b)

    class(rad_problem), intent(inout) :: this
    real(r8), intent(in)  :: time
    real(r8), intent(in)  :: temp(:)
    real(r8), intent(out) :: b(:)

    real(r8), dimension(this%nface_er) :: temp_er, rhs_er

    ASSERT(size(temp) == this%nface_hc)
    ASSERT(size(b) == this%nface_hc)

    !! Form the local temperature vector in the ER ordering.
    call reorder (this%perm_er_to_hc, temp_er, temp)

    select type (vf => this%vf)
    type is (moving_vf)
      call vf%set_time(time)
    end select

    call this%sol%rhs (time, temp_er, rhs_er)

    !! Form the local rhs vector in the HC ordering.
    call reorder (this%perm_hc_to_er, b, rhs_er)

  end subroutine rhs

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! WRITE_MESH_SURFACE
 !!
 !! This diagnostics routine writes a GMV graphics file of the surface from
 !! the heat conduction mesh that was identified as the enclosure surface.
 !! The displayed node numbers will be the external mesh node numbers (or by
 !! altering two line, the internal mesh node number).  The displayed face
 !! numbers will be the face numbers from the corresponding enclosure file.
 !! The per-process list of mesh face indices is passed in FACES, and the
 !! corresponding list of enclosure face indices is passed in EFACES.
 !!

!FIXME -- REWRITE TO OPERATE ON A UNSTR_MESH TYPE MESH
!  subroutine write_mesh_surface (file, mesh, faces, efaces)
!
!    use gmvwrite_c_binding
!    use parallel_communication
!    use unstr_mesh_type
!
!    character(len=*), intent(in) :: file
!    type(unstr_mesh), intent(in) :: mesh
!    integer, Intent(in) :: faces(:), efaces(:)
!
!    integer :: j, k, n, offset, bsize(nPE)
!    integer :: num_nodes, num_faces
!    logical,  allocatable :: tag(:)
!    integer,  allocatable :: nodes(:), map_l(:), fnode(:,:)
!    integer,  pointer :: map(:)
!    real(r8), pointer :: x(:), y(:), z(:)
!    character(len=8) :: name
!
!    ASSERT(size(faces) == size(efaces))
!    ASSERT(global_all(faces >= 1))
!    ASSERT(global_all(faces <= mesh%face_imap%onp_size))
!
!    !if (is_IOP) call gmvwrite_openfile_ir_f (file, 4, 8) ! bug with node ids
!    if (is_IOP) call gmvwrite_openfile_ir_ascii_f (file, 4, 8)
!
!    !! Tag the nodes belonging to the specified faces.
!    allocate(tag(mesh%nnode))
!    tag = .false.
!    do j = 1, size(faces)
!      tag(mesh%fnode(:,faces(j))) = .true.
!    end do
!    call mesh%node_imap%scatter_offp_or(tag)
!
!    !! Count the on-process surface nodes per process (BSIZE).
!    call collate (count(tag(:mesh%nnode_onP)), bsize)
!    call broadcast (bsize)
!
!    !! Create the local list of on-process surface node indices, and create
!    !! the local block of the inverse mapping from global mesh node indices
!    !! to global surface node indices: MAP_L(NODES(J)) = OFFSET + J where J
!    !! and OFFSET+J are the local and global surface node indices, resp.
!    allocate(nodes(bsize(this_PE)), map_l(mesh%nnode_onP))
!    map_l = 0
!    n = 0 ! local surface node index
!    offset = sum(bsize(:this_PE-1))
!    do j = 1, mesh%nnode_onP
!      if (tag(j)) then
!        n = n + 1
!        nodes(n) = j
!        map_l(j) = offset + n ! global surface node index
!      end if
!    end do
!    deallocate(tag)
!
!    !! Global mapping array.
!    call allocate_collated_array (map, mesh%node_imap%global_size)
!    call collate (map_l, map)
!    deallocate(map_l)
!
!    !! Write the node coordinate data.
!    num_nodes = global_sum(size(nodes))
!    call allocate_collated_array (x, num_nodes)
!    call allocate_collated_array (y, num_nodes)
!    call allocate_collated_array (z, num_nodes)
!    call collate (mesh%x(1,nodes), x)
!    call collate (mesh%x(2,nodes), y)
!    call collate (mesh%x(3,nodes), z)
!    if (is_IOP) call gmvwrite_node_data_f (num_nodes, x, y, z)
!    deallocate(x, y, z)
!
!    !! Collate the surface face node array, ...
!    num_faces = global_sum(size(faces))
!    allocate(fnode(size(mesh%fnode,dim=1),num_faces))
!    call collate (mesh%node_imap%global_index(mesh%fnode(:,faces)), fnode)
!    !! and remap mesh node numbers to surface node numbers.
!    if (is_IOP) then
!      do j = 1, size(fnode,dim=2)
!        do k = 1, size(fnode,dim=1)
!          fnode(k,j) = map(fnode(k,j))
!        end do
!      end do
!      INSIST(minval(fnode) >= 1 .and. maxval(fnode) <= num_nodes)
!    end if
!    deallocate(map)
!
!    !! Write the cell data.
!    if (is_IOP) then
!      call gmvwrite_cell_header_f (num_faces)
!      select case (size(fnode,dim=1))
!      case (3)
!        name = 'tri'
!      case (4)
!        name = 'quad'
!      case default
!        INSIST(.false.)
!      end select
!      do j = 1, num_faces
!        call gmvwrite_cell_type_f (name, size(fnode,dim=1), fnode(:,j))
!      end do
!    end if
!    deallocate (fnode)
!
!    !! Write mesh node numbers as the nodeids -- GMV uses these for display.
!    call allocate_collated_array (map, num_nodes)
!    !call collate (mesh%node_imap%global_index(nodes), map)  ! internal mesh node numbers
!    call collate (mesh%xnode(nodes), map) ! external mesh node numbers
!    if (is_IOP) call gmvwrite_nodeids_f (map)
!    deallocate (map, nodes)
!
!    !! Write the enclosure face indices as the cellids -- GMV uses these for display.
!    call allocate_collated_array (map, num_faces)
!    call collate (efaces, map)
!    if (is_IOP) call gmvwrite_cellids_f (map)
!    deallocate (map)
!
!    if (nPE > 1) then
!      !! Write the face partitioning info.
!      call allocate_collated_array (map, num_faces)
!      call collate (spread(this_PE, dim=1, ncopies=size(faces)), map)
!      if (is_IOP) then
!        call gmvwrite_flag_header_f ()
!        call gmvwrite_flag_name_f ('facepart', nPE, CELLDATA)
!        do j = 1, nPE
!          write(name,'(a,i0)') 'P', j
!          call gmvwrite_flag_subname_f (name)
!        end do
!        call gmvwrite_flag_data_f (CELLDATA, map)
!        call gmvwrite_flag_endflag_f ()
!      end if
!      deallocate(map)
!    end if
!
!    if (is_IOP) call gmvwrite_closefile_f
!
!  end subroutine write_mesh_surface

  subroutine write_encl_surface(file, encl)

    use rad_encl_gmv

    character(*), intent(in) :: file
    type(rad_encl), intent(in) :: encl

    call gmv_open(file)
    call gmv_write_enclosure(encl)
    call gmv_close

  end subroutine write_encl_surface

  subroutine update_moving_vf(this)
    use moving_vf_type
    class(rad_problem), intent(inout) :: this
    select type (vf => this%vf)
    type is (moving_vf)
      call vf%next_interval
    end select
  end subroutine

  subroutine add_moving_vf_events(this, eventq)
    use sim_event_queue_type
    class(rad_problem), intent(in) :: this
    type(sim_event_queue), intent(inout) :: eventq
    integer :: j
    select type (vf => this%vf)
    type is (moving_vf)
      do j = 1, size(vf%times)
        call eventq%add_event(vf%times(j), vf_event(DT_POLICY_FACTOR, 0.125_r8))
      end do
    end select
  end subroutine

  pure function vf_event_init_dt(this, dt_last, dt_next) result(dt)
    class(vf_event), intent(in) :: this
    real(r8), intent(in) :: dt_last, dt_next
    real(r8) :: dt
    select case (this%dt_policy)
    case (DT_POLICY_NONE, DT_POLICY_NEXT)
      dt = dt_next
    case (DT_POLICY_FACTOR)
      dt = this%c*dt_last
    case (DT_POLICY_VALUE)
      dt = this%c
    case default
      dt = 0.0_r8 ! should never be here
    end select
  end function

end module rad_problem_type
