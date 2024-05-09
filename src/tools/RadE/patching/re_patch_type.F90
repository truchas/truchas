!!
!! RE_PATCH_TYPE
!!
!! This module provides a derived type RE_PATCH, which encapsulates the data
!! describing enclosure radiation patches methods that operate on instances of
!! this type.
!!
!! This module also provides a procedure for parsing the PATCHES namelist. The
!! data is copied into a parameter list.
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 16 Aug 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  READ_PATCHES_NAMELIST(LUN, PARAMS, FOUND) reads the first instance of the
!!    PATCHES namelist from the file opened on logical unit LUN. The values
!!    read are returned in the parameter list PARAMS.
!!    FOUND is a logical scalar set to .TRUE. if the namelist is found, and
!!    .FALSE. otherwise. This is a collective procedure; the file is read on the
!!    I/O process and the data replicated to all other processes.  If an error
!!    occurs (I/O error or invalid data) a message is written and execution is
!!    halted.  The presence of this namelist serves to enable patching of the
!!    radiation enclosure, as well as defining the parameters of the selected
!!    patching algorithm.  If this namelist is not present, no patches will be
!!    generated.
!!
!!  The RE_PATCH type describes a clustering of the enclosure faces. It has the
!!  following type bound subroutines.
!!
!!  BCAST() broadcasts the contents of the object on process rank 1 to all other
!!    process ranks, replicating the patch data on all processes.
!!
!!  GENERATE_PATCHES(E, PARAMS) generates patches for the given enclosure with
!!    the given parameters.  E is the TYPE(ENCL) radiation enclosure that will
!!    be patched.  This is a collective procedure, but only for the purposes of
!!    error handling.  The patch data is stored on process rank 1; the object is
!!    ignored on all other processes.
!!
!!  READ_PATCH_DATA(PATH) initializes the object with the enclosure patch data
!!    read from the radiation enclosure dataset PATH.  If the dataset does not
!!    include patch data, then it is assumed that each face is a patch, and
!!    THIS%F2P_MAP will not be populated.  This is a collective procedure, but
!!    only for the purposes purposes of error handling.  The patch data is
!!    stored on process rank 1; the object is ignored on all other processes.
!!
!!  WRITE_PATCH_DATA(PATH) writes the patch data contained in the object to the
!!    radiation enclosure dataset PATH.  The dataset must already contain
!!    compatible enclosure data.  If the object has no patches (i.e. each face
!!    is a patch) then the face-to-patch map F2P_MAP will not be a field of the
!!    dataset.  This is a collective procedure, but only for the purposes
!!    purposes of error handling.  The patch data is read from process rank 1;
!!    the object is ignored on all other processes.
!!
!!  PATCH_COLORING(E) computes a coloring of the patches where adjacent patches
!!    have different colors.  E is the TYPE(ENCL) radiation enclosure whose
!!    patching is described by this object.  If the object has patch data, the
!!    returned array PCOLOR will be a rank-1 integer array of length NPATCH
!!    containing a valid patch coloring.  If the object has no patch data (i.e.
!!    each face is a patch) then PCOLOR will not be allocated.  This is a serial
!!    procedure (no communication) and must only be called on a process with an
!!    initialized object.
!!


#include "f90_assert.fpp"

module re_patch_type

  use,intrinsic :: iso_fortran_env, only: i8 => int64, r8 => real64
  use scl
  use re_encl_type
  use parameter_list_type
  implicit none
  private

  public :: read_patches_namelist

  !! Patch algorithms
  character(32), parameter :: PATCH_ALGORITHMS(6) = ['NONE ','VSA  ','VAC  ','PAVE ','METIS','FILE ']
  integer, parameter :: PATCH_ALG_NONE = 1
  integer, parameter :: PATCH_ALG_VSA = 2
  integer, parameter :: PATCH_ALG_VAC = 3
  integer, parameter :: PATCH_ALG_PAVE = 4
  integer, parameter :: PATCH_ALG_METIS = 5
  integer, parameter :: PATCH_ALG_FILE = 6

  !! Parameter defaults
  integer, parameter :: PATCH_ALGORITHM_DEFAULT = PATCH_ALG_METIS
  integer, parameter :: VERBOSITY_LEVEL_DEFAULT = 1
  real(r8), parameter :: MAX_ANGLE_DEFAULT = 20.0

  type, public :: re_patch
    integer :: npatch  ! Total number of enclosure patches
    integer :: nface   ! Total number of faces
    integer, allocatable :: f2p_map(:)     ! Face-to-patch map
    integer, allocatable :: global_ids(:)  ! Global patch IDs
    logical :: has_patches  ! False if each face is a patch, true otherwise.
  contains
    procedure, public :: generate_patches
    procedure, public :: read => re_patch_read
    procedure, public :: write => re_patch_write
    procedure, public :: bcast => re_patch_bcast
    procedure, public :: patch_coloring
    procedure, public :: patch_to_face_array
    procedure, private :: no_patches
    procedure, private :: file_patches
    procedure, private :: run_patches
  end type

contains

  subroutine read_patches_namelist(lun, params, found)

    use string_utilities, only: i_to_c, raise_case
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use re_utilities
    use vsa_patching_type, only: VSA_MAX_ITER_DEFAULT, VSA_MIN_DELTA_DEFAULT, &
      VSA_FACE_PATCH_RATIO_DEFAULT, VSA_MAX_PATCH_RADIUS_DEFAULT, VSA_NORMALIZE_DIST_DEFAULT
    use vac_patching_type, only: VAC_MERGE_LEVEL_DEFAULT, VAC_SPLIT_PATCH_SIZE_DEFAULT
    use metis_patching_type, only: METIS_FACE_PATCH_RATIO_DEFAULT, METIS_FACE_WEIGHT_DEFAULT, &
      METIS_EDGE_WEIGHT_DEFAULT

    integer, intent(in) :: lun
    type(parameter_list), intent(out) :: params
    logical, intent(out) :: found

    type(parameter_list), pointer :: plist
    logical :: is_IOP
    integer :: ios, stat, patch_alg
    character(9) :: string
    character(255) :: iom, patch_file

    !! The PATCHES namelist variables; user visible.
    character(32) :: patch_algorithm
    integer  :: verbosity_level
    real(r8) :: max_angle
    integer  :: vsa_max_iter
    real(r8) :: vsa_min_delta
    real(r8) :: vsa_face_patch_ratio
    real(r8) :: vsa_max_patch_radius
    logical  :: vsa_normalize_dist
    integer  :: vsa_random_seed
    integer  :: vac_merge_level
    integer  :: vac_split_patch_size
    integer  :: pave_merge_level
    integer  :: pave_split_patch_size
    integer  :: pave_random_seed
    real(r8) :: metis_face_patch_ratio
    logical :: metis_face_weight
    logical :: metis_edge_weight
    integer :: metis_ctype, metis_iptype, metis_objtype, metis_no2hop, &
      metis_contig, metis_minconn, metis_ufactor, metis_niter, &
      metis_ncuts, metis_seed, metis_dbglvl, metis_ptype
    namelist /patches/ patch_algorithm, verbosity_level, max_angle, &
      vsa_max_iter, vsa_min_delta, vsa_face_patch_ratio, vsa_max_patch_radius, &
        vsa_normalize_dist, vsa_random_seed, &
      vac_merge_level, vac_split_patch_size, &
      pave_merge_level, pave_split_patch_size, pave_random_seed, &
      metis_face_patch_ratio, metis_face_weight, metis_edge_weight, &
        metis_ctype, metis_iptype, metis_objtype, metis_no2hop, &
        metis_contig, metis_minconn, metis_ufactor, metis_niter, &
        metis_ncuts, metis_seed, metis_dbglvl, metis_ptype, &
      patch_file

    is_IOP = (scl_rank()==1)  ! process rank 1 does the reading
    call params%set_path('patches')  ! used when printing errors

    !! Seek to the first instance of the PATCH namelist.
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'PATCHES', found, iostat=ios)
    end if
    call scl_bcast(ios)
    if (ios /= 0) call re_halt('error reading input file: iostat=' // i_to_c(ios))

    !! If no namelist found, make no patches
    call scl_bcast(found)
    if (.not. found) then
      call params%set('patch-alg', PATCH_ALG_NONE)
      return
    end if

    !! Read the namelist, assigning default values first.
    call re_info('Reading PATCHES namelist ...')
    patch_algorithm = NULL_C
    verbosity_level = NULL_I
    max_angle = NULL_R

    vsa_max_iter = NULL_I
    vsa_min_delta = NULL_R
    vsa_face_patch_ratio = NULL_R
    vsa_max_patch_radius = NULL_R
    vsa_normalize_dist = VSA_NORMALIZE_DIST_DEFAULT
    vsa_random_seed = NULL_I

    vac_merge_level = NULL_I
    vac_split_patch_size = NULL_I
    pave_merge_level = NULL_I
    pave_split_patch_size = NULL_I
    pave_random_seed = NULL_I

    metis_face_patch_ratio = NULL_R
    metis_face_weight = METIS_FACE_WEIGHT_DEFAULT
    metis_edge_weight = METIS_EDGE_WEIGHT_DEFAULT
    metis_ctype = NULL_I
    metis_iptype = NULL_I
    metis_objtype = NULL_I
    metis_no2hop = NULL_I
    metis_contig = NULL_I
    metis_minconn = NULL_I
    metis_ufactor = NULL_I
    metis_niter = NULL_I
    metis_ncuts = NULL_I
    metis_seed = NULL_I
    metis_dbglvl = NULL_I
    metis_ptype = NULL_I

    patch_file = NULL_C

    if (is_IOP) read(lun,nml=patches,iostat=ios,iomsg=iom)
    call scl_bcast(ios)
    if (ios /= 0) call re_halt('error reading PATCHES namelist: ' // trim(iom))

    !! Replicate the namelist variables on all processes.
    call scl_bcast(patch_algorithm)
    call scl_bcast(verbosity_level)
    call scl_bcast(max_angle)

    call scl_bcast(vsa_max_iter)
    call scl_bcast(vsa_min_delta)
    call scl_bcast(vsa_face_patch_ratio)
    call scl_bcast(vsa_max_patch_radius)
    call scl_bcast(vsa_normalize_dist)
    call scl_bcast(vsa_random_seed)

    call scl_bcast(vac_merge_level)
    call scl_bcast(vac_split_patch_size)
    call scl_bcast(pave_merge_level)
    call scl_bcast(pave_split_patch_size)
    call scl_bcast(pave_random_seed)

    call scl_bcast(metis_face_patch_ratio)
    call scl_bcast(metis_face_weight)
    call scl_bcast(metis_edge_weight)
    call scl_bcast(metis_ctype)
    call scl_bcast(metis_iptype)
    call scl_bcast(metis_objtype)
    call scl_bcast(metis_no2hop)
    call scl_bcast(metis_contig)
    call scl_bcast(metis_minconn)
    call scl_bcast(metis_ufactor)
    call scl_bcast(metis_niter)
    call scl_bcast(metis_ncuts)
    call scl_bcast(metis_seed)
    call scl_bcast(metis_dbglvl)
    call scl_bcast(metis_ptype)

    call scl_bcast(patch_file)

    !! Check the values
    stat = 0

    !! Select patch algorithm
    patch_algorithm = raise_case(trim(patch_algorithm))
    select case (patch_algorithm)
      case (NULL_C)
        patch_alg = PATCH_ALGORITHM_DEFAULT
        patch_algorithm = PATCH_ALGORITHMS(PATCH_ALGORITHM_DEFAULT)
        call re_info('  using default PATCH_ALGORITHM='//trim(patch_algorithm))
      case ('NONE')
        patch_alg = PATCH_ALG_NONE
      case ('VSA')
        patch_alg = PATCH_ALG_VSA
      case ('VAC')
        patch_alg = PATCH_ALG_VAC
      case ('PAVE')
        patch_alg = PATCH_ALG_PAVE
      case ('METIS')
        patch_alg = PATCH_ALG_METIS
      case ('FILE')
        patch_alg = PATCH_ALG_FILE
      case default
        call data_err('unrecognized PATCH_ALGORITHM: '//trim(patch_algorithm))
    end select
    call params%set('patch-alg', patch_alg)

    !! General settings
    if (patch_alg /= PATCH_ALG_NONE) then
      if (verbosity_level == NULL_I) then
        verbosity_level = VERBOSITY_LEVEL_DEFAULT
        call re_info('  using default VERBOSITY_LEVEL='//i_to_c(verbosity_level))
      else if (verbosity_level < 0) then
        call data_err('VERBOSITY_LEVEL must be >= 0')
      end if
      call params%set('verbosity-level', verbosity_level)
      if (max_angle == NULL_R) then
        max_angle = MAX_ANGLE_DEFAULT
        write(string,fmt='(f6.2)') max_angle
        call re_info('  using default MAX_ANGLE='//trim(string))
      else if (max_angle < 0 .or. max_angle > 180) then
        call data_err('MAX_ANGLE must be >= 0 and <= 180')
      end if
      call params%set('max-angle', max_angle)
    end if

    !! VSA settings
    if (patch_alg == PATCH_ALG_VSA) then
      if (vsa_max_iter == NULL_I) then
        vsa_max_iter = VSA_MAX_ITER_DEFAULT
        call re_info('  using default VSA_MAX_ITER='//i_to_c(vsa_max_iter))
      else if (vsa_max_iter < 1) then
        call data_err('VSA_MAX_ITER must be >= 1')
      end if
      call params%set('max-iter', vsa_max_iter)
      if (vsa_min_delta == NULL_R) then
        vsa_min_delta = VSA_MIN_DELTA_DEFAULT
        write(string,fmt='(es9.2)') vsa_min_delta
        call re_info('  using default VSA_MIN_DELTA='//string)
      else if (vsa_min_delta < 0.0_r8) then
        call data_err('VSA_MIN_DELTA must be >= 0')
      end if
      call params%set('min-delta', vsa_min_delta)
      if (vsa_face_patch_ratio == NULL_R) then
        vsa_face_patch_ratio = VSA_FACE_PATCH_RATIO_DEFAULT
        write(string,fmt='(es9.2)') vsa_face_patch_ratio
        call re_info('  using default VSA_FACE_PATCH_RATIO='//string)
      else if (vsa_face_patch_ratio < 1.0_r8) then
        call data_err('VSA_FACE_PATCH_RATIO must be >= 1')
      end if
      call params%set('face-patch-ratio', vsa_face_patch_ratio)
      if (vsa_max_patch_radius == NULL_R) then
        vsa_max_patch_radius = VSA_MAX_PATCH_RADIUS_DEFAULT
        write(string,fmt='(es9.2)') vsa_max_patch_radius
        call re_info('  using default VSA_MAX_PATCH_RADIUS='//string)
      else if (vsa_max_patch_radius <= 0.0_r8) then
        call data_err('VSA_MAX_PATCH_RADIUS must be > 0')
      end if
      call params%set('max-patch-radius', vsa_max_patch_radius)
      if (vsa_normalize_dist .eqv. VSA_NORMALIZE_DIST_DEFAULT) then
        write(string,fmt='(l1)') vsa_normalize_dist
        call re_info('  using default VSA_NORMALIZE_DIST='//trim(string))
      end if
      call params%set('normalize-dist', vsa_normalize_dist)
      if (vsa_random_seed == NULL_I) then
        call re_info('  using default VSA_RANDOM_SEED')
      else if (vsa_random_seed < 0) then
        call data_err('VSA_RANDOM_SEED must be >= 0')
      else
        !! Optional, only set if specified by user and valid
        call params%set('random-seed', vsa_random_seed)
      end if
    end if

    !! VAC settings
    if (patch_alg == PATCH_ALG_VAC) then
      if (vac_merge_level == NULL_I) then
        vac_merge_level = VAC_MERGE_LEVEL_DEFAULT
        call re_info('  using default VAC_MERGE_LEVEL='//i_to_c(vac_merge_level))
      else if (vac_merge_level < 0) then
        call data_err('VAC_MERGE_LEVEL must be >= 0')
      end if
      call params%set('merge-level', vac_merge_level)
      if (vac_split_patch_size == NULL_I) then
        vac_split_patch_size = VAC_SPLIT_PATCH_SIZE_DEFAULT
        call re_info('  using default VAC_SPLIT_PATCH_SIZE='//i_to_c(vac_split_patch_size))
      else if (vac_split_patch_size < 0) then
        call data_err('VAC_SPLIT_PATCH_SIZE must be >= 0')
      end if
      call params%set('split-patch-size', vac_split_patch_size)
    end if

    !! PAVE settings
    if (patch_alg == PATCH_ALG_PAVE) then
      if (pave_merge_level == NULL_I) then
        pave_merge_level = VAC_MERGE_LEVEL_DEFAULT
        call re_info('  using default PAVE_MERGE_LEVEL='//i_to_c(pave_merge_level))
      else if (pave_merge_level < 0) then
        call data_err('PAVE_MERGE_LEVEL must be >= 0')
      end if
      call params%set('merge-level', pave_merge_level)
      if (pave_split_patch_size == NULL_I) then
        pave_split_patch_size = VAC_SPLIT_PATCH_SIZE_DEFAULT
        call re_info('  using default PAVE_SPLIT_PATCH_SIZE='//i_to_c(pave_split_patch_size))
      else if (pave_split_patch_size < 0) then
        call data_err('PAVE_SPLIT_PATCH_SIZE must be >= 0')
      end if
      call params%set('split-patch-size', pave_split_patch_size)
      if (pave_random_seed == NULL_I) then
        call re_info('  using default PAVE_RANDOM_SEED')
      else if (pave_random_seed < 0) then
        call data_err('PAVE_RANDOM_SEED must be >= 0')
      else
        !! Optional, only set if specified by user and valid
        call params%set('random-seed', pave_random_seed)
      end if
    end if

    !! METIS settings
    if (patch_alg == PATCH_ALG_METIS) then
      if (metis_face_patch_ratio == NULL_R) then
        metis_face_patch_ratio = METIS_FACE_PATCH_RATIO_DEFAULT
        write(string,fmt='(es9.2)') metis_face_patch_ratio
        call re_info('  using default METIS_FACE_PATCH_RATIO='//string)
      else if (metis_face_patch_ratio < 1.0_r8) then
        call data_err('METIS_FACE_PATCH_RATIO must be >= 1')
      end if
      call params%set('face-patch-ratio', metis_face_patch_ratio)
      if (metis_face_weight .eqv. METIS_FACE_WEIGHT_DEFAULT) then
        write(string,fmt='(l1)') metis_face_weight
        call re_info('  using default METIS_FACE_WEIGHT='//trim(string))
      end if
      call params%set('face-weight', metis_face_weight)
      if (metis_edge_weight .eqv. METIS_EDGE_WEIGHT_DEFAULT) then
        write(string,fmt='(l1)') metis_edge_weight
        call re_info('  using default METIS_EDGE_WEIGHT='//trim(string))
      end if
      call params%set('edge-weight', metis_edge_weight)
      !! METIS library options
      plist => params%sublist('metis-options')
      if (metis_ptype   /= NULL_I) call plist%set('ptype', metis_ptype)
      if (metis_ctype   /= NULL_I) call plist%set('ctype', metis_ctype)
      if (metis_iptype  /= NULL_I) call plist%set('iptype', metis_iptype)
      if (metis_objtype /= NULL_I) call plist%set('objtype', metis_objtype)
      if (metis_no2hop  /= NULL_I) call plist%set('no2hop', metis_no2hop)
      if (metis_contig  /= NULL_I) call plist%set('contig', metis_contig)
      if (metis_minconn /= NULL_I) call plist%set('minconn', metis_minconn)
      if (metis_ufactor /= NULL_I) call plist%set('ufactor', metis_ufactor)
      if (metis_niter   /= NULL_I) call plist%set('niter', metis_niter)
      if (metis_ncuts   /= NULL_I) call plist%set('ncuts', metis_ncuts)
      if (metis_seed    /= NULL_I) call plist%set('seed', metis_seed)
      if (metis_dbglvl  /= NULL_I) call plist%set('dbglvl', metis_dbglvl)
    end if

    !! FILE settings
    if (patch_alg == PATCH_ALG_FILE) then
      if (patch_file == NULL_C) call data_err('PATCH_FILE not specified')
      call params%set('patch-file', trim(patch_file))
    end if

    if (stat /= 0) call re_halt('errors found in PATCHES namelist variables')

  contains
    subroutine data_err (errmsg)
      character(len=*), intent(in) :: errmsg
      stat = 1
      call re_info ('  ERROR: ' // errmsg)
    end subroutine data_err
  end subroutine read_patches_namelist


  !! Select and run a patching algorithm
  subroutine generate_patches (this, e, params, stat, errmsg)

    class(re_patch), intent(out) :: this
    class(encl), intent(in)  :: e
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: patch_alg

    stat = 0
    if (scl_rank() == 1) then
      call params%get('patch-alg', patch_alg)
      select case (patch_alg)
      case (PATCH_ALG_NONE)
        write (*,'(a)') 'No patches will be generated'
        call this%no_patches(e)
      case (PATCH_ALG_FILE)
        write (*,'(a)') 'Reading patches from a disk file'
        call this%file_patches(e, params, stat, errmsg)
      case default
        call this%run_patches(e, params, stat, errmsg)
      end select
    end if
    call scl_bcast(stat)

  end subroutine generate_patches


  !! Generate the identity face-to-patch map
  subroutine no_patches (this, e)

    class(re_patch), intent(out) :: this
    class(encl), intent(in)  :: e

    integer :: i

    this%has_patches = .false.
    this%npatch = e%nface
    this%nface = e%nface

    allocate(this%f2p_map(e%nface), this%global_ids(e%nface))
    do i = 1, e%nface
      this%f2p_map(i) = i
      this%global_ids(i) = i
    end do

  end subroutine no_patches


  !! Read patches from a file
  subroutine file_patches (this, e, params, stat, errmsg)
    class(re_patch), intent(out) :: this
    class(encl), target, intent(in) :: e
    character(:), allocatable :: path
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    stat = 0
    call params%get('patch-file', path)
    call re_patch_read(this, path)
    if (this%nface /= e%nface) then
      stat = 1
      errmsg = "face count in patch file '"//trim(path)//"' does not match input mesh"
    end if
  end subroutine file_patches


  !! Generate patches using one of the patching algorithms
  subroutine run_patches (this, e, params, stat, errmsg)

    use patching_class
    use vsa_patching_type
    use vac_patching_type
    use metis_patching_type

    class(re_patch), intent(out) :: this
    class(encl), target, intent(in)  :: e
    type(parameter_list), intent(inout)  :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    class(patching), allocatable :: patch
    integer :: patch_alg

    !! Choose algorithm
    call params%get('patch-alg', patch_alg)
    select case (patch_alg)
    case (PATCH_ALG_VSA)
      write (*,'(a)') 'Generating patches using the VSA algorithm'
      allocate(vsa_patching :: patch)
    case (PATCH_ALG_VAC)
      write (*,'(a)') 'Generating patches using the VAC algorithm'
      allocate(vac_patching :: patch)
    case (PATCH_ALG_PAVE)
      write (*,'(a)') 'Generating patches using the PAVE algorithm'
      allocate(pave_patching :: patch)
    case (PATCH_ALG_METIS)
      write (*,'(a)') 'Generating patches using the METIS graph partitioner'
      allocate(metis_patching :: patch)
    case default
      stat = 1
      errmsg = 'no such patching algorithm'
    end select
    if (stat/=0) return

    !! Run patching algorithm
    this%has_patches = .true.
    this%nface = e%nface

    call patch%init(e, params, stat, errmsg)
    if (stat/=0) return

    call patch%run(stat, errmsg)
    if (stat/=0) return

    call patch%output(this%f2p_map, this%global_ids, this%npatch)

  end subroutine run_patches


  subroutine re_patch_bcast (this)
    class(re_patch), intent(inout) :: this
    call scl_bcast(this%npatch)
    call scl_bcast(this%nface)
    call scl_bcast_alloc(this%f2p_map)
    call scl_bcast_alloc(this%global_ids)
    call scl_bcast(this%has_patches)
  end subroutine re_patch_bcast


  subroutine re_patch_read (this, path)

    use rad_encl_file_type

    class(re_patch), intent(out) :: this
    character(*), intent(in) :: path

    type(rad_encl_file) :: file
    integer :: i, nface, npatch

    if (scl_rank() == 1) then
      call file%open_ro(path)
      call file%get_patch_dims(nface, npatch)
      this%npatch = npatch
      this%nface = nface
      this%has_patches = file%has_patches()

      allocate(this%global_ids(npatch))
      do i = 1, this%npatch
        this%global_ids(i) = i
      end do

      allocate(this%f2p_map(nface))
      if (this%has_patches) then
        call file%get_f2p_map(this%f2p_map)
      else
        ASSERT(nface==npatch)
        this%f2p_map = this%global_ids
      end if
    end if

  end subroutine re_patch_read


  subroutine re_patch_write (this, path)

    use rad_encl_file_type

    class(re_patch), intent(in) :: this
    character(*), intent(in) :: path

    type(rad_encl_file) :: file

    if (scl_rank() == 1) then
      call file%open_rw(path)
      call file%init_patch(INT(this%npatch,kind=i8), this%has_patches)
      if (this%has_patches) call file%put_f2p_map(this%f2p_map)
      call file%close
    end if

  end subroutine re_patch_write


  !! Color patches so that adjacent patches have different colors.
  !! The procedure does nothing if this%has_patches is false.
  subroutine patch_coloring (this, e, pcolor, stat, errmsg)

    use patching_tools, only: get_face_neighbor_array
    use graph_type

    class(re_patch), intent(in) :: this
    type(encl), intent(in) :: e
    integer, allocatable, intent(out) :: pcolor(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(graph) :: pgraph  ! Patch adjacency graph
    integer, allocatable :: xfnhbr(:), fnhbr(:)
    integer :: p1, p2
    integer :: i, f, n, ncolor

    if (.not. this%has_patches) return

    call pgraph%init(this%npatch)
    call get_face_neighbor_array(e%xface, e%fnode, xfnhbr, fnhbr, stat, errmsg)
    if (stat /= 0) then
      block
        integer :: n
        character(40) :: coord
        n = findloc(fnhbr, -1, dim=1) ! first side with invalid topology
        write(coord,'("(",es12.5,2(",",es12.5),")")') e%side_location(n)
        errmsg = errmsg // ' at ' // coord
      end block
      return
    end if

    !! Construct patch adjacency graph
    do f = 1, e%nface
      p1 = this%f2p_map(f)
      do i = xfnhbr(f), xfnhbr(f+1)-1
        n = fnhbr(i)  ! neighbor of f
        if (n <= 0) cycle
        p2 = this%f2p_map(n)
        call pgraph%add_edge(p1, p2)
      end do
    end do

    call pgraph%vertex_coloring(pcolor, ncolor)

  end subroutine patch_coloring


  !! Expands a patch-length array into a face-length array
  subroutine patch_to_face_array (this, pvec, fvec)

    class(re_patch), intent(in) :: this
    real, intent(in)  :: pvec(:)
    real, intent(out) :: fvec(:)

    ASSERT(size(pvec) == this%npatch)

    if (this%has_patches) then
      fvec = pvec(this%f2p_map)
    else
      fvec = pvec
    end if

  end subroutine patch_to_face_array


end module re_patch_type
