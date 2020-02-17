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
!! IMPORTANT: RE_PATCH methods are serial, and must only be run on the I/O
!!            process.  Adhering to convention, READ_PATCHES_NAMELIST is
!!            collective.  Each rank gets a full copy of the PATCH_PARAM data
!!            structure, even though it will only be used by the single process
!!            running the RE_PATCH type methods.
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
!!  READ_PATCHES_NAMELIST (LUN, PARAMS, FOUND) reads the first instance of the
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
!!  GENERATE_PATCHES(THIS, E, PARAMS) generates patches for the given enclosure
!!    with the given parameters.  THIS is a TYPE(RE_PATCH) object that will
!!    store the patch data internally.  E is the TYPE(ENCL) radiation enclosure
!!    that will be patched. This is a serial procedure and must only be called
!!    by the I/O process.
!!
!!  READ_PATCH_DATA (THIS, PATH) initializes the TYPE(RE_PATCH) object THIS with
!!    the enclosure patch data read from the radiation enclosure dataset PATH.
!!    If the dataset does not include patch data, then it is assumed that each
!!    face is a patch, and THIS%F2P_MAP will not be populated.  This is a serial
!!    procedure and must only be called by the I/O process.
!!
!!  WRITE_PATCH_DATA (THIS, PATH) writes the patch data contained in the
!!    TYPE(RE_PATCH) object THIS to the radiation enclosure dataset PATH.  The
!!    dataset must already contain compatible enclosure data.  If THIS has no
!!    patches (i.e. each face is a patch) then the face-to-patch map F2P_MAP
!!    will not be a field of the dataset.  This is a serial procedure and must
!!    only be called by the I/O process.
!!
!!  PATCH_COLORING (THIS, E) computes a coloring of the patches where adjacent
!!    patches have different colors.  THIS is the TYPE(RE_PATCH) object storing
!!    the patch data corresponding to the TYPE(ENCL) radiation enclosure E.  If
!!    THIS has patch data, the returned array PCOLOR will be a rank-1 integer
!!    array of length THIS%NPATCH containing a valid patch coloring.  If THIS
!!    has no patch data (i.e. each face is a patch) then PCOLOR will not be
!!    allocated.  This is a serial procedure and must only be called by the
!!    I/O process.
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
  character(32), parameter :: PATCH_ALGORITHMS(4) = ['NONE','VSA ','VAC ','PAVE']
  integer, parameter :: PATCH_ALG_NONE = 1
  integer, parameter :: PATCH_ALG_VSA = 2
  integer, parameter :: PATCH_ALG_VAC = 3
  integer, parameter :: PATCH_ALG_PAVE = 4

  !! Parameter defaults
  integer, parameter :: PATCH_ALGORITHM_DEFAULT = PATCH_ALG_PAVE
  integer, parameter :: VERBOSITY_LEVEL_DEFAULT = 1
  real(r8), parameter :: MAX_ANGLE_DEFAULT = 20.0

  type, public :: re_patch
    integer :: npatch  ! Total number of enclosure patches
    integer, allocatable :: f2p_map(:)     ! Face-to-patch map
    integer, allocatable :: global_ids(:)  ! Global patch IDs
    logical :: has_patches  ! False if each face is a patch, true otherwise.
  contains
    procedure, public :: generate_patches
    procedure, public :: read_patch_data
    procedure, public :: write_patch_data
    procedure, public :: patch_coloring
    procedure, public :: patch_to_face_array
    procedure, private :: no_patches
    procedure, private :: vsa_patches
    procedure, private :: vac_patches
    procedure, private :: pave_patches
  end type

contains

  subroutine read_patches_namelist(lun, params, found)

    use string_utilities, only: i_to_c, raise_case
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I, NULL_R
    use re_utilities
    use vsa_patching_type, only: VSA_MAX_ITER_DEFAULT, VSA_MIN_DELTA_DEFAULT, &
      VSA_AVG_FACES_PER_PATCH_DEFAULT
    use vac_patching_type, only: VAC_MERGE_LEVEL_DEFAULT, VAC_SPLIT_PATCH_SIZE_DEFAULT

    integer, intent(in) :: lun
    type(parameter_list), intent(out) :: params
    logical, intent(out) :: found

    logical :: is_IOP
    integer :: ios, stat, patch_alg
    character(9) :: string
    character(255) :: iom

    !! The PATCHES namelist variables; user visible.
    character(32) :: patch_algorithm
    integer  :: verbosity_level
    real(r8) :: max_angle
    integer  :: vsa_max_iter
    real(r8) :: vsa_min_delta
    real(r8) :: vsa_avg_faces_per_patch
    integer  :: vac_merge_level
    integer  :: vac_split_patch_size
    integer  :: pave_merge_level
    integer  :: pave_split_patch_size
    namelist /patches/ patch_algorithm, verbosity_level, max_angle, &
      vsa_max_iter, vsa_min_delta, vsa_avg_faces_per_patch, &
      vac_merge_level, vac_split_patch_size, &
      pave_merge_level, pave_split_patch_size

    is_IOP = (scl_rank()==1)  ! process rank 1 does the reading

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
    vsa_avg_faces_per_patch = NULL_R
    vac_merge_level = NULL_I
    vac_split_patch_size = NULL_I
    pave_merge_level = NULL_I
    pave_split_patch_size = NULL_I

    if (is_IOP) read(lun,nml=patches,iostat=ios,iomsg=iom)
    call scl_bcast(ios)
    if (ios /= 0) call re_halt('error reading PATCHES namelist: ' // trim(iom))

    !! Replicate the namelist variables on all processes.
    call scl_bcast(patch_algorithm)
    call scl_bcast(verbosity_level)
    call scl_bcast(max_angle)
    call scl_bcast(vsa_max_iter)
    call scl_bcast(vsa_min_delta)
    call scl_bcast(vsa_avg_faces_per_patch)
    call scl_bcast(vac_merge_level)
    call scl_bcast(vac_split_patch_size)
    call scl_bcast(pave_merge_level)
    call scl_bcast(pave_split_patch_size)

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
      call params%set('vsa-max-iter', vsa_max_iter)
      if (vsa_min_delta == NULL_R) then
        vsa_min_delta = VSA_MIN_DELTA_DEFAULT
        write(string,fmt='(es9.2)') vsa_min_delta
        call re_info('  using default VSA_MIN_DELTA='//string)
      else if (vsa_min_delta < 0.0_r8) then
        call data_err('VSA_MIN_DELTA must be >= 0')
      end if
      call params%set('vsa-min-delta', vsa_min_delta)
      if (vsa_avg_faces_per_patch == NULL_R) then
        vsa_avg_faces_per_patch = VSA_AVG_FACES_PER_PATCH_DEFAULT
        write(string,fmt='(es9.2)') vsa_avg_faces_per_patch
        call re_info('  using default VSA_AVG_FACES_PER_PATCH='//string)
      else if (vsa_avg_faces_per_patch < 1.0_r8) then
        call data_err('VSA_AVG_FACES_PER_PATCH must be >= 1')
      end if
      call params%set('vsa-avg-faces-per-patch', vsa_avg_faces_per_patch)
    end if

    !! VAC settings
    if (patch_alg == PATCH_ALG_VAC) then
      if (vac_merge_level == NULL_I) then
        vac_merge_level = VAC_MERGE_LEVEL_DEFAULT
        call re_info('  using default VAC_MERGE_LEVEL='//i_to_c(vac_merge_level))
      else if (vac_merge_level < 0) then
        call data_err('VAC_MERGE_LEVEL must be >= 0')
      end if
      call params%set('vac-merge-level', vac_merge_level)
      if (vac_split_patch_size == NULL_I) then
        vac_split_patch_size = VAC_SPLIT_PATCH_SIZE_DEFAULT
        call re_info('  using default VAC_SPLIT_PATCH_SIZE='//i_to_c(vac_split_patch_size))
      else if (vac_split_patch_size < 0) then
        call data_err('VAC_SPLIT_PATCH_SIZE must be >= 0')
      end if
      call params%set('vac-split-patch-size', vac_split_patch_size)
    end if

    !! PAVE settings
    if (patch_alg == PATCH_ALG_PAVE) then
      if (pave_merge_level == NULL_I) then
        pave_merge_level = VAC_MERGE_LEVEL_DEFAULT
        call re_info('  using default PAVE_MERGE_LEVEL='//i_to_c(pave_merge_level))
      else if (pave_merge_level < 0) then
        call data_err('PAVE_MERGE_LEVEL must be >= 0')
      end if
      call params%set('pave-merge-level', pave_merge_level)
      if (pave_split_patch_size == NULL_I) then
        pave_split_patch_size = VAC_SPLIT_PATCH_SIZE_DEFAULT
        call re_info('  using default PAVE_SPLIT_PATCH_SIZE='//i_to_c(pave_split_patch_size))
      else if (pave_split_patch_size < 0) then
        call data_err('PAVE_SPLIT_PATCH_SIZE must be >= 0')
      end if
      call params%set('pave-split-patch-size', pave_split_patch_size)
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
  subroutine generate_patches (this, e, params)

    class(re_patch), intent(out) :: this
    class(encl), intent(in)  :: e
    type(parameter_list), intent(inout) :: params

    integer :: patch_alg

    if (scl_rank() == 1) then
      call params%get('patch-alg', patch_alg)
      select case (patch_alg)
      case (PATCH_ALG_NONE)
        write (*,'(a)') 'No patches will be generated'
        call this%no_patches(e)
      case (PATCH_ALG_VSA)
        write (*,'(a)') 'Generating patches using the VSA algorithm'
        call this%vsa_patches(e, params)
      case (PATCH_ALG_VAC)
        write (*,'(a)') 'Generating patches using the VAC algorithm'
        call this%vac_patches(e, params)
      case (PATCH_ALG_PAVE)
        write (*,'(a)') 'Generating patches using the VAC PAVE algorithm'
        call this%pave_patches(e, params)
      case default
        !! Programming error, exit immediately.
        INSIST(.false.)
      end select
    end if

  end subroutine generate_patches


  !! Generate the identity face-to-patch map
  subroutine no_patches (this, e)

    class(re_patch), intent(out) :: this
    class(encl), intent(in)  :: e

    integer :: i

    this%has_patches = .false.
    this%npatch = e%nface

    allocate(this%f2p_map(e%nface), this%global_ids(e%nface))

    do i = 1, e%nface
      this%f2p_map(i) = i
      this%global_ids(i) = i
    end do

  end subroutine no_patches


  !! Generate patches using the VSA algorithm
  subroutine vsa_patches (this, e, params)

    use vsa_patching_type

    class(re_patch), intent(out) :: this
    class(encl), target, intent(in)  :: e
    type(parameter_list), intent(inout)  :: params

    type(vsa_patching) :: vsa
    integer :: verbosity, max_iter
    real(r8) :: max_angle, min_delta, avg_fpp

    call params%get('verbosity-level', verbosity)
    call params%get('vsa-avg-faces-per-patch', avg_fpp)
    call params%get('max-angle', max_angle)
    call params%get('vsa-min-delta', min_delta)
    call params%get('vsa-max-iter', max_iter)

    this%has_patches = .true.

    call vsa%init(e, avg_fpp, max_angle, verbosity)

    call vsa%run(min_delta, max_iter)

    call vsa%output(this%f2p_map, this%global_ids, this%npatch)

  end subroutine vsa_patches


  !! Generate patches using the VAC algorithm
  subroutine vac_patches (this, e, params)

    use vac_patching_type

    class(re_patch), intent(out) :: this
    class(encl), target, intent(in)  :: e
    type(parameter_list), intent(inout)  :: params

    type(vac_patching) :: vac
    integer :: verbosity, merge_level, split_patch_size
    real(r8) :: max_angle

    call params%get('verbosity-level', verbosity)
    call params%get('max-angle', max_angle)
    call params%get('vac-merge-level', merge_level)
    call params%get('vac-split-patch-size', split_patch_size)

    this%has_patches = .true.

    call vac%init(e, max_angle, merge_level, split_patch_size, verbosity)

    call vac%run()

    call vac%output(this%f2p_map, this%global_ids, this%npatch)

  end subroutine vac_patches


  !! Generate patches using the VAC PAVE algorithm
  subroutine pave_patches (this, e, params)

    use vac_patching_type

    class(re_patch), intent(out) :: this
    class(encl), target, intent(in)  :: e
    type(parameter_list), intent(inout)  :: params

    type(pave_patching) :: pave
    integer :: verbosity, merge_level, split_patch_size
    real(r8) :: max_angle

    call params%get('verbosity-level', verbosity)
    call params%get('max-angle', max_angle)
    call params%get('pave-merge-level', merge_level)
    call params%get('pave-split-patch-size', split_patch_size)

    this%has_patches = .true.

    call pave%init(e, max_angle, merge_level, split_patch_size, verbosity)

    call pave%run()

    call pave%output(this%f2p_map, this%global_ids, this%npatch)

  end subroutine pave_patches


  subroutine read_patch_data (this, path)

    use rad_encl_file_type

    class(re_patch), intent(out) :: this
    character(*), intent(in) :: path

    type(rad_encl_file) :: file
    integer :: nface_tot, npatch_tot

    call file%open_ro(path)
    call file%get_patch_dims(nface_tot, npatch_tot)
    this%npatch = npatch_tot
    this%has_patches = file%has_patches()
    if (.not. this%has_patches) return

    allocate(this%f2p_map(nface_tot))
    call file%get_f2p_map(this%f2p_map)

  end subroutine read_patch_data


  subroutine write_patch_data (this, path)

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

  end subroutine write_patch_data


  !! Color patches so that adjacent patches have different colors.
  !! The procedure does nothing if this%has_patches is false.
  function patch_coloring (this, e) result(pcolor)

    use patching_tools, only: get_face_neighbor_array
    use graph_type

    class(re_patch), intent(in) :: this
    type(encl), intent(in) :: e
    integer, allocatable :: pcolor(:)

    type(graph) :: pgraph  ! Patch adjacency graph
    integer, allocatable :: xfnhbr(:), fnhbr(:)
    integer :: p1, p2
    integer :: i, f, n, ncolor, stat

    if (.not. this%has_patches) return

    call pgraph%init(this%npatch)
    call get_face_neighbor_array(e%xface, e%fnode, xfnhbr, fnhbr, stat)
    ASSERT( stat == 0 )

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

  end function patch_coloring


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
