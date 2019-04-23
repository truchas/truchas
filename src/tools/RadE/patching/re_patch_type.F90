!!
!! RE_PATCH_TYPE
!!
!! This module provides a derived type RE_PATCH, which encapsulates the data
!! describing enclosure radiation patches methods that operate on instances of
!! this type.
!!
!! This module also provides a procedure for parsing the PATCHES namelist. The
!! data is copied to the parameter list PATCH_PARAM for later use.
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
!!  READ_PATCHES_NAMELIST (LUN, PPAR, FOUND) reads the first instance of the
!!    PATCHES namelist from the file opened on logical unit LUN.  PPAR is the
!!    TYPE(PATCH_PARAM) data structure in which to store the namelist data.
!!    FOUND is a logical scalar set to .TRUE. if the namelist is found, and
!!    .FALSE. otherwise. This is a collective procedure; the file is read on the
!!    I/O process and the data replicated to all other processes.  If an error
!!    occurs (I/O error or invalid data) a message is written and execution is
!!    halted.  The presence of this namelist serves to enable patching of the
!!    radiation enclosure, as well as defining the parameters of the selected
!!    patching algorithm.  If this namelist is not present, no patches will be
!!    generated.
!!
!!  GENERATE_PATCHES (THIS, E, PPAR) generates patches for the given enclosure
!!    with the given parameters.  THIS is a TYPE(RE_PATCH) object that will
!!    store the patch data internally.  E is the TYPE(ENCL) radiation enclosure
!!    that will be patched.  PPAR is the TYPE(PATCH_PARAM) patch parameters that
!!    determine which patching algorithm to use, and what parameters to pass to
!!    the patching algorithm.  This is a serial procedure and must only be
!!    called by the I/O process.
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


#include "f90_assert.fpp"

module re_patch_type

  use kinds, only: i8
  use scl
  use re_encl_type
  implicit none
  private

  public :: read_patches_namelist

  !! Patch algorithms
  character(32), parameter :: PATCH_ALGORITHMS(1) = ['NONE']
  integer, parameter :: PATCH_ALG_NONE = 1

  !! Parameter defaults
  integer, parameter :: PATCH_ALGORITHM_DEFAULT = PATCH_ALG_NONE
  integer, parameter :: VERBOSITY_LEVEL_DEFAULT = 1

  type, public :: patch_param
    private
    integer  :: verbosity  ! Verbosity level
    integer  :: patch_alg  ! Patching algorithm to run
  end type

  type, public :: re_patch
    integer :: npatch  ! Total number of enclosure patches
    integer, allocatable :: f2p_map(:)     ! Face-to-patch map
    integer, allocatable :: global_ids(:)  ! Global patch IDs
    logical :: has_patches  ! False if each face is a patch, true otherwise.
  contains
    procedure, public :: generate_patches
    procedure, public :: read_patch_data
    procedure, public :: write_patch_data
    procedure, public :: patch_to_face_array
    procedure, private :: no_patches
  end type


contains


  subroutine read_patches_namelist (lun, ppar, found)

    use string_utilities, only: i_to_c, raise_case
    use input_utilities, only: seek_to_namelist, NULL_C, NULL_I
    use re_utilities

    integer, intent(in) :: lun
    type(patch_param), intent(out) :: ppar
    logical, intent(out) :: found

    logical :: is_IOP
    integer :: ios, stat, patch_alg

    !! The PATCHES namelist variables; user visible.
    character(32) :: patch_algorithm
    integer  :: verbosity_level
    namelist /patches/ patch_algorithm, verbosity_level

    is_IOP = (scl_rank()==1)  ! process rank 1 does the reading

    !! Seek to the first instance of the PATCH namelist.
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist (lun, 'PATCHES', found, iostat=ios)
    end if
    call scl_bcast (ios)
    if (ios /= 0) call re_halt ('error reading file connected to unit ' // &
                                i_to_c(lun) // ': iostat=' // i_to_c(ios))

    !! If no namelist found, make no patches
    call scl_bcast (found)
    if (.not. found) then
      ppar%patch_alg = PATCH_ALG_NONE
      return
    end if

    !! Read the namelist, assigning default values first.
    call re_info ('Reading PATCHES namelist ...')
    if (is_IOP) then
      patch_algorithm = NULL_C
      verbosity_level = NULL_I
      read(lun,nml=patches,iostat=ios)
    end if
    call scl_bcast (ios)
    if (ios /= 0) call re_halt ('Error reading PATCHES namelist: iostat=' // i_to_c(ios))

    !! Replicate the namelist variables on all processes.
    call scl_bcast (patch_algorithm)
    call scl_bcast (verbosity_level)

    !! Check the values
    stat = 0

    !! Select patch algorithm
    patch_algorithm = raise_case(trim(patch_algorithm))
    select case (patch_algorithm)
      case (NULL_C)
        patch_alg = PATCH_ALGORITHM_DEFAULT
        patch_algorithm = PATCH_ALGORITHMS(PATCH_ALGORITHM_DEFAULT)
        call re_info ('  using default PATCH_ALGORITHM='//trim(patch_algorithm))
      case ('NONE')
        patch_alg = PATCH_ALG_NONE
      case default
        call data_err ('Unrecognized PATCH_ALGORITHM: '//trim(patch_algorithm))
    end select

    !! General settings
    if (verbosity_level == NULL_I) then
      verbosity_level = VERBOSITY_LEVEL_DEFAULT
      call re_info ('  using default VERBOSITY_LEVEL='//i_to_c(verbosity_level))
    else if (verbosity_level < 0) then
      call data_err ('VERBOSITY_LEVEL must be >= 0')
    end if

    if (stat /= 0) call re_halt ('errors found in PATCHES namelist variables')

    !! Everything checks out; write values into return data structure
    ppar%patch_alg = patch_alg
    ppar%verbosity = verbosity_level

  contains
    subroutine data_err (errmsg)
      character(len=*), intent(in) :: errmsg
      stat = 1
      call re_info ('  ERROR: ' // errmsg)
    end subroutine data_err
  end subroutine read_patches_namelist


  !! Select and run a patching algorithm
  subroutine generate_patches (this, e, ppar)

    class(re_patch), intent(out) :: this
    type(encl), intent(in)  :: e
    type(patch_param), intent(in) :: ppar

    select case (ppar%patch_alg)
      case (PATCH_ALG_NONE)
        write (*,'(a)') 'No patches will be generated'
        call this%no_patches(e)
      case default
        !! Programming error, exit immediately.
        INSIST(.false.)
    end select

  end subroutine generate_patches


  !! Generate the identity face-to-patch map
  subroutine no_patches (this, e)

    class(re_patch), intent(out) :: this
    type(encl), intent(in)  :: e

    integer :: i

    this%has_patches = .false.
    this%npatch = e%nface

    allocate(this%f2p_map(e%nface), this%global_ids(e%nface))

    do i = 1, e%nface
      this%f2p_map(i) = i
      this%global_ids(i) = i
    end do

  end subroutine no_patches


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

    call file%open_rw(path)
    call file%init_patch(INT(this%npatch,kind=i8), this%has_patches)
    if (this%has_patches) call file%put_f2p_map(this%f2p_map)
    call file%close

  end subroutine write_patch_data


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
