!!
!! RE_CHAPARRAL_VF
!!
!! This module provides a method for computing view factors for a radiation
!! enclosure using Sandia's CHAPARRAL package.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 4 Apr 2008; updated Feb 2020
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! 6 Aug 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  CALL READ_CHAPARRAL_NAMELIST(LUN, PARAMS, FOUND) reads the first occurrence
!!    of a CHAPARRAL namelist from the file opened on unit LUN. If the namelist
!!    is not found, the logical argument FOUND returns false.  Otherwise it
!!    returns true and the values read (and defaults for any others) are
!!    returned in the parameter list PARAMS.  The subroutine takes care
!!    of handling any IO errors and checks the namelist values for correctness,
!!    gracefully terminating execution of the program if any errors are
!!    encountered.  This is a collective procedure.  Input takes place on
!!    process rank 1, but the returned PARAMS is replicated on all processes.
!!
!!    The CHAPARRAL namelist contains the following variables:
!!      verbosity_level
!!      blocking_enclosure
!!      partial_enclosure
!!      partial_area
!!      BSP_max_tree_depth
!!      BSP_min_leaf_length
!!      spatial_tolerance
!!      hemicube_resolution
!!      max_subdivisions
!!      min_separation
!!      smoothing_weight
!!      smoothing_tolerance
!!      smoothing_max_iter
!!
!!  CALL CALCULATE_VF(E, PARAMS, EP, VF) uses the Chaparral package to calculate
!!    the view factors VF for the enclosure E whose faces are clustered as
!!    described by the enclosure patches object EP.  The control parameters for
!!    the Chaparral procedures are contained in the parameter list PARAMS that
!!    was initialized by a call to READ_CHAPARRAL_NAMELIST. This is a collective
!!    procedure. The enclosure E and parameters PARAMS must be replicated across
!!    all processes, and the enclosure patches EP must only be initialized on
!!    process rank 1. The view factors are returned in a distributed form in VF.
!!

#include "f90_assert.fpp"

module re_chaparral_vf

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use re_utilities
  use re_dist_vf_type
  use parameter_list_type
  use scl
  implicit none
  private

  public :: read_chaparral_namelist, calculate_vf

contains

  subroutine read_chaparral_namelist(lun, params, found)

    use string_utilities
    use input_utilities

    integer, intent(in) :: lun
    type(parameter_list), intent(out) :: params
    logical, intent(out) :: found

    logical :: is_IOP
    integer :: ios, stat
    character(9) :: string
    character(255) :: iom

    !! The CHAPARRAL namelist variables; user visible.
    logical  :: blocking_enclosure, partial_enclosure
    integer  :: verbosity_level, BSP_max_tree_depth, BSP_min_leaf_length
    real(r8) :: partial_area, spatial_tolerance
    integer  :: hemicube_resolution, max_subdivisions, smoothing_max_iter
    real(r8) :: min_separation, smoothing_weight, smoothing_tolerance
    namelist /chaparral/ verbosity_level, blocking_enclosure, partial_enclosure, &
        partial_area, BSP_max_tree_depth, BSP_min_leaf_length, spatial_tolerance, &
        hemicube_resolution, max_subdivisions, min_separation, smoothing_weight, &
        smoothing_tolerance, smoothing_max_iter

    is_IOP = (scl_rank()==1)  ! process rank 1 does the reading

    !! Seek to the first instance of the CHAPARRAL namelist.
    if (is_IOP) then
      rewind(lun)
      call seek_to_namelist(lun, 'CHAPARRAL', found, iostat=ios)
    end if
    call scl_bcast(ios)
    if (ios /= 0) call re_halt('error reading input file: iostat=' // i_to_c(ios))

    !! This is an optional namelist.
    call scl_bcast(found)
    if (.not.found) return

    !! Read the namelist, assigning default values first.
    call re_info('Reading CHAPARRAL namelist ...')
    verbosity_level     = NULL_I
    blocking_enclosure  = .true.
    partial_enclosure   = .false.
    partial_area        = NULL_R
    BSP_max_tree_depth  = NULL_I
    BSP_min_leaf_length = NULL_I
    spatial_tolerance   = NULL_R
    hemicube_resolution = NULL_I
    max_subdivisions    = NULL_I
    min_separation      = NULL_R
    smoothing_weight    = NULL_R
    smoothing_tolerance = NULL_R
    smoothing_max_iter  = NULL_I
    if (is_IOP) read(lun,nml=chaparral,iostat=ios,iomsg=iom)
    call scl_bcast(ios)
    if (ios /= 0) call re_halt('error reading CHAPARRAL namelist: ' // trim(iom))

    !! Replicate the namelist variables on all processes.
    call scl_bcast(verbosity_level)
    call scl_bcast(blocking_enclosure)
    call scl_bcast(partial_enclosure)
    call scl_bcast(partial_area)
    call scl_bcast(BSP_max_tree_depth)
    call scl_bcast(BSP_min_leaf_length)
    call scl_bcast(spatial_tolerance)
    call scl_bcast(hemicube_resolution)
    call scl_bcast(max_subdivisions)
    call scl_bcast(min_separation)
    call scl_bcast(smoothing_weight)
    call scl_bcast(smoothing_tolerance)
    call scl_bcast(smoothing_max_iter)

    !! Check the values !!

    stat = 0

    if (verbosity_level == NULL_I) then
      verbosity_level = 2
      call re_info('  using default VERBOSITY_LEVEL='//i_to_c(verbosity_level))
    else if (verbosity_level < 0) then
      call data_err('VERBOSITY_LEVEL must be >= 0')
    end if

    if (partial_enclosure) then
      if (partial_area == NULL_R) then
        call data_err('PARTIAL_AREA must be assigned a positive value')
      else if (partial_area <= 0.0_r8) then
        call data_err('PARTIAL_AREA must be > 0.0')
      end if
    else if (partial_area /= NULL_R) then
      call re_info('Ignoring PARTIAL_AREA; not a partial enclosure.')
    end if

    if (BSP_max_tree_depth == NULL_I) then
      BSP_max_tree_depth = 15
      call re_info('  using default BSP_MAX_TREE_DEPTH='//i_to_c(BSP_max_tree_depth))
    else if (BSP_max_tree_depth < 1) then
      call data_err('BSP_MAX_TREE_DEPTH must be >= 1')
    end if

    if (BSP_min_leaf_length == NULL_I) then
      BSP_min_leaf_length = 25
      call re_info('  using default BSP_MAX_LEAF_LENGTH='//i_to_c(BSP_min_leaf_length))
    else if (BSP_min_leaf_length < 1) then
      call data_err('BSP_MIN_LEAF_LENGTH must be >= 1')
    end if

    if (spatial_tolerance == NULL_R) then
      call data_err('SPATIAL_TOLERANCE must be assigned a value')
    else if (spatial_tolerance <= 0.0) then
      call data_err('SPATIAL_TOLERANCE must be > 0.0')
    end if

    if (hemicube_resolution == NULL_I) then
      call data_err('HEMICUBE_RESOLUTION must be assigned a value')
    else if (hemicube_resolution < 4) then
      call data_err('HEMICUBE_RESOLUTION must be >= 4')
    end if

    if (max_subdivisions == NULL_I) then
      call data_err('MAX_SUBDIVISIONS must be assigned a value')
    else if (max_subdivisions < 0) then
      call data_err('MAXSUBDIVISIONS must be >= 0')
    end if

    if (min_separation == NULL_R) then
      call data_err('MIN_SEPARATION must be assigned a value')
    else if (min_separation < 0.0) then
      call data_err('MIN_SEPARATION must be >= 0.0')
    end if

    if (smoothing_weight == NULL_R) then
      smoothing_weight = 2.0_r8
      write(string,fmt='(es9.2)') smoothing_weight
      call re_info('  using default SMOOTHING_WEIGHT='//string)
    else if (smoothing_weight <= 0.0) then
      call data_err('SMOOTHING_WEIGHT must be > 0.0')
    end if

    if (smoothing_tolerance == NULL_R) then
      call data_err('SMOOTHING_TOLERANCE must be assigned a value')
    else if (smoothing_tolerance <= 0.0) then
      call data_err('SMOOTHING_TOLERANCE must be > 0.0')
    end if

    if (smoothing_max_iter == NULL_I) then
      call data_err('SMOOTHING_MAX_ITER must be assigned a value')
    else if (smoothing_max_iter < 0) then
      call data_err('SMOOTHING_MAX_ITER must be >= 0')
    end if

    if (stat /= 0) call re_halt('errors found in CHAPARRAL namelist variables')

    !! Everything checks out; stuff the values into the return parameter list.
    call params%set('verbosity', verbosity_level)
    call params%set('nonblocking', merge(0, 1, blocking_enclosure))
    call params%set('partial', merge(1, 0, partial_enclosure))
    call params%set('asink', partial_area)
    call params%set('bsp-depth', BSP_max_tree_depth)
    call params%set('bsp-length', BSP_min_leaf_length)
    call params%set('spatial-tol', spatial_tolerance)
    call params%set('hc-resolution', hemicube_resolution)
    call params%set('hc-sub-divide', max_subdivisions)
    call params%set('hc-min-sep', min_separation)
    call params%set('smooth-wt', smoothing_weight)
    call params%set('smooth-tol', smoothing_tolerance)
    call params%set('smooth-max-iter', smoothing_max_iter)

  contains

    subroutine data_err (errmsg)
      character(*), intent(in) :: errmsg
      stat = 1
      call re_info('  ERROR: ' // errmsg)
    end subroutine data_err

  end subroutine read_chaparral_namelist


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! CALCULATE_VF
 !!
 !! This subroutine does all the serious work of setting up the calls to
 !! Chaparral functions to carry out the computation of the view factors.
 !! Some comments about what is going on here:
 !!
 !! 1) In Chaparral a patch corresponds to a DoF in the radiosity system
 !!   (e.g., a temperature or flux value), and each patch can be composed
 !!   of multiple faces.  The RE_PATCH object EP stores a face-to-patch
 !!   map that defines the patches of the radiosity system.
 !!
 !! 2) Rather than assume the patches are numbered consecutively from 1 (or 0),
 !!   Chaparral requires the user to specify the "global id" of each patch.
 !!   The EP object also provides these global IDs, although currently they are
 !!   just a 1-based consecutive numbering of the patches.
 !!
 !! 3) For a partial enclosure, Chaparral requires an additional "virtual
 !!   surface" patch that completes the enclosure.  Per requirements we
 !!   assign this patch the largest global id, and thus it is the last one.
 !!   Note that this patch will not be associated with any faces.
 !!
 !! 4) The EP object provides the face-to-patch map required by Chaparral.
 !!
 !! 5) Internally Chaparral reorders and redistributes the patches among the
 !!   processes, so for the calculation of the VF it doesn't really matter
 !!   how we supply the patches to Chaparral; we could supply them all on
 !!   process rank 1 if we wanted.  However, when we extract the VF matrix
 !!   from Chaparral we receive the row for a patch on the process the patch
 !!   was supplied.  For this reason we equidistribute the patches across
 !!   processes in order to extract a reasonably-well distributed VF matrix
 !!   which can be quite large.
 !!

  subroutine calculate_vf(e, params, ep, evf)

    use chaparral_c_binding
    use re_encl_type
    use re_patch_type

    class(encl), intent(in), target :: e
    type(parameter_list), intent(inout) :: params
    type(re_patch), intent(in) :: ep
    type(dist_vf), intent(out) :: evf

    integer :: j, n, handle, nfacets, nrotations, max_surfaces, nproc, my_rank
    integer :: npatch_tot       ! Total number of real patches
    integer :: npatch           ! Number of real patches owned by this rank
    integer :: npatch_tot_chap  ! Total number of patches, including virtual surface patch
    integer :: npatch_chap      ! Number of patches, including virtual surface patch, owned by this rank
    integer :: xmirror, ymirror, zmirror, face_offset, patch_offset
    integer, allocatable :: global_ids_g(:)  ! Collated global patch IDs. Owned by rank 1.
    integer, allocatable :: f2p_map_g(:)     ! Collated face-to-patch map. Owned by rank 1.
    integer, allocatable :: global_ids(:), f2p_map(:)  ! Distributed versions of the collated arrays
    integer, allocatable :: c(:,:), nsizes(:)
    real(r8), allocatable :: x(:), y(:), z(:)
    real(r8), allocatable :: parea(:)  ! Patch areas computed by Chaparral
    real(r8), allocatable :: farea(:)  ! Face areas computed from the enclosure
    integer :: partial, nonblocking, verbosity, bsp_depth, bsp_length, hc_sub_divide, hc_resolution, smooth_max_iter
    real(r8) :: asink, spatial_tol, hc_min_sep, smooth_wt, smooth_tol

    !! Hardwired Chaparral parameters
    integer, parameter :: GEOM_TYPE = 3     ! 3D geometry
    integer, parameter :: SYM_METHOD = 3    ! averaging

    nproc = scl_size()
    my_rank = scl_rank()

    !! Get patch data
    if (my_rank == 1) then
      npatch_tot = ep%npatch
      allocate(global_ids_g, source=ep%global_ids)
      allocate(f2p_map_g, source=ep%f2p_map)
    else
      allocate(global_ids_g(0),f2p_map_g(0))
    end if
    call scl_bcast(npatch_tot)

    !! Divvy up the faces and patches.
    nfacets = e%nface/nproc
    npatch  = npatch_tot/nproc
    if (my_rank <= modulo(e%nface,nproc)) nfacets = nfacets + 1
    if (my_rank <= modulo(npatch_tot,nproc)) npatch = npatch + 1
    INSIST( scl_global_sum(nfacets) == e%nface )
    INSIST( scl_global_sum(npatch) == npatch_tot )

    !! Assign virtual surface patch to last rank
    npatch_chap = npatch
    call params%get('partial', partial)
    if (partial==1 .and. my_rank==nproc) npatch_chap = npatch_chap + 1

    allocate(nsizes(nproc))
    call scl_allgather(nfacets, nsizes)
    face_offset = sum(nsizes(1:my_rank-1))
    call scl_allgather(npatch, nsizes)
    patch_offset = sum(nsizes(1:my_rank-1))
    deallocate(nsizes)

    !! Distribute global patch IDs and face-to-patch map.
    !!  Chaparral gathers them again internally.
    allocate(global_ids(npatch_chap), f2p_map(nfacets))
    call scl_scatter(global_ids_g, global_ids)
    call scl_scatter(f2p_map_g, f2p_map)

    !! Assign largest global ID to the virtual surface patch
    if (partial==1 .and. my_rank==nproc) global_ids(npatch_chap) = npatch_tot + 1

    !! Unpack our packed quad/tet connection array into a structured
    !! quad connection array, with dummy values for tets.  Really
    !! silly because after Chaparral gets this in converts it back
    !! to what we started with!
    allocate(c(4,nfacets))
    do j = 1, nfacets
      n = e%xface(face_offset+j+1) - e%xface(face_offset+j)
      c(:n,j) = e%fnode(e%xface(face_offset+j):e%xface(face_offset+j+1)-1)
      if (n == 3) c(4,j) = -1
    end do

    max_surfaces = e%nface
    if (partial==1) max_surfaces = max_surfaces + 1

    !! Chaparral is only capable of rotation about the y-axis, so we cyclically
    !! permute other rotation axes to the y-axis; view factors are invariant.
    !! This means we also have to transform the reflection axes accordingly.
    xmirror = 0
    ymirror = 0
    zmirror = 0
    select case (e%rot_axis)
    case (1)  ! n-fold rotation about x-axis
      x = e%x(3,:)
      y = e%x(1,:)
      z = e%x(2,:)
      nrotations = e%num_rot
      if (e%mirror(1)) ymirror = 1
      if (e%mirror(2)) zmirror = 1
      if (e%mirror(3)) xmirror = 1
    case (2)  ! n-fold rotation about y-axis
      x = e%x(1,:)
      y = e%x(2,:)
      z = e%x(3,:)
      nrotations = e%num_rot
      if (e%mirror(1)) xmirror = 1
      if (e%mirror(2)) ymirror = 1
      if (e%mirror(3)) zmirror = 1
    case (3)  ! n-fold rotation about z-axis
      x = e%x(2,:)
      y = e%x(3,:)
      z = e%x(1,:)
      nrotations = e%num_rot
      if (e%mirror(1)) zmirror = 1
      if (e%mirror(2)) xmirror = 1
      if (e%mirror(3)) ymirror = 1
    case default ! No rotations
      x = e%x(1,:)
      y = e%x(2,:)
      z = e%x(3,:)
      nrotations = 1
      if (e%mirror(1)) xmirror = 1
      if (e%mirror(2)) ymirror = 1
      if (e%mirror(3)) zmirror = 1
    end select

    !! Use Chaparral to compute the view factors.
    call VF_Setup
    call VF_SetNumEnclosures(1)
    call VF_SetMaxSurfaces(max_surfaces)
    !call VF_RandomizeSurfacesOff()
    !call VF_JitterOff()
    call params%get('nonblocking', nonblocking)
    call params%get('asink', asink)
    call params%get('verbosity', verbosity)
    handle = VF_DefineEnclosure(e%name, nonblocking, partial, asink, npatch_chap, global_ids, verbosity)
    call params%get('bsp-depth', bsp_depth)
    call params%get('bsp-length', bsp_length)
    call params%get('spatial-tol', spatial_tol)
    call VF_DefineTopology(handle, GEOM_TYPE, nfacets, e%nnode, x, y, z, &
                           c, 1, f2p_map, nrotations, xmirror, ymirror, zmirror, &
                           bsp_depth, bsp_length, spatial_tol, verbosity)
    call params%get('hc-sub-divide', hc_sub_divide)
    call params%get('hc-resolution', hc_resolution)
    call params%get('hc-min-sep', hc_min_sep)
    call VF_CalcHemicube(handle, hc_sub_divide, hc_resolution, hc_min_sep)
    call params%get('smooth-wt', smooth_wt)
    call params%get('smooth-tol', smooth_tol)
    call params%get('smooth-max-iter', smooth_max_iter)
    call VF_SmoothMatrix(handle, smooth_wt, smooth_tol, smooth_max_iter, SYM_METHOD, verbosity)
    call VF_OutputMatrixSummaryBanner

    !! Extract the VF from Chaparral and create the distributed VF structure.
    evf%npatch = npatch
    evf%offset = patch_offset
    evf%npatch_tot = npatch_tot
    allocate(evf%ia(npatch+1))
    call VF_GetRowCounts(handle, mode=1, count=evf%ia)

    !! Get VF matrix and ambient view factors
    n = sum(evf%ia(:npatch))
    allocate(evf%ja(n), evf%val(n))
    if (partial == 1) then
      evf%has_ambient = .true.
      allocate(evf%ambient(npatch))
      call VF_GetMatrix(handle, evf%ia(2:), evf%ja, evf%val, vf_virt=evf%ambient)
    else
      evf%has_ambient = .false.
      call VF_GetMatrix(handle, evf%ia(2:), evf%ja, evf%val)
    end if

    !! Convert the row counts into the local IA indexing array.
    evf%ia(1) = 1
    do j = 1, npatch
      evf%ia(j+1) = evf%ia(j) + evf%ia(j+1)
    end do

    !! Get patch areas
    evf%has_area = .true.
    npatch_tot_chap = npatch_tot + merge(1, 0, partial==1)
    n = merge(npatch_tot, 0, my_rank==1)
    allocate(parea(npatch_tot_chap), evf%area(n))
    call VF_GetMatrixAreas(parea)
    if (my_rank == 1) evf%area = parea(1:npatch_tot)

    !! Compute face weights
    evf%has_weight = ep%has_patches
    n = merge(e%nface, 0, my_rank==1 .and. evf%has_weight)
    allocate(farea(n), evf%w(n))
    if (n > 0) then
      call compute_face_area(e%xface, e%fnode, e%x, farea)
      !! Check face areas match patch areas computed by Chaparral.
      parea = 0.0_r8
      do j = 1, e%nface
        parea(f2p_map_g(j)) = parea(f2p_map_g(j)) + farea(j)
      end do
      do j = 1, npatch_tot
        INSIST( abs(parea(j)-evf%area(j)) <= 1.0D-14 )
      end do
      !! Compute face weights
      do j = 1, e%nface
        evf%w(j) = farea(j) / evf%area(f2p_map_g(j))
        !! Face weight cannot exceed 1
        INSIST( evf%w(j) <= 1.0_r8 + 1.0D-14 )
        if (evf%w(j) > 1.0_r8) evf%w(j) = 1.0_r8
      end do
      INSIST( all(0.0_r8 < evf%w) )
      INSIST( all(evf%w <= 1.0_r8) )
    end if

    !call VF_ResetTopology(handle)
    call VF_CleanUp

  end subroutine calculate_vf

end module re_chaparral_vf
