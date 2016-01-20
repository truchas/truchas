!!
!! CELL_IMPL
!!
!! Implements the legacy cell structure array and related data and procedures.
!!
!! NOTES
!!
!! The LINEAR_PROP_FACE procedure from LINEAR_MODULE have been incorporated
!! into this module.  Dependency issues forced this choice for the time being.
!! Initialization of the FACE_CENTROID data component uses LINEAR_MODULE which
!! in turns uses the FACE_CENTROID_L component. Besides solid mechanics, the
!! only user of FACE_CENTROID_L is LINEAR_MODULE so I considered moving that
!! data there, which makes a lot of sense, and have this module depend on
!! LINEAR_MODULE. However initialization of FACE_CENTROID_L requires other
!! data components from this module. I would not compromise on is having the
!! declaration of data in the same module as its initialization, which was
!! the approach of the original code.
!!
!! The remaining procedures from LINEAR_MODULE, LINEAR_PROP_VECTOR and
!! LINEAR_GRAD_VECTOR, were independent of FACE_CENTROID_L and used only by
!! solid_mechanics and so were moved there to NODE_OPERATOR_MODULE.
!! LINEAR_FACE_PROP_INT was entirely unused
!!

#include "f90_assert.fpp"

module cell_impl

  use kinds, only: r8
  use common_impl, only: ncells, new_mesh
  implicit none
  private

  public :: init_cell_impl
  public :: linear_prop

  !! A copy of the CELL_GEOMETRY type from MESH_MODULE
  type, public :: cell_geometry
    real(r8) :: centroid(3)          ! cell centroid (physical coordinates)
    real(r8) :: volume               ! cell volume
    real(r8) :: halfwidth(6)         ! cell halfwidth
    real(r8) :: face_centroid(3,6)   ! face centroid (physical coordinates)
    real(r8) :: face_centroid_l(3,6) ! face centroid (logical coordinates)
    real(r8) :: face_area(6)         ! face area
    real(r8) :: face_normal(3,6)     ! face unit normal
  end type cell_geometry

  !! Replacement for the CELL structure array from MESH_MODULE
  type(cell_geometry), allocatable, target, public :: cell(:)

  logical, protected, public :: orthogonal_mesh

  interface linear_prop
    module procedure linear_prop_face
  end interface

contains

  subroutine init_cell_impl

    use discrete_ops_data, only: use_ortho_face_gradient, discrete_ops_type

    integer :: j
    real(r8), allocatable :: centroid(:,:), face_area(:,:), face_normal(:,:,:)
    real(r8), allocatable :: face_centroid_l(:,:,:), face_centroid(:,:,:), halfwidth(:,:)

    allocate(cell(ncells))

    !! CELL%VOLUME
    call init_volume (cell%volume)
    !call test_init_volume

    !! CELL%CENTROID
    allocate(centroid(3,ncells))
    call init_centroid (centroid)
    do j = 1, ncells
      cell(j)%centroid = centroid(:,j)
    end do
    deallocate(centroid)
    !call test_init_centroid

    !! CELL%FACE_AREA, CELL%FACE_NORMAL
    allocate(face_area(6,ncells), face_normal(3,6,ncells))
    call init_face_normal (face_area, face_normal)
    do j = 1, ncells
      cell(j)%face_area = face_area(:,j)
      cell(j)%face_normal = face_normal(:,:,j)
    end do
    deallocate(face_area, face_normal)
    !call test_init_face_normal

    !! CELL%FACE_CENTROID_L
    allocate(face_centroid_l(3,6,ncells))
    call init_face_centroid_l (face_centroid_l)
    do j = 1, ncells
      cell(j)%face_centroid_l = face_centroid_l(:,:,j)
    end do
    deallocate(face_centroid_l)
    !call test_init_face_centroid_l

    !! CELL%FACE_CENTROID
    allocate(face_centroid(3,6,ncells))
    call init_face_centroid (face_centroid)
    do j = 1, ncells
      cell(j)%face_centroid = face_centroid(:,:,j)
    end do
    deallocate(face_centroid)
    !call test_init_face_centroid

    !! CELL%HALFWIDTH
    allocate(halfwidth(6,ncells))
    call init_halfwidth (halfwidth)
    do j = 1, ncells
      cell(j)%halfwidth = halfwidth(:,j)
    end do
    deallocate(halfwidth)
    !call test_init_halfwidth

    !! ORTHOGONAL_MESH
    call init_orthogonal_mesh (orthogonal_mesh)
    !call test_init_orthogonal_mesh

    !! The original cell_geometry_module::jacobian procedure included the
    !! following as a side effect of computing the orthogonal_mesh flag.
    !! This is an instance of "mesh" using "client" (not supposed to happen)
    !! which was not indentified during the initial code analysis.
    !! TODO: MOVE THIS OUTSIDE OF OLD_MESH_API -- IT DOES NOT BELONG HERE.
    if (discrete_ops_type == 'default') use_ortho_face_gradient = orthogonal_mesh

  end subroutine init_cell_impl

  !! Copy cell volumes from the new mesh.  Gap cell volumes are arbitrarily
  !! set to 2*alittle as was originally done in the old mesh.  Note that the
  !! volume formula used by the new mesh is very different and leads to
  !! relatively small volume differences (though as much as 5 or 6 least
  !! significant decimal digits). The original formula due to Dukowicz suffers
  !! from substantial cancellation error and is really quite poor.

  subroutine init_volume (volume)
    use common_impl, only: pcell_old_to_new
    use parallel_permutations, only: rearrange
    use cutoffs_module, only: alittle
    real(r8), intent(out) :: volume(:)
    ASSERT(size(volume) == ncells)
    call rearrange (pcell_old_to_new, volume, new_mesh%volume(:new_mesh%ncell_onP), default=2*alittle)
  end subroutine init_volume

  !! Use the legacy algorithm for centroids.  We use the new mesh volumes.
  !! If the old mesh volumes are used (commented out code) the centroid
  !! results are bit-for-bit the same, as expected.  But with the new
  !! volumes there are small differences; see test_init_centroid below.

  subroutine init_centroid (centroid)
    !use mesh_module, only: old_cell => cell
    use en_gather_impl, only: gather_vertex_coord
    use legacy_geometry, only: cell_centroid
    real(r8), intent(out) :: centroid(:,:)
    real(r8) :: x(3,8,ncells)
    integer :: j
    ASSERT(size(centroid,1) == 3)
    ASSERT(size(centroid,2) == ncells)
    call gather_vertex_coord (x)
    do j = 1, ncells
      !call cell_centroid (x(:,:,j), old_cell(j)%volume, centroid(:,j))
      call cell_centroid (x(:,:,j), cell(j)%volume, centroid(:,j))
    end do
  end subroutine init_centroid

  !! We could pull these values from the new mesh (and will eventually) but
  !! it is simpler for now to recompute them using the legacy algorithm.

  subroutine init_face_normal (area, normal)
    use en_gather_impl, only: gather_vertex_coord
    use legacy_geometry, only: face_area
    real(r8), intent(out) :: area(:,:), normal(:,:,:)
    real(r8) :: x(3,8,ncells)
    integer :: j
    ASSERT(size(area,1) == 6)
    ASSERT(size(area,2) == ncells)
    ASSERT(size(normal,1) == 3)
    ASSERT(size(normal,2) == 6)
    ASSERT(size(normal,3) == ncells)
    call gather_vertex_coord (x)
    do j = 1, ncells
      call face_area (x(:,:,j), area(:,j), normal(:,:,j))
    end do
  end subroutine init_face_normal

  subroutine init_face_centroid_l (centroid_l)
    use en_gather_impl, only: gather_vertex_coord
    use legacy_geometry, only: face_centroid_logical
    real(r8), intent(out) :: centroid_l(:,:,:)
    real(r8) :: x(3,8,ncells)
    integer :: j
    ASSERT(size(centroid_l,1) == 3)
    ASSERT(size(centroid_l,2) == 6)
    ASSERT(size(centroid_l,3) == ncells)
    call gather_vertex_coord (x)
    do j = 1, ncells
      call face_centroid_logical (x(:,:,j), cell(j)%face_area, cell(j)%face_normal, centroid_l(:,:,j))
    end do
  end subroutine init_face_centroid_l

  !! Essentially cribbed from cell_geometry_module::face_centroid_physical.

  subroutine init_face_centroid (centroid)
    use en_gather_impl, only: gather_vertex_coord
    real(r8), intent(out) :: centroid(:,:,:)
    integer :: i, j, n
    real(r8) :: coord(8,ncells)
    ASSERT(size(centroid,1) == 3)
    ASSERT(size(centroid,2) == 6)
    ASSERT(size(centroid,3) == ncells)
    do n = 1, 3
      call gather_vertex_coord (coord, dim=n)
      do j = 1, ncells
        do i = 1, 6
          call linear_prop_face_aux (cell(j)%face_centroid_l, i, coord(:,j), centroid(n,i,j))
        end do
      end do
    end do
  end subroutine init_face_centroid

  !! Essentially cribbed from cell_geometry_module::cell_width.

  subroutine init_halfwidth (halfwidth)
    real(r8), intent(out) :: halfwidth(:,:)
    integer :: j, k
    ASSERT(size(halfwidth,1) == 6)
    ASSERT(size(halfwidth,2) == ncells)
    do j = 1, ncells
      do k = 1, 6
        halfwidth(k,j) = norm2(cell(j)%face_centroid(:,k)-cell(j)%centroid)
      end do
    end do
  end subroutine init_halfwidth

  !! Essentially cribbed from cell_geometry_module::jacobian.

  subroutine init_orthogonal_mesh (orthogonal_mesh)
    use legacy_geometry, only: is_cell_orthog
    use parallel_communication, only: global_all
    logical, intent(out) :: orthogonal_mesh
    integer :: j
    orthogonal_mesh = .true.
    do j = 1, ncells
      orthogonal_mesh = is_cell_orthog(cell(j)%face_centroid)
      if (.not.orthogonal_mesh) exit
    end do
    orthogonal_mesh = global_all(orthogonal_mesh)
  end subroutine init_orthogonal_mesh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! LINEAR_PROP CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! This code is originally from LINEAR_MODULE.  It has been reorganized
  !! a bit but it is not essentially different.  The nesting of do-loops was
  !! changed and the part operating on a single cell called out as auxiliary
  !! procedure not depending on the CELL structure array.

  !=======================================================================
  ! Purpose(s):
  !   Property evaluation using 2-D or 3-D linear interpolation at
  !   the cell-face centroid vector location.
  !
  !               nvc
  !               ---
  !               \
  !   Prop_p(:) =  |  Coef_v(:) * Vrtx_v(:)
  !               /
  !               ---
  !               v=1
  !
  !               ndim
  !               ___
  !               | |
  !   Coef_v(:) = | | [LC1(n,v) + LC2(n,v)*Xi(n,:)]
  !               | |
  !               n=1
  !
  !=======================================================================

  pure subroutine linear_prop_face_aux (face_centroid_l, f, vrtx, prop)

    real(r8), intent(in)  :: face_centroid_l(:,:)
    integer,  intent(in)  :: f
    real(r8), intent(in)  :: vrtx(:)
    real(r8), intent(out) :: prop

    integer  :: n, v
    real(r8) :: Coef

    !! Originally from LINEAR_MODULE.  These assume a possibly degenerate hex.
    real(r8), parameter :: LC1(3,8) = reshape( &
        [ 0.0_r8,  1.0_r8,  1.0_r8, &
          0.0_r8,  0.0_r8,  1.0_r8, &
          1.0_r8,  0.0_r8,  1.0_r8, &
          1.0_r8,  1.0_r8,  1.0_r8, &
          0.0_r8,  1.0_r8,  0.0_r8, &
          0.0_r8,  0.0_r8,  0.0_r8, &
          1.0_r8,  0.0_r8,  0.0_r8, &
          1.0_r8,  1.0_r8,  0.0_r8 ], shape=[3,8])

    real(r8), parameter :: LC2(3,8) = reshape( &
        [ 1.0_r8, -1.0_r8, -1.0_r8, &
          1.0_r8,  1.0_r8, -1.0_r8, &
         -1.0_r8,  1.0_r8, -1.0_r8, &
         -1.0_r8, -1.0_r8, -1.0_r8, &
          1.0_r8, -1.0_r8,  1.0_r8, &
          1.0_r8,  1.0_r8,  1.0_r8, &
         -1.0_r8,  1.0_r8,  1.0_r8, &
         -1.0_r8, -1.0_r8,  1.0_r8 ], shape=[3,8])

    prop = 0.0_r8
    do v = 1, 8
      coef = 1.0_r8
      do n = 1, 3
        coef = coef * (lc1(n,v) + lc2(n,v)*face_centroid_l(n,f))
      end do
      prop = prop + coef*vrtx(v)
    end do

  end subroutine linear_prop_face_aux

  subroutine linear_prop_face (f, vrtx, prop)
    integer, intent(in) :: f
    real(r8), intent(in) :: vrtx(:,:)
    real(r8), intent(out) :: prop(:)
    integer :: j
    ASSERT(size(vrtx,1) == 8)
    ASSERT(size(vrtx,2) == ncells)
    ASSERT(size(prop) == ncells)
    do j = 1, ncells
      call linear_prop_face_aux (cell(j)%face_centroid_l, f, vrtx(:,j), prop(j))
    end do
  end subroutine linear_prop_face

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! TESTING CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! The legacy algorithm for computing cell volume suffers from substantial
  !! cancellation error and was a poor choice.  Note the relatively large
  !! error tolerance (5 or 6 decimal digits).  It is as small as it can be
  !! to have this test pass for all regression test meshes.  The regression
  !! tests themselves are not sensitive to the differences.

  subroutine test_init_volume
    use mesh_module, only: old_cell => cell
    use parallel_communication, only: global_any, global_maxval
    use truchas_logging_services
    integer :: j, unit
    logical :: error
    real(r8) :: volume(ncells)
    call init_volume (volume)
    unit = TLS_debug_unit()
    error = .false.
    do j = 1, ncells
      if (abs(volume(j) - old_cell(j)%volume) <= (2**16)*spacing(volume(j))) cycle
      if (.not. error) then
        write(unit,'(/,a)') 'CELL%VOLUME *********************************'
        error = .true.
      end if
      write(unit,'(a,i0,a,1x,g0)') 'old[', j, ']=', old_cell(j)%volume
      write(unit,'(a,i0,a,1x,g0)') 'new[', j, ']=', volume(j)
    end do
    if (global_any(error)) then
      write(unit,'(a,g0)') 'MAX REL ERROR = ',global_maxval(abs(volume - old_cell%volume)/volume)
      call TLS_fatal ('error validating VOLUME')
    end if
  end subroutine test_init_volume

  !! The legacy formula is still being used to compute centroids, but with the
  !! new volume.  Hence the error tolerance (as small as it can be and have the
  !! test pass for all regression test meshes) is somewhat larger than would
  !! normally be expected.  There is also a testing issue that arises when a
  !! centroid component is 0, or nearly so -- that is the reason for the
  !! inclusion of an absolute tolerance.

  subroutine test_init_centroid
    use mesh_module, only: old_cell => cell
    use parallel_communication, only: global_any
    use truchas_logging_services
    integer :: j, unit
    logical :: error
    real(r8) :: centroid(3,ncells)
    call init_centroid (centroid)
    unit = TLS_debug_unit()
    error = .false.
    do j = 1, ncells
      if (all(abs(centroid(:,j) - old_cell(j)%centroid) <= &
              max((2**13)*spacing(old_cell(j)%centroid),8*epsilon(1.0_r8)))) cycle
      if (.not. error) then
        write(unit,'(/,a)') 'CELL%CENTROID *********************************'
        error = .true.
      end if
      write(unit,'(a,i0,a,3(1x,g0))') 'old[', j, ']=', old_cell(j)%centroid
      write(unit,'(a,i0,a,3(1x,g0))') 'new[', j, ']=', centroid(:,j)
    end do
    if (global_any(error)) call TLS_fatal ('error validating CENTROID')
  end subroutine test_init_centroid

  !! As long as we are using the legacy algorithm for computing the face areas
  !! and normals we should get exactly the same results, and we are as reflected
  !! by the very small tolerance (2 least significant bits).

  subroutine test_init_face_normal
    use mesh_module, only: old_cell => cell
    use parallel_communication, only: global_any
    use truchas_logging_services
    integer :: j, k, unit
    logical :: error
    real(r8) :: face_area(6,ncells), face_normal(3,6,ncells)
    call init_face_normal (face_area, face_normal)
    unit = TLS_debug_unit()
    error = .false.
    do j = 1, ncells
      if (all(abs(face_area(:,j) - old_cell(j)%face_area) <= 2*spacing(old_cell(j)%face_area))) cycle
      if (.not. error) then
        write(unit,'(/,a)') 'CELL%FACE_AREA *********************************'
        error = .true.
      end if
      write(unit,'(a,i0,a,*(1x,g0))') 'old[', j, ']=', old_cell(j)%face_area
      write(unit,'(a,i0,a,*(1x,g0))') 'new[', j, ']=', face_area(:,j)
    end do
    if (global_any(error)) call TLS_fatal ('error validating FACE_AREA')
    error = .false.
    do j = 1, ncells
      do k = 1, 6
        if (all(abs(face_normal(:,k,j) - old_cell(j)%face_normal(:,k)) <= &
                2*spacing(old_cell(j)%face_normal(:,k)))) cycle
        if (.not. error) then
          write(unit,'(/,a)') 'CELL%FACE_NORMAL *********************************'
          error = .true.
        end if
        write(unit,'(2(a,i0),a,*(1x,g0))') 'old[', k, ',', j, ']=', old_cell(j)%face_normal(:,k)
        write(unit,'(2(a,i0),a,*(1x,g0))') 'new[', k, ',', j, ']=', face_normal(:,k,j)
      end do
    end do
    if (global_any(error)) call TLS_fatal ('error validating FACE_NORMAL')
  end subroutine test_init_face_normal

  subroutine test_init_face_centroid_l
    use mesh_module, only: old_cell => cell
    use parallel_communication, only: global_any
    use truchas_logging_services
    integer :: j, k, unit
    logical :: error
    real(r8) :: face_centroid_l(3,6,ncells)
    call init_face_centroid_l (face_centroid_l)
    unit = TLS_debug_unit()
    error = .false.
    do j = 1, ncells
      do k = 1, 6
        if (all(abs(face_centroid_l(:,k,j) - old_cell(j)%face_centroid_l(:,k)) <= &
                2*spacing(old_cell(j)%face_centroid_l(:,k)))) cycle
        if (.not. error) then
          write(unit,'(/,a)') 'CELL%FACE_CENTROID_L *********************************'
          error = .true.
        end if
        write(unit,'(2(a,i0),a,*(1x,g0))') 'old[', k, ',', j, ']=', old_cell(j)%face_centroid_l(:,k)
        write(unit,'(2(a,i0),a,*(1x,g0))') 'new[', k, ',', j, ']=', face_centroid_l(:,k,j)
      end do
    end do
    if (global_any(error)) call TLS_fatal ('error validating FACE_CENTROID_L')
  end subroutine test_init_face_centroid_l

  subroutine test_init_face_centroid
    use mesh_module, only: old_cell => cell
    use parallel_communication, only: global_any
    use truchas_logging_services
    integer :: j, k, unit
    logical :: error
    real(r8) :: face_centroid(3,6,ncells)
    call init_face_centroid (face_centroid)
    unit = TLS_debug_unit()
    error = .false.
    do j = 1, ncells
      do k = 1, 6
        if (all(abs(face_centroid(:,k,j) - old_cell(j)%face_centroid(:,k)) <= &
                max(8*spacing(old_cell(j)%face_centroid(:,k)),2*epsilon(1.0_r8)))) cycle
        if (.not. error) then
          write(unit,'(/,a)') 'CELL%FACE_CENTROID *********************************'
          error = .true.
        end if
        write(unit,'(2(a,i0),a,*(1x,g0))') 'old[', k, ',', j, ']=', old_cell(j)%face_centroid(:,k)
        write(unit,'(2(a,i0),a,*(1x,g0))') 'new[', k, ',', j, ']=', face_centroid(:,k,j)
      end do
    end do
    if (global_any(error)) call TLS_fatal ('error validating FACE_CENTROID')
  end subroutine test_init_face_centroid

  !! The error tolerance here needs to be surprisingly large in order for all
  !! the regression test meshes to pass.  The face centroids are very close to
  !! the old values, but the centroid above does differ much more significantly,
  !! so what might expect needing a similarly sized tolerance here.  But in fact
  !! the tolerance needs to be about 10 times larger.

  subroutine test_init_halfwidth
    use mesh_module, only: old_cell => cell
    use parallel_communication, only: global_any
    use truchas_logging_services
    integer :: j, unit
    logical :: error
    real(r8) :: halfwidth(6,ncells)
    call init_halfwidth (halfwidth)
    unit = TLS_debug_unit()
    error = .false.
    do j = 1, ncells
      if (all(abs(halfwidth(:,j) - old_cell(j)%halfwidth) <= &
              (2**17)*spacing(old_cell(j)%halfwidth))) cycle
      if (.not. error) then
        write(unit,'(/,a)') 'CELL%HALFWIDTH *********************************'
        error = .true.
      end if
      write(unit,'(a,i0,a,*(1x,g0))') 'old[', j, ']=', old_cell(j)%halfwidth
      write(unit,'(a,i0,a,*(1x,g0))') 'new[', j, ']=', halfwidth(:,j)
    end do
    if (global_any(error)) call TLS_fatal ('error validating HALFWIDTH')
  end subroutine test_init_halfwidth

  subroutine test_init_orthogonal_mesh
    use mesh_module, only: old_orthogonal_mesh => orthogonal_mesh
    use truchas_logging_services
    logical :: orthogonal_mesh
    call init_orthogonal_mesh (orthogonal_mesh)
    if (orthogonal_mesh .neqv. old_orthogonal_mesh) &
        call TLS_fatal ('error validating ORTHOGONAL_MESH')
  end subroutine test_init_orthogonal_mesh

end module cell_impl
