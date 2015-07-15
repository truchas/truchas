!!
!! CHACO_C_BINDING
!!
!! Fortran interface to the Chaco C library function.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! May 2015
!!

module chaco_c_binding

  use,intrinsic :: iso_c_binding, only: c_int, c_short, c_long, c_float, c_double, c_ptr
  implicit none
  private

  public :: interface
  
  !! In the true interface (commented out below) a null pointer may be passed
  !! to some of the arguments to indicate "no value".  This will be handled in
  !! F2015 by making those arguments optional. Until then we declare a special
  !! interface where those arguments we wish to omit are declared as C_PTR and
  !! pass by value, so that we can explicitly pass a C_NULL_PTR.

  interface
    function interface (nvtxs, start, adjacency, vwgts, ewgts, x, y, z, &
            outassignname, outfilename, assignment, architecture, ndims_tot, mesh_dims, &
            goal, global_method, local_method, rqi_flag, vmax, ndims, eigtol, seed) &
        result (ierr) bind(c,name='interface')
      import c_int, c_short, c_long, c_float, c_double, c_ptr
      integer(c_int),  value :: nvtxs             ! number of vertices in full graph
      integer(c_int), intent(in) :: start(*)      ! start of edge list for each vertex
      integer(c_int), intent(in) :: adjacency(*)  ! edge list data
      type(c_ptr), value         :: vwgts         ! weights for all vertices
      real(c_float),  intent(in) :: ewgts(*)      ! weights for all edges
      type(c_ptr), value         :: x, y, z       ! coordinates for inertial method
      type(c_ptr), value         :: outassignname ! name of assignment output file
      type(c_ptr), value         :: outfilename   ! output file name
      integer(c_short), intent(out) :: assignment(*)  ! set number of each vtx (length n)
      integer(c_int),  value :: architecture      ! 0 => hypercube, d => d-dimensional mesh
      integer(c_int),  value :: ndims_tot         ! total number of cube dimensions to divide
      integer(c_int)         :: mesh_dims(3)      ! dimensions of mesh of processors
      type(c_ptr),     value :: goal              ! desired set sizes for each set
      integer(c_int),  value :: global_method     ! global partitioning algorithm
      integer(c_int),  value :: local_method      ! local partitioning algorithm
      integer(c_int),  value :: rqi_flag          ! should I use RQI/Symmlq eigensolver?
      integer(c_int),  value :: vmax              ! how many vertices to coarsen down to?
      integer(c_int),  value :: ndims             ! number of eigenvectors (2^d sets)
      real(c_double),  value :: eigtol            ! tolerance on eigenvectors
      integer(c_long), value :: seed              ! for random graph mutations
      integer(c_int) :: ierr
    end function
  end interface

  !! This is the true interface.  A null pointer can be passed to some of the
  !! arguments to indicate "no value".  In Fortran 2015 these are interoperable
  !! with optional arguments; a null pointer is passed if the argument is absent

!  interface
!    function interface (nvtxs, start, adjacency, vwgts, ewgts, x, y, z, &
!            outassignname, outfilename, assignment, architecture, ndims_tot, mesh_dims, &
!            goal, global_method, local_method, rqi_flag, vmax, ndims, eigtol, seed) &
!        result (ierr) bind(c,name='interface')
!      import c_int, c_short, c_long, c_float, c_double, c_char
!      integer(c_int),  value :: nvtxs             ! number of vertices in full graph
!      integer(c_int), intent(in) :: start(*)      ! start of edge list for each vertex
!      integer(c_int), intent(in) :: adjacency(*)  ! edge list data
!      integer(c_int), intent(in) :: vwgts(*)      ! weights for all vertices
!      real(c_float),  intent(in) :: ewgts(*)      ! weights for all edges
!      real(c_float),  intent(in) :: x(*), y(*), z(*)          ! coordinates for inertial method
!      character(kind=c_char), intent(in) :: outassignname(*)  ! name of assignment output file
!      character(kind=c_char), intent(in) :: outfilename(*)    ! output file name
!      integer(c_short), intent(out) :: assignment(*)          ! set number of each vtx (length n)
!      integer(c_int),  value :: architecture      ! 0 => hypercube, d => d-dimensional mesh
!      integer(c_int),  value :: ndims_tot         ! total number of cube dimensions to divide
!      integer(c_int)         :: mesh_dims(3)      ! dimensions of mesh of processors
!      real(c_double), intent(in) :: goal(*)       ! desired set sizes for each set
!      integer(c_int),  value :: global_method     ! global partitioning algorithm
!      integer(c_int),  value :: local_method      ! local partitioning algorithm
!      integer(c_int),  value :: rqi_flag          ! should I use RQI/Symmlq eigensolver?
!      integer(c_int),  value :: vmax              ! how many vertices to coarsen down to?
!      integer(c_int),  value :: ndims             ! number of eigenvectors (2^d sets)
!      real(c_double),  value :: eigtol            ! tolerance on eigenvectors
!      integer(c_long), value :: seed              ! for random graph mutations
!      integer(c_int) :: ierr
!      optional :: vwgts, ewgts, x, y, z, outassignname, outfilename, goal
!    end function
!  end interface

end module chaco_c_binding
