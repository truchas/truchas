!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vof_init
   !----------------------------------------------------------------------------
   ! purpose:
   !
   !    Calculate the volume of each body in the cells, for volume fraction
   !    initialization.
   !
   !    There are two methods available, selected by the VOF_METHOD
   !    input flag in the INTERFACES namelist.  Valid values are
   !    'points' and 'divide'.
   !
   !    The points method is what's been in the code for a long time.  It
   !    identifies interface cells in a first pass, by sampling at a
   !    user-specified number of points in each cell (and adding some extra
   !    cells), then makes a second pass over the interface cells, sampling at a
   !    user-specified number of points in each interface cell.
   !
   !    The divide method is new.  It checks each cell for the presence
   !    of multiple bodies in a cell (an interface).  If it finds an interface,
   !    it subdivides the cell into subcells and checks each subcell for an
   !    interface.  It continues this subdivision recursively until the subcells
   !    are contained entirely in one body, or are smaller than a user-specified
   !    size.  An approximation to the partial volume of of each body in the
   !    subcell is then made and the partial volumes are returned to the caller,
   !    to be accumulated into the original cells.
   !
   ! contacts:
   !    Bryan Lally, lally@lanl.gov
   !    Marianne Francois, mmfran@lanl.gov
   !
   ! to do:
   !    detect degenerate vertices and skip testing them multiple times
   !    improve creation of subcells when degenerate vertices are involved
   !----------------------------------------------------------------------------

   use kinds, only: r8
   use legacy_mesh_api, only: ndim, nvc, ncells, nfc, nec
   use interfaces_module, only: nbody, vof_method, vof_tolerance, vof_max_recursion
   use parameter_module, only: nmat
   use truchas_logging_services
   implicit none
   private
   
   public :: vof_initialize

   ! local development/testing control flags
   logical, parameter :: test_vertex_once     = .true.
   logical, parameter :: degenerate_vertices  = .false.
   logical, parameter :: make_gmv_picture     = .false.

   ! symbolic constants
   integer, parameter :: BODY_ID_UNKNOWN = -1

   ! vertex_info type is a collection of a vertex's location and the number of
   ! the body that the vertex is in, maybe other information

   type vertex_info
      integer :: body_id
      real(r8), dimension(ndim) :: location
   end type vertex_info

contains
   function vof_initialize()
      !-------------------------------------------------------------------------
      ! top-level drive with method switch
      !-------------------------------------------------------------------------

      use tally_module,      only: vof_init_tally_points
      use legacy_mesh_api,   only: cell

      ! return value
      real(r8), dimension(nbody,ncells) :: vof_initialize

      !-------------------------------------------------------------------------

      ! If there is only have one body (nbody=1) and one material (nmat=1) then
      ! there is no need to compute the partial volumes.  The entire
      ! domain is made of the same stuff.  Return the entire volume of each cell.
      if (nbody == 1 .and. nmat == 1) then
         vof_initialize(1,:) = Cell%Volume 
         vof_initialize(2:,:) = 0.0d0
         return
      end if

      ! select a method and calculate the partial volumes
      select case (trim(vof_method))
      case ('points')
         vof_initialize = vof_init_tally_points()
      case ('divide')
         vof_initialize = vof_init_divide_conquer()
      case default
         call TLS_fatal ('VOF_INITIALIZE: invalid VOF initialization method: '// trim(vof_method))
      end select

!     get some results out for a quick check
!      write (*,*) 'vof_initialize:'
!      write (*,*) vof_initialize

   end function vof_initialize

   !----------------------------------------------------------------------------

   function vof_init_divide_conquer()
      !-------------------------------------------------------------------------
      ! purpose:
      !    Call the volume of a body in a mesh cell the partial volume of that
      !    body in a cell.  This procedure calculates the partial volume of each
      !    body in each cell, using a recursive divide and conquer algorithm.
      !
      ! arguments:
      !    none, requires a properly initialized mesh and body descriptions
      !
      ! return:
      !    nbody x ncells array of partial volumes
      !
      ! algorithm:
      !    for cell in cells:
      !        partial_volumes(nbody,cell) = partial_volumes_in_cell(recursion_depth, tolerance, cell_vertices)
      !-------------------------------------------------------------------------

      use legacy_mesh_api,  only: mesh, cell, gather_vertex_coord
      use PGSLib_module,    only: PGSLib_GLOBAL_MINVAL
      use input_utilities,  only: NULL_R, NULL_I

      ! return value
      real(r8), dimension(nbody,ncells) :: vof_init_divide_conquer

      ! local variables
      integer :: i
      integer :: v                               ! vertex loop index
      integer :: c                               ! cell loop index
      integer :: vtx                             ! temporary, holds a vertex number
      real(r8) :: tolerance                       ! volume at which recursive algorithm terminates
      type (vertex_info), dimension(nvc) :: vertices        ! nvc array of vertex_info structures - location + body ids
      real(r8) :: vtx_coord(ndim,nvc,ncells)

      !-------------------------------------------------------------------------

      if (ndim /= 3) then
         call TLS_fatal ('VOF_INIT_DIVIDE_CONQUER: Vof_Method divide requires ndim = 3')
      end if

      ! if user hasn't set vof_tolerance nor vof_max_recursion, apply some defaults
      if (vof_tolerance == NULL_R) then
         vof_tolerance = 0.001d0
      end if

      if (vof_max_recursion == NULL_I) then
         vof_max_recursion = 100
      end if

      ! Look for the smallest cell in the mesh, and set the tolerance by
      ! multiplying the smallest cell volume by an input parameter.
      ! tolerance = input_parameter * smallest_cell_volume
      tolerance = vof_tolerance * PGSLIB_GLOBAL_MINVAL(cell%volume)
      
      !! NNC: I'm working around a parallel bug in the next loop, and I want
      !! to muck with it as little as possible.  With more time and effort
      !! I could probably merge the two and eliminate this extra storage but ...
      call gather_vertex_coord (vtx_coord)

      ! loop over the cells, calculating the volume of each body in the cell
      do c = 1, ncells

         ! We need the locations of the nvc vertices of a cell - store them in
         ! an nvc array of vertex_info's.  Also initialize the body identifier
         ! for each vertex to an impossible value (BODY_ID_UNKNOWN - the bodies
         ! are numbered with positive integers starting at 1), which means we
         ! haven't figured out what body this vertex is in.
         do v = 1, nvc
            ! Find the vertex number of each vertex (mesh(c)%ngbr_vrtx(v)), use
            ! that to look up the coordinates of that vertex.  Use a temporary
            ! (vtx) so we can see the indirection values with a debugger.
            ! NNC: the next commented lines don't work in parallel!
            !vtx = mesh(c)%ngbr_vrtx(v)
            !vertices(v)%location(:) = vertex(vtx)%coord(:)
            vertices(v)%location(:) = vtx_coord(:,v,c)
            vertices(v)%body_id = BODY_ID_UNKNOWN
         end do

         ! If we have degenerate hexahedral cells, we have multiple vertices at
         ! the same coordinates.  We can avoid some work by detecting them now
         ! and only finding the body that they are contained in once.  We can do
         ! that here.  We'll have to add some logic when we subdivide cells into
         ! subcells as well.  That's where the real savings will occur.
         if (degenerate_vertices) then
            continue
         end if

         ! From the vertex location information, get the body partial volumes.
         vof_init_divide_conquer(:,c) = partial_volumes(1, tolerance, vertices(:))

      end do
      
   end function vof_init_divide_conquer

   !----------------------------------------------------------------------------

   recursive function partial_volumes (recursion_depth, tolerance, vertices) result (pv)
      !-------------------------------------------------------------------------
      ! purpose:
      !    Given the vertices of a cell, calculate the partial volume of each
      !    body in that cell.
      !
      ! arguments:
      !    recursion_depth    current number of partial_volumes calls on the stack
      !    tolerance          volume below which partial volumes are approximated
      !    vertices           nvc array of vertex_info structures - locations and body ids
      !
      ! return:
      !    an nbody array of the partial volumes of each body in the specified cell
      !
      ! algorithm:
      !    if all vertices and the centroid are in the same body
      !       return partial volumes of the bodies in the cell
      !          (all partial volumes are zero, except for the body #, which
      !          is cell volume)
      !    if cell volume < tolerance
      !       approximate partial volumes inside the body
      !       return partial volumes
      !    else
      !       subdivide into octants
      !       partial_volumes_sum = 0
      !       for cell in octant
      !          partial_volumes_sum += volume_of_body_in_cell(body, tolerance, cell vertices)
      !       return partial_volumes_sum
      !-------------------------------------------------------------------------

      ! arguments
      integer,  intent(in) :: recursion_depth
      real(r8), intent(in) :: tolerance
      type(vertex_info), dimension(nvc), intent(inout) :: vertices

      ! return value
      real(r8), dimension(nbody) :: pv

      ! local variables
      logical :: all_the_same
      integer :: v
      integer :: n
      real(r8) :: volume
      real(r8), dimension(ndim) :: cell_centroid

      ! nvc vertices for each cell, nvc subcells per cell
      type (vertex_info), dimension(nvc,nvc) :: subcell_vertices

      !-------------------------------------------------------------------------

      ! calculate the cell volume
      volume = cell_volume(vertices(:))

      ! find the body each vertex is in, store it in the vertex's body_id
      do v = 1, nvc
         if (test_vertex_once) then
            ! check the body id only if we haven't already figured it out
            if (vertices(v)%body_id == BODY_ID_UNKNOWN) then
               vertices(v)%body_id = body_id_from_vertex(vertices(v)%location(:))
            end if
         else
            ! brute force - check every time
            vertices(v)%body_id = body_id_from_vertex(vertices(v)%location(:))
         end if
      end do

      ! If all vertices and the cell centroid are in the same body,
      ! return the cell volume for that body's partial volume, zero
      ! for the rest.
      all_the_same = ALL(vertices(2:nvc)%body_id == vertices(1)%body_id)

      ! (maybe) make an additional check at cell centroid
      if (all_the_same) then
         do n=1, ndim
            cell_centroid(n) = 0.0_r8
            do v = 1, nvc
               cell_centroid(n) = cell_centroid(n) + vertices(v)%location(n)
            end do
            cell_centroid(n) = cell_centroid(n) / nvc
         end do
         if (body_id_from_vertex(cell_centroid) == vertices(1)%body_id) then
            pv(:) = 0.0d0
            pv(vertices(1)%body_id) = volume        ! all the body ids are the same
            goto 999
         end if
      end if

100   continue

      ! If cell volume has become small enough (from
      ! subdivision/recursion), or if maximum recursion has been
      ! reached, make an approximation to the partial volumes.  This is
      ! the termination of the recursive algorithm.

      if (volume < tolerance .or. recursion_depth > vof_max_recursion) then
         pv(:) = approximation(vertices, volume)
         goto 999
      end if

      ! Subdivide cell into quadrants (2D) or octants (3D) and recursively
      ! calculate the partial volume contributions from each subcell.
      subcell_vertices = subdivide_cell(vertices)
      pv(:) = 0.0d0
      do v = 1, nvc
         pv(:) = pv(:) + partial_volumes(recursion_depth+1, tolerance, subcell_vertices(:,v))

         ! After each subcell is calculated, we know the body identifications
         ! for some of the vertices in the subcells still to do.  If we're using
         ! the test_vertex_once optimizations, we need to update the subcell
         ! data.  This requires knowledge of the cell division being done in
         ! subdivide_cell().
         !
         ! The numbering of subcell_vertices(:,:) is (vertex, subcell#)

         if (test_vertex_once) then
            select case (v)
            case (1)
               subcell_vertices(1,2)%body_id = subcell_vertices(2,1)%body_id
               subcell_vertices(4,2)%body_id = subcell_vertices(3,1)%body_id
               subcell_vertices(5,2)%body_id = subcell_vertices(6,1)%body_id
               subcell_vertices(8,2)%body_id = subcell_vertices(7,1)%body_id
               subcell_vertices(1,3)%body_id = subcell_vertices(3,1)%body_id
               subcell_vertices(5,3)%body_id = subcell_vertices(7,1)%body_id
               subcell_vertices(2,4)%body_id = subcell_vertices(3,1)%body_id
               subcell_vertices(1,4)%body_id = subcell_vertices(4,1)%body_id
               subcell_vertices(6,4)%body_id = subcell_vertices(7,1)%body_id
               subcell_vertices(5,4)%body_id = subcell_vertices(8,1)%body_id
               subcell_vertices(1,5)%body_id = subcell_vertices(5,1)%body_id
               subcell_vertices(2,5)%body_id = subcell_vertices(6,1)%body_id
               subcell_vertices(3,5)%body_id = subcell_vertices(7,1)%body_id
               subcell_vertices(4,5)%body_id = subcell_vertices(8,1)%body_id
               subcell_vertices(1,6)%body_id = subcell_vertices(6,1)%body_id
               subcell_vertices(4,6)%body_id = subcell_vertices(7,1)%body_id
               subcell_vertices(1,7)%body_id = subcell_vertices(7,1)%body_id
               subcell_vertices(2,8)%body_id = subcell_vertices(7,1)%body_id
               subcell_vertices(1,8)%body_id = subcell_vertices(8,1)%body_id
            case (2)
               subcell_vertices(2,3)%body_id = subcell_vertices(3,2)%body_id
               subcell_vertices(6,3)%body_id = subcell_vertices(7,2)%body_id
               subcell_vertices(2,6)%body_id = subcell_vertices(6,2)%body_id
               subcell_vertices(3,6)%body_id = subcell_vertices(7,2)%body_id
               subcell_vertices(2,7)%body_id = subcell_vertices(7,2)%body_id
            case (3)
               subcell_vertices(3,4)%body_id = subcell_vertices(4,3)%body_id
               subcell_vertices(7,4)%body_id = subcell_vertices(8,3)%body_id
               subcell_vertices(3,7)%body_id = subcell_vertices(7,3)%body_id
               subcell_vertices(4,7)%body_id = subcell_vertices(8,3)%body_id
               subcell_vertices(3,8)%body_id = subcell_vertices(8,3)%body_id
            case (4)
               subcell_vertices(4,8)%body_id = subcell_vertices(8,4)%body_id
            case (5)
               subcell_vertices(5,6)%body_id = subcell_vertices(6,5)%body_id
               subcell_vertices(8,6)%body_id = subcell_vertices(7,5)%body_id
               subcell_vertices(5,7)%body_id = subcell_vertices(7,5)%body_id
               subcell_vertices(6,8)%body_id = subcell_vertices(7,5)%body_id 
               subcell_vertices(5,8)%body_id = subcell_vertices(8,5)%body_id 
            case (6)
               subcell_vertices(6,7)%body_id = subcell_vertices(7,6)%body_id 
            case (7)
               subcell_vertices(7,8)%body_id = subcell_vertices(8,7)%body_id 
            case (8)
            case default
            end select
         end if

      end do

999   continue

      if (make_gmv_picture) then
         call add_cell_to_gmv_file (vertices)
      end if

   end function partial_volumes

   !----------------------------------------------------------------------------

   function subdivide_cell (vertices) result (subcell_vertices)
      !-------------------------------------------------------------------------
      ! purpose:
      !    Divide a cell (defined by its vertices, in standard cell order) into
      !    2^ndim (nvc) subcells.  Only works in 3d.
      !
      ! arguments:
      !    vertices, an nvc array of vertex_info structures
      !
      ! return:
      !    nvc x nvc array of vertex_info structures
      !
      ! algorithm:
      !    calculate subcells by finding the bisection points of cell edges and
      !    the centroids of cell faces and volumes
      !
      !    if test_vertex_once flag is set, update the body_id value for any
      !    vertices that we create where we know the answer
      !-------------------------------------------------------------------------

      ! arguments 
      type (vertex_info), dimension(:), intent(in) :: vertices

      ! return value
      type (vertex_info), dimension(nvc,nvc) :: subcell_vertices

      ! local variables
      integer :: v1
      integer :: v2
      integer :: v3
      integer :: v4
      integer :: e
      integer :: f
      integer :: n
      integer :: v
      real(r8), dimension(ndim,nec) :: Xec ! Xec: coordinates of edge centers
      real(r8), dimension(ndim,nfc) :: Xfc ! Xfc: coordinates of face centers (face centroid)
      real(r8), dimension(ndim)     :: Xcc ! Xcc: coordinates of cell center  (cell centroid)

      !-------------------------------------------------------------------------

      ! Find the coordinates of edge centers Xec.  Do this for nec edges.  For
      ! each edge e, pick the vertex numbers that define the edge, then
      ! calculate the midpoint of that edge.
      do e=1, nec
         select case (e)
         case (1)                       ! edge 1
            v1=1; v2=2
         case (2)                       ! edge 2
            v1=2; v2=3
         case (3)                       ! edge 3
            v1=3; v2=4
         case (4)                       ! edge 4
            v1=4; v2=1
         case (5)                       ! edge 5
            v1=2; v2=6
         case (6)                       ! edge 6
            v1=3; v2=7
         case (7)                       ! edge 7
            v1=4; v2=8
         case (8)                       ! edge 8
            v1=1; v2=5
         case (9)                       ! edge 9
            v1=5; v2=6
         case (10)                      ! edge 10
            v1=6; v2=7
         case (11)                      ! edge 11
            v1=7; v2=8
         case (12)                      ! edge 12
            v1=8; v2=5
         end select

         ! midpoint
         do n=1, ndim        
            Xec(n,e) = (vertices(v1)%location(n) &
                     +  vertices(v2)%location(n)) / 2.0d0
         end do
      end do

      ! Find the coordinates of face centers Xfc.  Do this for nfc faces.  For
      ! each face f, pick the vertex numbers that define the face, then
      ! calculate the centroid of that face.
      do f=1, nfc
         select case (f)
         case (1)                       ! face 1 (x = 0)
            v1=4; v2=8; v3=7; v4=3
         case (2)                       ! face 2 (x = 1)
            v1=5; v2=1; v3=2; v4=6
         case (3)                       ! face 3 (y = 0)
            v1=5; v2=8; v3=4; v4=1
         case (4)                       ! face 4 (y = 1)
            v1=6; v2=2; v3=3; v4=7
         case (5)                       ! face 5 (z = 0)
            v1=3; v2=2; v3=1; v4=4
         case (6)                       ! face 6 (z = 1)
            v1=7; v2=8; v3=5; v4=6
         end select

         ! centroid
         do n=1, ndim
            Xfc(n,f) = (vertices(v1)%location(n) &
                     +  vertices(v2)%location(n) &
                     +  vertices(v3)%location(n) &
                     +  vertices(v4)%location(n)) / 4.0d0
         end do
      end do

      ! Find the coordinates of the cell centroid Xcc.  Use all the vertices to
      ! define the cell centroid.  This is an odd weighting in cells that are degenerate hexes.
      do n=1, ndim
         Xcc(n) = (vertices(1)%location(n) &
                +  vertices(2)%location(n) &
                +  vertices(3)%location(n) &
                +  vertices(4)%location(n) &
                +  vertices(5)%location(n) &
                +  vertices(6)%location(n) &
                +  vertices(7)%location(n) &
                +  vertices(8)%location(n)) / 8.0d0
      end do

      ! Assign the coordinates of the subcell vertices.
      do n=1, ndim

         ! subcell 1 (near vertex 1)
         subcell_vertices(1,1)%location(n) = vertices(1)%location(n)
         subcell_vertices(2,1)%location(n) = Xec(n,1)
         subcell_vertices(3,1)%location(n) = Xfc(n,5)
         subcell_vertices(4,1)%location(n) = Xec(n,4)
         subcell_vertices(5,1)%location(n) = Xec(n,8)
         subcell_vertices(6,1)%location(n) = Xfc(n,2)
         subcell_vertices(7,1)%location(n) = Xcc(n)
         subcell_vertices(8,1)%location(n) = Xfc(n,3)

         ! subcell 2 (near vertex 2)
         subcell_vertices(1,2)%location(n) = Xec(n,1)
         subcell_vertices(2,2)%location(n) = vertices(2)%location(n)
         subcell_vertices(3,2)%location(n) = Xec(n,2)
         subcell_vertices(4,2)%location(n) = Xfc(n,5)
         subcell_vertices(5,2)%location(n) = Xfc(n,2)
         subcell_vertices(6,2)%location(n) = Xec(n,5)
         subcell_vertices(7,2)%location(n) = Xfc(n,4)
         subcell_vertices(8,2)%location(n) = Xcc(n)

         ! subcell 3 (near vertex 3)
         subcell_vertices(1,3)%location(n) = Xfc(n,5)
         subcell_vertices(2,3)%location(n) = Xec(n,2)
         subcell_vertices(3,3)%location(n) = vertices(3)%location(n)
         subcell_vertices(4,3)%location(n) = Xec(n,3)
         subcell_vertices(5,3)%location(n) = Xcc(n)
         subcell_vertices(6,3)%location(n) = Xfc(n,4)
         subcell_vertices(7,3)%location(n) = Xec(n,6)
         subcell_vertices(8,3)%location(n) = Xfc(n,1)

         ! subcell 4 (near vertex 4)
         subcell_vertices(1,4)%location(n) = Xec(n,4)
         subcell_vertices(2,4)%location(n) = Xfc(n,5)
         subcell_vertices(3,4)%location(n) = Xec(n,3)
         subcell_vertices(4,4)%location(n) = vertices(4)%location(n)
         subcell_vertices(5,4)%location(n) = Xfc(n,3)
         subcell_vertices(6,4)%location(n) = Xcc(n)
         subcell_vertices(7,4)%location(n) = Xfc(n,1)
         subcell_vertices(8,4)%location(n) = Xec(n,7)

         ! subcell 5 (near vertex 5)
         subcell_vertices(1,5)%location(n) = Xec(n,8)
         subcell_vertices(2,5)%location(n) = Xfc(n,2)
         subcell_vertices(3,5)%location(n) = Xcc(n)
         subcell_vertices(4,5)%location(n) = Xfc(n,3)
         subcell_vertices(5,5)%location(n) = vertices(5)%location(n)
         subcell_vertices(6,5)%location(n) = Xec(n,9)
         subcell_vertices(7,5)%location(n) = Xfc(n,6)
         subcell_vertices(8,5)%location(n) = Xec(n,12)

         ! subcell 6 (near vertex 6)
         subcell_vertices(1,6)%location(n) = Xfc(n,2)
         subcell_vertices(2,6)%location(n) = Xec(n,5)
         subcell_vertices(3,6)%location(n) = Xfc(n,4)
         subcell_vertices(4,6)%location(n) = Xcc(n)
         subcell_vertices(5,6)%location(n) = Xec(n,9)
         subcell_vertices(6,6)%location(n) = vertices(6)%location(n)
         subcell_vertices(7,6)%location(n) = Xec(n,10)
         subcell_vertices(8,6)%location(n) = Xfc(n,6)

         ! subcell 7 (near vertex 7)
         subcell_vertices(1,7)%location(n) = Xcc(n)
         subcell_vertices(2,7)%location(n) = Xfc(n,4)
         subcell_vertices(3,7)%location(n) = Xec(n,6)
         subcell_vertices(4,7)%location(n) = Xfc(n,1)
         subcell_vertices(5,7)%location(n) = Xfc(n,6)
         subcell_vertices(6,7)%location(n) = Xec(n,10)
         subcell_vertices(7,7)%location(n) = vertices(7)%location(n)
         subcell_vertices(8,7)%location(n) = Xec(n,11)

         ! subcell 8 (near vertex 8)
         subcell_vertices(1,8)%location(n) = Xfc(n,3)
         subcell_vertices(2,8)%location(n) = Xcc(n)
         subcell_vertices(3,8)%location(n) = Xfc(n,1)
         subcell_vertices(4,8)%location(n) = Xec(n,7)
         subcell_vertices(5,8)%location(n) = Xec(n,12)
         subcell_vertices(6,8)%location(n) = Xfc(n,6)
         subcell_vertices(7,8)%location(n) = Xec(n,11)
         subcell_vertices(8,8)%location(n) = vertices(8)%location(n)

      end do

      ! Initialize the body identifier of each vertex to an impossible value
      ! (BODY_ID_UNKNOWN), which means we haven't figured out what body these
      ! vertices are in.

      subcell_vertices(:,:)%body_id = BODY_ID_UNKNOWN

      ! We already know the body identifiers for some of the new vertices, we
      ! can insert that data now and save the work of looking them up again.

      if (test_vertex_once) then
         do v = 1, nvc
            subcell_vertices(v,v)%body_id = vertices(v)%body_id
         end do
      end if
 
   end function subdivide_cell

   !----------------------------------------------------------------------------

   function approximation (vertices, volume)
      !-------------------------------------------------------------------------
      ! purpose:
      !    Approximate the partial volumes of a cell that we're not going to
      !    further subdivide.
      !
      ! arguments:
      !    vertices, an nvc array of vertex_info structures
      !    volume, the total volume of the cell
      !
      ! return:
      !    nbody array of partial volumes
      !
      ! algorithm:
      !    1) assign the subcell's volume to the body that the subcell's
      !    centroid belongs to
      !
      !    2) divide the volume of the subcell into partial volumes by looking
      !    at what body each vertex is in and adding 1/nvc * volume to that
      !    body's volume - this "leaks" material if there is an interface on a
      !    cell face
      !-------------------------------------------------------------------------

      ! arguments
      type(vertex_info), dimension(nvc), intent(in) :: vertices
      real(r8), intent(in) :: volume

      ! return value
      real, dimension(nbody) :: approximation

      ! local variables
      integer :: n
      real(r8), dimension(ndim) :: cell_centroid

      !-------------------------------------------------------------------------

      ! initialize the partial volumes to zero, overwrite as needed according to
      ! the approximation used
      approximation(:) = 0.0d0

      ! method 1) assign the contribution to the partial volumes to the body the
      ! subcell centroid is in

      ! find the coordinates of the cell centroid
      do n=1, ndim
         cell_centroid(n) = (vertices(1)%location(n) &
                          +  vertices(2)%location(n) &
                          +  vertices(3)%location(n) &
                          +  vertices(4)%location(n) &
                          +  vertices(5)%location(n) &
                          +  vertices(6)%location(n) &
                          +  vertices(7)%location(n) &
                          +  vertices(8)%location(n)) / 8.0d0
      end do

      ! assign the volume to the body the centroid is in
      approximation(body_id_from_vertex (cell_centroid)) = volume

      ! method 2) distribute the contribution to the partial volumes according
      ! to the number of vertices in each body

      ! do v = 1, nvc
      !    approximation(vertices(v)%body_id) = approximation(vertices(v)%body_id) + volume / nvc
      ! end do

      if (make_gmv_picture) then
         call add_cell_to_gmv_file (vertices)
      end if

   end function approximation

   !----------------------------------------------------------------------------
   
   function cell_volume (vertices)

      !-------------------------------------------------------------------------
      ! purpose:
      !    Calculate the volume of a cell defined by its vertices.  Only works
      !    in 3d
      !
      ! arguments:
      !    vertices, an nvc array of vertex_info structures
      !
      ! return:
      !    the volume of the cell
      !
      ! algorithm:
      !-------------------------------------------------------------------------

      ! arguments
      type (vertex_info), dimension(nvc), intent(in) :: vertices

      ! return value
      real(r8) :: cell_volume 

      ! local variables
      integer :: n
      integer :: f
      integer :: v1
      integer :: v2
      integer :: v3
      integer :: v4
      integer :: v5
      integer :: v6
      real(r8), dimension(ndim) :: X1
      real(r8), dimension(ndim) :: X2
      real(r8), dimension(ndim) :: X3

      !-------------------------------------------------------------------------

      ! cell_volume will be calculated by an accumulation technique, so we
      ! initialize cell_volume to zero.

      cell_volume = 0.0d0

      ! Loop over faces.  Select some vertices depending on we don't know what,
      ! in a particular order, calculate something from their positions,
      ! accumulating the cell volume.

      do f = 1, nfc
         select case (f)
         case (1)
            v1=8; v2=4; v3=7; v4=8; v5=3; v6=4
         case (2)
            v1=6; v2=2; v3=5; v4=6; v5=1; v6=2
         case (3)
            v1=5; v2=1; v3=8; v4=5; v5=4; v6=1
         case (4)
            v1=7; v2=3; v3=6; v4=7; v5=2; v6=3
         case (5)
            v1=3; v2=4; v3=2; v4=3; v5=1; v6=4
         case (6)
            v1=6; v2=5; v3=7; v4=6; v5=8; v6=5
         end select

         do n = 1, ndim
            X1(n) = vertices(v1)%location(n) + vertices(v2)%location(n)
            X2(n) = vertices(v3)%location(n) + vertices(v4)%location(n)
            X3(n) = vertices(v5)%location(n) + vertices(v6)%location(n)
         end do

         do n = 1, ndim
            select case (n)
            case (1)
               v1=1; v2=2; v3=3
            case (2)
               v1=2; v2=3; v3=1
            case (3)
               v1=3; v2=1; v3=2
            end select

            cell_volume = cell_volume + X1(v1) * (X2(v2) * X3(v3) - X3(v2) * X2(v3))

         end do
      end do

      cell_volume = cell_volume / nec

   end function cell_volume

   !----------------------------------------------------------------------------

   function body_id_from_vertex (X)
      !-------------------------------------------------------------------------
      ! given a vertex (represented as (x,y,z), return the body that
      ! the vertex is in
      !
      ! this code is ugly, undocumented, and based on what was already there
      !-------------------------------------------------------------------------

      use interfaces_module, only: Ar, Cosa, Isrftype, nbody, Nsurf, Offset, &
                                   Rotangl, Rotpt, Sgeom, Sina
      use parameter_module,  only: nrot
      use tally_module,      only: PRE_TALLY
 
      ! arguments
      real(r8), dimension(ndim), intent(in) :: X

      ! return value
      integer :: body_id_from_vertex

      ! local variables
      logical :: Lsf                             ! logical surface flag
      logical :: Lbf                             ! logical body flag
      integer :: ib
      integer :: is
      integer :: isrf
      integer :: n
      integer :: n1
      integer :: n2
      integer :: na
      real(r8) :: Total
      real(r8) :: coeff
      real(r8), dimension(ndim) :: Xloc

      logical, save :: first_time = .true.

      !-------------------------------------------------------------------------

      ! get the geometric information of the bodies once
      if (first_time) then
        call PRE_TALLY() 
        first_time = .false.
      endif 

      ! assume that everything is in the background body unless we
      ! find some other body
      body_id_from_vertex = nbody

      do ib = 1, nbody
         ! Loop over the bodies in the domain and determine
         ! whether the point falls inside any of the bodies.
         Lbf = .true.

         do is = 1, Nsurf(ib)
            ! Loop over each surface that defines the body and
            ! determine whether the point lies on the correct side
            ! of the surface.

            ! Initialize the coordinates local to the surface
            Xloc = X

            do n = 1, nrot

               ! Select the rotation axis
               na = n
               if (ndim == 2) then
                  na = 3
               end if

               ! Select the rotation plane axes (two of them).
               select case (na)
               case (1)
                  n1=2; n2=3; coeff =  1.0d0
               case (2)
                  n1=3; n2=1; coeff = -1.0d0
               case (3)
                  n1=1; n2=2; coeff =  1.0d0
               end select

               ! Rotate the surface coordinate system.
               if (Rotangl(n,is,ib) /= 0.0_r8) then
                  Total    = Cosa(n,is,ib)*(Xloc(n1)-Rotpt(n1,is,ib)) &
                           - coeff*Sina(n,is,ib)*(Xloc(n2)-Rotpt(n2,is,ib)) + Rotpt(n1,is,ib)
                  Xloc(n2) = Cosa(n,is,ib)*(Xloc(n2)-Rotpt(n2,is,ib)) &
                           + coeff*Sina(n,is,ib)*(Xloc(n1)-Rotpt(n1,is,ib)) + Rotpt(n2,is,ib)
                  Xloc(n1) = Total
               end if

            end do

            ! Translate the surface coordinate system
            do n = 1, ndim
               Xloc(n) = Xloc(n) - Offset(n,is,ib)
            end do

            ! Jump to specified surface type and compute whether
            ! the point lies on the body side of the surface.
            isrf = ABS(Isrftype(is,ib))

            select case (isrf)

            case (1)
               ! Surface is a plane which, after translation and rotation
               ! contains the local origin and is normal to one of the
               ! coordinate axes.  Parameter Sgeom(1,is,ib) specifies the
               ! normal direction to the plane, with absolute values of
               ! 1., 2., or 3. implying the local normal to the plane is
               ! the x-, y-, or z-axis, respectively.  Positive values of
               ! Isrftype(is,ib) imply the body lies on the positive side
               ! of the plane, and negative values imply the body is on the
               ! negative side.

               do n = 1, ndim
                  if (ABS(Sgeom(1,is,ib)) == REAL(n)) then
                     Lsf = Xloc(n) > 0.0_r8
                  end if
               end do
               if (Isrftype(is,ib) < 0) then
                  Lsf = .not. Lsf
               end if
               Lbf = Lbf .and. Lsf

            case (2)
               ! Surface is a right parallelopiped whose local origin
               ! is its centroid.  Its x-y-z half-lengths are half of
               ! Sgeom(1,is,ib), Sgeom(2,is,ib), and Sgeom(3,is,ib),
               ! respectively.

               Lsf = .true.
               do n = 1, ndim
                  Lsf = Lsf .and. ABS(Xloc(n)) < Ar(n,is,ib)
               end do
               if (Isrftype(is,ib) < 0) then
                  Lsf = .not. Lsf
               end if
               Lbf = Lbf .and. Lsf

            case (3:4)
               ! Surface is a sphere or ellipsoid

               Total = 0.0_r8
               do n = 1, ndim
                  Total = Total + (Ar(n,is,ib)*Xloc(n))**2
               end do
               Lsf = Total <= 1.0_r8
               if (Isrftype(is,ib) < 0) then
                  Lsf = .not. Lsf
               end if
               Lbf = Lbf .and. Lsf
            
            case (5)
               ! Surface is a right circular cylinder.
               ! Parameter Sgeom(1,is,ib) specifies the axis along
               ! which the cylinder is aligned.  Values of 1., 2., or 3.
               ! specify the x-, y-, or z-axis, respectively.  Parameter
               ! Sgeom(2,is,ib) is the cylinder radius and Sgeom(3,is,ib)
               ! is the cylinder height.

               Total = 0.0_r8
               do n = 1, ndim
                  Total = Total + (Ar(n,is,ib)*Xloc(n))**2
               end do
               Lsf = Total <= 1.0_r8
               do n = 1, ndim
                  if (ABS(Sgeom(1,is,ib)) == REAL(n)) then
                     Lsf = Lsf .and. (Xloc(n) - Sgeom(3,is,ib))*Xloc(n) <= 0.0_r8
                  end if
               end do
               if (Isrftype(is,ib) < 0) then
                  Lsf = .not. Lsf
               end if
               Lbf = Lbf .and. Lsf

            case (6)
               ! Surface is a right circular cone.  Parameter Sgeom(1,is,
               ! ib) specifies the axis along which the cone is aligned.
               ! Values of 1., 2., or 3. specify the x-, y-, or z-axis,
               ! respectively.  Parameter Sgeom(2,is,ib) is the cone height
               ! relative to the cone origin at the base.  A positive
               ! height implies the cone apex points in the positive axis
               ! direction, and while a negative height implies the apex
               ! points in the negative axis direction.  Parameter
               ! Sgeom(3,is,ib) is the radius of the cone at its base.

               do n = 1,ndim

                  select case (n)
                  case (1)
                     n1=2; n2=3
                     if (ndim == 2) then
                        n2 = n1
                     end if
                  case (2)
                     n1=3; n2=1
                     if (ndim ==2) then
                        n1 = n2
                     end if
                  case (3)
                     n1=1; n2=2
                  end select

                  if (ABS(Sgeom(1,is,ib)) == REAL(n)) then
                     Total = Sgeom(3,is,ib)/Sgeom(2,is,ib)*(Sgeom(2,is,ib) - Xloc(n))
                     Lsf = Xloc(n1)**2 + Xloc(n2)**2 <= Total**2 &
                         .and. Sgeom(2,is,ib)*(Sgeom(2,is,ib) - Xloc(n)) >= 0.0_r8 &
                         .and. Sgeom(2,is,ib)*Xloc(n) >= 0.0_r8
                  end if

               end do

               if (Isrftype(is,ib) < 0) then
                  Lsf = .not. Lsf
               end if
               Lbf = Lbf .and. Lsf

            case (7)

               ! Body is a tabular n-sided polygon, to be rotated if
               ! parameter Sgeom(1,is,ib) = 0, or translated if = 1,
               ! about the axis (1=x, 2=y, 3=z) specified by the
               ! parameter Sgeom(2,is,ib). The integer part
               ! of parameter Sgeom(3,is,ib) is npoly, the # of vertices
               ! in the polygon, given by Rtab(0:npoly,ib) and
               ! Ztab(0:npoly,ib).

               call TLS_fatal ('BODY_ID_FROM_VERTEX: Vof_Method divide with tabular bodies not yet implemented')

!                npoly = NINT(Sgeom(3,is,ib))
!                Rtab(0,ib) = Rtab(npoly,ib)
!                Ztab(0,ib) = Ztab(npoly,ib)
! 
!                ! Now rotate or translate
!                do n = 1, ndim
!                   select case (n)
!                   case (1)
!                      n1=2; n2=3
!                      if (ndim == 2) then
!                         n2 = n1
!                      end if
! 
!                   case (2)
!                      n1=3; n2=1
!                      if (ndim == 2) then
!                         n1 = n2
!                      end if
! 
!                   case (3)
!                      n1=1; n2=2
! 
!                   end select
! 
!                   if (Sgeom(1,is,ib) == 0.0_r8) then
!                      ! Rotate
!                      Xloc(n1) = SQRT(Xloc(n1)**2 + Xloc(n2)**2)
!                      na = n
!                   else
!                      ! Translate
!                      na = n2
!                   end if
! 
!                   ! Select axis of rotation or translation
!                   if (ABS(Sgeom(2,is,ib)) == REAL(n)) then
!                      call IN_OUT_CRSE (Xloc(n1), Xloc(na), Total, npoly, ib)
!                   end if
!                end do
! 
!                Lsf = Total == 1
!                if (Isrftype(is,ib) < 0) then
!                   Lsf = .not. Lsf
!                end if
!                Lbf = Lbf .and. Lsf

            case (8)
               ! Get material distribution from mesh file; whereever the Mesh_material
               ! array matches the material number of this body, we have a hit.

               call TLS_fatal ('BODY_ID_FROM_VERTEX: Vof_Method divide with mesh materials not yet implemented')

!                Lsf = .false.
!                if (ASSOCIATED(Mesh_Material)) then
!                   Lsf = Mesh_Material == Mesh_Matnum(ib)
!                   if (Isrftype(is,ib) < 0) then
!                      Lsf = .not. Lsf
!                   end if
! 
!                   Lbf = Lbf .and. Lsf
!                end if

            end select

         end do

         if (lbf) then
            body_id_from_vertex = ib
            exit 
         endif

       end do

   end function body_id_from_vertex

   !----------------------------------------------------------------------------

   subroutine add_cell_to_gmv_file (vertices)

      ! arguments
      type (vertex_info), dimension(nvc), intent(in) :: vertices

      ! local variables
      logical, save :: first_time = .true.
      integer, save :: lun
      integer, save :: n
      integer :: v

      !-------------------------------------------------------------------------

      if (first_time) then
         open (newunit=lun, file='gmv.data')
         first_time = .false.
         n = 1
      end if

      do v = 1, nvc
         write (lun, 100) vertices(v)%location(1), vertices(v)%location(2), vertices(v)%location(3)
      end do
      write (lun, 110) n, n+1, n+2, n+3, n+4, n+5, n+6, n+7
      n = n + 8

100   format ('vrtx:',3(' ',e12.3))
110   format ('cell:',8(' ',i7))

   end subroutine add_cell_to_gmv_file

   !----------------------------------------------------------------------------

end module vof_init
