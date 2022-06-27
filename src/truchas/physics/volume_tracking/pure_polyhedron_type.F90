!!
!! polyhedron_type
!!
!! This module defines an arbitrary polyhedron type, along with routines for
!! calculating volume, splitting polyhedra, locating intersections, etc.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! October 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! References:
!!     1. Hopcroft and Kahn. A Paradigm for Robust Geometric Algorithms. Algorithmica, 1992
!!

#include "f90_assert.fpp"

module pure_polyhedron_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use polygon_type
  implicit none
  private

  integer, parameter :: ndim = 3

  type, public :: pure_polyhedron
    !private
    integer               :: nVerts, nEdges, nFaces       ! number of vertices, edges, and faces
    real(r8), allocatable :: x(:,:),face_normal(:,:)      ! vertex positions and face outward normals
    integer, allocatable :: face_vid(:,:), edge_vid(:,:) ! face and edge IDs
    integer, allocatable :: vertex_faces(:,:), edge_faces(:,:), face_eid(:,:)
    real(r8)              :: vol ! should be private, but need inheritance
  contains
    procedure, private :: init_polyhedron
    procedure, private :: init_polyhedron_null
    procedure, private :: init_tet
    generic :: init => init_polyhedron, init_polyhedron_null, init_tet
    procedure :: volume
    procedure :: intersection_verts
    procedure :: split
    procedure :: volume_behind_plane
    procedure :: is_inside
    procedure :: print_data
    procedure :: tesselated
    procedure :: has_nonplanar_face
    procedure, private :: face_is_nonplanar
    procedure, private :: calculate_edge_faces
    procedure, private :: calculate_vertex_faces
    !procedure, private :: edge_containing_vertices
    procedure, private :: compute_scaled_polyhedron
    procedure, private :: polyhedron_on_side_of_plane
    procedure, private :: update_face_normals
    procedure, private :: is_valid
    !procedure, private :: remove_dangling_vertices
  end type pure_polyhedron

contains

  subroutine init_polyhedron (this, ierr, x, face_v, edge_v, face_normal, vol)

    class(pure_polyhedron),  intent(out) :: this
    integer,            intent(out) :: ierr
    real(r8),           intent(in)  :: x(:,:)
    integer,            intent(in)  :: face_v(:,:), edge_v(:,:)
    real(r8), optional, intent(in)  :: face_normal(:,:), vol

    this%nVerts = size(x,     dim=2)
    this%nEdges = size(edge_v,dim=2)
    this%nFaces = size(face_v,dim=2)
    ierr = 0

    this%x = x
    this%edge_vid = edge_v
    this%face_vid = face_v

    this%vol = merge(vol, 0.0_r8, present(vol))

    if (present(face_normal)) then
      this%face_normal = face_normal
    else
      allocate(this%face_normal(ndim,this%nFaces))
      call this%update_face_normals()
    end if

    !call this%remove_dangling_vertices ()
    ierr = this%is_valid()
    if (ierr/=0) then
      write(*,*) 'tried to create invalid polyhedron!'
      call this%print_data()
      write(*,*)
    end if

    ! if the faces are of type polygon
    ! do f = 1,this%nFaces
    !   nV = count(face_v(:,f) /= 0) ! number of vertices on this face
    !   call this%face(f)%init (x(:,face_v(1:nV,f)))
    ! end do
    ! note that by taking the cross product between edges described in a
    ! counter-clockwise manner, we guarantee the normal to be outward facing

    call this%calculate_edge_faces(ierr)
    if (ierr /= 0) return
    call this%calculate_vertex_faces(ierr)
    if (ierr /= 0) return
    !if (this%tesselated) call this%tesselate_nonplanar_faces()

  end subroutine init_polyhedron

  subroutine init_tet (this, ierr, x, face_normal, vol, set_face_normals)

    use cell_geometry, only: tet_volume
    use cell_topology, only: TET4_FACES, TET4_EDGES

    class(pure_polyhedron),  intent(out) :: this
    integer,            intent(out) :: ierr
    real(r8),           intent(in)  :: x(:,:)
    real(r8), optional, intent(in)  :: face_normal(:,:), vol
    logical, optional, intent(in) :: set_face_normals

    logical :: set_face_normalsh

    ierr = 0
    ASSERT(size(x, dim=2) == 4)

    this%nVerts = 4
    this%nEdges = 6
    this%nFaces = 4

    this%x = x
    this%edge_vid = reshape(TET4_EDGES, [2,6])
    this%face_vid = reshape(TET4_FACES, [3,4])
    this%face_eid = reshape([1,5,3, 4,6,5, 3,6,2, 2,4,1], [3,4])
    this%edge_faces = reshape([1,4, 3,4, 1,3, 2,4, 1,2, 2,3], [2,6])
    this%vertex_faces = reshape([1,3,4, 1,2,4, 2,3,4, 1,2,3], [3,4])

    if (present(face_normal)) then
      this%face_normal = face_normal
    else
      if (present(set_face_normals)) then
        set_face_normalsh = set_face_normals
      else
        set_face_normalsh = .true.
      end if

      if (set_face_normalsh) then
        allocate(this%face_normal(3,this%nFaces))
        call this%update_face_normals()
      end if
    end if

    this%vol = merge(vol, tet_volume(this%x), present(vol))

    if (this%vol < 0) then
      call this%print_data()
    end if

    ASSERT(this%vol >= 0)

  end subroutine init_tet

  function is_valid (this) result(ierr)

    class(pure_polyhedron), intent(in) :: this
    integer :: ierr

    ierr = 0

    if (this%nVerts < 4) ierr = 1 ! must have >= 4 vertices
    if (this%nFaces < 3) ierr = 1 ! must have >= 3 faces
    if (any(count(this%face_vid /= 0, dim=1) < 3)) ierr = 1 ! each face must have >= 3 vertices

    ! ! no dangling vertices
    ! ! every vertex must be connected to at least three other vertices
    ! do v = 1,this%nVerts
    !   if (count(this%edge_vid==v)<3) then
    !     write(*,*) 'dangling vertex ',v
    !     ierr = 1
    !   end if
    ! end do

  end function is_valid

  subroutine update_face_normals (this)

    use cell_geometry, only: cross_product, normalized

    class(pure_polyhedron), intent(inout) :: this

    integer :: f

    if (.not.allocated(this%face_normal)) allocate(this%face_normal(3,this%nFaces))

    do f = 1,this%nFaces
      this%face_normal(:,f) = normalized(cross_product (&
          this%x(:,this%face_vid(2,f)) - this%x(:,this%face_vid(1,f)), &
          this%x(:,this%face_vid(3,f)) - this%x(:,this%face_vid(1,f))))
    end do

  end subroutine update_face_normals

  ! subroutine remove_dangling_vertices (this)

  !   class(pure_polyhedron), intent(inout) :: this

  !   integer :: v,vv

  !   ! no dangling vertices
  !   ! every vertex must be connected to at least three other vertices
  !   ! vertices connected to only two other vertices may be removed, and
  !   ! its connected vertices joined
  !   ! WARNING: this assumes the dangling vertex is still attached to two other vertices,
  !   !          and doesn't consider vertices with only one edge
  !   v = 1
  !   do while (v <= this%nVerts)
  !     if (count(this%edge_vid==v)<3) then
  !       ! delete this vertex
  !       do vv = v+1,this%nVerts
  !         this%x(:,vv-1) = this%x(:,vv)
  !       end do
  !       this%nVerts = this%nVerts - 1

  !       ! delete the two edges it is connected to, and add a new one between those two vertices
  !       call delete_edges_containing_vertex (this,v)

  !       ! remove the vertex from the face it is connected to
  !       call delete_faces_containing_vertex (this,v)
  !     else
  !       v = v+1
  !     end if
  !   end do

  ! contains

  !   subroutine delete_edges_containing_vertex (this,v)

  !     class(pure_polyhedron), intent(inout) :: this
  !     integer, intent(in) :: v

  !     integer :: e,ee,vn(2),j

  !     e = 1; j=1
  !     do while (e <= this%nEdges)
  !       if (any(this%edge_vid(:,e)==v)) then
  !         ! store the other vertex connected to this edge
  !         vn(j) = this%edge_vid(1,e)
  !         if (vn(j)==v) vn(j) = this%edge_vid(2,e)
  !         j = j+1

  !         ! delete edges containing this vertex
  !         do ee = e+1,this%nEdges
  !           this%edge_vid(:,ee-1) = this%edge_vid(:,ee)
  !         end do
  !         this%nEdges = this%nEdges - 1
  !       else
  !         e = e+1
  !       end if
  !     end do

  !     ! add the vertex connecting [v1,v2]
  !     this%nEdges = this%nEdges + 1
  !     this%edge_vid(:,this%nEdges) = vn

  !     ! update the edge structure
  !     where (this%edge_vid>v) this%edge_vid = this%edge_vid - 1

  !   end subroutine delete_edges_containing_vertex

  !   subroutine delete_faces_containing_vertex (this,v)

  !     use array_utils, only: first_true_loc

  !     class(pure_polyhedron), intent(inout) :: this
  !     integer, intent(in) :: v

  !     integer :: f,ff,fv,nV

  !     f = 1
  !     do while (f <= this%nFaces)
  !       if (any(this%face_vid(:,f)==v)) then
  !         if (count(this%face_vid(:,f)/=0)==3) then
  !           ! if the face only has three vertices, delete the entire face
  !           do ff = f+1,this%nFaces
  !             this%face_vid(:,ff-1) = this%face_vid(:,ff)
  !             this%face_normal(:,ff-1) = this%face_normal(:,ff)
  !           end do
  !           this%nFaces = this%nFaces - 1
  !         else
  !           nV = count(this%face_vid(:,f)/=0)
  !           fv = first_true_loc (this%face_vid(:,f)==v)
  !           do ff = fv+1,nV
  !             this%face_vid(ff-1,f) = this%face_vid(ff,f)
  !           end do
  !           this%face_vid(nV,f) = 0
  !           f = f+1
  !         end if
  !       else
  !         f = f+1
  !       end if
  !     end do

  !     ! update the face structure
  !     where (this%face_vid>v) this%face_vid = this%face_vid - 1

  !   end subroutine delete_faces_containing_vertex

  ! end subroutine remove_dangling_vertices

  subroutine init_polyhedron_null (this)

    class(pure_polyhedron), intent(out) :: this

    this%nVerts = 0
    this%nEdges = 0
    this%nFaces = 0
    this%vol = 0

  end subroutine init_polyhedron_null

  subroutine compute_scaled_polyhedron (this, scaled, volume_factor)

    class(pure_polyhedron), intent(in) :: this
    type(pure_polyhedron), intent(out) :: scaled
    real(r8), intent(out) :: volume_factor

    integer :: v
    real(r8) :: x0(3), xl(3)

    scaled = this

    x0 = minval(scaled%x,dim=2)
    !x0 = [0.0_r8, 0.0_r8, 0.0_r8] ! WARN: no offsetting
    xl = maxval(scaled%x,dim=2) - x0
    xl = [1.0_r8, 1.0_r8, 1.0_r8] ! WARN: no scaling
    volume_factor = product(xl)
    do v = 1,scaled%nVerts
      scaled%x(:,v) = (scaled%x(:,v) - x0) / xl
    end do
    call scaled%update_face_normals()

  end subroutine compute_scaled_polyhedron

  ! This function calculates the volume of a polehedron following the
  ! algorithm presented by [1]. It calculates the sum of the surface
  ! integral of x over all faces, which is equal to the volume by the
  ! divergence theorem.
  !
  ! note 1: this is done so we can find the volume of very tiny ones without hitting precision limits
  !         in some cases, very tiny polyhedra (without this trick) would produce negative volumes
  !         because their volume was on the order of floating point errors. Scaling the polyhedron
  !         up to some normalized size before calculating the volume, then scaling back, seems to
  !         counter this rather well, but seems a bit tricky. It might instead be worth just setting
  !         the volume of a polyhedron to zero if it is calculated as below zero and very close to
  !         zero. Note this occurs even though the polyhedra splitting routine is "robust" and will
  !         not allow vertices to be within some distance alpha (=1e-9). These tiny distances can
  !         still produce volumes on the order of e-27, far below the double precision limit of e-16.
  real(r8) function volume (this)

    use cell_geometry,   only: cross_product
    use near_zero_function

    class(pure_polyhedron), intent(inout) :: this
    !integer, intent(out) :: ierr

    real(r8), parameter :: alittle = 1e-9_r8
    integer :: f,nV,v,i,n
    real(r8) :: volume_factor, xn(3), xi(3), xc(3), vt, area, t
    type(pure_polyhedron) :: scaled

    !ierr = 0
    volume = this%vol
    if (.not.allocated(this%x) .or. volume > 0.0_r8) return

    ! scale the polyhedron (see note 1)
    ! scaled = this
    ! volume_factor = 1
    call this%compute_scaled_polyhedron(scaled, volume_factor)

    ! sum up the integral of n_x*x over all faces (could just as easily be any other direction)
    volume = 0
    do f = 1,scaled%nFaces
      nV = count(scaled%face_vid(:,f) /= 0) ! number of vertices on this face
      xc = sum(scaled%x(:,scaled%face_vid(:nV,f)), dim=2) / nV
      area = 0
      vt = 0
      do v = 1,nV
        i = scaled%face_vid(v,f)
        n = scaled%face_vid(modulo(v,nV)+1,f)
        xi = scaled%x(:,i) - xc
        xn = scaled%x(:,n) - xc
        t = dot_product(cross_product(xi, xn), scaled%face_normal(:,f))
        vt = vt + dot_product(xi + xn, scaled%face_normal(:,f)) * t
        area = area + t
      end do
      area = area / 2
      volume = volume + vt / 18 + dot_product(xc, scaled%face_normal(:,f)) * area / 3
    end do
    volume = volume * volume_factor
    this%vol = volume

    if (this%vol < 0) then
      if (near_zero(this%vol, 1e5_r8*alittle)) then
        ! if this polyhedron has a volume of almost zero, make it zero
        ! this seems to be necessary for very tiny polyhedrons,
        ! where floating point errors may make the volume calculation drop below zero
        volume = 0.0_r8
        this%vol = 0.0_r8
        deallocate(this%x)
        this%nVerts = 0
      else
        volume = 0.0_r8
        this%vol = 0.0_r8
        deallocate(this%x)
        this%nVerts = 0

        print *, 'vf: ',volume_factor
        call this%print_data ()
        print *

        !ierr = 1
        !call tls_fatal ("calculated negative polyhedron volume!")
      end if
    end if

    ! if (this%vol == 0) then
    !   print *, 'vf: ',volume_factor

    !   do f = 1,scaled%nFaces
    !     nV = count(scaled%face_vid(:,f) /= 0) ! number of vertices on this face
    !     tmp = 0
    !     do v = 1,nV
    !       tmp = tmp + cross_product(&
    !           scaled%x(:,scaled%face_vid(v,f)), &
    !           scaled%x(:,scaled%face_vid(modulo(v,nV)+1,f)))
    !     end do
    !     print *, dot_product(scaled%face_normal(:,f),scaled%x(:,scaled%face_vid(1,f))) &
    !         * dot_product(scaled%face_normal(:,f),tmp)
    !     print *, 'n: ', scaled%face_normal(:,f)
    !     print *, 'x: ', scaled%x(:,scaled%face_vid(1,f))
    !     print *, dot_product(scaled%face_normal(:,f),scaled%x(:,scaled%face_vid(1,f)))
    !     print *, dot_product(scaled%face_normal(:,f),tmp)
    !     print *
    !   end do
    ! end if

  end function volume

  ! Given an equation of a plane and a polyhedron, return a polygon from the
  ! points where the plane intersects the polyhedron edges
  type(polygon) function intersection_verts (this,P,v_assoc_pe)

    use plane_type
    !use set_type

    class(pure_polyhedron), intent(in)  :: this
    class(plane),      intent(in)  :: P
    integer, optional, intent(out) :: v_assoc_pe(:)

    integer  :: e,Nintersections, on_point
    real(r8) :: x(ndim,this%nEdges),intx(ndim)
    !type(set_integer) :: vertex_faces(this%nEdges) !vertex_faces(size(this%vertex_faces, dim=1), this%nEdges)

    if (present(v_assoc_pe)) v_assoc_pe = -1

    ! loop through all edges
    Nintersections = 0
    do e = 1,this%nEdges
      ! check if the P intersects this edge
      if (P%intersects(this%x(:,this%edge_vid(:,e)))) then
        ! if it does, find the point where they intersect
        call P%intersection_point(intx, on_point, this%x(:,this%edge_vid(:,e)))
        ! print *, on_point, norm2(intx - this%x(:,this%edge_vid(1,e))), norm2(intx - this%x(:,this%edge_vid(2,e)))
        ! print *, 'edge: ',this%edge_vid(:,e), norm2(this%x(:,2) - this%x(:,3))
        ! print *

        if (on_point > 0 .and. Nintersections>0 .and. containsPoint(intx,x(:,1:Nintersections))) then
          ! if the point was already found, note that this edge intersects with that found point
          ! this particularly comes into effect when we intersect a vertex
          if (present(v_assoc_pe)) v_assoc_pe(e) = pointIndex(intx, x(:,1:Nintersections))
        else
          Nintersections = Nintersections + 1
          x(:,Nintersections) = intx
          if (present(v_assoc_pe)) v_assoc_pe(e) = Nintersections

          ! if (on_point > 0) then
          !   nf = count(this%vertex_faces(:,this%edge_vid(on_point,e)) > 0)
          !   call vertex_faces(Nintersections)%add(this%vertex_faces(:nf,this%edge_vid(on_point,e)))
          ! else
          !   call vertex_faces(Nintersections)%add(this%edge_faces(:,e))
          ! end if

        end if
      end if
    end do

    ! pass the intersection points to the polygon constructor
    if (Nintersections>2) then
      !print *, norm2(x(:,1) - x(:,2))
      call intersection_verts%init (x(:,1:Nintersections))

      ! this probably doesn't need to be called every time this function is used
      ! index_sort = intersection_verts%sort_order()

      ! !call convex_corrections()

      ! ! update arrays with new sorting
      ! !print *, 'sort0'
      ! intersection_verts%x = intersection_verts%x(:,index_sort)
      ! !print *, 'sort1'

      ! index_sort = invert (index_sort)
      ! do i = 1,size(v_assoc_pe)
      !   if (v_assoc_pe(i)>0) v_assoc_pe(i) = index_sort(v_assoc_pe(i))
      ! end do
      ! !print *, 'sort2', v_assoc_pe

      call intersection_verts%order(v_assoc_pe)

      call intersection_verts%update_plane_normal()

      ! make sure the vertices are ordered counter-clockwise
      if (dot_product(intersection_verts%norm,P%normal) < 0) then
        ! reverse both x and v_assoc_pe
        intersection_verts%x = intersection_verts%x(:,size(intersection_verts%x, dim=2):1:-1)

        if (present(v_assoc_pe)) then
          do e = 1,size(v_assoc_pe)
            if (v_assoc_pe(e)>0) v_assoc_pe(e) = Nintersections - v_assoc_pe(e)+1
          end do
        end if

        call intersection_verts%update_plane_normal ()
      end if
      !print *, 'sort3'
    ! else
    !   print *, Nintersections, " intersection points"
    !   call this%print_data()
    !   call P%print_data()
    !   call tls_fatal ("not enough intersection points")
    end if

  ! contains

  !   subroutine convex_corrections()

  !     integer :: i, in, v, p, s, n, j
  !     logical :: intersects(20)

  !     ! print *, 'nv: ', intersection_verts%nVerts
  !     ! do j = 1,intersection_verts%nVerts
  !     !   v = index_sort(j)
  !     !   p = index_sort(modulo(j-2,intersection_verts%nVerts)+1)
  !     !   intersects(j) = .not.vertex_faces(v).intersects(vertex_faces(p))
  !     ! end do
  !     ! print *, intersects(:intersection_verts%nVerts)

  !     ! ! sort the vertices (correct the previous sorting)
  !     ! !print '(a,12i6)', 'sort0', index_sort
  !     ! do i = 1,intersection_verts%nVerts
  !     !   v = index_sort(i)
  !     !   p = index_sort(modulo(i-2,intersection_verts%nVerts)+1)

  !     !   ! if this vertex shares no faces with the
  !     !   ! previous vertex, there is an order mismatch
  !     !   ! => swap this vertex with one of the next ones
  !     !   if (.not.vertex_faces(v).intersects(vertex_faces(p))) then
  !     !     print *, i
  !     !     do s = 0,intersection_verts%nVerts - i - 1
  !     !       n = index_sort(modulo(i+s,intersection_verts%nVerts)+1)
  !     !       print *, s, vertex_faces(n).intersects(vertex_faces(p))
  !     !       if (vertex_faces(n).intersects(vertex_faces(p))) exit
  !     !     end do
  !     !     !s = 0
  !     !     in = modulo(i+s,intersection_verts%nVerts)+1
  !     !     n = index_sort(in)
  !     !     index_sort(i) = n
  !     !     index_sort(in) = v

  !     !     do j = 1,intersection_verts%nVerts
  !     !       v = index_sort(j)
  !     !       p = index_sort(modulo(j-2,intersection_verts%nVerts)+1)
  !     !       intersects(j) = .not.vertex_faces(v).intersects(vertex_faces(p))
  !     !     end do
  !     !     print *, intersects(:intersection_verts%nVerts)
  !     !     print *
  !     !   end if
  !     ! end do

  !     ! sort the vertices (correct the previous sorting)
  !     !print '(a,12i6)', 'sort0', index_sort
  !     do i = 1,intersection_verts%nVerts-2
  !       v = index_sort(i)
  !       p = index_sort(modulo(i-2,intersection_verts%nVerts)+1)

  !       ! if this vertex shares no faces with the
  !       ! previous vertex, there is an order mismatch
  !       ! => swap this vertex with one of the next ones
  !       if (.not.vertex_faces(v)%intersects(vertex_faces(p))) then
  !         !print *, i
  !         do s = 1,intersection_verts%nVerts - i - 1
  !           n = index_sort(i+s)
  !           !print *, s, vertex_faces(n)%intersects(vertex_faces(p))
  !           if (vertex_faces(n)%intersects(vertex_faces(p))) exit
  !         end do
  !         !s = 0
  !         in = i+s
  !         n = index_sort(in)
  !         index_sort(i) = n
  !         index_sort(in) = v

  !         do j = 1,intersection_verts%nVerts
  !           v = index_sort(j)
  !           p = index_sort(modulo(j-2,intersection_verts%nVerts)+1)
  !           intersects(j) = .not.vertex_faces(v)%intersects(vertex_faces(p))
  !         end do
  !         !print *, intersects(:intersection_verts%nVerts)
  !         !print *
  !       end if
  !     end do


  !     ! do j = 1,intersection_verts%nVerts
  !     !   v = index_sort(j)
  !     !   p = index_sort(modulo(j-2,intersection_verts%nVerts)+1)
  !     !   intersects(j) = .not.vertex_faces(v).intersects(vertex_faces(p))
  !     ! end do
  !     ! print *, intersects(:intersection_verts%nVerts)
  !     ! print *
  !     ! print *
  !     !stop

  !     ! print '(a,12i6)', 'sort5', index_sort
  !     ! print *, 'lens: ', size(intersection_verts%x, dim=1), size(v_assoc_pe)

  !   end subroutine convex_corrections

  end function intersection_verts

  function tesselated(this)

    class(pure_polyhedron), intent(in) :: this
    type(pure_polyhedron), allocatable :: tesselated(:)

    real(r8) :: xtet(3,4)
    integer :: f, v, nV, t, ierr

    ! calculate cell center
    xtet(:,1) = sum(this%x, dim=2) / this%nVerts

    ! for each vertex on each edge of each face, we create a tet
    allocate(tesselated(count(this%face_vid>0)))
    t = 0
    do f = 1,this%nFaces
      nV = count(this%face_vid(:,f) > 0)
      xtet(:,2) = sum(this%x(:,this%face_vid(:nV,f)), dim=2) / nV
      do v = 1,nV
        ! xtet(:,3) = this%x(:,this%face_vid(v,f))
        ! xtet(:,4) = sum(this%x(:,this%face_vid([v,modulo(v,nV)+1],f)), dim=2) / 2
        ! t = t+1
        ! call tesselated(t)%init_tet(ierr, xtet)

        ! xtet(:,3) = xtet(:,4)
        ! xtet(:,4) = this%x(:,this%face_vid(modulo(v,nV)+1,f))
        ! t = t+1
        ! call tesselated(t)%init_tet(ierr, xtet)

        xtet(:,3) = this%x(:,this%face_vid(v,f))
        xtet(:,4) = this%x(:,this%face_vid(modulo(v,nV)+1,f))
        t = t+1
        call tesselated(t)%init_tet(ierr, xtet)
      end do
    end do

    ASSERT(t == size(tesselated))

  end function tesselated

  ! calculate a list of faces touching each vertex
  subroutine calculate_vertex_faces(this, ierr)

    class(pure_polyhedron), intent(inout) :: this
    integer, intent(out) :: ierr

    integer :: v, f, j(this%nVerts), nV, vid
    integer, allocatable :: tmp(:,:)

    ierr = 0
    allocate(tmp(10, this%nVerts))

    j = 0
    do f = 1,this%nFaces
      nV = count(this%face_vid(:,f) > 0)
      do v = 1,nV
        vid = this%face_vid(v,f)
        j(vid) = j(vid) + 1
        tmp(j(vid),vid) = f
      end do
    end do

    if (any(j>10)) call Tls_Fatal("not large enough initial allocation of vertex_faces")
    this%vertex_faces = tmp(:maxval(j),:)

  end subroutine calculate_vertex_faces

  subroutine calculate_edge_faces(this, ierr)

    class(pure_polyhedron), intent(inout) :: this
    integer, intent(out) :: ierr

    integer :: e, f, nV, v, j(this%nEdges)

    ierr = 0
    allocate(this%face_eid(maxval(count(this%face_vid > 0, dim=2)), this%nFaces), &
        this%edge_faces(2, this%nEdges))
    this%edge_faces = 0; this%face_eid = 0

    j = 0
    faces: do f = 1,this%nFaces
      nV = count(this%face_vid(:,f) > 0)
      do v = 1,nV
        e = pair_index(this%face_vid([v, modulo(v,nV)+1],f), this%edge_vid)
        if (e == -1) then
          ierr = 1
          exit faces
        end if

        j(e) = j(e) + 1
        if (j(e) > 2) then
          ierr = 1
          exit faces
        end if

        this%edge_faces(j(e), e) = f
        this%face_eid(v,f) = e
      end do
    end do faces

    if (ierr==1 .or. any(this%edge_faces==0)) then
      ierr = 1
      print *, 'edge_faces failure'
      ! print *, this%face_vid([v, modulo(v,nV)+1],f)
      call this%print_data()
      ! call Tls_Fatal("could not find pair in list")
    end if

  end subroutine calculate_edge_faces

  pure integer function pair_index(pair, list)

    integer, intent(in) :: pair(:), list(:,:)

    integer :: revpair(2)

    revpair = pair(size(pair):1:-1)
    do pair_index = 1,size(list, dim=2)
      if (all(pair==list(:,pair_index)) .or. all(revpair==list(:,pair_index))) return
    end do
    pair_index = -1

  end function pair_index

  ! return a list of the edge ids for edges intersected by the plane
  function intersected_edges (this,P) result(inte)
    use plane_type

    class(pure_polyhedron), intent(in) :: this
    class(plane),      intent(in) :: P

    integer :: e
    logical :: inte(this%nEdges)

    do e = 1,this%nEdges
      inte(e) = P%intersects(this%x(:,this%edge_vid(:,e)))
    end do

  end function intersected_edges

  ! return two polyhedrons produced by dividing a given polyhedron with a plane
  ! the first element returned is in front of the plane
  ! the second element returned is behind the plane
  !
  ! note 1: this part is not part of the Hopcroft and Kahn algorithm
  !         When splitting a very thin polyhedron, the intersection
  !         vertices may end up within alpha of each other. In this
  !         case, say both polyhedra are null, since we are within
  !         the cutvof anyways.
  subroutine split (this,P,split_poly,ierr)

    use plane_type

    class(pure_polyhedron), intent(in) :: this
    class(plane),      intent(in) :: P
    type(pure_polyhedron), intent(inout) :: split_poly(:)
    !type(polygon_box), intent(out) :: interface_polygons ! interface polygon
    integer,          intent(out) :: ierr

    integer       :: v, v_assoc_pe(this%nEdges),side(this%nVerts)
    type(polygon) :: intpoly
    real(r8)      :: dist, tmp1, tmp2

    ASSERT(size(split_poly)==2)
    ierr = 0
    !interface_polygons%n_elements = 0

    ! check which side of the plane each vertex lies
    ! vertices within distance alpha of the plane are considered to lie on it
    !  dist  >  alpha => side =  1
    !  dist  < -alpha => side = -1
    ! |dist| <  alpha => side =  0
    do v = 1,this%nVerts
      dist = P%signed_distance (this%x(:,v)) ! calculate the signed distance from the plane
      side(v) = merge(int(sign(1.0_r8, dist)), 0, abs(dist) > alpha) ! decide where it lies
    end do

    if (.not.any(side<0)) then
      split_poly(1) = this
      call split_poly(2)%init ()
    else if (.not.any(side>0)) then
      split_poly(2) = this
      call split_poly(1)%init ()
    else
      intpoly = this%intersection_verts (P,v_assoc_pe)

      if (intpoly%nVerts > 2) then
        split_poly(1) = this%polyhedron_on_side_of_plane (P,  1, side, intpoly, v_assoc_pe, ierr)
        if (ierr/=0) call tls_fatal ("polyhedron split failed: invalid child")
        split_poly(2) = this%polyhedron_on_side_of_plane (P, -1, side, intpoly, v_assoc_pe, ierr)
        if (ierr/=0) call tls_fatal ("polyhedron split failed: invalid child")

        ! interface_polygons%n_elements = 1
        ! interface_polygons%elements = [intpoly]
      else ! see note 1
        call split_poly(1)%init ()
        call split_poly(2)%init ()
      end if
    end if

    ! if any of the polyhedrons have a face with less than 3 vertices, throw an error
    do v = 1,2
      if (allocated(split_poly(v)%face_vid)) then
        if (any(count(split_poly(v)%face_vid /= 0, dim=1) < 3)) then
          call split_poly(v)%print_data (normalized=.true.)
          write(*,*)
          ierr = 1
        end if
      end if
    end do
    if (ierr/=0) call tls_fatal ("polyhedron split failed--one of the children has an invalid face")

    ! if either of the volumes are less than zero, throw an error
    ! TODO: this is a more expensive check. might be worth skipping when sufficiently confident.
    tmp1 = split_poly(1)%volume()
    tmp2 = split_poly(2)%volume()
    if (tmp1 < 0.0_r8 .or. tmp2 < 0.0_r8) then
      write(*,*)
      write(*,*) 'parent:'
      call this%print_data ()
      write(*,*)

      write(*,*) 'child1:'
      call split_poly(1)%print_data ()
      write(*,*)

      write(*,*) 'child2:'
      call split_poly(2)%print_data ()
      write(*,*)

      write(*,*) 'other data:'
      write(*,'(15i3)') side
      call P%print_data()

      write(*,*)
      call intpoly%print_data ()
      write(*,*)

      write(*,*) 'problematic vols: ',tmp1,tmp2
      ierr = 1
      !call tls_fatal ('polyhedron split failed: invalid volume')
    end if

  end subroutine split

  ! return the volume behind (opposite normal) a plane and inside the polyhedron
  real(r8) function volume_behind_plane (this,P,ierr)

    use ieee_arithmetic, only: ieee_is_nan
    use plane_type

    class(pure_polyhedron), intent(inout) :: this
    class(plane),      intent(in) :: P
    integer,           intent(out) :: ierr

    real(r8)         :: dist
    integer          :: v, side(this%nVerts), v_assoc_pe(this%nEdges)
    type(pure_polyhedron) :: behind
    type(polygon)    :: intpoly

    ierr = 0

    ! check which side of the plane each vertex lies
    ! vertices within distance alpha of the plane are considered to lie on it
    !  dist  >  alpha -> side =  1
    !  dist  < -alpha -> side = -1
    ! |dist| <  alpha -> side =  0
    do v = 1,this%nVerts
      dist = P%signed_distance (this%x(:,v)) ! calculate the signed distance from the plane
      side(v) = merge(int(sign(1.0_r8, dist)), 0, abs(dist) > alpha) ! decide where it lies
    end do

    if (.not.any(side>0)) then
      volume_behind_plane = this%volume()
    else if (any(side>0) .and. any(side<0)) then
      intpoly = this%intersection_verts (P,v_assoc_pe)
      if (intpoly%nVerts > 2) then
        behind = this%polyhedron_on_side_of_plane (P, -1, side, intpoly, v_assoc_pe,ierr)
        if (ierr /= 0) then
          call dumpData ()
          write(*,*) 'volume_behind_plane: polyhedron split failed'
          volume_behind_plane = 0
          return
        end if

        ! if any of the polyhedrons have a face with less than 3 vertices, throw an error
        if (allocated(behind%face_vid)) then
          if (any(count(behind%face_vid /= 0, dim=1) < 3)) then
            call dumpData ()
            call tls_fatal ("polyhedron split failed: invalid face")
          end if
        end if

        ! calculate the volume of the polyhedron behind the plane
        volume_behind_plane = behind%volume ()

        if (ieee_is_nan(volume_behind_plane) .or. volume_behind_plane < 0) then
          call dumpData ()
          write(*,*) 'problematic vol: ',volume_behind_plane
          call tls_fatal ('polyhedron split failed: invalid volume')
        end if
      else
        volume_behind_plane = 0
      end if
    else
      volume_behind_plane = 0
    end if

  contains

    subroutine dumpData ()
      write(*,*)
      write(*,*) 'parent:'
      call this%print_data ()
      write(*,*)

      write(*,*) 'child:'
      call behind%print_data ()
      write(*,*)

      write(*,*) 'other data:'
      write(*,'(15i3)') side
      call P%print_data()

      write(*,*)
      call intpoly%print_data ()
      write(*,*)
    end subroutine dumpData

  end function volume_behind_plane

  ! WARNING: need to update this routine to pass face normals down to the child polyhedron
  ! Reference [1]
  ! note 1: In this case, the polyhedron may not have been clearly all within the half-space.
  !         Some vertices on the face may have landed on each side of the plane, but if
  !         more than three landed on the plane itself, we consider this equivalent to
  !         the entire face landing on the plane.
  type(pure_polyhedron) function polyhedron_on_side_of_plane (this,P,valid_side,side,intpoly,v_assoc_pe,ierr)

    use plane_type

    class(pure_polyhedron), intent(in) :: this
    type(plane),       intent(in) :: P
    integer,           intent(in) :: side(:)       ! gives which side of the plane vertices lie on
    integer,           intent(in) :: valid_side    ! +/- 1, indicating the side of the plane we want
    type(polygon),     intent(in) :: intpoly       ! polygon of intersection coordinates
    integer,           intent(in) :: v_assoc_pe(:) ! intersection polygon vertex id for a given parent edge id
    integer,           intent(out) :: ierr

    integer :: pcf, nVerts, nParVerts, nEdges, nFaces, tmp, &
         p2c_vid(this%nVerts), & ! parent to child vertex id conversion table for cases they correspond
         edge_vid(2,this%nEdges+intpoly%nVerts), & ! can have intpoly%nVerts more edges than parent
         face_vid(max(size(this%face_vid,dim=1)+2,intpoly%nVerts),this%nFaces+1) ! could have 1 more face and faces could be intpoly%nverts longer than parent

    ! note: an updated planar face can only include 1 more node than the original,
    !       but I'm not sure if there is a limit to how many nodes the entirely new face can have.
    !       For cubes the number is 2.
    real(r8)      :: x(3,this%nVerts+intpoly%nVerts), face_normal(3, this%nFaces+1)

    call polyhedron_on_side_of_plane%init ()

    if (.not.any(side==-valid_side)) then ! the entire polyhedron is within the half-space
      polyhedron_on_side_of_plane = this
    else if (any(side==-valid_side) .and. any(side==valid_side)) then ! the polyhedron is split
      pcf = plane_coinciding_face(this,side,valid_side)
      if (pcf > 0) then
        ! the polyhedron really landed entirely inside or entirely outside the half-space (see note 1)
        if (dot_product(this%face_normal(:,pcf), P%normal) > 0.0_r8) &
            polyhedron_on_side_of_plane = this
      else
        ! make a list of vertices for the new polyhedron
        !print *, 'split0'
        call generate_new_verts (x,p2c_vid,nVerts,nParVerts, this,side,v_assoc_pe,valid_side)

        ! find new edges, knowing the original structure and the interface polygon
        !print *, 'split1'
        call find_edges (edge_vid,nEdges, this,side,valid_side,p2c_vid,nParVerts,nVerts,v_assoc_pe)

        ! construct a set of faces from the edge information
        ! note these will not be in a particular order, like pececillo would expect for hexes
        !print *, 'split2'
        call find_faces (face_vid,face_normal,nFaces,ierr, &
            this,side,valid_side,p2c_vid,nParVerts,nVerts,intpoly%nVerts,v_assoc_pe, P%normal)
        if (ierr /= 0) call fatal()
        !print *, 'split3'

        ! initialize final polyhedron
        tmp = maxval(count(face_vid(:,:) /= 0,dim=1)) ! the maximum number of vertices on a face
        if (nVerts < 3) then
          call this%print_data()
          write(*,*)
          call P%print_data()
          call tls_fatal ("not enough vertices for a polygon!")
        end if
        call polyhedron_on_side_of_plane%init (ierr, x(:,1:nVerts), face_vid(1:tmp,1:nFaces), &
            edge_vid(:,1:nEdges), face_normal(:,1:nFaces))
        !print *, 'split4'

        !write(*,*) 'poly', nVerts, tmp, nFaces, nEdges
        if (ierr /= 0) call fatal()
      end if
    end if

  contains

    subroutine fatal ()
      print *, 'parent: '
      call this%print_data()
      write(*,*)

      write(*,*) 'other data:'
      write(*,'(15i3)') side
      call P%print_data()

      write(*,*)
      call intpoly%print_data ()
      write(*,*)

      call tls_fatal ("polyhedron_on_side_of_plane failed: invalid child polyhedron")
    end subroutine fatal

    ! finds the polyhedron face which coincides with the face, if one exists
    integer function plane_coinciding_face (this,side,valid_side)

      class(pure_polyhedron), intent(in) :: this
      integer, intent(in) :: side(:), valid_side

      integer :: nV

      do plane_coinciding_face = 1,this%nFaces
        nV = count(this%face_vid(:,plane_coinciding_face) > 0)
        ! the face coincides with the plane if at least 3 of its vertices lie on the plane
        if (count(side(this%face_vid(1:nV,plane_coinciding_face))==0)>=3) return
      end do
      plane_coinciding_face = -1

    end function plane_coinciding_face

    ! note 1: these are added entirely at the end, rather than including the parent vertices
    !         that lie on the plane in the above loop because having all these vertices together
    !         makes it easy for constructing the new face later.
    subroutine generate_new_verts (x,p2c_vid,nVerts,nParVerts, this,side,v_assoc_pe,valid_side)

      use near_zero_function

      real(r8),         intent(out) :: x(:,:)
      integer,          intent(out) :: p2c_vid(:), nVerts, nParVerts
      type(pure_polyhedron), intent(in)  :: this
      integer,          intent(in)  :: side(:), v_assoc_pe(:), valid_side

      integer :: v,iv

      nVerts = 0; p2c_vid = 0

      ! first get the vertices fully on the valid side of the intersection plane
      do v = 1,this%nVerts
        if (side(v)==valid_side) then
          nVerts = nVerts+1
          x(:,nVerts) = this%x(:,v)
          p2c_vid(v) = nVerts ! parent to child vertex id
        end if
      end do
      nParVerts = nVerts ! number of vertices acquired from parent

      ! then get the vertices from the plane-polyhedron intersection (see note 1)
      x(:,nParVerts+1:nParVerts+intpoly%nVerts) = intpoly%x
      nVerts = nVerts + intpoly%nVerts

      if (nVerts < 3) call tls_fatal("not enough vertices to make a polyhedron")

      ! update the parent to child vertex id table with vertices that lie on the plane
      do v = 1,this%nVerts
        if (side(v)==0) then
          ! search for the interface polygon vertex which coincides with this vertex
          do iv = 1,intpoly%nVerts
            !if (near_zero(magnitude(this%x(:,v) - intpoly%x(:,iv)))) then
            if (all(near_zero(this%x(:,v) - intpoly%x(:,iv)))) then
              p2c_vid(v) = nParVerts + iv
              exit
            end if
          end do
          if (p2c_vid(v)==0) call tls_fatal ("parent point not found in intersecting polygon")
        end if
      end do

    end subroutine generate_new_verts

    subroutine find_edges (edge_vid,nEdges, this,side,valid_side,p2c_vid,nParVerts,nVerts,v_assoc_pe)

      integer,          intent(out) :: edge_vid(:,:), nEdges
      type(pure_polyhedron), intent(in)  :: this
      integer,          intent(in)  :: side(:), valid_side, p2c_vid(:), nParVerts, nVerts, &
          v_assoc_pe(:)

      integer :: e, v, n_on_valid_side, new_edge(2)

      nEdges = 0; edge_vid = 0

      ! add edges that were a part of the parent
      do e = 1,this%nEdges
        ! how many of the edge vertices are on the valid side of the plane?
        n_on_valid_side = count(side(this%edge_vid(:,e))==valid_side)

        if (n_on_valid_side>0) then
          if (n_on_valid_side==2) then ! this entire edge is on the valid side of the plane
            new_edge = p2c_vid(this%edge_vid(:,e))
          else if (n_on_valid_side==1) then ! this edge was divided (or one of the vertices lies on the plane)
            ! the edge consists of the vertex on this side of the plane and the intersection point
            if (side(this%edge_vid(1,e))==valid_side) then
              new_edge = [p2c_vid(this%edge_vid(1,e)),nParVerts+v_assoc_pe(e)]
            else
              new_edge = [p2c_vid(this%edge_vid(2,e)),nParVerts+v_assoc_pe(e)]
            end if
          end if

          ! if this edge isn't already listed, add it
          ! the edge might already be listed in cases where we have very close vertices (O(alpha)),
          ! which are then combined in the new polyhedron.
          if (.not.containsPair(new_edge, edge_vid(:,1:nEdges))) then
            nEdges = nEdges + 1
            edge_vid(:,nEdges) = new_edge
          end if
        end if
      end do

      ! points on the plane make up edges with each other
      do v = nParVerts+1,nVerts-1
        nEdges = nEdges + 1
        edge_vid(:,nEdges) = [v,v+1]
      end do
      nEdges = nEdges + 1
      edge_vid(:,nEdges) = [nVerts,nParVerts+1] ! complete the loop

    end subroutine find_edges

    ! note 1: There may not be a valid edge between the previously found vertex and this one.
    !         This particularly happens when there are multiple vertices on this face which
    !         also lie on the plane. In that case, we loop through those points, adding them
    !         until we find an edge between one and the next vertex.
    subroutine find_faces (face_vid,face_normal,nFaces,ierr, &
        this,side,valid_side,p2c_vid,nParVerts,nVerts,nPolyVerts,v_assoc_pe, plane_normal)

      integer,          intent(out) :: face_vid(:,:),nFaces,ierr
      real(r8), intent(out) :: face_normal(:,:)
      type(pure_polyhedron), intent(in)  :: this
      integer,          intent(in)  :: side(:), valid_side, p2c_vid(:), nParVerts, nVerts, &
          nPolyVerts, v_assoc_pe(:)
      real(r8), intent(in) :: plane_normal(:)

      integer :: v,nV,f, e, edge_cont_verts(this%nVerts,this%nVerts)

      ! first make a lookup table for finding the edge associated with a pair of vertices
      !print *, 'findfaces0'
      edge_cont_verts = 0
      do e = 1,this%nEdges
        edge_cont_verts(this%edge_vid(1,e),this%edge_vid(2,e)) = e
        edge_cont_verts(this%edge_vid(2,e),this%edge_vid(1,e)) = e
      end do
      !print *, 'findfaces1'
      nFaces = 0; face_vid = 0; ierr = 0

      ! loop over all of the parent's faces, adding, modifying, or ignoring it's faces as needed
      do f = 1,this%nFaces
        !print *, 'findfaces2',f
        ! if any vertices for this original face are on the valid side of the plane,
        ! then the face structure can be preserved or modified.
        ! Otherwise, it is thrown out since the entire face is behind the plane.
        nV = count(this%face_vid(:,f) /= 0) ! number of vertices on this face
        if (any(side(this%face_vid(1:nV,f))==valid_side)) then
          nFaces = nFaces + 1
          face_normal(:,nFaces) = this%face_normal(:,f)
          if (all(side(this%face_vid(1:nV,f))/=-valid_side)) then
            ! if all vertices are on the valid side of the face, this face is preserved exactly
            face_vid(1:nV,nFaces) = p2c_vid(this%face_vid(1:nV,f))
          else
            ! if some vertices are on the valid side of the face, this face is modified
            ! write(*,*)
            ! write(*,*) 'f: ',nFaces
            ! write(*,*) this%face_vid(1:nV,f)
            ! write(*,*) side(this%face_vid(1:nV,f))
            ! write(*,'(20i3)') v_assoc_pe
            call modified_parent_face (face_vid(:,nFaces),ierr, this%face_vid(1:nV,f), p2c_vid, &
                edge_vid(:,1:nEdges), side, valid_side, nParVerts, nPolyVerts, v_assoc_pe, &
                edge_cont_verts)

            ! check that the face is valid by making sure each pair of points corresponds to an edge
            nV = count(face_vid(:,nFaces)/=0)
            do v = 1,nV
              if (.not.containsPair(face_vid([v,modulo(v,nV)+1],nFaces),edge_vid(:,1:nEdges))) then
                write(*,*) "find_faces: face contains edge that doesn't exist"
                write(*,*) "This often happens when the parent polyhedron has a non-planar face."
                write(*,*) 'vertices: ',face_vid([v,modulo(v,nV)+1],nFaces)
                ierr = 1
                exit
              end if
            end do

            if (ierr /= 0) then
              nV = count(this%face_vid(:,f) /= 0)
              write(*,*) 'invalid face'
              write(*,'(a,12i3)') 'face: ',face_vid(:,nFaces)
              write(*,'(a,  i3)') 'valid_side: ',valid_side
              write(*,'(a,15i3)') 'parent vertex sides ',side
              write(*,'(a,15i3)') 'parent face vertices ',this%face_vid(1:nV,f)
              write(*,'(a,15i3)') 'parent face vertex sides ',side(this%face_vid(1:nV,f))
              write(*,'(a,15i3)') 'p2c_vid ',p2c_vid(this%face_vid(1:nV,f))
              write(*,'(a,15i3)') 'p2c_vid ',p2c_vid
              do e = 1,nEdges
                write(*,'(a,i4,a,2i4)') 'edge ',e,': ',edge_vid(:,e)
              end do
              do e = 1,nVerts
                write(*,'(a,i3,a,3es35.25)') 'x ',e,':  ',x(:,e)
              end do
              write(*,*)
              return
            end if

            ! ! if this face is invalid, delete it.
            ! ! this line is here to counter a common, but easily fixable problem:
            ! ! invalid (and unnecessary) faces generated when two edges are very
            ! ! close together (<2alpha).
            ! ! WARNING: this is also a prime target for subtle errors, and should
            ! !          really be figured out and handled in a more elegant manner
            ! if (count(face_vid(:,nFaces)>0) < 3) nFaces = nFaces - 1
          end if
        end if
      end do

      ! add the new face from points lying on the plane
      nFaces = nFaces + 1
      if (valid_side>0) then ! needs to be opposite the input order
        face_vid(1:nVerts-nParVerts,nFaces) = [(f, f=nVerts,nParVerts+1, sign(1,nParVerts+1-nVerts))]
      else
        face_vid(1:nVerts-nParVerts,nFaces) = [(f, f=nParVerts+1,nVerts, sign(1,nVerts-nParVerts-1))]
      end if

      face_normal(:,nFaces) = plane_normal

    end subroutine find_faces

    ! calculate the face by modifying the parent face structure with the new edge from the plane
    !
    ! TODO: Can the "next" point on the intersecting plane be found by going backwards on the
    !       intersecting polygon? i.e., if v on the polygon is on this face, is v-1 (with
    !       appropriate modulos) the necessary next point on the face? This might prevent
    !       the more expensive searching with edge_cont_verts.
    !
    ! note 1: in cases where several vertices on the parent face lie on the plane,
    !         this vertex may not be connected to the previous one found, but instead
    !         one of the other previous ones. If this edge does not exist and the
    !         last point was also on the plane, delete the last found point
    subroutine modified_parent_face (face_vid, ierr, par_face_vid, p2c_vid, edge_vid, side, &
        valid_side, nParVerts, nPolyVerts, v_assoc_pe, edge_cont_verts)

      integer, intent(out) :: face_vid(:), ierr
      integer, intent(in) :: par_face_vid(:), p2c_vid(:), edge_vid(:,:), side(:), valid_side, &
          nParVerts, nPolyVerts, v_assoc_pe(:), edge_cont_verts(:,:)

      integer :: v, cvid, ecv, pv, nV

      face_vid = 0; ierr = 0
      nV = size(par_face_vid)
      cvid = 1

      !print *, 'modface0'
      ! start with the first valid vertex after the section intersected by the plane
      v = findloc(side(par_face_vid), valid_side, dim=1)
      if (v==1) v = modulo(findloc(side(par_face_vid)/=valid_side, .true., dim=1, back=.true.),nV)+1
      !print *
      !write(*,*) 'v', v, nV, modulo(v-2,nV)+1, modulo(v-1,nV)+1

      ! the edge [v-1,v] *must* intersect with the plane
      ecv = edge_cont_verts(par_face_vid(modulo(v-2,nV)+1),par_face_vid(modulo(v-1,nV)+1))
      if (ecv < 1) then
        write(*,*) 'find faces failed: edge missing'
        write(*,*) 'vertices: ',par_face_vid(modulo(v-2,nV)+1),par_face_vid(modulo(v-1,nV)+1)
        ierr = 1
        return
      end if
      face_vid(cvid) = nParVerts + v_assoc_pe(ecv)
      cvid = cvid+1
      !write(*,*) '1',face_vid(cvid-1), ecv, nParVerts, v_assoc_pe(ecv)

      ! add all the parent's valid vertices in sequence, stop when we hit one that is invalid
      ! this stopping point indicates the next edge [v-1,v] which *must* be intersected by the plane
      do while (side(par_face_vid(modulo(v-1,nV)+1))==valid_side)
        face_vid(cvid) = p2c_vid(par_face_vid(modulo(v-1,nV)+1))
        v = v+1; cvid = cvid+1
        !write(*,*) '2',face_vid(cvid-1)
      end do

      ! the second new vertex is on this edge [v-1,v]
      ecv = edge_cont_verts(par_face_vid(modulo(v-2,nV)+1),par_face_vid(modulo(v-1,nV)+1))
      if (ecv < 1) then
        write(*,*) 'find faces failed: edge missing'
        write(*,*) 'vertices: ',par_face_vid(modulo(v-2,nV)+1),par_face_vid(modulo(v-1,nV)+1)
        ierr = 1
        return
      end if
      face_vid(cvid) = nParVerts + v_assoc_pe(ecv)
      cvid = cvid+1
      !write(*,*) '3',face_vid(cvid-1), ecv, nParVerts, v_assoc_pe(ecv)

      ! add additional vertices on this face which lie on the plane
      ! stop when a complete face is formed
      !write(*,*) 'edge: ',[face_vid(1),face_vid(cvid-1)]
      v = v+1; pv = 1
      do while (.not.containsPair([face_vid(1),face_vid(cvid-1)],edge_vid) .and. pv<=nV)
        !print *, 'try', pv
        if (side(par_face_vid(modulo(v-1,nV)+1))==0) then
          face_vid(cvid) = p2c_vid(par_face_vid(modulo(v-1,nV)+1))
          cvid = cvid+1
          !write(*,*) '4',face_vid(cvid-1)
        end if
        v = v+1; pv = pv+1
      end do
      if (pv>nV .and. .not.containsPair([face_vid(1),face_vid(cvid-1)],edge_vid)) then
        write(*,*) 'final pair not found: ',[face_vid(1),face_vid(cvid-1)]
        ierr = 1
      end if

    end subroutine modified_parent_face

  end function polyhedron_on_side_of_plane

  ! subroutine tesselate_nonplanar_faces(this)

  !   class(pure_polyhedron), intent(inout) :: this

  !   integer :: f

  !   do f = 1,this%nFaces
  !     if (this%face_is_nonplanar(f)) call this%tesselate_face(f)
  !   end do

  ! end subroutine tesselate_nonplanar_faces

  ! check if the given point is inside the polyhedron
  logical function is_inside(this, x)

    class(pure_polyhedron), intent(in) :: this
    real(r8), intent(in) :: x(:)

    integer :: f

    ! check that the point is behind each face
    is_inside = .true.
    do f = 1,this%nFaces
      is_inside = is_inside .and. &
          dot_product(this%face_normal(:,f), x - this%x(:,this%face_vid(1,f))) <= 0
      if (.not.is_inside) exit
    end do

  end function is_inside

  logical function has_nonplanar_face(this)

    class(pure_polyhedron), intent(in) :: this

    integer :: f

    has_nonplanar_face = .false.
    do f = 1, this%nFaces
      has_nonplanar_face = this%face_is_nonplanar(f)
      if (has_nonplanar_face) exit
    end do

  end function has_nonplanar_face

  logical function face_is_nonplanar(this, f)

    use plane_type

    class(pure_polyhedron), intent(in) :: this
    integer, intent(in) :: f

    integer :: v, nV
    type(plane) :: p

    face_is_nonplanar = .false.

    nV = count(this%face_vid(:,f) /= 0)
    if (nV==3) return

    ! generate plane from first three points
    p%normal = this%face_normal(:,f)
    p%rho = - dot_product(p%normal, this%x(:,this%face_vid(1,f)))

    ! check if remaining points are more than a distance alpha from the plane
    do v = 4,nV
      if (abs(p%signed_distance(this%x(:,this%face_vid(v,f)))) > alpha) then
        face_is_nonplanar = .true.
        exit
      end if
    end do

  end function face_is_nonplanar

  ! subroutine tesselate_face(this, f)

  !   use cell_geometry, only: cross_product, normalized

  !   class(pure_polyhedron), intent(inout) :: this
  !   integer, intent(in) :: f

  !   integer :: nV, n_new, ff, i, j
  !   integer, allocatable :: tmpi(:,:)
  !   real(r8), allocatable :: tmpr(:,:)

  !   nV = count(this%face_vid(:,f) /= 0)
  !   if (nV == 3) return

  !   ! new vertex is at face center
  !   this%nVerts = this%nVerts + 1
  !   allocate(tmpr(3,this%nVerts))
  !   tmpr(:,:this%nVerts-1) = this%x
  !   tmpr(:,this%nVerts) = sum(this%x(:,this%face_vid(:nV,f)), dim=2) / nV
  !   call move_alloc(tmpr, this%x)

  !   ! all vertices on this face have a new edge to the new vertex
  !   n_new = this%nEdges + nV
  !   allocate(tmpi(2,n_new))
  !   tmpi(:,:this%nEdges) = this%edge_vid
  !   tmpi(1,this%nEdges+1:) = this%face_vid(:nV,f)
  !   tmpi(2,this%nEdges+1:) = this%nVerts
  !   call move_alloc(tmpi, this%edge_vid)
  !   this%nEdges = n_new

  !   ! add new faces
  !   n_new = this%nFaces + nV - 1
  !   allocate(tmpi(size(this%face_vid, dim=1), n_new))
  !   tmpi(:,:this%nFaces) = this%face_vid
  !   tmpi(:3,f) = [this%face_vid(1,f), this%face_vid(2,f), this%nVerts] ! first new face replaces old
  !   tmpi(4:,f) = 0
  !   do ff = this%nFaces + 1, n_new
  !     i = ff - this%nFaces + 1
  !     j = modulo(i, nV) + 1
  !     tmpi(:3,ff) = [this%face_vid(i,f), this%face_vid(j,f), this%nVerts]
  !     tmpi(4:,ff) = 0
  !   end do
  !   call move_alloc(tmpi, this%face_vid)

  !   ! calculate new face normal vectors
  !   allocate(tmpr(3,n_new))
  !   tmpr(:,:this%nFaces) = this%face_normal
  !   tmpr(:,f) = normalized(cross_product (&
  !         this%x(:,this%face_vid(2,f)) - this%x(:,this%face_vid(1,f)), &
  !         this%x(:,this%face_vid(3,f)) - this%x(:,this%face_vid(1,f))))
  !   do ff = this%nFaces+1, n_new
  !     tmpr(:,ff) = normalized(cross_product (&
  !         this%x(:,this%face_vid(2,ff)) - this%x(:,this%face_vid(1,ff)), &
  !         this%x(:,this%face_vid(3,ff)) - this%x(:,this%face_vid(1,ff))))
  !   end do
  !   call move_alloc(tmpr, this%face_normal)
  !   this%nFaces = n_new

  ! end subroutine tesselate_face

  subroutine print_data (this,normalized)

    class(pure_polyhedron), intent(in) :: this
    logical, optional, intent(in) :: normalized

    integer :: v,e,f
    real(r8) :: x0(3),xl(3)
    logical :: normalizedh

    normalizedh = merge(normalized, .false., present(normalized))

    write(*,*) 'POLYHEDRON DATA:'

    if (allocated(this%x)) then
      if (normalizedh) then
        x0 = minval(this%x,dim=2)
        xl = maxval(this%x,dim=2) - x0
        do v = 1,this%nVerts
          write(*,'(a,i3,a,3es25.15)') 'x ',v,':  ',(this%x(:,v)-x0)/xl
        end do
      else
        do v = 1,this%nVerts
          write(*,'(a,i3,a,3es35.25)') 'x ',v,':  ',this%x(:,v)
        end do
      end if
      write(*,*)
    end if

    if (allocated(this%edge_vid)) then
      do e = 1,this%nEdges
        write(*,'(a,i3,a,2i4)') 'edge ',e,':  ',this%edge_vid(:,e)
      end do
      write(*,*)
    end if

    if (allocated(this%face_vid)) then
      do f = 1,this%nFaces
        write(*,'(a,i3,a,10i4)') 'face ',f,':  ',this%face_vid(:,f)
      end do
      write(*,*)
    end if

    if (allocated(this%face_normal)) then
      do f = 1,this%nFaces
        write(*,'(a,i3,a,3es35.25)') 'norm ',f,':  ',this%face_normal(:,f)
      end do
      write(*,*)
    end if

    write(*,'(a,es20.10)') 'volume ',this%vol

  end subroutine print_data

  pure logical function containsPair (pair,list)

    integer, intent(in) :: pair(:), list(:,:)

    integer :: i

    containsPair = .false.
    do i = 1,size(list, dim=2)
      containsPair = containsPair .or. &
          all(list(:,i)==pair) .or. all(list(:,i)==pair(size(pair):1:-1))
    end do

  end function containsPair

  pure logical function containsPoint (x,list)

    use near_zero_function

    real(r8), intent(in) :: x(:), list(:,:)

    integer :: i

    containsPoint = .false.
    do i = 1,size(list, dim=2)
      containsPoint = containsPoint .or. all(near_zero(x - list(:,i)))
      if (containsPoint) return
    end do

  end function containsPoint

  pure integer function pointIndex(x,list)

    use near_zero_function

    real(r8), intent(in) :: x(:), list(:,:)

    do pointIndex = 1,size(list,dim=2)
      if (all(near_zero(x - list(:,pointIndex)))) return
    end do
    pointIndex = -1

  end function pointIndex

  function invert (x)
    integer, intent(in) :: x(:)
    integer             :: invert(size(x))

    integer             :: i

    do i = 1,size(x)
      invert(x(i)) = i
    end do

  end function invert

end module pure_polyhedron_type
