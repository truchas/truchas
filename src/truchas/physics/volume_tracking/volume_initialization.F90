!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module volume_initialization

  use kinds, only: r8
  use unstr_mesh_type
  use cell_topology
  use cell_geometry
  use truchas_logging_services
  implicit none
  private

  public :: compute_initial_volumes

  integer, parameter :: tet4_ntets = 8
  integer, parameter :: pyr5_ntets = 4
  integer, parameter :: wed6_ntets = 17
  integer, parameter :: hex8_ntets = 24
  integer, parameter :: nrecurse = 5
  integer, parameter :: nverts_max = 14

  integer, parameter :: pyr5_to_tet(4,pyr5_ntets) = &
      reshape(source=[1,2,5,6, 2,3,5,6, 3,4,5,6, 1,4,5,6], shape=[4,pyr5_ntets])

  integer, parameter :: wed6_to_tet(4,wed6_ntets) = &
      reshape(source=[1,4,7,9, 2,5,7,8, 3,6,8,9, 7,8,9,10, 7,8,9,11, 1,7,9,10, 2,7,8,10, &
      3,8,9,10, 1,2,7,10, 1,3,9,10, 2,3,8,10, 5,7,8,11, 6,8,9,11, 4,7,9,11, &
      5,6,8,11, 4,5,7,11, 4,6,9,11], shape=[4,wed6_ntets])

  integer, parameter :: hex8_to_tet(4,hex8_ntets) = &
      reshape(source=[1,5,9,12, 2,6,9,10, 3,7,10,11, 4,8,11,12, 1,2,9,13, 2,3,10,13, 3,4,11,13, &
      1,4,12,13, 5,6,9,14, 6,7,10,14, 7,8,11,14, 5,8,12,14, 9,10,13,14, &
      10,11,13,14, 11,12,13,14, 9,12,13,14, 5,9,12,14, 6,9,10,14, 7,10,11,14, &
      8,11,12,14, 1,9,12,13, 2,9,10,13, 3,10,11,13, 4,11,12,13], shape=[4,hex8_ntets])

  integer, parameter :: tet4_to_tet(4,tet4_ntets) = &
      reshape(source=[1,5,6,7, 2,5,8,9, 3,6,8,10, 4,7,9,10, 5,6,7,10, 5,6,8,10, 5,7,9,10, 5,8,9,10],&
      shape=[4,tet4_ntets])

contains

  subroutine compute_initial_volumes(m, vols)
    type(unstr_mesh), intent(in) :: m
    real(r8), intent(out) :: vols(:,:)

    integer :: i, j, id(nverts_max), nnodes, nfaces
    real(r8) :: v_tot, v(size(vols,dim=1)), verts(3,nverts_max)

    do i = 1, m%ncell_onp
      v_tot = m%volume(i)
      if (size(vols,dim=1) == 1) then
        vols(1,i) = v_tot
        cycle
      end if

      v(:) = 0.0_r8

      associate(cn => m%cnode(m%xcnode(i):m%xcnode(i+1)-1))
        select case (size(cn))
        case (4) ! tet
          nnodes = size(cn)
          verts(:,1:nnodes) = m%x(:,cn)
          call vertex_ids(verts, id, nnodes)
          call split_tet_recursive(verts, id, v, 0)

        case (5) ! pyramid
          nnodes = size(cn)
          nfaces = 5
          verts(:,1:nnodes) = m%x(:,cn)
          call vertex_ids(verts, id, nnodes)

          if (all(id(:nnodes).eq.id(1))) then
            v(id(1)) = v_tot
          else
            ! split the only potentially non-planar face (face 5)
            verts(:,nnodes+1) = face_centroid(m, cn, nfaces)
            call vertex_ids(verts(:,nnodes+1:), id(nnodes+1:), 1)
            do j = 1,pyr5_ntets
              call split_tet_recursive(verts(:,pyr5_to_tet(:,j)), id(pyr5_to_tet(:,j)), v, 0)
            end do
          end if
        case (6) ! wedge
          nnodes = size(cn)
          nfaces = 5
          verts(:,1:nnodes) = m%x(:,cn)
          call vertex_ids(verts, id, nnodes)

          if (all(id(:nnodes).eq.id(1))) then
            v(id(1)) = v_tot
          else
            do j = 1,nfaces
              verts(:,nnodes+j) = face_centroid(m, cn, j)
            end do
            call vertex_ids(verts(:,nnodes+1:), id(nnodes+1:), nfaces)
            do j = 1,wed6_ntets
              call split_tet_recursive(verts(:,wed6_to_tet(:,j)), id(wed6_to_tet(:,j)), v, 0)
            end do
          end if
        case (8) ! hex
          nnodes = size(cn)
          nfaces = 6
          verts(:,1:nnodes) = m%x(:,cn)
          call vertex_ids(verts, id, nnodes)

          if (all(id(:nnodes).eq.id(1))) then
            v(id(1)) = v_tot
          else
            do j = 1,nfaces
              verts(:,nnodes+j) = face_centroid(m, cn, j)
            end do
            call vertex_ids(verts(:,nnodes+1:), id(nnodes+1:), nfaces)
            do j = 1,hex8_ntets
              call split_tet_recursive(verts(:,hex8_to_tet(:,j)), id(hex8_to_tet(:,j)), v, 0)
            end do
          end if
        case default
          call tls_fatal("pc load letter")
        end select
      end associate

      ! correct for mistakes in planar approximations
      vols(:,i) = v * v_tot/sum(v)

    end do

  contains

    subroutine vertex_ids(x, id, n)
      real(r8), intent(in) :: x(:, :)
      integer, intent(out) :: id(:)
      integer, intent(in) :: n

      integer :: i

      do i=1,n
        id(i) = body_id_from_vertex(x(:,i))
      end do

    end subroutine vertex_ids

    function face_centroid(m, cn, fi) result(v)
      type(unstr_mesh), intent(in) :: m
      integer, intent(in) :: cn(:), fi
      real(r8) :: v(3)

      select case (size(cn))
      case (5)
        v = sum(m%x(:,cn(pyr5_faces(pyr5_xface(fi):pyr5_xface(fi+1)-1))),dim=2)&
            /real(pyr5_fsize(fi), r8)
      case (6)
        v = sum(m%x(:,cn(wed6_faces(wed6_xface(fi):wed6_xface(fi+1)-1))),dim=2)&
            /real(wed6_fsize(fi), r8)
      case (8)
        v = sum(m%x(:,cn(hex8_faces(hex8_xface(fi):hex8_xface(fi+1)-1))),dim=2)&
            /real(hex8_fsize(fi), r8)
      end select
    end function face_centroid

    recursive subroutine split_tet_recursive(verts, id, vol, lvl)
      real(r8), intent(in) :: verts(:,:)
      integer, intent(in) :: id(:), lvl
      real(r8), intent(inout) :: vol(:)
      !-
      real(r8) :: nv(3,10)
      integer :: nid(10), i

      if (all(id(:4).eq.id(1))) then
        vol(id(1)) = vol(id(1)) + abs(tet_volume(verts(:,:4)))
        return
      end if

      if (lvl.ge.nrecurse) then
        ! get id of centroid and ascibe tet volume to it
        nv(:,1) = 0.25_r8*sum(verts(:,:4),dim=2)
        nid(1) = body_id_from_vertex(nv(:,1))
        vol(nid(1)) = vol(nid(1)) + abs(tet_volume(verts(:,:4)))
        return
      end if

      ! build node + edge midpoint vertex collection
      nv(:,:4) = verts(:,:4)
      nid(:4) = id(:4)

      do i = 1,6
        nv(:,4+i) = 0.5_r8*sum(verts(:,tet4_edges(:,i)),dim=2)
      end do
      call vertex_ids(nv(:,5:), nid(5:), 6)

      do i = 1,tet4_ntets
        call split_tet_recursive(nv(:,tet4_to_tet(:,i)), &
            nid(tet4_to_tet(:,i)), vol, lvl+1)
      end do

    end subroutine split_tet_recursive

  end subroutine compute_initial_volumes

  !==
  ! copied from vof_init.f90
  !==
  function body_id_from_vertex (x)
    !-------------------------------------------------------------------------
    ! given a vertex (represented as (x,y,z), return the body that
    ! the vertex is in
    !
    ! this code is ugly, undocumented, and based on what was already there
    !-------------------------------------------------------------------------

    use interfaces_module, only: ar, cosa, isrftype, nbody, nsurf, offset, &
        rotangl, rotpt, sgeom, sina
    use parameter_module,  only: nrot
    use tally_module,      only: pre_tally
    implicit none
    ! arguments
    real(r8), dimension(3), intent(in) :: x

    ! return value
    integer :: body_id_from_vertex

    ! local variables
    logical :: lsf                             ! logical surface flag
    logical :: lbf                             ! logical body flag
    integer :: ib
    integer :: is
    integer :: isrf
    integer :: n
    integer :: n1
    integer :: n2
    integer :: na
    real(r8) :: total
    real(r8) :: coeff
    real(r8), dimension(3) :: xloc

    logical, save :: first_time = .true.

    !-------------------------------------------------------------------------

    ! get the geometric information of the bodies once
    if (first_time) then
      call pre_tally()
      first_time = .false.
    endif

    ! assume that everything is in the background body unless we
    ! find some other body
    body_id_from_vertex = nbody

    do ib = 1, nbody
      ! loop over the bodies in the domain and determine
      ! whether the point falls inside any of the bodies.
      lbf = .true.

      do is = 1, nsurf(ib)
        ! loop over each surface that defines the body and
        ! determine whether the point lies on the correct side
        ! of the surface.

        ! initialize the coordinates local to the surface
        xloc = x

        do n = 1, nrot

          ! select the rotation axis
          na = n
          if (3 == 2) then
            na = 3
          end if

          ! select the rotation plane axes (two of them).
          select case (na)
          case (1)
            n1=2; n2=3; coeff =  1.0d0
          case (2)
            n1=3; n2=1; coeff = -1.0d0
          case (3)
            n1=1; n2=2; coeff =  1.0d0
          end select

          ! rotate the surface coordinate system.
          if (rotangl(n,is,ib) /= 0.0_r8) then
            total    = cosa(n,is,ib)*(xloc(n1)-rotpt(n1,is,ib)) &
                - coeff*sina(n,is,ib)*(xloc(n2)-rotpt(n2,is,ib)) + rotpt(n1,is,ib)
            xloc(n2) = cosa(n,is,ib)*(xloc(n2)-rotpt(n2,is,ib)) &
                + coeff*sina(n,is,ib)*(xloc(n1)-rotpt(n1,is,ib)) + rotpt(n2,is,ib)
            xloc(n1) = total
          end if

        end do

        ! translate the surface coordinate system
        do n = 1, 3
          xloc(n) = xloc(n) - offset(n,is,ib)
        end do

        ! jump to specified surface type and compute whether
        ! the point lies on the body side of the surface.
        isrf = abs(isrftype(is,ib))

        select case (isrf)

        case (1)
          ! surface is a plane which, after translation and rotation
          ! contains the local origin and is normal to one of the
          ! coordinate axes.  parameter sgeom(1,is,ib) specifies the
          ! normal direction to the plane, with absolute values of
          ! 1., 2., or 3. implying the local normal to the plane is
          ! the x-, y-, or z-axis, respectively.  positive values of
          ! isrftype(is,ib) imply the body lies on the positive side
          ! of the plane, and negative values imply the body is on the
          ! negative side.

          do n = 1, 3
            if (abs(sgeom(1,is,ib)) == real(n)) then
              lsf = xloc(n) > 0.0_r8
            end if
          end do
          if (isrftype(is,ib) < 0) then
            lsf = .not. lsf
          end if
          lbf = lbf .and. lsf

        case (2)
          ! surface is a right parallelopiped whose local origin
          ! is its centroid.  its x-y-z half-lengths are half of
          ! sgeom(1,is,ib), sgeom(2,is,ib), and sgeom(3,is,ib),
          ! respectively.

          lsf = .true.
          do n = 1, 3
            lsf = lsf .and. abs(xloc(n)) < ar(n,is,ib)
          end do
          if (isrftype(is,ib) < 0) then
            lsf = .not. lsf
          end if
          lbf = lbf .and. lsf

        case (3:4)
          ! surface is a sphere or ellipsoid

          total = 0.0_r8
          do n = 1, 3
            total = total + (ar(n,is,ib)*xloc(n))**2
          end do
          lsf = total <= 1.0_r8
          if (isrftype(is,ib) < 0) then
            lsf = .not. lsf
          end if
          lbf = lbf .and. lsf

        case (5)
          ! surface is a right circular cylinder.
          ! parameter sgeom(1,is,ib) specifies the axis along
          ! which the cylinder is aligned.  values of 1., 2., or 3.
          ! specify the x-, y-, or z-axis, respectively.  parameter
          ! sgeom(2,is,ib) is the cylinder radius and sgeom(3,is,ib)
          ! is the cylinder height.

          total = 0.0_r8
          do n = 1, 3
            total = total + (ar(n,is,ib)*xloc(n))**2
          end do
          lsf = total <= 1.0_r8
          do n = 1, 3
            if (abs(sgeom(1,is,ib)) == real(n)) then
              lsf = lsf .and. (xloc(n) - sgeom(3,is,ib))*xloc(n) <= 0.0_r8
            end if
          end do
          if (isrftype(is,ib) < 0) then
            lsf = .not. lsf
          end if
          lbf = lbf .and. lsf

        case (6)
          ! surface is a right circular cone.  parameter sgeom(1,is,
          ! ib) specifies the axis along which the cone is aligned.
          ! values of 1., 2., or 3. specify the x-, y-, or z-axis,
          ! respectively.  parameter sgeom(2,is,ib) is the cone height
          ! relative to the cone origin at the base.  a positive
          ! height implies the cone apex points in the positive axis
          ! direction, and while a negative height implies the apex
          ! points in the negative axis direction.  parameter
          ! sgeom(3,is,ib) is the radius of the cone at its base.

          do n = 1,3

            select case (n)
            case (1)
              n1=2; n2=3
              if (3 == 2) then
                n2 = n1
              end if
            case (2)
              n1=3; n2=1
              if (3 ==2) then
                n1 = n2
              end if
            case (3)
              n1=1; n2=2
            end select

            if (abs(sgeom(1,is,ib)) == real(n)) then
              total = sgeom(3,is,ib)/sgeom(2,is,ib)*(sgeom(2,is,ib) - xloc(n))
              lsf = xloc(n1)**2 + xloc(n2)**2 <= total**2 &
                  .and. sgeom(2,is,ib)*(sgeom(2,is,ib) - xloc(n)) >= 0.0_r8 &
                  .and. sgeom(2,is,ib)*xloc(n) >= 0.0_r8
            end if

          end do

          if (isrftype(is,ib) < 0) then
            lsf = .not. lsf
          end if
          lbf = lbf .and. lsf

        case (7)

          ! body is a tabular n-sided polygon, to be rotated if
          ! parameter sgeom(1,is,ib) = 0, or translated if = 1,
          ! about the axis (1=x, 2=y, 3=z) specified by the
          ! parameter sgeom(2,is,ib). the integer part
          ! of parameter sgeom(3,is,ib) is npoly, the # of vertices
          ! in the polygon, given by rtab(0:npoly,ib) and
          ! ztab(0:npoly,ib).

          call TLS_fatal ('body_id_from_vertex: vof_method divide with tabular bodies not yet implemented')

          !                npoly = nint(sgeom(3,is,ib))
          !                rtab(0,ib) = rtab(npoly,ib)
          !                ztab(0,ib) = ztab(npoly,ib)
          !
          !                ! now rotate or translate
          !                do n = 1, 3
          !                   select case (n)
          !                   case (1)
          !                      n1=2; n2=3
          !                      if (3 == 2) then
          !                         n2 = n1
          !                      end if
          !
          !                   case (2)
          !                      n1=3; n2=1
          !                      if (3 == 2) then
          !                         n1 = n2
          !                      end if
          !
          !                   case (3)
          !                      n1=1; n2=2
          !
          !                   end select
          !
          !                   if (sgeom(1,is,ib) == 0.0_r8) then
          !                      ! rotate
          !                      xloc(n1) = sqrt(xloc(n1)**2 + xloc(n2)**2)
          !                      na = n
          !                   else
          !                      ! translate
          !                      na = n2
          !                   end if
          !
          !                   ! select axis of rotation or translation
          !                   if (abs(sgeom(2,is,ib)) == real(n)) then
          !                      call in_out_crse (xloc(n1), xloc(na), total, npoly, ib)
          !                   end if
          !                end do
          !
          !                lsf = total == 1
          !                if (isrftype(is,ib) < 0) then
          !                   lsf = .not. lsf
          !                end if
          !                lbf = lbf .and. lsf

        case (8)
          ! get material distribution from mesh file; whereever the mesh_material
          ! array matches the material number of this body, we have a hit.

          call TLS_fatal ('body_id_from_vertex: vof_method divide with mesh materials not yet implemented')

          !                lsf = .false.
          !                if (associated(mesh_material)) then
          !                   lsf = mesh_material == mesh_matnum(ib)
          !                   if (isrftype(is,ib) < 0) then
          !                      lsf = .not. lsf
          !                   end if
          !
          !                   lbf = lbf .and. lsf
          !                end if

        end select

      end do

      if (lbf) then
        body_id_from_vertex = ib
        exit
      endif

    end do

  end function body_id_from_vertex
end module volume_initialization
