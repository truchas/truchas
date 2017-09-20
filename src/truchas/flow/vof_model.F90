
#include "f90_assert.fpp"

module vof_model

  use kinds, only: r8
  use unstr_mesh_type
  use parameter_list_type
  use truchas_logging_services
  use cell_topology
  use cell_geometry
  implicit none
  private

  type, public :: vof_model_t
   contains
     procedure :: initial_state_cell
     procedure :: initial_state_all
  end type vof_model_t

  integer, parameter :: tet4_ntets = 8
  integer, parameter :: pyr5_ntets = 4
  integer, parameter :: wed6_ntets = 17
  integer, parameter :: hex8_ntets = 24
  integer, parameter :: nrecurse = 5
  integer, parameter :: nverts_max = 14

  integer, parameter :: pyr5_to_tet(4,pyr5_ntets) = &
       reshape(source=[1,2,5,6,  2,3,5,6,  3,4,5,6,  1,4,5,6],&
               shape=[4,pyr5_ntets])

  integer, parameter :: wed6_to_tet(4,wed6_ntets) = &
       reshape(source=[1,4,7,9,  2,5,7,8,  3,6,8,9,  7,8,9,10,  7,8,9,11, &
                       1,7,9,10,  2,7,8,10,  3,8,9,10,  1,2,7,10,  1,3,9,10, &
                       2,3,8,10,  5,7,8,11,  6,8,9,11,  4,7,9,11,  5,6,8,11, &
                       4,5,7,11,  4,6,9,11],&
               shape=[4,wed6_ntets])

  integer, parameter :: hex8_to_tet(4,hex8_ntets) = &
       reshape(source=[1,5,9,12,  2,6,9,10,  3,7,10,11,  4,8,11,12,  1,2,9,13,&
                       2,3,10,13,  3,4,11,13,  1,4,12,13,  5,6,9,14, 6,7,10,14,&
                       7,8,11,14, 5,8,12,14,  9,10,13,14,  10,11,13,14,&
                       11,12,13,14,  9,12,13,14, 5,9,12,14,  6,9,10,14,&
                       7,10,11,14,  8,11,12,14,  1,9,12,13, 2,9,10,13,&
                       3,10,11,13,  4,11,12,13], &
               shape=[4,hex8_ntets])

  integer, parameter :: tet4_to_tet(4,tet4_ntets) = &
       reshape(source=[1,5,6,7,  2,5,8,9,  3,6,8,10,  4,7,9,10,  5,6,7,10,&
                       5,6,8,10,  5,7,9,10,  5,8,9,10], shape=[4,tet4_ntets])


contains

  ! return array of (nbody) floats indicitating vofs for this cell
  function initial_state_cell(this, m, nbody, ci) result (v)
    implicit none
    class(vof_model_t), intent(inout) :: this
    type(unstr_mesh), intent(in) :: m
    integer, intent(in) :: nbody, ci
    real(r8) :: v(nbody)
    !-
    real(r8) :: verts(3, nverts_max), v_tot, v_local, dv
    integer :: id(nverts_max)
    integer :: i, fi, nnodes, nfaces

    v_tot = m%volume(ci)
    v = 0.0_r8
    v_local = 0.0_r8

    if (nbody.eq.1) then
       v = v_tot
       return
    end if

    associate (cell_nodes => m%cnode(m%xcnode(ci):m%xcnode(ci+1)-1))
      select case (size(cell_nodes))
      case (4) ! tet
         nnodes = size(cell_nodes)
         verts(:,1:nnodes) = m%x(:,cell_nodes)

         call vertex_ids(verts, id, nnodes)
         call split_tet_recursive(verts, id, v, 0)

      case (5) ! pyramid 5 nodes and 5 faces
         nnodes = 5
         nfaces = 5
         verts(:,1:nnodes) = m%x(:,cell_nodes)

         call vertex_ids(verts, id, nnodes)

         if (all(id(:nnodes).eq.id(1))) then
            v(id(1)) = v_tot
         else
            ! split the only potentially non-planar face (face 5)
            verts(:,nnodes+1) = face_centroid(m, cell_nodes, nfaces)
            call vertex_ids(verts(:,nnodes+1:), id(nnodes+1:), 1)

            do i = 1,pyr5_ntets
               call split_tet_recursive(verts(:,pyr5_to_tet(:,i)), &
                    id(pyr5_to_tet(:,i)), v, 0)
            end do
         end if
      case (6) ! wedge 6 nodes and 5 faces
         nnodes = 6
         nfaces = 5
         verts(:,1:nnodes) = m%x(:,cell_nodes)

         call vertex_ids(verts, id, nnodes)

         if (all(id(:nnodes).eq.id(1))) then
            v(id(1)) = v_tot
            !print *, 'hex body id:', id(1)
         else
            do i = 1,nfaces
               verts(:,nnodes+i) = face_centroid(m, cell_nodes, i)
            end do
            call vertex_ids(verts(:,nnodes+1:), id(nnodes+1:), nfaces)

            do i = 1,wed6_ntets
               call split_tet_recursive(verts(:,wed6_to_tet(:,i)), &
                    id(wed6_to_tet(:,i)), v, 0)
            end do
         end if
      case (8) ! hex
         !print *, '<<< HEX <<< '
         nnodes = 8
         nfaces = 6
         verts(:,1:nnodes) = m%x(:,cell_nodes)

         call vertex_ids(verts, id, nnodes)
         if (all(id(:nnodes).eq.id(1))) then
            v(id(1)) = v_tot
            return
         else
            do i = 1,nfaces
               verts(:,nnodes+i) = face_centroid(m, cell_nodes, i)
            end do
            call vertex_ids(verts(:,nnodes+1:), id(nnodes+1:), nfaces)

!!$            write(*,'(a,/,8(3es15.5, /))') "HEX VERTS", verts(:,:nnodes)
!!$            write(*,'(a, 8i3)') "HEX NODES", id(:nnodes)
            do i = 1,hex8_ntets
               ! local volume sum
               v_local = v_local + abs(tet_volume(verts(:,hex8_to_tet(:,i))))
!!$               dv = sum(v)
!!$               write(*,'(a,i3,4es15.5)') '**PRE:  ', i, &
!!$                    abs(tet_volume(verts(:,hex8_to_tet(:,i)))), v
               call split_tet_recursive(verts(:,hex8_to_tet(:,i)), &
                    id(hex8_to_tet(:,i)), v, 0)
!!$               dv = sum(v)-dv
!!$               write(*,'(a,i3,4es15.5)') '**POST: ', i, &
!!$                    dv, v
            end do
         end if
      case default
         call TLS_fatal ("PC LOAD LETTER")
      end select
    end associate

    write(*,'(a,4(es12.5))') 'cell_vol, local, body_vol, err:', &
         v_tot, v_local, sum(v), abs(v_tot-sum(v))
    ! correct for any mistakes in intitial planar approximation
    v = v * v_tot/sum(v)

  contains

    subroutine vertex_ids(x, id, n)
      implicit none
      real(r8), intent(in) :: x(:, :)
      integer, intent(out) :: id(:)
      integer, intent(in) :: n
      integer :: i
      !-
      do i=1,n
         id(i) = body_id_from_vertex(x(:,i))
      end do

    end subroutine vertex_ids

    function face_centroid(m, cn, fi) result(v)
      implicit none
      type(unstr_mesh), intent(in) :: m
      integer, intent(in) :: cn(:), fi
      real(r8) :: v(3)

      select case (size(cn))
      case (5)
         v = sum(m%x(:,cn(PYR5_FACES(PYR5_XFACE(fi):PYR5_XFACE(fi+1)-1))),dim=2)&
              /real(PYR5_FSIZE(fi), r8)
      case (6)
         v = sum(m%x(:,cn(WED6_FACES(WED6_XFACE(fi):WED6_XFACE(fi+1)-1))),dim=2)&
              /real(WED6_FSIZE(fi), r8)

      case (8)
         v = sum(m%x(:,cn(HEX8_FACES(HEX8_XFACE(fi):HEX8_XFACE(fi+1)-1))),dim=2)&
              /real(HEX8_FSIZE(fi), r8)
      end select
    end function face_centroid

    recursive subroutine split_tet_recursive(verts, id, vol, lvl)
      implicit none
      real(r8), intent(in) :: verts(:,:)
      integer, intent(in) :: id(:), lvl
      real(r8), intent(inout) :: vol(:)
      !-
      real(r8) :: nv(3,10)
      integer :: nid(10), i

      if (all(id(:4).eq.id(1))) then
         vol(id(1)) = vol(id(1)) + abs(tet_volume(verts(:,:4)))
         !print *, '>>>tet id, vol, lvl:', id(1), abs(tet_volume(verts)), lvl
         return
      end if

      if (lvl.ge.nrecurse) then
         ! get id of centroid and ascibe tet volume to it
         nv(:,1) = 0.25_r8*sum(verts(:,:4),dim=2)
         nid(1) = body_id_from_vertex(nv(:,1))
         vol(nid(1)) = vol(nid(1)) + abs(tet_volume(verts(:,:4)))
!!$         write(*,'(a,2i4, es15.5)') '***lvl, cent id, vol:', lvl, nid(1), &
!!$              abs(tet_volume(verts(:,:4)))
         return
      end if

      ! build node + edge midpoint vertex collection
      nv(:,:4) = verts(:,:4)
      nid(:4) = id(:4)

      do i = 1,6
         nv(:,4+i) = 0.5_r8*sum(verts(:,TET4_EDGES(:,i)),dim=2)
      end do
      call vertex_ids(nv(:,5:), nid(5:), 6)

      do i = 1,tet4_ntets
         call split_tet_recursive(nv(:,tet4_to_tet(:,i)), &
              nid(tet4_to_tet(:,i)), vol, lvl+1)
      end do

    end subroutine split_tet_recursive

  end function initial_state_cell

  ! return array of (nbody) floats indicitating vofs for this cell
  subroutine initial_state_all(this, mesh, nbody, vof)
    implicit none
    class(vof_model_t), intent(inout) :: this
    type(unstr_mesh), intent(in) :: mesh
    integer, intent(in) :: nbody
    real(r8), intent(out) :: vof(:, :)
    !-
    integer :: i

    print*, 'nbody: ', nbody

    do i=1,mesh%ncell_onP
       print *, 'cell:', i
       vof(:,i) = this%initial_state_cell(mesh, nbody, i)
    end do

  end subroutine initial_state_all

  !==
  ! Copied from vof_init.F90
  !==
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
      implicit none
      ! arguments
      real(r8), dimension(3), intent(in) :: X

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
      real(r8), dimension(3) :: Xloc

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
               if (3 == 2) then
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
            do n = 1, 3
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

               do n = 1, 3
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
               do n = 1, 3
                  Lsf = Lsf .and. ABS(Xloc(n)) < Ar(n,is,ib)
               end do
               if (Isrftype(is,ib) < 0) then
                  Lsf = .not. Lsf
               end if
               Lbf = Lbf .and. Lsf

            case (3:4)
               ! Surface is a sphere or ellipsoid

               Total = 0.0_r8
               do n = 1, 3
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
               do n = 1, 3
                  Total = Total + (Ar(n,is,ib)*Xloc(n))**2
               end do
               Lsf = Total <= 1.0_r8
               do n = 1, 3
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

end module vof_model
