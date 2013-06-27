MODULE TALLY_MODULE
   !----------------------------------------------------------------------------
   ! purpose:
   !
   !    Define various procedures needed to generate and place particles in
   !    each cell, counting the number of particles that fall on either side of
   !    a user-specified interface
   !
   ! public interface:
   !
   !    vf = VOF_INIT_TALLY_POINTS ()
   !
   !       return (in an nbody,ncells array) the volume of each body in each
   !       cell
   !
   ! contains:
   !    TALLY
   !    PRE_TALLY
   !    HEX_VOLUME
   !    IN_OUT_CRSE
   !    IN_OUT_FINE
   !    PARTICLE_VOLUME
   !    TALLY_CRSE
   !    TALLY_FINE
   !
   ! authors:
   !    Douglas B. Kothe (LANL, dbk@lanl.gov)
   !    S. Jay Mosso (LANL, sjm@lanl.gov)
   !----------------------------------------------------------------------------
   use kinds, only: r8
   use truchas_logging_services
   implicit none
   private

   ! public procedures
   public :: VOF_INIT_TALLY_POINTS, PRE_TALLY

   ! public data
   public :: Mesh_Material, Mesh_Material_Tot

   ! material type associated with each cell (read in from the mesh file) 
   ! this shouldn't be here ...
   integer, dimension(:), pointer :: Mesh_Material
   integer, dimension(:), pointer :: Mesh_Material_Tot

CONTAINS

   !--- public routines --------------------------------------------------------

   function VOF_INIT_TALLY_POINTS ()
      !-------------------------------------------------------------------------
      ! purpose:
      !
      !    Generate and place particles in each cell, counting the number of
      !    particles that fall on either side of a user-specified interface
      !-------------------------------------------------------------------------

      use cutoffs_module,       only: cutvof
      use interfaces_module,    only: int_particles, nbody, vof_particles
      use gs_module,            only: EN_MIN_GATHER, EN_MIN_SCATTER, EN_OR_GATHER, EN_OR_SCATTER
      use mesh_module,          only: Cell
      use parameter_module,     only: nicells, nicells_tot, nnodes, ncells, ndim
      use pgslib_module,        only: PGSLib_GLOBAL_SUM

      ! return value
      real(r8), dimension(nbody,ncells) :: VOF_INIT_TALLY_POINTS

      ! local variables
      logical, dimension(ncells) :: Interface_Cell
      logical, dimension(ncells) :: Interface_Cell_old
      logical, dimension(nnodes) :: Node_mask
      integer :: ib
      integer :: iterpack
      integer :: iterpackmax
      integer :: max_particles
      integer :: particlesf
      real(r8) :: hitsmax
      real(r8) :: hitsmaxc
      real(r8), dimension(nnodes) :: Tmp1
      real(r8), dimension(ncells) :: Tmp3
      real(r8), dimension(nbody,ncells) :: Hits_Vol
      character(128) :: message

      !-------------------------------------------------------------------------

      ! Initialize the Hits_Vol array to contain nothing - this is a local
      ! temporary array
      Hits_Vol = 0.0d0

      ! Compute the geometric constants that define the rotated and translated
      ! bodies that were specified in the user input.
      call PRE_TALLY()

      ! Set the maximum allowed number of particles (based on 0.1*cutvof)
      max_particles = (10/cutvof)**(1.0_r8/ndim)

      ! Coarse particle mapping using approximately (particles/4)**ndim
      ! random points in each cell.  Count how many points lie within
      ! each body.
      particlesf    = vof_particles
      vof_particles = int_particles
      nicells       = ncells
      nicells_tot = PGSLib_GLOBAL_SUM(nicells)

      ! announce what's about to happen
      if (particlesf > vof_particles) then
         write (message, 80) nicells_tot, vof_particles, ndim
80       format(4x,'Interface cells will be identified from ', &
            i10,' total cells with ',i3,'**',i1,' particles/cell')
      else
         write (message, 90) nicells_tot, vof_particles, ndim
90       format(4x,'Volume fractions will be computed in ',    &
            i10,' total cells with ',i3,'**',i1,' particles/cell')
      end if
      call TLS_info ('')
      call TLS_info (message)

      ! Do the coarse particle tally

      ! This identifies interface cells, using only a few particles in each
      ! cell.  Using many particles to identify interfaces would be too
      ! expensive.  We could probably look at just the cell vertices and skip
      ! this step.
      call TALLY_CRSE (Hits_Vol)

      ! Fine particle mapping using particles**3 random
      ! points in designated subset of cells.

      ! this whole algorithm needs to be documented - right now it's magic

      FINE_PARTICLE_TALLY: if (particlesf > vof_particles) then

         iterpack      = 1
         iterpackmax   = ncells**(1.0_r8/ndim)
         hitsmaxc      = vof_particles**ndim
         vof_particles = particlesf
         hitsmax       = vof_particles**ndim
         hitsmaxc      = hitsmaxc - 0.01_r8
         Interface_Cell_old = .false.

100      continue

         Interface_Cell = .false.

         do ib = 1,nbody

            where (Hits_Vol(ib,:) > 0                                  &
               .and. (Cell%Volume - Hits_Vol(ib,:)) > cutvof*Cell%Volume) &
               Interface_Cell = .true.

            ! Find the minimum number of hits in the neighborhood of each cell.
            ! First find the min of cells around a vertex.  Then the min of the
            ! result of vertices around a cell.

            Tmp3 = Hits_Vol(ib,:)
            call EN_MIN_SCATTER (Tmp1, Tmp3)
            call EN_MIN_GATHER (Tmp3, Tmp1)

            ! If a cell has some of body IB in it and it has a neighbor that
            ! had NO hits, then fine zone it.
            where (Hits_Vol(ib,:) > 0 .and. Tmp3 == 0) Interface_Cell = .true.

         end do

         ! Find the logical OR of Interface_Cell in the neighborhood of each
         ! cell.  First form the logical OR of cells around a vertex, then the
         ! OR of the result of vertices around a cell.

         call EN_OR_SCATTER (Node_mask, Interface_Cell)
         call EN_OR_GATHER (Interface_Cell, Node_mask)

         where (Interface_Cell_old) Interface_Cell = .false.
         nicells = COUNT(Interface_Cell)
         nicells_tot = PGSLib_GLOBAL_SUM(nicells)

         if (nicells_tot > 0) then

            if (iterpack == 1) then
               do ib = 1,nbody
                  where (Hits_Vol(ib,:) > 0) Hits_Vol(ib,:) = Cell%Volume
               end do
            end if

            do ib = 1,nbody
               where (Interface_Cell) Hits_Vol(ib,:) = 0
            end do

            write (message, 110) nicells_tot, vof_particles, ndim
110         format('  Computing volume fractions in ',i6, &
               ' interface cells with ',i3, '**',i1,' particles/cell')
            call TLS_info ('')
            call TLS_info (message)

            ! Do the fine particle tally
            call TALLY_FINE (Interface_Cell, Hits_Vol)

            Interface_Cell_old = Interface_Cell .or. Interface_Cell_old
            iterpack = iterpack + 1
            if (iterpack <= iterpackmax) go to 100

         end if

      end if FINE_PARTICLE_TALLY

      VOF_INIT_TALLY_POINTS = Hits_Vol

   end function VOF_INIT_TALLY_POINTS

   !----------------------------------------------------------------------------

  SUBROUTINE PRE_TALLY ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the needed geometric constants that define the input
    !   bodies.
    !
    !=======================================================================
    use constants_module,  only: pi
    use interfaces_module, only: Ar, Cosa, Isrftype, nbody, Nsurf, &
                                 Offset, Rotangl, Sina, Sgeom
    use parameter_module,  only: ndim, nrot

    ! Local Variables
    integer :: ib, is, n, n1, n2, na
    real(r8) :: deg, tmp_s

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize relevant variables
    deg = pi/180.

    ! Compute the needed geometric constants.
    BODY_GEOM: do ib = 1,nbody

       SURFACE_GEOM: do is = 1,Nsurf(ib)

          SURFACE_TYPE: if (Isrftype(is,ib) /= 0) then

             do n = 1,nrot
                Cosa(n,is,ib) = COS(-deg*Rotangl(n,is,ib))
                Sina(n,is,ib) = SIN(-deg*Rotangl(n,is,ib))
             end do

             select case (ABS(Isrftype(is,ib)))

                case default

                   ! Do nothing.

                case (2)

                   ! Surface is a box of lengths Sgeom(i,is,ib),i=1,3
                   do n = 1,ndim
                      Ar(n,is,ib)=0.5_r8*ABS(Sgeom(n,is,ib))
                   end do

                case (3)

                   ! Surface is a sphere of radius Sgeom(1,is,ib).
                   Ar(1,is,ib) = 1.0_r8/ABS(Sgeom(1,is,ib))
                   do n = 2,ndim
                      Ar(n,is,ib) = Ar(1,is,ib)
                   end do

                case (4)

                   ! Surface is an ellipsoid.
                   do n = 1,ndim
                      Ar(n,is,ib) = 1.0_r8/ABS(Sgeom(n,is,ib))
                   end do

                case (5)

                   do n = 1,ndim
                      Ar(n,is,ib) = 1.0_r8/ABS(Sgeom(2,is,ib))
                   end do

                   do n = 1,ndim
                      if (Sgeom(1,is,ib) == n) Ar(n,is,ib) = 0
                   end do

             end select

          end if SURFACE_TYPE

       end do SURFACE_GEOM

    end do BODY_GEOM

    ! Define the translations relative to the rotated coordinate system.
    BODY_TRANS: do ib = 1,nbody

       SURFACE_TRANS: do is = 1,Nsurf(ib)

          do n = 1,nrot

             ! Select the rotation axis
             na = n
             if (ndim == 2) na = 3

             ! Select the rotation plane axes (two of them).
             select case (na)
             case (1)
                n1 = 2; n2 = 3
             case (2)
                n1 = 1; n2 = 3
             case (3)
                n1 = 1; n2 = 2
             end select

             ! Rotate the offset coordinates in the rotation plane.
             if (Rotangl(n,is,ib) /= 0) then
                tmp_s            = Cosa(n,is,ib)*Offset(n1,is,ib) - Sina(n,is,ib)*Offset(n2,is,ib)
                Offset(n2,is,ib) = Sina(n,is,ib)*Offset(n1,is,ib) + Cosa(n,is,ib)*Offset(n2,is,ib)
                Offset(n1,is,ib) = tmp_s
             end if

          end do

       end do SURFACE_TRANS

    end do BODY_TRANS

  END SUBROUTINE PRE_TALLY

  !--- private routines --------------------------------------------------------

  SUBROUTINE HEX_VOLUME (nicells, Xv, Volume)
    !=======================================================================
    ! Purpose(s): 
    !
    !   Compute the volume V of nicells hexahedral cells defined by
    !   vertices Xv
    !
    !=======================================================================
    use parameter_module, only: ndim, nvc, nfc, nec
    use pgslib_module,    only: PGSLib_GLOBAL_ANY

    ! Arguments
    integer, intent(IN)  :: nicells
    real(r8), dimension(ndim,nvc,nicells), intent(IN)  :: Xv
    real(r8), dimension(nicells), intent(OUT) :: Volume

    ! Local Variables
    integer :: n, f, v1, v2, v3, v4, v5, v6
    real(r8), dimension(ndim,nicells) :: X1, X2, X3

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize relevant quantities
    Volume = 0

    ! Loop over faces, accumulating the cell volume.
    FACE_LOOP: do f = 1,nfc

       FACES: select case (f)

       case (1)

          ! Side 1
          v1 = 8; v2 = 4; v3 = 7; v4 = 8; v5 = 3; v6 = 4

       case (2)

          ! Side 2
          v1 = 6; v2 = 2; v3 = 5; v4 = 6; v5 = 1; v6 = 2

       case (3)

          ! Side 3
          v1 = 5; v2 = 1; v3 = 8; v4 = 5; v5 = 4; v6 = 1

       case (4)

          ! Side 4
          v1 = 7; v2 = 3; v3 = 6; v4 = 7; v5 = 2; v6 = 3

       case (5)

          ! Side 5
          v1 = 3; v2 = 4; v3 = 2; v4 = 3; v5 = 1; v6 = 4

       case (6)

          ! Side 6
          v1 = 6; v2 = 5; v3 = 7; v4 = 6; v5 = 8; v6 = 5

       end select FACES

       DIMENSION: select case (ndim)

       case (2)
             ! 2-D
          X1 = 0.5_r8
          v1 = v2
          v3 = v5
          v4 = v6
          do n = 1,ndim
             X2(n,:) = Xv(n,v3,:) + Xv(n,v5,:)
             X3(n,:) = Xv(n,v2,:) + Xv(n,v4,:)
          end do
          v1 = 1; v2 = 1; v3 = 2
          Volume = Volume + X1(v1,:)*(X2(v2,:)*X3(v3,:) - X3(v2,:)*X2(v3,:))

       case (3)
          ! 3-D
          do n = 1,ndim
             X1(n,:) = Xv(n,v1,:) + Xv(n,v2,:)
             X2(n,:) = Xv(n,v3,:) + Xv(n,v4,:)
             X3(n,:) = Xv(n,v5,:) + Xv(n,v6,:)
          end do
          do n = 1,ndim
             select case (n)
             case (1)
                v1 = 1; v2 = 2; v3 = 3
             case (2)
                v1 = 2; v2 = 3; v3 = 1
             case (3)
                v1 = 3; v2 = 1; v3 = 2
             end select
             Volume = Volume + X1(v1,:)*(X2(v2,:)*X3(v3,:) - X3(v2,:)*X2(v3,:))
          end do

       end select DIMENSION

    end do FACE_LOOP

    Volume = Volume/REAL(nec)

    ! Make sure volumes are OK.
    if (PGSLib_GLOBAL_ANY(Volume < 0)) then
       call TLS_fatal ('HEX_VOLUME: Negative particle volumes found during volume fraction tallying')
    end if

  END SUBROUTINE HEX_VOLUME

  SUBROUTINE IN_OUT_CRSE (R0, Z0, Xin, vertices, body)
    !=======================================================================
    ! Purpose(s): 
    !
    !   Determine whether an array of points (R0,Z0) lie within an n-sided
    !   polygon, whose vertices are sequential (either clockwise or
    !   counter-clockwise) pairs (r,z) in vectors spanning 0 to n. By
    !   periodic boundary conditions, r(0) = r(n) and z(0) = z(n).
    !
    !=======================================================================
    use interfaces_module, only: Rtab, Ztab
    use parameter_module,  only: ncells

    ! Arguments
    integer, intent(IN)  :: vertices
    integer, intent(IN)  :: body
    real(r8), dimension(ncells), intent(IN)  :: R0, Z0
    real(r8), dimension(ncells), intent(OUT) :: Xin

    ! Local Variables
    integer :: i

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! For a vertical line passing thru R0, count up # of
    ! crossings on sides of the polygon that lie below Z0.
    Xin = 0

    ! Loop over n pairs of vertices (r,z)
    do i = 1,vertices

       ! Test for a vector between vertices (i-1) -> (i)
       ! Cross product (0,i-1) x (0,i) > 0 ...
       where (Rtab(i-1,body) <= R0 .and. R0 < Rtab(i,body)) &
            Xin = Xin +  &
            MAX( 0.0_r8, SIGN(1.0_r8,(Rtab(i-1,body)-R0)*(Ztab(i,body)-Z0) - &
            (Ztab(i-1,body)-Z0)*(Rtab(i,body)-R0)) )

       ! Test for a vector between vertices (i) <- (i-1)
       ! Cross product (0,i-1) x (0,i) < 0 ...
       where (Rtab(i-1,body) > R0 .and. R0 >= Rtab(i,body)) &
            Xin = Xin - &
            MIN( 0.0_r8, SIGN(1.0_r8,(Rtab(i-1,body)-R0)*(Ztab(i,body)-Z0) - &
            (Ztab(i-1,body)-Z0)*(Rtab(i,body)-R0)) )

    end do

    ! # of crossings == ODD  => point(X0,Y0) INSIDE  polygon (in = 1).
    !                == EVEN => point(X0,Y0) OUTSIDE polygon (in = 0).
    Xin = MOD(NINT(Xin),2)

  END SUBROUTINE IN_OUT_CRSE

  SUBROUTINE IN_OUT_FINE (R0, Z0, Xin, vertices, body, nicells)
    !=======================================================================
    ! Purpose(s): 
    !
    !   As part of a fine particle mapping on a subset of cells, determine
    !   whether an array of points (R0,Z0) lie within an n-sided
    !   po1lygon, whose vertices are sequential (either clockwise or
    !   counter-clockwise) pairs (r,z) in vectors spanning 0 to n. By
    !   periodic boundary conditions, r(0) = r(n) and z(0) = z(n).
    !
    !=======================================================================
    use interfaces_module, only: Rtab, Ztab

    ! Arguments
    integer, intent(IN)  :: body
    integer, intent(IN)  :: vertices
    integer, intent(IN)  :: nicells
    real(r8), dimension(nicells), intent(IN)  :: R0, Z0
    real(r8), dimension(nicells), intent(OUT) :: Xin

    ! Local Variables
    integer :: i

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! For a vertical line passing thru R0, count up # of
    !    crossings of sides of the polygon that lie below Z0.
    Xin = 0.0_r8

    ! Loop over n pairs of vertices (r,z)
    do i = 1,vertices

       ! Test for a vector between vertices (i-1) -> (i)
       ! Cross product (0,i-1) x (0,i) > 0 ...
       where (Rtab(i-1,body) <= R0 .and. R0 < Rtab(i,body)) &
            Xin = Xin + &
            MAX( 0.0_r8, SIGN(1.0_r8,(Rtab(i-1,body)-R0)*(Ztab(i,body)-Z0) - &
            (Ztab(i-1,body)-Z0)*(Rtab(i,body)-R0)) )

       ! Test for a vector between vertices (i) <- (i-1)
       ! Cross product (0,i-1) x (0,i) < 0 ...
       where (Rtab(i-1,body) > R0 .and. R0 >= Rtab(i,body)) &
            Xin = Xin - &
            MIN( 0.0_r8, SIGN(1.0_r8,(Rtab(i-1,body)-R0)*(Ztab(i,body)-Z0) - &
            (Ztab(i-1,body)-Z0)*(Rtab(i,body)-R0)) )

    end do

    ! # of crossings: ODD  => point(X0,Y0) INSIDE  polygon (in = 1).
    !                 EVEN => point(X0,Y0) OUTSIDE polygon (in = 0).
    Xin = MOD(NINT(Xin),2)

  END SUBROUTINE IN_OUT_FINE

  SUBROUTINE PARTICLE_VOLUME (Xv, i, j, k, nicells, delta, Volume)
    !=======================================================================
    ! Purpose(s): 
    !
    !   Compute particle subvolumes
    !
    !=======================================================================
    use parameter_module, only: ndim, nvc

    ! Arguments
    integer, intent(IN)    :: i
    integer, intent(IN)    :: j
    integer, intent(IN)    :: k
    integer, intent(IN)    :: nicells
    real(r8), intent(IN)    :: delta
    real(r8), dimension(ndim,nvc,nicells), intent(INOUT) :: Xv
    real(r8), dimension(nicells),          intent(OUT)   :: Volume

    ! Local Variables
    integer :: v, vp
    real(r8) :: coeff, xi, eta, zeta
    real(r8), dimension(ndim,nvc,nicells) :: Xvp

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Interpolation coefficients
    Xvp  = 0.0_r8

    ! Compute physical positions of vertices bounding the particle subvolumes
    do vp = 1,nvc
       select case (vp)
       case (1)
          xi   =    i   *delta
          eta  = (j - 1)*delta
          zeta = (k - 1)*delta
       case (2)
          xi   =    i   *delta
          eta  =    j   *delta
          zeta = (k - 1)*delta
       case (3)
          xi   = (i - 1)*delta
          eta  =    j   *delta
          zeta = (k - 1)*delta
       case (4)
          xi   = (i - 1)*delta
          eta  = (j - 1)*delta
          zeta = (k - 1)*delta
       case (5)
          xi   =    i   *delta
          eta  = (j - 1)*delta
          zeta =    k   *delta
       case (6)
          xi   =    i   *delta
          eta  =    j   *delta
          zeta =    k   *delta
       case (7)
          xi   = (i - 1)*delta
          eta  =    j   *delta
          zeta =    k   *delta
       case (8)
          xi   = (i - 1)*delta
          eta  = (j - 1)*delta
          zeta =    k   *delta
       end select
       do v = 1,nvc
          select case (v)
          case (1)
             coeff = xi*(1.0_r8 - eta)*(1.0_r8 - zeta)
          case (2)
             coeff = xi*eta*(1.0_r8 - zeta)
          case (3)
             coeff = (1.0_r8 - xi)*eta*(1.0_r8 - zeta)
          case (4)
             coeff = (1.0_r8 - xi)*(1.0_r8 - eta)*(1.0_r8 - zeta)
          case (5)
             coeff = xi*(1.0_r8 - eta)*zeta
          case (6)
             coeff = xi*eta*zeta
          case (7)
             coeff = (1.0_r8 - xi)*eta*zeta
          case (8)
             coeff = (1.0_r8 - xi)*(1.0_r8 - eta)*zeta
          end select
          Xvp(:,vp,:) = Xvp(:,vp,:) + coeff*Xv(:,v,:)
       end do
    end do

    ! Now compute the volume of the hexahedra surrounding each particle
    call HEX_VOLUME (nicells, Xvp, Volume)

  END SUBROUTINE PARTICLE_VOLUME

  SUBROUTINE TALLY_CRSE (Hits_Vol)
    !=======================================================================
    ! Purpose(s): 
    !
    !   Do a coarse particle tally from body definitions to mesh.
    !   The primary purpose of this tally is to find the interface
    !   cells.
    !
    !=======================================================================
    use cutoffs_module,       only: cutvof
    use interfaces_module,    only: Ar, Cosa, Isrftype, nbody, Nsurf, Offset, &
                                    Rotangl, Rotpt, Rtab, Sgeom, Sina,        &
                                    vof_particles, Ztab, Mesh_Matnum
    use gs_module,            only: EN_GATHER
    use mesh_module,          only: Cell, Vertex, Vrtx_Bdy, Mesh, mesh_has_cblockid_data
    use parameter_module,     only: ndim, ncells, nvc, nrot, string_len
    use random_module,        only: GENERATE_RANDOM
    use utilities_module,     only: TIMESTAMP

    ! Arguments
    real(r8), dimension(nbody,ncells), intent(OUT) :: Hits_Vol

    ! Local Variables
    logical, dimension(ncells) :: Lsf, Lbf
    integer :: i, j, k, ib, n, overlap, is, isrf, npoly, back_body, &
               Overlap_Body(nbody), v, total_particles, iend, jend, kend, &
               n1, n2, na
    real(r8) :: xi, coeff, eta, zeta, ovn
    real(r8), dimension(ndim) :: Rand
    real(r8), dimension(ncells) :: Subvol, Total, Hitscurrent
    real(r8), dimension(ndim,ncells)     :: X, Xloc
    real(r8), dimension(ndim,nvc,ncells) :: X_v
    character(LEN = string_len) :: run_date
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize particle data.
    total_particles = vof_particles**ndim

    ! Gather the cell vertices
    do n = 1,ndim
       call EN_GATHER (X_v(n,:,:), Vertex%Coord(n), BOUNDARY=Vrtx_Bdy(n)%Data)
    end do

    ! Initialize relevant quantities
    Hits_Vol = 0.0_r8
    ovn      = 1.0_r8/vof_particles
    Subvol   = Cell%Volume/total_particles
    Overlap_Body = 0

    ! Inform the user what's about to happen
    call TIMESTAMP (run_date)
    write (message, 7) vof_particles, ndim, total_particles
7   format(4x,'Tallying ',i8,'**',i1,' =',i7,' interface particles')
    call TLS_info ('')
    call TLS_info (message)
    call TLS_info ('')

    ! Generate total_particles random points in each cell;
        iend = vof_particles
        jend = vof_particles
        kend = vof_particles

    ! Count how many points lie within each body.
    ZETA_LOOP: do k = 1,kend
       ETA_LOOP: do j = 1,jend
          XI_LOOP: do i = 1,iend
             xi = 0.0_r8; eta = 0.0_r8; zeta = 0.0_r8
             call GENERATE_RANDOM(0.0_r8,1.0_r8,ndim,Rand)

             ! Coordinates (x,y,z) are for each particle,
             ! placed in a random rectangular pattern
             do n = 1,ndim
                select case (n)
                   case (1)
                      xi = MIN(1.0_r8,MAX((i - 1.0_r8 + Rand(n))*ovn,0.0_r8))
                   case (2)
                      eta = MIN(1.0_r8,MAX((j - 1.0_r8 + Rand(n))*ovn,0.0_r8))
                   case (3)
                      zeta = MIN(1.0_r8,MAX((k - 1.0_r8 + Rand(n))*ovn,0.0_r8))
                end select
             end do

             ! Compute particle physical positions via linear interpolation
             X = 0.0_r8
             do v = 1,nvc
                select case (v)
                   case (1)
                     coeff = xi*(1.0_r8 - eta)*(1.0_r8 - zeta)
                   case (2)
                     coeff = xi*eta*(1.0_r8 - zeta)
                   case (3)
                     coeff = (1.0_r8 - xi)*eta*(1.0_r8 - zeta)
                   case (4)
                     coeff = (1.0_r8 - xi)*(1.0_r8 - eta)*(1.0_r8 - zeta)
                   case (5)
                     coeff = xi*(1.0_r8 - eta)*zeta
                   case (6)
                     coeff = xi*eta*zeta
                   case (7)
                     coeff = (1.0_r8 - xi)*eta*zeta
                   case (8)
                     coeff = (1.0_r8 - xi)*(1.0_r8 - eta)*zeta
                end select
                X = X + coeff*X_v(:,v,:)
             end do

             Hitscurrent = 0.0_r8

             BODY_LOOP: do ib = 1,nbody

                ! Loop over the bodies in the domain and determine
                ! whether the point falls inside any of the bodies.
                Lbf = .true.

                SURFACE_LOOP: do is = 1,Nsurf(ib)

                   ! Loop over each surface that defines the body and
                   ! determine whether the point lies on the correct side
                   ! of the surface.
                   ! [for fill material (ib=mmat), Nsurf = 0 => skip over]

                   ! Initialize the coordinates local to the surface
                   Xloc = X

                   do n = 1,nrot

                      ! Select the rotation axis
                      na = n
                      if (ndim == 2) na = 3

                      ! Select the rotation plane axes (two of them).
                      select case (na)
                      case (1)
                         n1 = 2; n2 = 3; coeff = 1.0_r8
                      case (2)
                         n1 = 3; n2 = 1; coeff = -1.0_r8
                      case (3)
                         n1 = 1; n2 = 2; coeff = 1.0_r8
                      end select

                      ! Rotate the surface coordinate system.
                      if (Rotangl(n,is,ib) /= 0.0_r8) then
                         Total      = Cosa(n,is,ib)*(Xloc(n1,:)-Rotpt(n1,is,ib)) - &
                                      coeff*Sina(n,is,ib)*(Xloc(n2,:)-Rotpt(n2,is,ib)) + Rotpt(n1,is,ib)
                         Xloc(n2,:) = Cosa(n,is,ib)*(Xloc(n2,:)-Rotpt(n2,is,ib)) + &
                                      coeff*Sina(n,is,ib)*(Xloc(n1,:)-Rotpt(n1,is,ib)) + Rotpt(n2,is,ib)
                         Xloc(n1,:) = Total
                      end if

                   end do

                   ! Translate the surface coordinate system
                   do n = 1,ndim
                      Xloc(n,:) = Xloc(n,:) - Offset(n,is,ib)
                   end do

                   ! Jump to specified surface type and compute whether
                   ! the point lies on the body side of the surface.
                   isrf = ABS(Isrftype(is,ib))

                   select case (isrf)

                   case default

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
                         do n = 1,ndim
                            if (ABS(Sgeom(1,is,ib)) == REAL(n)) Lsf = Xloc(n,:) > 0.0_r8
                         end do
                         if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
                         Lbf = Lbf .and. Lsf

                      case (2)

                         ! Surface is a right parallelopiped whose local origin
                         ! is its centroid.  Its x-y-z half-lengths are half of
                         ! Sgeom(1,is,ib), Sgeom(2,is,ib), and Sgeom(3,is,ib),
                         ! respectively.
                         Lsf = .true.
                         do n = 1,ndim
                            Lsf = Lsf .and. ABS(Xloc(n,:)) < Ar(n,is,ib)
                         end do
                         if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
                         Lbf = Lbf .and. Lsf

                      case (3:4)

                         ! Surface is a sphere or ellipsoid
                         Total = 0.0_r8
                         do n = 1,ndim
                            Total = Total + (Ar(n,is,ib)*Xloc(n,:))**2
                         end do
                         Lsf = Total <= 1.0_r8
                         if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
                         Lbf = Lbf .and. Lsf

                      case (5)

                         ! Surface is a right circular cylinder:
                         ! Parameter Sgeom(1,is,ib) specifies the axis along
                         ! which the cylinder is aligned.  Values of 1., 2., or 3.
                         ! specify the x-, y-, or z-axis, respectively.  Parameter
                         ! Sgeom(2,is,ib) is the cylinder radius and Sgeom(3,is,ib)
                         ! is the cylinder height.
                         Total = 0.0_r8
                         do n = 1,ndim
                            Total = Total + (Ar(n,is,ib)*Xloc(n,:))**2
                         end do
                         Lsf = Total <= 1.0_r8

                         do n = 1,ndim
                            if (ABS(Sgeom(1,is,ib)) == REAL(n)) then
                               Lsf = Lsf .and. (Xloc(n,:) - Sgeom(3,is,ib))*Xloc(n,:) <= 0.0_r8
                            end if
                         end do
                         if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
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
                                  n1 = 2; n2 = 3
                                  if (ndim ==2) n2 = n1
                               case (2)
                                  n1 = 3; n2 = 1
                                  if (ndim ==2) n1 = n2
                               case (3)
                                  n1 = 1; n2 = 2
                            end select
                            if (ABS(Sgeom(1,is,ib)) == REAL(n)) then

                               Total = Sgeom(3,is,ib)/Sgeom(2,is,ib)*(Sgeom(2,is,ib) - Xloc(n,:))
                               Lsf = Xloc(n1,:)**2 + Xloc(n2,:)**2 <= Total**2 .and. &
                                    Sgeom(2,is,ib)*(Sgeom(2,is,ib) - Xloc(n,:)) >= 0.0_r8 .and. &
                                    Sgeom(2,is,ib)*Xloc(n,:) >= 0.0_r8
                            end if
                         end do
                         if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
                         Lbf = Lbf .and. Lsf

                      case (7)

                         ! Body is a tabular n-sided polygon, to be rotated if
                         ! parameter Sgeom(1,is,ib) = 0, or translated if = 1,
                         ! about the axis (1=x, 2=y, 3=z) specified by the
                         ! parameter Sgeom(2,is,ib). The integer part
                         ! of parameter Sgeom(3,is,ib) is npoly, the # of vertices
                         ! in the polygon, given by Rtab(0:npoly,ib) and
                         ! Ztab(0:npoly,ib).
                         npoly = NINT(Sgeom(3,is,ib))
                         Rtab(0,ib) = Rtab(npoly,ib)
                         Ztab(0,ib) = Ztab(npoly,ib)

                         ! Now rotate or translate
                         do n = 1,ndim
                            select case (n)
                            case (1)
                               n1 = 2; n2 = 3
                               if (ndim ==2) n2 = n1
                            case (2)
                               n1 = 3; n2 = 1
                               if (ndim ==2) n1 = n2
                            case (3)
                               n1 = 1; n2 = 2
                            end select
                            if (Sgeom(1,is,ib) == 0.0_r8) then
                               ! Rotate
                               Xloc(n1,:) = SQRT(Xloc(n1,:)**2 + Xloc(n2,:)**2)
                               na = n
                            else
                               ! Translate
                               na = n2
                            end if
                            ! Select axis of rotation or translation
                            if (ABS(Sgeom(2,is,ib)) == REAL(n)) then
                               call IN_OUT_CRSE (Xloc(n1,:), Xloc(na,:), Total, npoly, ib)
                            end if
                         end do
                         Lsf = Total == 1
                         if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
                         Lbf = Lbf .and. Lsf

                      case (8)

                         ! Get material distribution from mesh file; whereever the cell block ID
                         ! array matches the material number of this body, we have a hit.
                         Lsf = .false.
                         if (mesh_has_cblockid_data) then
                            Lsf = (Mesh%CBlockID == Mesh_Matnum(ib))
                            if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
                            Lbf = Lbf .and. Lsf
                         end if

                   end select

                end do SURFACE_LOOP

                if (Nsurf(ib) > 0) then     ! skip background bodies

                   ! Clip current body if overlapped by previously
                   ! defined body by ensuring that particle can
                   ! be in only one body.  Note that Hitscurrent will
                   ! exceed 1.0 if particle is in both the current
                   ! body and a previously defined body.
                   where (Lbf) Hitscurrent = Hitscurrent + 1.0_r8

                   overlap = COUNT(Hitscurrent > 1.0_r8)

                   if (overlap > 0) Overlap_Body(ib) = overlap

                   where (Lbf .and. Hitscurrent < 2.0_r8) Hits_Vol(ib,:) = Hits_Vol(ib,:) + Subvol

                   Hitscurrent = MIN(Hitscurrent,1.0_r8)

                end if
                
             end do BODY_LOOP

          end do XI_LOOP

          ! Done tallying this line of particles
          if (ndim == 2) then
             call TIMESTAMP (run_date)
             write (message, 92) j, vof_particles, j*iend
92           format(4x,'Tallied ',i3,' x ',i3,' =',i7, ' interface particles')
             call TLS_info (message)
          end if

       end do ETA_LOOP

       ! Done tallying this plane of particles
       call TIMESTAMP (run_date)
       write (message, 93) k, vof_particles, k*iend*jend
93     format(4x,'Tallied ',i3,' x ',i3,'**2 =',i7, ' interface particles')
       call TLS_info (message)

    end do ZETA_LOOP

    ! Ensure that cummulative hits in cell <= maximum,
    !    write out body clipping statistics,
    !    and set background body number, if any
    do ib = 1,nbody

       if (Overlap_Body(ib) > 0) then
          write (message, 110) ib, Overlap_Body(ib)
          call TLS_info ('')
          call TLS_info (message)
       end if

    end do

    Total = 0.0_r8
    back_body = 0

    do ib = 1,nbody

       if (Nsurf(ib) == 0) then
          back_body = ib
          cycle
       end if

       Xloc(1,:) = MAX(Total + Hits_Vol(ib,:) - Cell%Volume, 0.0_r8)
       overlap = COUNT(Xloc(1,:) > cutvof*Cell%Volume)

       if(overlap > 0) then
          write (message, 100) ib, overlap
100       format (4x,'Body #',i2,' was clipped in',i7,' cells because of excessive hits')
          call TLS_info ('')
          call TLS_info (message)
       end if

       Hits_Vol(ib,:) = Hits_Vol(ib,:) - Xloc(1,:)
       Total = Total + Hits_Vol(ib,:)

    end do

    ! Process the background body last if we have one
    if (back_body > 0) then

       ib = back_body
       Hits_Vol(ib,:) = Cell%Volume - Total
       overlap = COUNT(Total > 0.0_r8)

       if (overlap > 0) then
          write (message, 110) ib, overlap
110       format (4x,'Body #',i2,' was clipped in',i7, &
                  ' cells because of previously defined bodies')
          call TLS_info ('')
          call TLS_info (message)
       end if

       Total = Total + Hits_Vol(ib,:)

    end if

  END SUBROUTINE TALLY_CRSE

  SUBROUTINE TALLY_FINE (Mask, Hits_Vol)
    !=======================================================================
    ! Purpose(s): 
    !
    !   Do a fine particle tally from body definitions to interface
    !   cells, which are a subset of mesh cells.
    !
    !=======================================================================
    use cutoffs_module,       only: cutvof
    use interfaces_module,    only: Ar, Cosa, Isrftype, nbody, Nsurf, Offset, &
                                    Rotangl, Rotpt, Rtab, Sgeom, Sina,        &
                                    vof_particles, Ztab, Mesh_Matnum
    use gs_module,            only: EN_GATHER
    use mesh_module,          only: Cell, Vertex, Mesh, mesh_has_cblockid_data
    use parameter_module,     only: ndim, nicells, ncells, nvc, nrot, string_len
    use PGSlib_module,        only: PGSlib_GLOBAL_SUM
    use utilities_module,     only: TIMESTAMP

    ! Arguments
    logical, dimension(ncells), intent(IN) :: Mask
    real(r8), dimension(nbody,ncells), intent(INOUT) :: Hits_Vol

    ! Local Variables
    integer :: i, iend, j, jend, k, kend, ib, overlap, is, isrf, n, npoly, p, &
                back_body, v, n1, n2, na
    Character(LEN = string_len) :: run_date
    real(r8) :: xi, eta, zeta, coeff, ovn
    integer,  dimension(nbody) :: Overlap_Body
    real(r8), dimension(ndim) :: Rand, Seed
    logical,  dimension(nicells) :: Lsf, Lbf
    real(r8), dimension(nvc,ncells) :: Xv
    real(r8), dimension(ndim,nvc,nicells) :: Xv_int
    real(r8), dimension(nbody,nicells) :: Hits_Vol_Packed
    real(r8), dimension(ndim,nicells) :: X, Xloc
    real(r8), dimension(nicells) :: Total, Subvol, Hitscurrent, Cell_Volume
    real(r8), allocatable :: Mesh_Matl(:)
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Put the cell volume into a temporary array
    Cell_Volume = PACK (Cell%Volume, Mask)

    ! Gather cell vertices into a packed, interface-cell array
    Xv_Int = 0.0_r8
    do n = 1,ndim
       call EN_GATHER (Xv, Vertex%Coord(n))
       do v = 1,nvc
          Xv_Int(n,v,:) = PACK (Xv(v,:), Mask)
       end do
    end do

    ! If the mesh has cell block ID data, then
    ! allocate the smaller (packed) Mesh_Matl array.
    if (mesh_has_cblockid_data) then
       ALLOCATE (Mesh_Matl(nicells))
       Mesh_Matl = PACK (Mesh%CblockID, Mask)
    end if

    ! Generate total_particles random points in each cell;
    iend = 1; jend = 1; kend = 1
    do n = 1,ndim
       select case (n)
          case (1)
             iend    = vof_particles
             Seed(n) = MOD(SQRT(2.0_r8),1.0_r8)
          case (2)
             jend    = vof_particles
             Seed(n) = MOD(0.5_r8*(SQRT(5.0_r8) + 1.0_r8),1.0_r8)
          case (3)
             kend    = vof_particles
             Seed(n) = MOD(SQRT(3.0_r8),1.0_r8)
       end select
    end do

    ! Initialize relevant variables
    Hits_Vol_Packed = 0.0_r8
    ovn = 1.0_r8/vof_particles
    Overlap_Body = 0

    ! Inform the user what's about to happen
    call TIMESTAMP (run_date)
    write (message, 7) vof_particles, ndim, vof_particles**ndim
7   format(4x,'Tallying ',i8,'**',i1,' =',i7,' VOF particles')
    call TLS_info ('')
    call TLS_info (message)
    call TLS_info ('')

    ! Generate vof_particles**ndim random points in each cell;
    !    count how many points lie within each body.
    p = 0
    ZETA_LOOP: do k = 1,kend

       ETA_LOOP: do j = 1,jend

          XI_LOOP: do i = 1,iend

             ! Increment particle count.
             p = p + 1

             ! Compute particle volumee.s
             call PARTICLE_VOLUME (Xv_Int, i, j, k, nicells, ovn, Subvol)

             ! Set the ramdom numbers.
             xi = 0.0_r8; eta = 0.0_r8; zeta = 0.0_r8
             do n = 1,ndim
                Rand(n) = MOD(p*MOD(p*Seed(n),1.0_r8),1.0_r8)
             end do

             ! Coordinates (x,y,z) are for each particle,
             ! placed in a random rectangular pattern
             do n = 1,ndim
                select case (n)
                   case (1)
                      xi = MIN(1.0_r8,MAX((i - 1.0_r8 + Rand(n))*ovn,0.0_r8))
                   case (2)
                      eta = MIN(1.0_r8,MAX((j - 1.0_r8 + Rand(n))*ovn,0.0_r8))
                   case (3)
                      zeta = MIN(1.0_r8,MAX((k - 1.0_r8 + Rand(n))*ovn,0.0_r8))
                end select
             end do

             ! Compute particle physical positions via linear interpolation
             X = 0.0_r8
             do v = 1,nvc
                select case (v)
                   case (1)
                     coeff = xi*(1.0_r8 - eta)*(1.0_r8 - zeta)
                   case (2)
                     coeff = xi*eta*(1.0_r8 - zeta)
                   case (3)
                     coeff = (1.0_r8 - xi)*eta*(1.0_r8 - zeta)
                   case (4)
                     coeff = (1.0_r8 - xi)*(1.0_r8 - eta)*(1.0_r8 - zeta)
                   case (5)
                     coeff = xi*(1.0_r8 - eta)*zeta
                   case (6)
                     coeff = xi*eta*zeta
                   case (7)
                     coeff = (1.0_r8 - xi)*eta*zeta
                   case (8)
                     coeff = (1.0_r8 - xi)*(1.0_r8 - eta)*zeta
                end select
                do n = 1,ndim
                   X(n,:) = X(n,:) + coeff*Xv_Int(n,v,:)
                end do
             end do

             Hitscurrent = 0.0_r8

             BODY_LOOP: do ib = 1,nbody

                ! Loop over the bodies in the domain and determine
                ! whether the point falls inside any of the bodies.
                Lbf = .true.

                SURFACE_LOOP: do is = 1,Nsurf(ib)

                   ! Loop over each surface that defines the body and
                   ! determine whether the point lies on the correct side
                   ! of the surface.
                   ! [for fill material (ib=mmat), Nsurf = 0 => skip over]

                   ! Initialize the coordinates local to the surface
                   Xloc = X

                   do n = 1,nrot

                      ! Select the rotation axis
                      na = n
                      if (ndim == 2) na = 3

                      ! Select the rotation plane axes (two of them).
                      select case (na)
                      case (1)
                         n1 = 2; n2 = 3; coeff = 1.0_r8
                      case (2)
                         n1 = 3; n2 = 1; coeff = -1.0_r8
                      case (3)
                         n1 = 1; n2 = 2; coeff = 1.0_r8
                      end select

                      ! Rotate the surface coordinate system.
                      if (Rotangl(n,is,ib) /= 0.0_r8) then
                         Total      = Cosa(n,is,ib)*(Xloc(n1,:)-Rotpt(n1,is,ib)) - &
                                      coeff*Sina(n,is,ib)*(Xloc(n2,:)-Rotpt(n2,is,ib)) + Rotpt(n1,is,ib)
                         Xloc(n2,:) = Cosa(n,is,ib)*(Xloc(n2,:)-Rotpt(n2,is,ib)) + &
                                      coeff*Sina(n,is,ib)*(Xloc(n1,:)-Rotpt(n1,is,ib)) + Rotpt(n2,is,ib)
                         Xloc(n1,:) = Total
                      end if

                   end do

                   ! Translate the surface coordinate system
                   do n = 1,ndim
                      Xloc(n,:) = Xloc(n,:) - Offset(n,is,ib)
                   end do

                   ! Jump to specified surface type and compute whether
                   ! the point lies on the body side of the surface.
                   isrf = ABS(Isrftype(is,ib))

                   select case (isrf)

                   case default

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
                         do n = 1,ndim
                            if (ABS(Sgeom(1,is,ib)) == REAL(n)) Lsf = Xloc(n,:) > 0.0_r8
                         end do
                         if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
                         Lbf = Lbf .and. Lsf

                      case (2)

                         ! Surface is a right parallelopiped whose local origin
                         ! is its centroid.  Its x-y-z half-lengths are half of
                         ! Sgeom(1,is,ib), Sgeom(2,is,ib), and Sgeom(3,is,ib),
                         ! respectively.
                         Lsf = .true.
                         do n = 1,ndim
                            Lsf = Lsf .and. ABS(Xloc(n,:)) < Ar(n,is,ib)
                         end do
                         if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
                         Lbf = Lbf .and. Lsf

                      case (3:4)

                         ! Surface is a sphere or ellipsoid
                         Total = 0.0_r8
                         do n = 1,ndim
                            Total = Total + (Ar(n,is,ib)*Xloc(n,:))**2
                         end do
                         Lsf = Total <= 1.0_r8
                         if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
                         Lbf = Lbf .and. Lsf

                      case (5)

                         ! Surface is a right circular cylinder:
                         ! Parameter Sgeom(1,is,ib) specifies the axis along
                         ! which the cylinder is aligned.  Values of 1., 2., or 3.
                         ! specify the x-, y-, or z-axis, respectively.  Parameter
                         ! Sgeom(2,is,ib) is the cylinder radius and Sgeom(3,is,ib)
                         ! is the cylinder height.
                         Total = 0.0_r8
                         do n = 1,ndim
                            Total = Total + (Ar(n,is,ib)*Xloc(n,:))**2
                         end do
                         Lsf = Total <= 1.0_r8

                         do n = 1,ndim
                            if (ABS(Sgeom(1,is,ib)) == REAL(n)) then
                               Lsf = Lsf .and. (Xloc(n,:) - Sgeom(3,is,ib))*Xloc(n,:) <= 0.0_r8
                            end if
                         end do
                         if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
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
                                  n1 = 2; n2 = 3
                                  if (ndim ==2) n2 = n1
                               case (2)
                                  n1 = 3; n2 = 1
                                  if (ndim ==2) n1 = n2
                               case (3)
                                  n1 = 1; n2 = 2
                            end select
                            if (ABS(Sgeom(1,is,ib)) == REAL(n)) then

                               Total = Sgeom(3,is,ib)/Sgeom(2,is,ib)*(Sgeom(2,is,ib) - Xloc(n,:))
                               Lsf = Xloc(n1,:)**2 + Xloc(n2,:)**2 <= Total**2 .and. &
                                    Sgeom(2,is,ib)*(Sgeom(2,is,ib) - Xloc(n,:)) >= 0.0_r8 .and. &
                                    Sgeom(2,is,ib)*Xloc(n,:) >= 0.0_r8
                            end if
                         end do
                         if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
                         Lbf = Lbf .and. Lsf

                      case (7)

                         ! Body is a tabular n-sided polygon, to be rotated if
                         ! parameter Sgeom(1,is,ib) = 0, or translated if = 1,
                         ! about the axis (1=x, 2=y, 3=z) specified by the
                         ! parameter Sgeom(2,is,ib). The integer part
                         ! of parameter Sgeom(3,is,ib) is npoly, the # of vertices
                         ! in the polygon, given by Rtab(0:npoly,ib) and
                         ! Ztab(0:npoly,ib).
                         npoly = NINT(Sgeom(3,is,ib))
                         Rtab(0,ib) = Rtab(npoly,ib)
                         Ztab(0,ib) = Ztab(npoly,ib)

                         ! Now rotate or translate
                         do n = 1,ndim
                            select case (n)
                            case (1)
                               n1 = 2; n2 = 3
                               if (ndim ==2) n2 = n1
                            case (2)
                               n1 = 3; n2 = 1
                               if (ndim ==2) n1 = n2
                            case (3)
                               n1 = 1; n2 = 2
                            end select
                            if (Sgeom(1,is,ib) == 0.0_r8) then
                               ! Rotate
                               Xloc(n1,:) = SQRT(Xloc(n1,:)**2 + Xloc(n2,:)**2)
                               na = n
                            else
                               ! Translate
                               na = n2
                            end if
                            ! Select axis of rotation or translation
                            if (ABS(Sgeom(2,is,ib)) == REAL(n)) then
                               call IN_OUT_CRSE (Xloc(n1,:), Xloc(na,:), Total, npoly, ib)
                            end if
                         end do
                         Lsf = Total == 1
                         if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
                         Lbf = Lbf .and. Lsf

                      case (8)

                         ! Get material distribution from mesh file; whereever the cell block ID
                         ! array matches the material number of this body, we have a hit.
                         Lsf = .false.
                         if (ALLOCATED(Mesh_Matl)) then
                            Lsf = (Mesh_Matl == Mesh_Matnum(ib))
                            if (Isrftype(is,ib) < 0) Lsf = .not. Lsf
                            Lbf = Lbf .and. Lsf
                         end if

                   end select

                end do SURFACE_LOOP

                ! loop over surfaces
                if (Nsurf(ib) > 0) then     ! skip background bodies

                   ! Clip current body if overlapped by previously
                   ! defined body by ensuring that particle can
                   ! be in only one body.  Note that Hitscurrent will
                   ! exceed 1.0 if particle is in both the current
                   ! body and a previously defined body.
                   where (Lbf) Hitscurrent = Hitscurrent + 1.0_r8

                   overlap = COUNT(Hitscurrent > 1.0_r8)

                   if (overlap > 0) Overlap_Body(ib) = overlap

                   where (Lbf .and. Hitscurrent < 2.0_r8) &
                        Hits_Vol_Packed(ib,:) = Hits_Vol_Packed(ib,:) + Subvol

                   Hitscurrent = MIN(Hitscurrent,1.0_r8)

                end if

             end do BODY_LOOP

          end do XI_LOOP

          ! Done tallying this line of particles
          if (ndim == 2) then
             call TIMESTAMP (run_date)
             write (message, 92) j, vof_particles, j*iend
92           format(4x,'Tallied ',i3,' x ',i3,' =',i7, ' interface particles')
             call TLS_info (message)
          end if

       end do ETA_LOOP

       ! Done tallying this plane of particles
       call TIMESTAMP (run_date)
       write (message, 93) k, vof_particles, k*iend*jend
93     format(4x,'Tallied ',i3,' x ',i3,'**2 =',i7, ' VOF particles')
       call TLS_info (message)

    end do ZETA_LOOP

    ! Ensure that cummulative hits in cell <= maximum,
    !    write out body clipping statistics,
    !    and set background body number, if any
    do ib = 1,nbody

       if (Overlap_Body(ib) > 0) then
          write (message, 110) ib, Overlap_Body(ib)
          call TLS_info ('')
          call TLS_info (message)
       end if

    end do

    Total = 0.0_r8
    back_body = 0
    do ib = 1,nbody

       if (Nsurf(ib) == 0) then
          back_body = ib
          cycle
       end if

       Xloc(1,:) = MAX(Total + Hits_Vol_Packed(ib,:) - Cell_Volume, 0.0_r8)
       overlap = COUNT(Xloc(1,:) > cutvof*Cell_Volume)

       if (overlap > 0) then

          write (message, 100) ib, overlap
100       format (4x,'Body #',i2,' was clipped in',i7,' cells because of excessive hits')
          call TLS_info ('')
          call TLS_info (message)

       end if

       Hits_Vol_Packed(ib,:) = Hits_Vol_Packed(ib,:) - Xloc(1,:)
       Total = Total + Hits_Vol_Packed(ib,:)

    end do

    ! Process the background body last if we have one
    if (back_body > 0) then

       ib = back_body
       Hits_Vol_Packed(ib,:) = Cell_Volume - Total
       overlap = COUNT(Total > 0.0_r8)
       overlap = PGSlib_GLOBAL_SUM( overlap )

       if (overlap > 0) then

          write (message, 110) ib, overlap
110       format (4x,'Body #',i2,' was clipped in',i7, &
               ' cells because of previously defined bodies')
          call TLS_info ('')
          call TLS_info (message)

       end if

       Total = Total + Hits_Vol_Packed(ib,:)

    end if

    ! Scatter hits to full array
    do ib = 1,nbody
       Total = Hits_Vol_Packed(ib,:)
       Hits_Vol(ib,:) = UNPACK (Total, Mask, Hits_Vol(ib,:))
    end do

    ! Deallocate the smaller (packed) Mesh_Matl array if necessary.
    if (ALLOCATED(Mesh_Matl)) deallocate(Mesh_Matl)

  END SUBROUTINE TALLY_FINE

END MODULE TALLY_MODULE




