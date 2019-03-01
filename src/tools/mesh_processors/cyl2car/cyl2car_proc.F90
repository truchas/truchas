!!
!! CYL2CAR_PROC
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! 19 November 2004
!!
!! This module provides the mesh transformation procedures for
!! the utility program cyl2car.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module cyl2car_proc

  use exodus
  use string_utilities
  implicit none
  private

  public :: transform_mesh

  integer, parameter :: r8 = kind(1.0d0)
  real(kind=r8), parameter :: RAD_PER_DEG = 1.745329251994330d-2

  !! Signatures of the 6 faces of a hexahedral element.
  integer, save :: FACE_SIG(6)
  data FACE_SIG /b'00110011', &
                 b'01100110', &
                 b'11001100', &
                 b'10011001', &
                 b'00001111', &
                 b'11110000'/

  !! EVTX(:,E,V,F) are the local endpoints of an edge (parametrized by E and V)
  !! on face F of a hexahedron.  A face may be collapsed by collapsing opposite
  !! edges on a face (E=1,2), and there are two choices of opposing edges (V=1,2).

  integer, parameter :: EVTX(2,2,2,6) = reshape( (/ 1,2, 5,6,  1,5, 2,6, &
                                                    2,3, 6,7,  2,6, 3,7, &
                                                    3,4, 7,8,  3,7, 4,8, &
                                                    1,4, 5,8,  1,5, 4,8, &
                                                    1,2, 3,4,  1,4, 2,3, &
                                                    5,6, 7,8,  5,8, 6,7 /), shape=(/2,2,2,6/) )

  !! WVERT(:,V,F) is the list of 6 local hex vertices that define the wedge
  !! obtained by collapsing face F through the collapse of opposite edges
  !! specified by V (=1,2).

  integer, parameter :: WVERT(6,2,6) = reshape( (/ 2,3,4,6,7,8,  1,4,8,2,3,7, &
                                                   3,4,1,7,8,5,  2,1,5,3,4,8, &
                                                   4,1,2,8,5,6,  3,2,6,4,1,5, &
                                                   1,2,3,5,6,7,  5,6,2,8,7,3, &
                                                   1,5,6,4,8,7,  4,8,5,3,7,6, &
                                                   6,2,1,7,3,4,  5,1,4,6,2,3 /), shape=(/6,2,6/))

  !! WFACE(:,V,F) is the mapping from local hex faces to local wedge faces
  !! when the wedge is formed by collapsing face F in the manner V.

  integer, parameter :: WFACE(6,2,6) = reshape( (/ 0,1,2,3,4,5,  0,5,2,4,1,3, &
                                                   3,0,1,2,4,5,  4,0,5,2,1,3, &
                                                   2,3,0,1,4,5,  2,4,0,5,1,3, &
                                                   1,2,3,0,4,5,  4,2,5,0,3,1, &
                                                   4,3,5,1,0,2,  3,5,1,4,0,2, &
                                                   4,1,5,3,2,0,  1,5,3,4,2,0 /), shape=(/6,2,6/))
contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! TRANSFORM_MESH
 !!
 !! This subroutine takes INMESH with cylindrical-coordinate data and maps
 !! it to OUTMESH with Cartesian-coordinate data, transforming element types
 !! and side set data as necessary.
 !!

  subroutine transform_mesh (inmesh, outmesh)

    type(exodus_mesh), intent(in)  :: inmesh
    type(exodus_mesh), intent(out) :: outmesh

    integer :: j, n, m, l, bi, bo, offset
    logical :: nmask(inmesh%num_node)
    integer :: nmap(inmesh%num_node), tspec(inmesh%num_elem), emap(inmesh%num_elem)
    type(elem_blk) :: hblk, wblk
    type(side_set) :: sset

    ASSERT( defined(inmesh) )

    if (inmesh%num_dim /= 3)  call error ('coordinate data is not three-dimensional!')

    call examine_mesh (inmesh, nmask, nmap, tspec)

    outmesh%title    = inmesh%title
    outmesh%num_dim  = inmesh%num_dim
    outmesh%num_node = maxval(nmap)
    outmesh%num_elem = inmesh%num_elem

   !!!
   !!! COORDINATES

    allocate(outmesh%coord(outmesh%num_dim,outmesh%num_node))

    !! Map the coordinate data (r,theta,z) -> (x,y,z) for nodes that survive.
    m = 0
    do j = 1, inmesh%num_node
      if (nmask(j)) then
        m = m + 1
        outmesh%coord(1,m) = inmesh%coord(1,j) * cos(RAD_PER_DEG*inmesh%coord(2,j))
        outmesh%coord(2,m) = inmesh%coord(1,j) * sin(RAD_PER_DEG*inmesh%coord(2,j))
        outmesh%coord(3,m) = inmesh%coord(3,j)
      end if
    end do

    !! Ensure that nodes mapped to the axis are positioned exactly on the axis.
    do j = 1, inmesh%num_node
      if (.not.nmask(j)) outmesh%coord(1:2,nmap(j)) = 0.0_r8
    end do

   !!!
   !!! ELEMENT BLOCKS

    !! Count the number of element blocks in the transformed mesh.
    offset = 0
    outmesh%num_eblk = inmesh%num_eblk
    do n = 1, inmesh%num_eblk
      m = count(tspec(offset+1:offset+inmesh%eblk(n)%num_elem) == 0)
      if (m > 0 .and. m < inmesh%eblk(n)%num_elem) outmesh%num_eblk = outmesh%num_eblk + 1
      offset = offset + inmesh%eblk(n)%num_elem
    end do
    allocate(outmesh%eblk(outmesh%num_eblk))

    m = 0
    bo = 0
    offset = 0
    do bi = 1, inmesh%num_eblk

      n = inmesh%eblk(bi)%num_elem

      call transform_eblk (inmesh%eblk(bi), tspec(offset+1:offset+n), nmap, hblk, wblk)

      if (hblk%num_elem > 0) then
        bo = bo + 1
	outmesh%eblk(bo) = hblk
      end if

      if (wblk%num_elem > 0) then
        bo = bo + 1
	outmesh%eblk(bo) = wblk
	if (hblk%num_elem > 0) then
	  outmesh%eblk(bo)%ID = hblk%ID + maxval(inmesh%eblk%ID)
	  call info (i_to_c(wblk%num_elem) // ' of ' // i_to_c(n) // ' hexes in block ID ' // &
	    i_to_c(hblk%ID) // ' transformed to wedges and moved to new block ID ' //  &
	    i_to_c(outmesh%eblk(bo)%ID) // '.')
	else
	  call info ('All ' // i_to_c(n) // ' hexes in block ID ' // i_to_c(wblk%ID) // &
	    ' transformed into wedge elements.')
	end if
      end if

      !! Generate the element number mapping.
      do j = 1, n
        if (tspec(offset+j) == 0) then
          m = m + 1
          emap(offset+j) = m
        end if
      end do

      do j = 1, n
        if (tspec(offset+j) /= 0) then
          m = m + 1
          emap(offset+j) = m
        end if
      end do

      offset = offset + n

    end do

   !!!
   !!! SIDE SETS

    !! Count the number of side sets in the transformed mesh.
    outmesh%num_sset = 0
    do n = 1, inmesh%num_sset
      if (any(abs(tspec(inmesh%sset(n)%elem)) /= inmesh%sset(n)%face)) then
        outmesh%num_sset = outmesh%num_sset + 1
      end if
    end do

    if (outmesh%num_sset > 0) allocate(outmesh%sset(outmesh%num_sset))

    !! Transform the side sets.
    m = 0
    do n = 1, inmesh%num_sset
      call transform_sset (inmesh%sset(n), tspec, emap, sset)
      if (sset%num_side > 0) then
        m = m + 1
        outmesh%sset(m) = sset
      end if
    end do

   !!!
   !!! NODE SETS

    outmesh%num_nset = 0
    if (inmesh%num_nset > 0) then
      call warn ('Unable to transform node sets; omitting them from the transformed mesh.')
    end if

  end subroutine transform_mesh

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! EXAMINE_MESH
 !!
 !! This auxillary subroutine takes a mesh with (r,\theta,z) coordinate data
 !! and identifies the topological changes that would occur in transforming
 !! the mesh to cartesian coordinates.  The returned element-based array TSPEC
 !! specifies elements whose type would be changed into a wedge by the
 !! transformation: the absolute value specifes the face that is collapsed,
 !! and its sign one of two possible ways the collapse occurs.  While the mesh
 !! may be comprised of arbitrary types of elements, we only know how to
 !! transform hexahedra, and then only if it has one {r=0} face, and that face
 !! is collapsed by the collapse of opposite edges on the face, to form a
 !! wedge.  Nodes that would be coincident after the transformation are
 !! replaced by a single node.  The true values in NMASK indicate nodes that
 !! are retained.  NMAP(j) gives the new node number of node j.  For the nodes
 !! with true values in NMASK this is a order-preserving, 1-1 mapping onto
 !! {1,2,...,M} where M is the number of nodes in the transformed mesh.  For
 !! the other nodes, their NMAP value is one of these preceding values.
 !!

  subroutine examine_mesh (mesh, nmask, nmap, tspec)

    type(exodus_mesh), intent(in)  :: mesh
    logical,           intent(out) :: nmask(:)
    integer,           intent(out) :: nmap(:), tspec(:)

    integer :: sig, offset, n, j, l, l1, l2, f, v, e
    integer :: equiv(size(nmap))
    integer, pointer :: elem(:)
    logical :: hex_block
    real(kind(mesh%coord)) :: eps, dz(2), min_dz_A, max_dz_A, min_dz_B, max_dz_B
    real(kind(mesh%coord)), pointer :: r(:), z(:)

    ASSERT( defined(mesh) )
    ASSERT( size(nmask) == mesh%num_node )
    ASSERT( size(nmask) == size(nmap) )
    ASSERT( size(tspec) == mesh%num_elem )

    !! Convenient aliases to the coordinate data.
    r => mesh%coord(1,:)
    z => mesh%coord(3,:)

    !! Choose a distance epsilon.
    eps = 1.0d-5 * min(maxval(r) - minval(r), maxval(z) - minval(z))

    if (minval(r) <= -eps) then ! not a valid mesh
      call error ('input mesh has negative R coordinates!')
    end if

    !! Mark all nodes that would be mapped to the z-axis, {r=0}.
    nmask = (abs(r) < eps)

    tspec = 0
    equiv = 0
    offset = 0
    do n = 1, mesh%num_eblk

      !! We only know how to transform 8-node hexahedra.
      hex_block = ((raise_case(mesh%eblk(n)%elem_type(1:3)) == 'HEX') .and. &
                   (size(mesh%eblk(n)%connect,dim=1) == 8))

      EXAMINE_ELEMENT: do j = 1, mesh%eblk(n)%num_elem

        elem => mesh%eblk(n)%connect(:,j)

        if (any(nmask(elem))) then  ! this element needs to be collapsed

          if (.not.hex_block) then  ! we don't know how to do it :-(
            call error ('unable to collapse non-HEX8 element in element block ID ' // &
              i_to_c(mesh%eblk(n)%ID) // '.')
          end if

          !! Identify the face that will be collapsed by the transformation.
          sig = signature(nmask(elem))
          f = 1
          do while (sig /= FACE_SIG(f))
            f = f + 1
            if (f > size(FACE_SIG)) then  ! we don't know how to do it :-(
              call error ('unable to collapse hex ' // i_to_c(j) // ' in element block ID ' // &
                i_to_c(mesh%eblk(n)%ID) // '.')
            end if
          end do

          !! Identify the manner in which the transformation collapses the face.
                dz = abs(z(elem(evtx(2,:,1,f))) - z(elem(evtx(1,:,1,f))))
          min_dz_A = minval(dz)
          max_dz_A = maxval(dz)
                dz = abs(z(elem(evtx(2,:,2,f))) - z(elem(evtx(1,:,2,f))))
          min_dz_B = minval(dz)
          max_dz_B = maxval(dz)

          if (max_dz_A < eps .and. min_dz_B > 2*eps) then
            v = 1
            tspec(offset+j) = f
          else if (max_dz_B < eps .and. min_dz_A > 2*eps) then
            v = 2
            tspec(offset+j) = -f
          else  !! we don't know how to do it :-(
            call error ('unable to collapse hex ' // i_to_c(j) // ' in element block ID ' // &
              i_to_c(mesh%eblk(n)%ID) // '.')
          end if

          !! Condense the nodes on each of the edges collapsed by the transformation.
          do e = 1, 2
            l1 = elem(EVTX(1,e,v,f))
            do while (equiv(l1) /= 0)
              l1 = equiv(l1)
            end do
            l2 = elem(EVTX(2,e,v,f))
            do while (equiv(l2) /= 0)
              l2 = equiv(l2)
            end do
            if (l2 /= l1) equiv(l1) = l2
          end do

        end if

      end do EXAMINE_ELEMENT

      offset = offset + mesh%eblk(n)%num_elem
    end do

   !!!
   !!! GENERATE THE RENUMBERING OF THE NODES.

    nmap = 0
    nmask = (equiv == 0)  ! these are the representative nodes.

    !! Number the represetative nodes.
    n = 0
    do j = 1, size(nmap)
      if (nmask(j)) then
        n = n + 1
        nmap(j) = n
      end if
    end do

    !! Number the remaining nodes.
    do j = 1, size(nmap)
      if (.not.nmask(j)) then ! chase down the representative
        l = equiv(j)
        do while (equiv(l) /= 0)
          l = equiv(l)
        end do
        nmap(j) = nmap(l)
      end if
    end do

  contains   

    integer function signature (mask)
      logical, intent(in) :: mask(:)
      integer :: j
      ASSERT( size(mask) == 8 )
      ASSERT( size(mask) <= bit_size(signature) )
      signature = 0
      do j = 1, size(mask)
        if(mask(j)) signature = ibset(signature, j-1)
      end do
    end function signature

  end subroutine examine_mesh

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! TRANSFORM_EBLK
 !!
 !! This auxillary subroutine transforms the element block data EBLK according
 !! to the node map array NMAP and the element transformation specification
 !! array TSPEC.  Elements whose type is unchanged (TSPEC == 0) are returned
 !! in the element block HBLK, and the remaining elements (necessarily hexes),
 !! which are transformed into wedge elements, are returned in the element
 !! block WBLK.
 !!
 !! The list of nodes that define each element are remapped according to NMAP.
 !! For the hexahedral elements that are transformed into wedge elements, the
 !! parameter array WVERT specifies the 6 hex vertices which define the wedge
 !! element for each of the 12 possible ways a hexahedron transforms into a
 !! wedge (6 possible faces to collapse, 2 ways per face).  The nonzero values
 !! in TSPEC encode which way a particular element is transformed.
 !!
 !! Because an element block can only contain a single type of element, any
 !! transformed wedge elements must be placed in separate element block.
 !! Both transformed element blocks are assigned the ID of the input element
 !! block; the caller will need to adjust these in order to ensure that each
 !! element block has a unique ID.  The ordering of the elements withing the
 !! transformed element blocks preserve the relative ordering within the input
 !! block .
 !!

  subroutine transform_eblk (eblk, tspec, nmap, hblk, wblk)

    type(elem_blk), intent(in)  :: eblk
    integer,        intent(in)  :: tspec(:), nmap(:)
    type(elem_blk), intent(out) :: hblk, wblk

    integer :: j, l, f, v
    logical :: mask(size(tspec))

    ASSERT( defined(eblk) )
    ASSERT( size(tspec) == eblk%num_elem )
    ASSERT( minval(eblk%connect) > 0 )
    ASSERT( maxval(eblk%connect) <= size(nmap) )
    ASSERT( all(abs(tspec) <= 6) )

    !! These cells persist as is.
    mask = (tspec == 0)

    hblk%num_elem = count(mask)
    wblk%num_elem = eblk%num_elem - hblk%num_elem

    if (hblk%num_elem > 0) then
      hblk%ID = eblk%ID
      hblk%elem_type = eblk%elem_type
      allocate(hblk%connect(size(eblk%connect,dim=1),hblk%num_elem))

      !! Pack the elements that remain hexes into the element block.
      l = 0
      do j = 1, eblk%num_elem
        if (mask(j)) then
          l = l + 1
          hblk%connect(:,l) = nmap(eblk%connect(:,j))
        end if
      end do

    end if

    if (wblk%num_elem > 0) then

      wblk%ID = eblk%ID
      wblk%elem_type = 'WEDGE'
      allocate(wblk%connect(6,wblk%num_elem))

      l = 0
      do j = 1, eblk%num_elem
        if (.not.mask(j)) then
          l = l + 1
          if (tspec(j) > 0) then
            f = tspec(j)
            v = 1
          else
            f = -tspec(j)
            v = 2
          end if
          wblk%connect(:,l) = nmap(eblk%connect(WVERT(:,v,f),j))
        end if
      end do

    end if

  end subroutine transform_eblk

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! TRANSFORM_SSET
 !!
 !! This auxillary subroutine transforms the side set data SS_IN according to
 !! the element map array EMAP, and the element transformation array TSPEC.
 !! The transformed side set is returned in SS_OUT.
 !!
 !! Sides sets are a list of (element, face) pairs.  The element numbers must
 !! be remapped according to EMAP.  The face number is unchanged for elements
 !! whose type is unchanged (TSPEC == 0), but for the remaining elements
 !! (necessarily hexahedra) which are transformed into wedge elements, the face
 !! numbers must be mapped appropriately.  The parameter array WFACE defines
 !! this mapping for each of the 12 possible ways a hexahedra transforms into
 !! a wedge (6 possible faces to collapse, 2 ways per face). The nonzero values
 !! in TSPEC encode which way a particular element is transformed.  Note that
 !! some (element, face) pairs in the side set may disappear.
 !!

  subroutine transform_sset (ss_in, tspec, emap, ss_out)

    type(side_set), intent(in)  :: ss_in
    integer,        intent(in)  :: tspec(:), emap(:)
    type(side_set), intent(out) :: ss_out

    integer :: n, j, f
    logical :: mask(ss_in%num_side)

    ASSERT( defined(ss_in) )
    ASSERT( size(tspec) == size(emap) )
    ASSERT( maxval(ss_in%elem) <= size(tspec) )

    mask = (abs(tspec(ss_in%elem)) /= ss_in%face)
    ss_out%num_side = count(mask)

    if (ss_out%num_side > 0) then

      ss_out%ID = ss_in%ID
      allocate(ss_out%elem(ss_out%num_side), ss_out%face(ss_out%num_side))

      n = 0
      do j = 1, ss_in%num_side
        if (mask(j)) then
          n = n + 1
          ss_out%elem(n) = emap(ss_in%elem(j))
          f = tspec(ss_in%elem(j))
          select case (f)
          case (0)
            ss_out%face(n) = ss_in%face(j)
          case (1:)
            ss_out%face(n) = WFACE(ss_in%face(j),1,f)
          case (:-1)
            ss_out%face(n) = WFACE(ss_in%face(j),2,-f)
          end select
        end if
      end do

      n = ss_in%num_side - ss_out%num_side
      if (n > 0) call info (i_to_c(n) // ' of ' // i_to_c(ss_in%num_side) // &
        ' sides removed from side set ID ' // i_to_c(ss_in%ID) // '.')

    else

      call info ('No sides remain in side set ID ' // i_to_c(ss_in%ID) // &
        '; side set will be dropped.')

    end if

  end subroutine transform_sset

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! INFO, WARN, ERROR
 !!
 !! Informational and error message output procedures.
 !!

  subroutine info (string)
    character(len=*), intent(in) :: string
    write(unit=*,fmt='(a)') string
  end subroutine info

  subroutine warn (string)
    character(len=*), intent(in) :: string
    write(unit=0,fmt='(a)') 'WARNING: ' // string
  end subroutine warn

  subroutine error (string)
    character(len=*), intent(in) :: string
    write(unit=0,fmt='(a)') 'ERROR: ' // string
    stop
  end subroutine error

end module cyl2car_proc
