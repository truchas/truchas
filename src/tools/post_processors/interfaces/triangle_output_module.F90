!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE TRIANGLE_OUTPUT_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Write interface triangle data to files in a form compatible
  !   with visualization software.
  !
  ! Public Interface(s):
  !
  !   * call AVS_OUTPUT (nicells, Xv)
  !
  !     Write interface triangle information out in AVS UCD format.
  !
  !   * call GMV_OUTPUT (nicells, Xv)
  !
  !     Write interface triangle information out in GMV format.
  !
  ! Contains: AVS_OUTPUT
  !           GMV_OUTPUT
  !
  ! Author(s): Douglas B. Kothe (LANL Group T-3, dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only: string_len
  implicit none
  private

  ! Public Procedures
  public :: AVS_OUTPUT, GMV_OUTPUT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  character(string_len), public, save :: avs_file, gmv_file
  logical,    public, save :: write_mesh, write_avs, &
                                               write_gmv

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE AVS_OUTPUT (nicells, Xv, N_Order, Nvrt, Pv)
    !=======================================================================
    ! Purpose(s):
    !
    !   Write surface triangle information out in AVS UCD format.
    !
    !=======================================================================
    use parameter_module,  only: nvc, ndim, nec

    ! Arguments
    integer,                              intent(IN)  :: nicells
    real(r8),   dimension(nvc,ndim,nicells), intent(IN)  :: Xv
    integer, dimension(nec,nicells),      intent(IN)  :: N_Order
    integer, dimension(nicells),          intent(IN)  :: Nvrt
    real(r8),   dimension(ndim,nec,nicells), intent(IN)  :: Pv

    ! Local Variables
    character(10)      :: comp_label, comp_unit
    integer :: l, n1, n2, node_data, ncell_data, &
                                nmod_data, len_trim, num_data_comps, &
                                nsize_comp, nodes, ntri, index, f, &
                                edge, material_id, num_vert, maxtri, col

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize variables
    node_data = 1
    ncell_data = 0
    nmod_data = 0

    num_data_comps = 1
    nsize_comp = 1

    comp_label = 'zax,'
    comp_unit = 'scd'

    ! Compute number of triangles and triangle vertices
    nodes = 0
    do n1 = 1, nicells
       nodes = nodes + Nvrt(n1)
    end do
    ntri = nodes - 2*nicells
    l = LEN_TRIM(avs_file)

    ! Open the output file
    open(2, file = avs_file(:l), status = 'unknown')

    ! Write number of vertices and triangles; inform user
    write (*,105) avs_file(:l)
105 format(/,3x,'AVS UCD file ',a,' written:')
    write (*,100) ntri, nodes
100 format(3x,i5,' interface triangles (',i5,' vertices)')

    ! Alter node and element number if mesh
    !    vertices are to be written out
    if (write_mesh) then
       nodes = nodes + nvc*nicells
       ntri  = ntri + 12*nicells
       write (*,110) 12*nicells, nvc*nicells
110    format(3x,i5,' mesh lines',5x,'(',i5,' vertices)')
    end if

    ! Write out the AVS UCD header line
    write(2,10) nodes, ntri, node_data, ncell_data, nmod_data
10  format (i6,i6,i2,i2,i2)

    ! Write out mesh vertices if requested
    if (write_mesh) then
       index = 0
       do n1 = 1,nicells
          do n2 = 1,nvc
             index = index + 1
             write (2,15) index, Xv(n2,1,n1), Xv(n2,2,n1), Xv(n2,3,n1)
15           format(i5,3f8.3)
          end do
       end do
    end if

    ! Write out triangle vertex coordinates
    if (.not. write_mesh) index = 0
    do n1 = 1, nicells
       num_vert = Nvrt(n1)
       do n2 = 1, num_vert
          edge = N_order(n2, n1)
          index = index + 1
          write (2,20) index, Pv(1,edge,n1), Pv(2,edge,n1), Pv(3,edge,n1)
20        format(i5,3f8.3)
       end do
    end do

    ! Write out cell information: the vertex
    !    numbers that make up each cell
    if (write_mesh) then
       material_id = 2
       index = 0
       f = 1
       do n1 = 1,nicells
          do n2 = 1,12
             index = index + 1
             write(2,25) index, material_id, 'line', &
                  f + (n2 - 1)/3, f + 4 + mod(n2 - 1,4)
25           format(i5, i2, 1x, a4, i6, i6)
          end do
          f = f + nvc
       end do
    end if

    ! Write out triangle information: the vertex
    !    numbers that make up each triangle
    if (write_mesh) then
       index = 12*nicells
       f = nvc*nicells + 1
    else
       index = 0
       f = 1
    end if
    material_id = 1
    do n1 = 1, nicells
       maxtri = Nvrt(n1) - 2
       do n2 = 1, maxtri
          index = index + 1
          write(2,30) index, material_id, 'tri', &
               f, f + n2, f + n2 + 1
30        format(i5, i2, 1x, a3, i6, i6, i6)
       end do
       f = f + Nvrt(n1)
    end do

    ! Write out nodal header information
    write(2,40) num_data_comps, nsize_comp
40  format(2i2)

    write(2,50) comp_label, comp_unit
50  format(2a10)

    ! Write out nodal data
    if (write_mesh) then
       col = 150
       do n1 = 1, nvc*nicells
          write(2,60) n1, col
60        format(i5,i4)
       end do
       col = 50
       do n1 = nvc*nicells + 1,nodes
          write(2,70) n1, col
70        format(i5,i4)
       end do
    else
       col = 50
       do n1 = 1, nodes
          write(2,80) n1, col
80        format(i5,i4)
       end do
    end if

    ! Close the output file
    close (2)

  END SUBROUTINE AVS_OUTPUT

  SUBROUTINE GMV_OUTPUT (nicells, Xv,Normal, Nvrt, Pv, Perm, Mdpnt)
    !=======================================================================
    ! Purpose(s):
    !
    !   Write surface triangle information out in GMV format.
    !
    !=======================================================================
    use parameter_module,  only: ndim, nvc, nec

    ! Arguments
    integer,                              intent(IN)  :: nicells
    real(r8),   dimension(nvc,ndim,nicells), intent(IN)  :: Xv
    real(r8),   dimension(ndim,nicells),     intent(IN)  :: Normal
    integer, dimension(nicells),          intent(IN)  :: Nvrt
    integer, dimension(nicells,nec),      intent(IN)  :: Perm
    real(r8),   dimension(ndim,nec,nicells), intent(IN)  :: Pv
    real(r8),   dimension(nicells,ndim),     intent(IN)  :: Mdpnt

    ! Local Variables
    integer :: l, n1, n2, nodes, ntri, f, material_id, &
                                maxtri, nmats, mix_end, mix_beg, edge1, &
                                edge2, edge3

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize variables
    nmats = 1
    mix_beg = 1
    mix_end = nicells 

    ! Compute number of triangles and triangle vertices
    nodes = 0
    do n1 = mix_beg, mix_end
       nodes = nodes + Nvrt(n1)
    end do
    ntri = nodes - 2*(mix_end - mix_beg + 1)
    l = LEN_TRIM(gmv_file)

    ! Open the output file
    open(2, file = gmv_file(:l), status = 'unknown')

    ! Write number of vertices and triangles; inform user
    write (*,105) gmv_file(:l)
105 format(/,3x,'GMV input file ',a,' written:')
    write (*,100) ntri, nodes
100 format(3x,i5,' interface triangles (',i5,' vertices)')

    ! Write number of vertices and hexes; inform user
    nodes = nvc*(mix_end - mix_beg + 1)
    write (*,110) mix_end - mix_beg + 1, nodes
110 format(3x,i5,' hexes',5x,'(',i5,' vertices)')

    ! Write out the GMV header line
    write (2,4)
4   format ('gmvinput ascii')
    write (2,5) nodes
5   format ('nodes ',i6)

    ! Write out mesh vertices
    do n1 = mix_beg,mix_end
       do n2 = 1,nvc
          write (2,10) Xv(n2,1,n1)
10        format(f8.3)
       end do
    end do
    do n1 = mix_beg,mix_end
       do n2 = 1,nvc
          write (2,10) Xv(n2,2,n1)
       end do
    end do
    do n1 = mix_beg,mix_end
       do n2 = 1,nvc
          write (2,10) Xv(n2,3,n1)
       end do
    end do

    ! Write out cell information: the vertex
    !    numbers that make up each cell
    f = 1
    write (2,15) mix_end - mix_beg + 1
15  format ('cells ',i6)
    do n1 = mix_beg,mix_end
       write (2,20) f,f+1,f+2,f+3,f+4,f+5,f+6,f+7
20     format ('hex 8 ',8(i6,1x))
       f = f + nvc
    end do

    ! Write out cell-based material data
    material_id = 1
    write (2,25) nmats, 0
25  format ('material ',i2,1x,i2)
    do n1 = 1,nmats
       write (2,'(''copper'')')
    enddo
    do n1 = mix_beg,mix_end
       write (2,30) material_id
30     format (i2)
    end do

    ! Write out triangle vertex coordinates
    write (2,35)
35  format ('polygons')

    material_id = 1
    do n1 = mix_beg, mix_end

       maxtri = Nvrt(n1) - 1
       do n2 = 1, maxtri
          edge1 = Perm(n1,n2)
          edge2 = Perm(n1,n2+1)
          write (2,40) material_id, 3, Mdpnt(n1,1), Pv(1,edge1,n1), Pv(1,edge2,n1)
40        format(i2,1x,i1,3(1x,f8.3))
          write (2,41) Mdpnt(n1,2), Pv(2,edge1,n1), Pv(2,edge2,n1)
41        format(4x,3(1x,f8.3))
          write (2,42) Mdpnt(n1,3), Pv(3,edge1,n1), Pv(3,edge2,n1)
42        format(4x,3(1x,f8.3))
       end do

          edge1 = Perm(n1,Nvrt(n1))
          edge2 = Perm(n1,1)
          write (2,40) material_id, 3, Mdpnt(n1,1), Pv(1,edge1,n1), Pv(1,edge2,n1)
          write (2,41) Mdpnt(n1,2), Pv(2,edge1,n1), Pv(2,edge2,n1)
          write (2,42) Mdpnt(n1,3), Pv(3,edge1,n1), Pv(3,edge2,n1)
    end do

    write (2,45)
45  format ('endpoly')

    write (2,48)
48  format ('variables')

    write (2,49)
49  format ('Nx 0')
    write(2,50)Normal(1,:)

50  format(F10.6)

    write (2,51)
51  format ('Ny 0')
    write(2,50)Normal(2,:)

    write (2,52)
52  format ('Nz 0')
    write(2,50)Normal(3,:)

    write (2,53)
53  format ('endvars')

    ! End of GMV input
    write (2,58)
58  format ('endgmv')

    ! Close the output file
    close (2)

  END SUBROUTINE GMV_OUTPUT

END MODULE TRIANGLE_OUTPUT_MODULE
