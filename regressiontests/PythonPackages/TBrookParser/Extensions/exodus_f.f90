subroutine parser_read_exodus_mesh_size (path, nnodes, ncells)

   use exodus
   implicit none

   character(len=*), intent(in) :: path
   integer, intent(out) :: nnodes, ncells
   
   call read_exodus_mesh_size (path, nnodes, ncells)

end subroutine parser_read_exodus_mesh_size


subroutine parser_read_exodus_mesh (path, connect, coord, cblock)

  use exodus
  implicit none

  character(len=*), intent(in) :: path
  integer, intent(out) :: connect(8,*), cblock(*)
  double precision, intent(out) :: coord(3,*)

  !! Truchas degenerate hex node numbering to ExodusII element node numbering.
  integer, parameter ::   TET_NODE_MAP(8) = (/1,1,2,3,4,4,4,4/)
  integer, parameter :: WEDGE_NODE_MAP(8) = (/1,4,5,2,3,6,6,3/)
  
  !! ExodusII element types identified by their number of nodes.
  integer, parameter :: TET=4, WEDGE=6, HEX=8
  
  integer :: n, i, j, nodes_per_elem
  logical :: non_hex_mesh
  type(exodus_mesh) :: mesh 

  call read_exodus_mesh (path, mesh)
  
  if (mesh%num_dim /= 3) then
    write(0,'(a,i0)') 'PARSER_READ_EXODUS_MESH: FATAL: not a 3-D mesh: num_dim=', mesh%num_dim
    stop
  end if
  
  do j = 1, mesh%num_node
    coord(:,j) = mesh%coord(:,j)
  end do
  
  !! Translate Exodus element connectivity to Truchas convention.
  n = 0
  non_hex_mesh = .false.
  do i = 1, mesh%num_eblk
    nodes_per_elem = size(mesh%eblk(i)%connect,dim=1)
    !! Translate Exodus element connectivity to Truchas convention.  See the ExodusII
    !! section in MESH_READ from MESH_INPUT_MODULE and the later section that converts
    !! non-hex elements into degenerate hexes as the reference for what must be done here.
    select case (nodes_per_elem)
    case (HEX)
      do j = 1, mesh%eblk(i)%num_elem
        n = n + 1
        cblock(n) = mesh%eblk(i)%ID
        connect(:,n) = mesh%eblk(i)%connect(:,j)
      end do
    case (TET)
      non_hex_mesh = .true.
      do j = 1, mesh%eblk(i)%num_elem
        n = n + 1
        cblock(n) = mesh%eblk(i)%ID
        connect(:,n) = mesh%eblk(i)%connect(TET_NODE_MAP,j)
      end do
    case (WEDGE)
      non_hex_mesh = .true.
      do j = 1, mesh%eblk(i)%num_elem
        n = n + 1
        cblock(n) = mesh%eblk(i)%ID
        connect(:,n) = mesh%eblk(i)%connect(WEDGE_NODE_MAP,j)
      end do
    case default
      write(0,'(a,i0,a)') 'PARSER_READ_EXODUS_MESH: FATAL: unknown element type: ', &
                          nodes_per_elem, '-node element'
      stop
    end select
  end do
  
  if (non_hex_mesh) then
    write(0,*) 'PARSER_READ_EXODUS_MESH: WARNING: mesh contains degenerate hex elements.'
  end if
  
  call destroy (mesh)

end subroutine parser_read_exodus_mesh


subroutine getcellfields (ndim, ncells, nnodes, nvc, connectivity, coords, centroids, numneighbors)

   integer, parameter                                       :: double_kind = KIND(1.0d0)

   integer,                                    INTENT(IN)  :: ndim, ncells, nnodes, nvc
   integer,            dimension(nvc,ncells),  INTENT(IN)  :: connectivity
   real (double_kind), dimension(ndim,nnodes), INTENT(IN)  :: coords
   real (double_kind), dimension(ndim,ncells), INTENT(OUT) :: centroids
   integer,            dimension(ncells),      INTENT(OUT) :: numneighbors

   integer                                                  :: i, j, count, thisnode
   real (double_kind)                                       :: sumx, sumy, sumz
      
   !calculates cell centroid positions

   do i=1, ncells
      sumx  = 0.0
      sumy  = 0.0
      sumz  = 0.0
      count = 0
      do j=1, nvc
         if (connectivity(j,i) > 0) then
            count    = count + 1
            thisnode = connectivity(j,i)
            sumx     = sumx + coords(1,thisnode)
            sumy     = sumy + coords(2,thisnode)
            if (ndim > 2) then
               sumz  = sumz + coords(3,thisnode)
            end if
         end if
      end do
      centroids(1,i) = sumx/count
      centroids(2,i) = sumy/count
      if (ndim > 2) then
         centroids(3,i) = sumz/count
      end if

      !sjc to do: calculate accurate cell num neighbors
      numneighbors(i)  = 18
   end do

   return
end subroutine getcellfields
