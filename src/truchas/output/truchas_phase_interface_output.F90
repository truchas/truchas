!!
!! TRUCHAS_PHASE_INTERFACE_OUTPUT
!!
!! Robert Chiodi  <robertchiodi@lanl.gov> June 2019
!!
!!
!! Output of an unstructured polygon mesh that represents the phase interface
!! to an HDF5 file. This can then be processed with a python utility to
!! generate an XDMF descriptor file which can be visualized in paraview.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module truchas_phase_interface_output

  use truchas_danu_output_data, only : io_group_size
  use,intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR, c_associated
  use parallel_communication
  use truchas_logging_services
  use truchas_h5_outfile, only: th5_file, th5_mesh_group
  use kinds, only: r8
  implicit none
  private

  public :: TPIO_open, TPIO_close
  public :: TPIO_write_mesh

  type(th5_file) :: interface_outfile
  integer, private, save :: written_times = 0

contains

  subroutine TPIO_open ()
    use truchas_env, only: output_file_name
    call interface_outfile%open (output_file_name('interface.h5'), io_group_size, is_IOP)
  end subroutine TPIO_open

  subroutine TPIO_close ()
    call interface_outfile%close()
  end subroutine TPIO_close

  subroutine TPIO_write_mesh(a_polygon_array)

    use irl_fortran_interface
    use kinds, only: r8
    use time_step_module, only: t, dt, cycle_number    
    use legacy_mesh_api, only: ndim, nvc, ncells, ncells_tot, nnodes, nnodes_tot
    use legacy_mesh_api, only: vertex, mesh, unpermute_mesh_vector, unpermute_vertex_vector, mesh_has_cblockid_data
    use output_control,  only: part
    use truchas_logging_services

    type(Poly_type), intent(in) :: a_polygon_array(:)
    
    integer :: n,j, curr_vert, curr_poly
    integer :: total_verts, max_verts, current_verts, npoly
    integer :: global_totalverts, global_maxverts, global_npoly
    character(10) :: name
    type(th5_mesh_group) :: out_mesh
    real(r8), allocatable :: x(:,:)
    integer, allocatable :: cnode(:,:)

    written_times = written_times + 1

    !! Create the mesh entry.
    write(name, '(a,i0.4)') 'MESH', written_times
    print*,name
    print*,' '
    print*,' '
    print*,' '
    print*,' '     
    call interface_outfile%add_interface_mesh_group(trim(name), nvc, ndim, out_mesh)
    call out_mesh%write_attr('cycle', cycle_number)    
    call out_mesh%write_attr('sequence_number', written_times)
    call out_mesh%write_attr('time', t)    
    call out_mesh%write_attr('time step', dt)
    

    print*,'Wrote attributes'
    print*,' '
    print*,' '
    print*,' '
    print*,' '     

    total_verts = 0
    npoly = 0    
    max_verts = -100000
    do j = 1, ncells      
      current_verts = getNumberOfVertices(a_polygon_array(j))      
      total_verts = total_verts + current_verts
      max_verts = max(max_verts, current_verts)
      if(current_verts > 0) then
        npoly = npoly + 1
      end if
    end do

    global_totalverts = global_sum(total_verts)
    global_npoly  = global_sum(npoly)
    global_maxverts = global_maxval(max_verts)    

    print*,global_totalverts, global_npoly, global_maxverts
    
    !! Write the node coordinates.
    allocate(x(3,total_verts))
    curr_vert = 0
    do j = 1, ncells
      do n = 1, getNumberOfVertices(a_polygon_array(j))
        curr_vert = curr_vert + 1
        x(:,curr_vert) = getPt(a_polygon_array(j),n-1)
      end do
    end do
    call out_mesh%write_coordinates(global_totalverts, x)
    print*,'Wrote mesh'
    print*,' '
    print*,' '
    print*,' '
    print*,' '      
    deallocate(x)

    !! Write the polygon connectivity.
    allocate(cnode(global_maxverts+1, npoly))
    curr_vert = 0
    curr_poly = 0
    do j = 1, ncells
      if(getNumberOfVertices(a_polygon_array(j)) > 0 ) then
        curr_poly = curr_poly + 1
        cnode(1, curr_poly) = global_maxverts
        do n = 1, getNumberOfVertices(a_polygon_array(j))
          cnode(n+1, curr_poly) = curr_vert
          curr_vert = curr_vert + 1                    
        end do
        cnode(getNumberOfVertices(a_polygon_array(j))+1:global_maxverts+1, curr_poly) = curr_vert-1
      end if
    end do
    call out_mesh%write_connectivity(global_npoly, cnode)
    print*,'Wrote connectivity'
    print*,' '
    print*,' '
    print*,' '
    print*,' '      
    deallocate(cnode)


  end subroutine TPIO_write_mesh

end module truchas_phase_interface_output
