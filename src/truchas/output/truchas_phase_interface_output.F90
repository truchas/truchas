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
  use truchas_h5_outfile, only: th5_file, th5_mesh_group
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
    use legacy_mesh_api, only: ncells, ncells_tot, ndim, nvc
    use output_control,  only: part
    use truchas_logging_services

    type(Poly_type), intent(in) :: a_polygon_array(:)
    
    integer :: n,j, curr_vert, ind
    integer :: total_verts, current_verts, npoly
    integer :: global_totalverts, global_npoly
    character(10) :: name
    type(th5_mesh_group) :: out_mesh
    real(r8), allocatable :: x(:,:)
    integer, allocatable :: cnode(:)
    integer, allocatable :: number_of_verts(:)
    integer :: conn_offset

    ! NOTE: This probably won't work for restarts.
    ! See how Danu output handles it.
    written_times = written_times + 1

    !! Create the mesh entry.
    write(name, '(a,i0.4)') 'MESH', written_times
    call interface_outfile%add_interface_mesh_group(trim(name), nvc, ndim, out_mesh)
    call out_mesh%write_attr('cycle', 0)!cycle_number)    
    call out_mesh%write_attr('sequence_number', written_times)
    call out_mesh%write_attr('time', t)    
    call out_mesh%write_attr('time step', dt)
    
    total_verts = 0
    npoly = 0    
    do j = 1, ncells      
      current_verts = getNumberOfVertices(a_polygon_array(j))      
      total_verts = total_verts + current_verts
      if(current_verts > 0) then
        npoly = npoly + 1
      end if
    end do

    global_totalverts = global_sum(total_verts)
    global_npoly  = global_sum(npoly)

    !! Write the node coordinates.
    if(global_npoly >1) then
      allocate(x(3,total_verts))
      curr_vert = 0
      do j = 1, ncells
        do n = 1, getNumberOfVertices(a_polygon_array(j))
          curr_vert = curr_vert + 1
          x(:,curr_vert) = getPt(a_polygon_array(j),n-1)
        end do
      end do
    else
      allocate(x(3,3))
      x = 0.0_r8
      global_totalverts = 3*nPE
    end if
    call out_mesh%write_coordinates(global_totalverts, x)
    deallocate(x)

    !! Write the polygon connectivity.
    if(global_npoly > 1) then
      allocate(cnode(2*npoly+total_verts))

      ! This really should be an ALLGATHER
      ! Does PGSLib have such a thing?
      ! Also, here I am assuming that the surface
      ! nodes are written in order of proccessor
      ! this_PE, with order preserved during
      ! writing in Scorpio.
      allocate(number_of_verts(nPE))
      number_of_verts = 0
      number_of_verts(this_PE) = total_verts
      do j = 1, nPE
         number_of_verts(j) = global_sum(number_of_verts(j))
      end do
      conn_offset = sum(number_of_verts(1:this_PE-1))        
      deallocate(number_of_verts)

      ! NOTE: It appears mixed polygon topologies with XDMF
      ! do not respect setting BaseOffset (to 1). So here,
      ! the connectivity is written out as 0 referenced.
      ! The write-xdmf utility also ignores the Offset
      ! attribute written, which then defaults the XDMF
      ! to 0-based as well.
      curr_vert = conn_offset
      ind = 0
      do j = 1, ncells
        if(getNumberOfVertices(a_polygon_array(j)) > 0 ) then
          ind = ind + 1
          cnode(ind) = 3 ! Polygon type identifier in XDMF
          ind = ind + 1
          cnode(ind) = getNumberOfVertices(a_polygon_array(j))
          do n = 1, getNumberOfVertices(a_polygon_array(j))
            ind = ind + 1
            cnode(ind) = curr_vert
            curr_vert = curr_vert + 1            
          end do
        end if
      end do
    else
      allocate(cnode(5))
      global_npoly = nPE
      cnode(1:2) = 3
      cnode(3:5) = [0, 1, 2]
    end if
    call out_mesh%write_connectivity_1DConnect(global_npoly, 2*global_npoly+global_totalverts, cnode)
    deallocate(cnode)


  end subroutine TPIO_write_mesh

end module truchas_phase_interface_output
