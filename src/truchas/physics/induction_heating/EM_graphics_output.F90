!!
!!  The EM_graphics_output Module
!!
!!    Neil N. Carlson <nnc@newmexico.com>
!!    Last revised 2 Apr 2004
!!
!!  Provides procedures for producing DX graphics output for the induction
!!  heating simulation.  This code originated with the standalone EM solver,
!!  and has been temporarily carried over until the graphics output can be
!!  done through Truchas.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module EM_graphics_output

  use kinds, only: rk => r8
  use string_utilities, only: i_to_c
  use parallel_communication
  use simpl_mesh_type
  use mimetic_discretization
  use data_explorer
!NNC!  use field_probes

  implicit none
  private
  
  public :: export_mesh, initialize_field_output, export_fields, finalize_field_output
!NNC!  public :: initialize_probes, update_probes, write_probes
  
  real(kind=rk), allocatable, public :: probe_point(:,:)
  
  !! DX variables
  type(dx_file),   save, target :: dxf_m, dxf_e, dxf_b, dxf_q
  type(dx_object), save :: dxcon, dxpos
  type(dx_series), save :: eseries, bseries, qseries
  
  !! Probe variables
!NNC!  type(ProbeArray), save :: probes

  !! Counters used to generate sequential data file names.
  integer, save :: sim_num = 0, cycle_num
  
contains
  
  subroutine export_mesh (mesh, eps, mu, sigma)

    use truchas_env, only: output_file_name

    type(simpl_mesh),  intent(in) :: mesh
    real(kind=rk), intent(in) :: eps(:), mu(:), sigma(:)
    
    type(dx_object) :: dxfld
    integer, allocatable :: cnode(:,:), cblock(:)
    integer, pointer :: pdata(:)
    real(kind=rk), allocatable :: x(:,:)
    character(len=256) :: file

    sim_num = sim_num + 1
    cycle_num = 0
    file = output_file_name('EM-mesh-' // i_to_c(sim_num))

    !! Extract the bare collated mesh.
    call mesh%get_global_cnode_array (cnode)
    call mesh%get_global_x_array (x)
    call mesh%get_global_cblock_array (cblock)
    
    allocate(pdata(merge(mesh%cell_ip%global_size,0,is_iop)))
    call collate (spread(this_PE, dim=1, ncopies=mesh%cell_ip%onp_size), pdata)
    
    if (is_IOP) then
      !! Open the DX output file.
      call dx_open (dxf_m, file, append=.false.)

      !! Export the mesh.
      call dx_export_connections (dxf_m, dxcon, cnode, type='tetrahedra')
      call dx_export_array (dxf_m, dxpos, real(x))
      
      call dx_export_field (dxf_m, dxfld, dxcon, dxpos, real(pdata), cc=.true., name='partition')
      call dx_export_field (dxf_m, dxfld, dxcon, dxpos, real(cblock), cc=.true., name='region')
    end if
    
    call export_scalar_field (dxf_m, eps(:mesh%ncell_onP), 'epsilon')
    call export_scalar_field (dxf_m, mu(:mesh%ncell_onP), 'mu')
    call export_scalar_field (dxf_m, sigma(:mesh%ncell_onP), 'sigma')
    
    !! NB: this leaves the file name component intact, which is necessary since
    !! it is referenced by DXCON and DXPOS which are continually used.
    if (is_IOP) call dx_close (dxf_m)

    deallocate(cnode, x, cblock, pdata)

  end subroutine export_mesh


  subroutine initialize_field_output ()

    use truchas_env, only: output_file_name

    character(len=8) :: suffix
    character(len=256) :: file

    cycle_num = cycle_num + 1
    suffix = i_to_c(sim_num) // '.' // i_to_c(cycle_num)

    if (is_IOP) then
      file = output_file_name('EM-Efield-' // trim(suffix))
      call dx_open (dxf_e, file, append=.false.)
      call dx_new_series (eseries)

      file = output_file_name('EM-Bfield-' // trim(suffix))
      call dx_open (dxf_b, file, append=.false.)
      call dx_new_series (bseries)

      file = output_file_name('-EM-Qfield-' // trim(suffix))
      call dx_open (dxf_q, file, append=.false.)
      call dx_new_series (qseries)
    end if

  end subroutine initialize_field_output


  subroutine export_fields (mesh, t, efield, bfield, qfield)
  
    type(simpl_mesh), intent(in) :: mesh
    real(kind=rk),   intent(in) :: t
    real(kind=rk),   intent(in) :: efield(:)
    real(kind=rk),   intent(in) :: bfield(:)
    real(kind=rk),   intent(in) :: qfield(:)

    real :: v(3,mesh%ncell)
    real, pointer :: g_v(:,:), g_q(:)
    type(dx_object) :: dxfld
    
    allocate(g_v(3,merge(mesh%cell_ip%global_size,0,is_iop)))
    allocate(g_q(merge(mesh%cell_ip%global_size,0,is_iop)))
    
    v = w1_vector_on_cells(mesh, efield)
    call collate (v(:,:mesh%ncell_onP), g_v)
    if (is_IOP) then
      call dx_export_field (dxf_e, dxfld, dxcon, dxpos, g_v, cc=.true.)
      call dx_append_to_series (eseries, real(t), dxfld)
    end if
    
    v = w2_vector_on_cells(mesh, bfield)
    call collate (v(:,:mesh%ncell_onP), g_v)
    if (is_IOP) then
      call dx_export_field (dxf_b, dxfld, dxcon, dxpos, g_v, cc=.true.)
      call dx_append_to_series (bseries, real(t), dxfld)
    end if
    
    call collate (real(qfield(:mesh%ncell_onP)), g_q)
    if (is_IOP) then
      call dx_export_field (dxf_q, dxfld, dxcon, dxpos, g_q, cc=.true.)
      call dx_append_to_series (qseries, real(t), dxfld)
    end if
    
    deallocate(g_v, g_q)

  end subroutine export_fields


  subroutine finalize_field_output (q_avg)

    real(kind=rk), intent(in) :: q_avg(:)

    type(dx_object) :: dxo

    call export_scalar_field (dxf_q, q_avg, 'Qavg')

    if (is_IOP) then
      call dx_export_series (dxf_e, dxo, eseries, name='E')
      call dx_delete_series (eseries)
      call dx_close (dxf_e)

      call dx_export_series (dxf_b, dxo, bseries, name='B')
      call dx_delete_series (bseries)
      call dx_close (dxf_b)

      call dx_export_series (dxf_q, dxo, qseries, name='Q')
      call dx_delete_series (qseries)
      call dx_close (dxf_q)
    end if

  end subroutine finalize_field_output


  subroutine export_scalar_field (dxf, field, name)
    type(dx_file), intent(inout), target :: dxf
    real(kind=rk), intent(in) :: field(:)
    character(len=*), intent(in) :: name
    type(dx_object) :: dxfld
    integer :: n
    real, pointer :: g_field(:)
    n = global_sum(size(field))
    allocate(g_field(merge(n,0,is_iop)))
    call collate (real(field), g_field)
    if (is_IOP) then
      call dx_export_field (dxf, dxfld, dxcon, dxpos, g_field, cc=.true., name=trim(name))
    end if
    deallocate(g_field)
  end subroutine export_scalar_field


!NNC!  subroutine initialize_probes (mesh)
!NNC!    type(simpl_mesh), intent(in) :: mesh
!NNC!    integer :: j
!NNC!    if (.not.allocated(probe_point)) return
!NNC!    call new_probe_array (probes, size(probe_point,dim=2))
!NNC!    do j = 1, size(probe_point,dim=2)
!NNC!      call set_probe_point (probes, j, probe_point(:,j))
!NNC!    end do
!NNC!    call initialize_probe_array (probes, mesh, mask)
!NNC!  end subroutine initialize_probes
!NNC!  
!NNC!  
!NNC!  subroutine update_probes (t, efield, bfield, qfield)
!NNC!    real(kind=rk), intent(in) :: t
!NNC!    real(kind=rk), intent(in) :: efield(:)
!NNC!    real(kind=rk), intent(in) :: bfield(:)
!NNC!    real(kind=rk), intent(in) :: qfield(:)
!NNC!    if (.not.allocated(probe_point)) return
!NNC!    call update_probe_array (probes, real(t), w1_field=efield, w2_field=bfield, w3_field=qfield)
!NNC!  end subroutine update_probes
!NNC!  
!NNC!  
!NNC!  subroutine write_probes (name)
!NNC!    character(len=*), intent(in) :: name
!NNC!    real, pointer :: times(:), values(:,:), svalues(:)
!NNC!    integer :: i, j, lun
!NNC!    call new_unit (lun)
!NNC!    call get_probe_times (probes, times)
!NNC!    !! Write the E-probes
!NNC!    do j = 1, num_w1_probes(probes)
!NNC!      open(unit=lun,file=trim(name)//'-E-probe-'//i_to_c(j)//'.dat', status='replace',action='write')
!NNC!      call get_w1_probe_values (probes, j, values)
!NNC!      write(unit=lun,fmt='(4es13.5)') (times(i), values(:,i), i = 1, size(times))
!NNC!      close(unit=lun)
!NNC!    end do
!NNC!    !! Write the B-probes
!NNC!    do j = 1, num_w2_probes(probes)
!NNC!      open(unit=lun,file=trim(name)//'-B-probe-'//i_to_c(j)//'.dat', status='replace',action='write')
!NNC!      call get_w2_probe_values (probes, j, values)
!NNC!      write(unit=lun,fmt='(4es13.5)') (times(i), values(:,i), i = 1, size(times))
!NNC!      close(unit=lun)
!NNC!    end do
!NNC!    !! Write the Q-probes
!NNC!    do j = 1, num_w2_probes(probes)
!NNC!      open(unit=lun,file=trim(name)//'-Q-probe-'//i_to_c(j)//'.dat', status='replace',action='write')
!NNC!      call get_w3_probe_values (probes, j, svalues)
!NNC!      write(unit=lun,fmt='(4es13.5)') (times(i), svalues(i), i = 1, size(times))
!NNC!      close(unit=lun)
!NNC!    end do
!NNC!  end subroutine write_probes

end module EM_graphics_output
