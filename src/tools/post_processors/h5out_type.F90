!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module h5out_type

  use danu_module
  use kinds, only: r8
  use,intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR, c_associated
  implicit none
  private
  
  type, public :: h5out
    type(c_ptr) :: file_id  ! Danu HDF file handle
    type(c_ptr) :: mesh_id  ! Danu HDF mesh group handle
    type(c_ptr) :: sim_id   ! Danu HDF simulation group handle
    integer :: nnodes, ncells ! number of cells and nodes in the mesh
    integer, pointer :: cmap(:) => null() ! cell number map (file to original)
    integer, pointer :: nmap(:) => null() ! node number map (file to original)
    integer, pointer :: cycle_list(:) => null() ! list of available cycles
    real(r8), pointer :: time_list(:) => null() ! corresponding list of times
    type(probe), pointer :: probe_list(:) => null()
  end type h5out
  
  public :: h5out_init, h5out_delete
  public :: h5out_seq_num, h5out_cycle_num, h5out_last_seq_num, h5out_last_cycle_num
  public :: h5out_write_cycle_list, h5out_write_restart_file
  public :: h5out_write_probe_list, h5out_write_probe_data, h5out_num_probes
  public :: h5out_num_nodes, h5out_num_cells
  public :: h5out_get_coordinates, h5out_get_connectivity, h5out_get_block_ids
  public :: has_feature ! expose this for mapped restarts

  type :: probe
    character(32) :: dsname, name
    type(c_ptr) :: pid
  end type
  
contains

  subroutine h5out_init (this, h5file, stat, errmsg)
  
    type(h5out), intent(out) :: this
    character(*), intent(in) :: h5file
    integer, intent(out) :: stat
    character(*), intent(out) :: errmsg
    
    integer :: n
    type(c_ptr) :: seq_id
    
    call danu_file_open_rdonly (h5file, this%file_id, stat)
    if (stat /= DANU_SUCCESS) then
      stat = -1
      errmsg = 'unable to open file ' // trim(h5file)
      return
    end if
  
    call mesh_open (this%file_id, 'DEFAULT', this%mesh_id, stat)
    INSIST(stat == DANU_SUCCESS)
    call simulation_open (this%file_id, 'MAIN', this%sim_id, stat)
    INSIST(stat == DANU_SUCCESS)
    
    !! Read the cell and node maps
    call attribute_read (this%mesh_id, 'Number of Elements', this%ncells, stat)
    INSIST(stat == DANU_SUCCESS)
    call attribute_read (this%mesh_id, 'Number of Nodes', this%nnodes, stat)
    INSIST(stat == DANU_SUCCESS)
    allocate(this%cmap(this%ncells))
    call data_read (this%sim_id, 'CELLMAP', this%cmap, stat)
    INSIST(stat == DANU_SUCCESS)
    allocate(this%nmap(this%nnodes))
    call data_read (this%sim_id, 'NODEMAP', this%nmap, stat)
    INSIST(stat == DANU_SUCCESS)
    
    !! Store a list of the available cycles and times.
    call sequence_count (this%sim_id, n, stat)
    INSIST(stat == DANU_SUCCESS)
    allocate(this%cycle_list(n), this%time_list(n))
    do n = 1, size(this%cycle_list)
      call sequence_get_id (this%sim_id, n, seq_id, stat)
      INSIST(stat == DANU_SUCCESS)
      call attribute_read (seq_id, 'cycle', this%cycle_list(n))
      INSIST(stat == DANU_SUCCESS)
      call attribute_read (seq_id, 'time',  this%time_list(n))
      INSIST(stat == DANU_SUCCESS)
    end do
    
    stat = 0
    errmsg = ''

  end subroutine h5out_init
  
  subroutine h5out_delete (this)
    type(h5out), intent(inout) :: this
    if (c_associated(this%file_id)) call danu_file_close (this%file_id)
    this%file_id = C_NULL_PTR
    this%mesh_id = C_NULL_PTR
    this%sim_id  = C_NULL_PTR
    if (associated(this%cmap)) deallocate(this%cmap)
    if (associated(this%nmap)) deallocate(this%nmap)
    if (associated(this%cycle_list)) deallocate(this%cycle_list)
    if (associated(this%time_list)) deallocate(this%time_list)
    if (associated(this%probe_list)) deallocate(this%probe_list)
  end subroutine h5out_delete

  subroutine h5out_write_cycle_list (this, unit)
    type(h5out), intent(in) :: this
    integer, intent(in) :: unit
    integer :: n
    if (size(this%cycle_list) > 0) then
      write(unit,'(a)') ' cycle  time'
      write(unit,'(a)') repeat('-',29)
      do n = 1, size(this%cycle_list)
        write(unit,'(i6,1x,es22.15)') this%cycle_list(n), this%time_list(n)
      end do
    else
      write(unit,'(a)') 'no cycle output exists'
    end if
  end subroutine h5out_write_cycle_list
  
  subroutine h5out_write_probe_list (this, unit)
    type(h5out), intent(inout) :: this
    integer, intent(in) :: unit
    integer :: n
    call init_probe_list (this)
    if (size(this%probe_list) > 0) then
      write(unit,'(a)') 'Index  Probe name'
      write(unit,'(a)') repeat('-',7+len(this%probe_list%name))
      do n = 1, size(this%probe_list)
        write(unit,'(i4,3x,a)') n, trim(this%probe_list(n)%name)
      end do
    else
      write(unit,'(a)') 'no probe data exists'
    end if
  end subroutine h5out_write_probe_list
  
  integer function h5out_seq_num (this, cycle_num)
    type(h5out), intent(in) :: this
    integer, intent(in) :: cycle_num
    do h5out_seq_num = size(this%cycle_list), 1, -1
      if (this%cycle_list(h5out_seq_num) == cycle_num) return
    end do
  end function h5out_seq_num
  
  integer function h5out_cycle_num (this, seq_num)
    type(h5out), intent(in) :: this
    integer, intent(in) :: seq_num
    if (seq_num > 0 .and. seq_num <= size(this%cycle_list)) then
      h5out_cycle_num = this%cycle_list(seq_num)
    else
      h5out_cycle_num = -1
    end if
  end function h5out_cycle_num
  
  integer function h5out_last_seq_num (this)
    type(h5out), intent(in) :: this
    h5out_last_seq_num = size(this%cycle_list)
  end function h5out_last_seq_num
  
  integer function h5out_last_cycle_num (this)
    type(h5out), intent(in) :: this
    if (size(this%cycle_list) > 0) then
      h5out_last_cycle_num = this%cycle_list(size(this%cycle_list))
    else
      h5out_last_cycle_num = -1
    end if
  end function h5out_last_cycle_num
  
  integer function h5out_num_probes (this)
    type(h5out), intent(inout) :: this
    call init_probe_list (this)
    h5out_num_probes = size(this%probe_list)
  end function h5out_num_probes
  
  integer function h5out_num_nodes (this)
    type(h5out), intent(in) :: this
    h5out_num_nodes = this%nnodes
  end function h5out_num_nodes
  
  integer function h5out_num_cells (this)
    type(h5out), intent(in) :: this
    h5out_num_cells = this%ncells
  end function h5out_num_cells
  
  
  subroutine h5out_write_restart_file (this, seq_num, unit)

    use permutations, only: reorder

    type(h5out), intent(inout) :: this
    integer, intent(in) :: seq_num
    integer, intent(in) :: unit

    integer :: stat
    type(c_ptr) :: seq_id

    call sequence_get_id (this%sim_id, seq_num, seq_id, stat)
    INSIST(stat == DANU_SUCCESS)

    call write_header_segment
    call write_mesh_segment
    call write_core_data_segment
    if (has_feature(seq_id,'solid_mechanics')) call write_solid_mech_data_segment
    if (has_feature(seq_id,'species'))         call write_species_data_segment
    if (has_feature(seq_id,'joule_heat'))      call write_joule_heat_data_segment

  contains
  
    subroutine write_header_segment

      real(r8) :: t, dt
      integer  :: n, stat, cycle_num
      character(8), parameter :: fileid = 'TRF-3'
      character(32), pointer :: feat_list(:)

      !! Item 1 -- file format magic number
      write(unit) fileid

      !! Item 2 -- file features
      feat_list => feature_list(seq_id)
      write(unit) size(feat_list)
      do n = 1, size(feat_list)
        write(unit) feat_list(n)
      end do
      deallocate(feat_list)

      !! Item 3 -- simulation specification (free use -- not using any now)
      write(unit) 0

      !! Item 4 -- global data
      call attribute_read (seq_id, 'time', t, stat)
      INSIST(stat == DANU_SUCCESS)
      call attribute_read (seq_id, 'time step', dt, stat)
      INSIST(stat == DANU_SUCCESS)
      call attribute_read (seq_id, 'cycle', cycle_num, stat)
      INSIST(stat == DANU_SUCCESS)
      !call attribute_read (this%mesh_id, 'Number of Elements', this%ncells, stat)
      !INSIST(stat == DANU_SUCCESS)
      !call attribute_read (this%mesh_id, 'Number of Nodes', this%nnodes, stat)
      !INSIST(stat == DANU_SUCCESS)

      write(unit) t
      write(unit) dt
      write(unit) cycle_num
      write(unit) this%ncells
      write(unit) this%nnodes

    end subroutine write_header_segment

    subroutine write_mesh_segment

      integer :: stat, j, k
      integer, allocatable :: cnode(:,:), blockid(:)
      real(r8), allocatable :: x(:), y(:), z(:)

      allocate(cnode(8,this%ncells))
      call mesh_read_connectivity (this%mesh_id, cnode, stat)
      INSIST(stat == DANU_SUCCESS)
      !! Map the internal node numbers to their original values
      do j = 1, size(cnode,2)
        do k = 1, size(cnode,1)
          cnode(k,j) = this%nmap(cnode(k,j))
        end do
      end do
      !! Re-order the cell-based array to the original cell ordering.
      call reorder (cnode, this%cmap, forward=.true.)
      do k = 1, size(cnode,1)
        write(unit) cnode(k,:)
      end do
      deallocate(cnode)

      if (data_exists(this%sim_id, 'BLOCKID')) then
        allocate(blockid(this%ncells))
        call data_read (this%sim_id, 'BLOCKID', blockid, stat)
        INSIST(stat == DANU_SUCCESS)
        !! Re-order the cell-based array to the original cell ordering.
        call reorder (blockid, this%cmap, forward=.true.)
        write(unit) 1
        write(unit) blockid
        deallocate(blockid)
      else
        write(unit) 0
      end if

      allocate(x(this%nnodes), y(this%nnodes), z(this%nnodes))
      call mesh_read_coordinates (this%mesh_id, x, y, z, stat)
      INSIST(stat == DANU_SUCCESS)
      !! Re-order the node-based arrays to the original node ordering.
      call reorder (x, this%nmap, forward=.true.)
      call reorder (y, this%nmap, forward=.true.)
      call reorder (z, this%nmap, forward=.true.)
      write(unit) x
      write(unit) y
      write(unit) z
      deallocate(x, y, z)

    end subroutine write_mesh_segment

    subroutine write_core_data_segment

      integer :: stat, k, nmat, vof_shape(2)
      real(r8) :: array(this%ncells)
      real(r8), allocatable :: array2(:,:)

      !! RHO -- cell densities
      call simulation_data_read (seq_id, 'Z_RHO', array, stat)
      INSIST(stat == DANU_SUCCESS)
      call reorder (array, this%cmap, forward=.true.)
      write(unit) array

      !! TEMP -- cell temperatures
      call simulation_data_read (seq_id, 'Z_TEMP', array, stat)
      INSIST(stat == DANU_SUCCESS)
      call reorder (array, this%cmap, forward=.true.)
      write(unit) array

      !! ENTHALPY -- cell enthalpies
      call simulation_data_read (seq_id, 'Z_ENTHALPY', array, stat)
      INSIST(stat == DANU_SUCCESS)
      call reorder (array, this%cmap, forward=.true.)
      write(unit) array

      if (has_feature(seq_id, 'fluid_flow')) then
        !! PRESSURE -- cell pressures
        call simulation_data_read (seq_id, 'Z_P', array, stat)
        INSIST(stat == DANU_SUCCESS)
        call reorder (array, this%cmap, forward=.true.)
        write(unit) array
        !! VC_X,Y,Z -- cell velocities
        allocate(array2(3,this%ncells))
        call simulation_data_read (seq_id, 'Z_VC', array2, stat)
        INSIST(stat == DANU_SUCCESS)
        call reorder (array2, this%cmap, forward=.true.)
        do k = 1, 3
          write(unit) array2(k,:)
        end do
        deallocate(array2)
        !! VF_1,2,3,4,5,6 -- face fluxing velocities
        allocate(array2(6,this%ncells))
        INSIST(stat == DANU_SUCCESS)
        call simulation_data_read (seq_id, 'Face_Vel', array2, stat)
        call reorder (array2, this%cmap, forward=.true.)
        do k = 1, size(array2,1)
          write(unit) array2(k,:)
        end do
        deallocate(array2)
      else  ! write dummy data
        array = 0.0_r8
        do k = 1, 10
          write(unit) array
        end do
      end if

      !! VF_* -- phase volume fractions
      if (simulation_data_exists(seq_id, 'VOF')) then
        call simulation_data_dimensions (seq_id, 'VOF', vof_shape, stat)
        INSIST(stat == DANU_SUCCESS)
        nmat = vof_shape(1)
        allocate(array2(nmat,this%ncells))
        call simulation_data_read (seq_id, 'VOF', array2, stat)
        call reorder (array2, this%cmap, forward=.true.)
        INSIST(stat == DANU_SUCCESS)
        write(unit) nmat
        do k = 1, nmat
          write(unit) array2(k,:)
        end do
        deallocate(array2)
      else  ! single phase problem -- volume fraction is implicit
        array = 1.0_r8
        write(unit) 1
        write(unit) array
      end if

    end subroutine write_core_data_segment

    subroutine write_solid_mech_data_segment

      integer :: n, k, stat
      real(r8) :: array(this%ncells), tensor(6,this%ncells), vector(3,this%nnodes)
      character(32) :: name

      do n = 1, 12
        !! TOTAL_STRAIN components
        write(name,'(a,i2.2)') 'TOTAL_STRAIN_', n
        call simulation_data_read (seq_id, name, tensor, stat)
        INSIST(stat == DANU_SUCCESS)
        call reorder (tensor, this%cmap, forward=.true.)
        do k = 1, 6
          write(unit) tensor(k,:)
        end do
        !! ELASTIC_STRESS components
        write(name,'(a,i2.2)') 'ELASTIC_STRESS_', n
        call simulation_data_read (seq_id, name, tensor, stat)
        INSIST(stat == DANU_SUCCESS)
        call reorder (tensor, this%cmap, forward=.true.)
        do k = 1, 6
          write(unit) tensor(k,:)
        end do
        !! PLASTIC_STRAIN components
        write(name,'(a,i2.2)') 'PLASTIC_STRAIN_', n
        call simulation_data_read (seq_id, name, tensor, stat)
        INSIST(stat == DANU_SUCCESS)
        call reorder (tensor, this%cmap, forward=.true.)
        do k = 1, 6
          write(unit) tensor(k,:)
        end do
        !! PLASTIC_STRAIN_RATE
        write(name,'(a,i2.2)') 'PLASTIC_STRAIN_RATE_', n
        call simulation_data_read (seq_id, name, array, stat)
        INSIST(stat == DANU_SUCCESS)
        call reorder (array, this%cmap, forward=.true.)
        write(unit) array
      end do

      !! TOTAL_STRAIN components
      call simulation_data_read (seq_id, 'epsilon', tensor, stat)
      INSIST(stat == DANU_SUCCESS)
      call reorder (tensor, this%cmap, forward=.true.)
      do k = 1, 6
        write(unit) tensor(k,:)
      end do
      !! ELASTIC_STRESS components
      call simulation_data_read (seq_id, 'sigma', tensor, stat)
      INSIST(stat == DANU_SUCCESS)
      call reorder (tensor, this%cmap, forward=.true.)
      do k = 1, 6
        write(unit) tensor(k,:)
      end do
      !! PLASTIC_STRAIN components
      call simulation_data_read (seq_id, 'e_plastic', tensor, stat)
      INSIST(stat == DANU_SUCCESS)
      call reorder (tensor, this%cmap, forward=.true.)
      do k = 1, 6
        write(unit) tensor(k,:)
      end do
      !! PLASTIC_STRAIN_RATE
      call simulation_data_read (seq_id, 'epsdot', array, stat)
      INSIST(stat == DANU_SUCCESS)
      call reorder (array, this%cmap, forward=.true.)
      write(unit) array

      !! RHS
      call simulation_data_read (seq_id, 'RHS', vector, stat)
      INSIST(stat == DANU_SUCCESS)
      call reorder (vector, this%nmap, forward=.true.)
      do k = 1, 3
        write(unit) vector(k,:)
      end do

      !! THERMAL_STRAIN components
      call simulation_data_read (seq_id, 'epstherm', tensor, stat)
      INSIST(stat == DANU_SUCCESS)
      call reorder (tensor, this%cmap, forward=.true.)
      do k = 1, 6
        write(unit) tensor(k,:)
      end do

      !! PC_STRAIN components
      call simulation_data_read (seq_id, 'epspc', tensor, stat)
      INSIST(stat == DANU_SUCCESS)
      call reorder (tensor, this%cmap, forward=.true.)
      do k = 1, 6
        write(unit) tensor(k,:)
      end do

      !! DISP -- displacement components
      call simulation_data_read (seq_id, 'Displacement', vector, stat)
      INSIST(stat == DANU_SUCCESS)
      call reorder (vector, this%nmap, forward=.true.)
      do k = 1, 3
        write(unit) vector(k,:)
      end do

    end subroutine write_solid_mech_data_segment

    subroutine write_species_data_segment

      integer :: stat, n, num_species
      character(8) :: name
      real(r8) :: array(this%ncells)

      call attribute_read (this%sim_id, 'NUM_SPECIES', num_species, stat)
      INSIST(stat == DANU_SUCCESS)
      INSIST(num_species > 0)
      write(unit) num_species

      do n = 1, num_species
        write(name,'(a,i0)') 'phi', n
        call simulation_data_read (seq_id, name, array, stat)
        INSIST(stat == DANU_SUCCESS)
        call reorder (array, this%cmap, forward=.true.)
        write(unit) array
      end do

    end subroutine write_species_data_segment

    subroutine write_joule_heat_data_segment

      use,intrinsic :: iso_fortran_env, only: output_unit

      integer  :: stat, n
      real(r8) :: rvalue, t, t_em
      integer  :: shape1(1), shape2(2)
      real(r8), allocatable :: array1(:), array2(:,:)
      type(c_ptr) :: em_id, sim_id
      character(8) :: name

      ! scan through the EMnnn simulations, in order, looking for
      ! the last one whose TIME attribute value is <= the time
      ! attribute for the requested series sequence number.

      call attribute_read (seq_id, 'time', t, stat)
      INSIST(stat == DANU_SUCCESS)

      n = 0
      em_id = C_NULL_PTR
      do
        write(name,'(a,i3.3)') 'EM', n+1
        if (.not.simulation_exists(this%file_id, name)) exit
        call simulation_open (this%file_id, name, sim_id, stat)
        INSIST(stat == DANU_SUCCESS)
        call attribute_read (sim_id, 'TIME', t_em, stat)
        INSIST(stat == DANU_SUCCESS)
        if (t_em > t) exit
        n = n + 1
        em_id = sim_id
      end do
      INSIST(c_associated(em_id))

      !! Diagnostic output concerning which EM simulation is being used.
      write(name,'(a,i3.3)') 'EM', n
      call attribute_read (em_id, 'TIME', t_em, stat)
      INSIST(stat == DANU_SUCCESS)
      write(output_unit,'(a,es13.6,a)') 'Using EM simulation ' // trim(name) // ' (t =', t_em, ')'

      !! FREQ -- source frequency
      call data_read (em_id, 'FREQ', rvalue, stat)
      INSIST(stat == DANU_SUCCESS)
      write(unit) rvalue

      !! UHFS -- uniform field strength
      call data_read (em_id, 'UHFS', rvalue, stat)
      INSIST(stat == DANU_SUCCESS)
      write(unit) rvalue

      !! Coil data
      call data_dimensions (em_id, 'COILS', shape2, stat)
      INSIST(stat == DANU_SUCCESS)
      INSIST(shape2(1) == 7)
      allocate(array2(shape2(1),shape2(2)))
      call data_read (em_id, 'COILS', array2, stat)
      write(unit) size(array2,2)  ! NCOILS
      do n = 1, size(array2,2)
        write(unit) array2(1,n)   ! CURRENT
        write(unit) array2(2:4,n) ! CENTER
        write(unit) array2(5,n)   ! LENGTH
        write(unit) array2(6,n)   ! RADIUS
        write(unit) nint(array2(7,n)) ! NTURNS
      end do
      deallocate(array2)

      !! NMU and MU -- electric permeabilities (on tet mesh)
      call data_dimensions (em_id, 'MU', shape1, stat)
      INSIST(stat == DANU_SUCCESS)
      allocate(array1(shape1(1)))
      call data_read (em_id, 'MU', array1, stat)
      INSIST(stat == DANU_SUCCESS)
      write(unit) size(array1)
      write(unit) array1

      !! NSIGMA and SIGMA -- electric conductivities (on tet mesh)
      call data_read (em_id, 'SIGMA', array1, stat)
      INSIST(stat == DANU_SUCCESS)
      write(unit) size(array1)
      write(unit) array1
      deallocate(array1)

      !! NSIGMA and SIGMA -- electric conductivities (on hex mesh)
      allocate(array1(this%ncells))
      call data_read (em_id, 'JOULE', array1, stat)
      INSIST(stat == DANU_SUCCESS)
      call reorder (array1, this%cmap, forward=.true.)
      write(unit) this%ncells
      write(unit) array1
      deallocate(array1)

    end subroutine write_joule_heat_data_segment

  end subroutine h5out_write_restart_file
  
  subroutine h5out_write_probe_data (this, n, unit)
    type(h5out), intent(inout) :: this
    integer, intent(in) :: n, unit
    
    integer :: stat, index, n1, n2
    character(128) :: string
    real(r8) :: coord(3)
    real(r8), allocatable :: array(:,:)
    type(c_ptr) :: pid
    
    call init_probe_list (this)
    INSIST(n > 0 .and. n <= size(this%probe_list))
    
    pid = this%probe_list(n)%pid
    
    !! Write probe metadata
    call attribute_read (pid, 'NAME', string, stat)
    INSIST(stat == DANU_SUCCESS)
    write(unit,'(2a)') '# probe name: ', trim(string)
    call attribute_read (pid, 'DESCRIPTION', string, stat)
    INSIST(stat == DANU_SUCCESS)
    write(unit,'(2a)') '# probe description: ', trim(string)
    call attribute_read (pid, 'X', coord(1), stat)
    INSIST(stat == DANU_SUCCESS)
    call attribute_read (pid, 'Y', coord(2), stat)
    INSIST(stat == DANU_SUCCESS)
    call attribute_read (pid, 'Z', coord(3), stat)
    INSIST(stat == DANU_SUCCESS)
    write(unit,'(a,3es14.6,a)') '# probe location: [', coord, ']'
    call attribute_read (pid, 'CELL_INDEX', index, stat)
    INSIST(stat == DANU_SUCCESS)
    write(unit,'(a,i0)') '# nearest cell index: ', index
    call attribute_read (pid, 'CELL_X', coord(1), stat)
    INSIST(stat == DANU_SUCCESS)
    call attribute_read (pid, 'CELL_Y', coord(2), stat)
    INSIST(stat == DANU_SUCCESS)
    call attribute_read (pid, 'CELL_Z', coord(3), stat)
    INSIST(stat == DANU_SUCCESS)
    write(unit,'(a,3es14.6,a)') '# nearest cell location: [', coord, ']'
    call attribute_read (pid, 'NODE_INDEX', index, stat)
    INSIST(stat == DANU_SUCCESS)
    write(unit,'(a,i0)') '# nearest node index: ', index
    call attribute_read (pid, 'NODE_X', coord(1), stat)
    INSIST(stat == DANU_SUCCESS)
    call attribute_read (pid, 'NODE_Y', coord(2), stat)
    INSIST(stat == DANU_SUCCESS)
    call attribute_read (pid, 'NODE_Z', coord(3), stat)
    INSIST(stat == DANU_SUCCESS)
    write(unit,'(a,3es14.6,a)') '# nearest node location: [', coord, ']'
    write(unit,'(a)') '# time, value(s)'
    
    !! Write the probe data
    call probe_data_dimensions (this%sim_id, this%probe_list(n)%dsname, n1, n2, stat)
    INSIST(stat == DANU_SUCCESS)
    allocate(array(n1,n2))
    call probe_data_read (pid, array, stat)
    INSIST(stat == DANU_SUCCESS)
    write(string,'(a,i0,a)') '(', size(array,dim=1) - 1, 'es14.6)'
    write(unit,fmt=trim(string)) array(2:,:)
    deallocate(array)
    
  end subroutine h5out_write_probe_data
  
  logical function has_feature (seq_id, feature)
    type(c_ptr), intent(in) :: seq_id
    character(*), intent(in) :: feature
    select case (feature)
    case ('fluid_flow')
      has_feature = simulation_data_exists(seq_id, 'Z_VC')
    case ('solid_mechanics')
      has_feature = simulation_data_exists(seq_id, 'epsilon')
    case ('species')
      has_feature = simulation_data_exists(seq_id, 'phi1')
    case ('joule_heat')
      has_feature = simulation_data_exists(seq_id, 'Joule_P')
    case default
      INSIST(.false.)
    end select
  end function has_feature

  function feature_list (seq_id)
    type(c_ptr), intent(in) :: seq_id
    character(32), pointer :: feature_list(:)
    integer :: n
    character(32) :: name(4)
    logical :: exists(4)
    name(1) = 'fluid_flow'
    name(2) = 'solid_mechanics'
    name(3) = 'species'
    name(4) = 'joule_heat'
    do n = 1, 4
      exists(n) = has_feature(seq_id, trim(name(n)))
    end do
    allocate(feature_list(count(exists)))
    feature_list = pack(name, exists)
  end function feature_list
  
  subroutine init_probe_list (this)
    type(h5out), intent(inout) :: this
    integer :: n, stat
    character(32) :: name
    if (associated(this%probe_list)) return
    call probe_data_count (this%sim_id, n, stat)
    INSIST(stat == DANU_SUCCESS)
    allocate(this%probe_list(n))
    call probe_data_list (this%sim_id, this%probe_list%dsname, stat)
    INSIST(stat == DANU_SUCCESS)
    do n = 1, size(this%probe_list)
      call probe_data_open (this%sim_id, this%probe_list(n)%dsname, this%probe_list(n)%pid, stat)
      INSIST(stat == DANU_SUCCESS)
      call attribute_read (this%probe_list(n)%pid, 'NAME', name, stat)
      INSIST(stat == DANU_SUCCESS)
      this%probe_list(n)%name = trim(name) // ': ' // trim(this%probe_list(n)%dsname(6:))
    end do
  end subroutine init_probe_list
  
  subroutine h5out_get_coordinates (this, coord)
    type(h5out), intent(in) :: this
    real(r8), intent(out) :: coord(:,:)
    integer :: stat
    INSIST(size(coord,1) == 3)
    INSIST(size(coord,2) == this%nnodes)
    call mesh_read_coordinates (this%mesh_id, coord(1,:), coord(2,:), coord(3,:), stat)
    INSIST(stat == DANU_SUCCESS)
  end subroutine h5out_get_coordinates
  
  subroutine h5out_get_connectivity (this, connect)
    type(h5out), intent(in) :: this
    integer, intent(out) :: connect(:,:)
    integer :: stat
    INSIST(size(connect,1) == 8)
    INSIST(size(connect,2) == this%ncells)
    call mesh_read_connectivity (this%mesh_id, connect, stat)
    INSIST(stat == DANU_SUCCESS)
  end subroutine h5out_get_connectivity
  
  subroutine h5out_get_block_ids (this, blockid)
    type(h5out), intent(in) :: this
    integer, intent(out) :: blockid(:)
    integer :: stat
    INSIST(size(blockid) == this%ncells)
    call data_read (this%sim_id, 'BLOCKID', blockid, stat)
    INSIST(stat == DANU_SUCCESS)
  end subroutine h5out_get_block_ids
    
end module h5out_type
