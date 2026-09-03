!!
!! HT_2D_VTKHDF_WRITER_TYPE
!!
!! This module defines the VTKHDF writer for the two-dimensional heat-transfer
!! simulation. It writes the fixed unstructured mesh, mesh-associated global
!! and pedigree identifiers, ghost-cell metadata, and temporal cell enthalpy
!! and temperature data. For multimaterial simulations it also writes one
!! temporal scalar volume-fraction dataset per material.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module ht_2d_vtkhdf_writer_type

  use,intrinsic :: iso_fortran_env, only: int8, int32, int64, r8 => real64
  use simulation_environment_type
  use unstr_2d_mesh_type
  use material_model_type
  use parameter_list_type
  use vtkhdf_ug_file_type, only: vtkhdf_ug_file, vtkhdf_cell_data_handle, vtkhdf_field_data_handle, UG_FIXED_MESH
  implicit none
  private

  integer, parameter :: int32_field = 1, int64_field = 2, real64_field = 3

  type :: temporal_field
    character(:), allocatable :: name
    integer :: value_type
    type(vtkhdf_field_data_handle) :: handle
  end type

  type, public :: ht_2d_vtkhdf_writer
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(vtkhdf_ug_file) :: file
    type(vtkhdf_cell_data_handle) :: enthalpy
    type(vtkhdf_cell_data_handle) :: temperature
    type(vtkhdf_cell_data_handle), allocatable :: volume_fraction(:)
    type(temporal_field), allocatable :: temporal_fields(:)
    logical :: is_open = .false.
  contains
    procedure :: open
    procedure :: write_solution
    procedure :: close
  end type ht_2d_vtkhdf_writer

contains

  subroutine open(this, env, mesh, matl_model, temporal_output, stat, errmsg)

    use vtkhdf_vtk_cell_types, only: VTK_TRIANGLE, VTK_QUAD

    class(ht_2d_vtkhdf_writer), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(material_model), intent(in) :: matl_model
    type(parameter_list), target, intent(in) :: temporal_output
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: x(:,:)
    integer, allocatable :: xcnode(:), cnode(:), global_cell_ids(:), global_node_ids(:)
    integer :: j, nnode
    character(:), allocatable :: name
    integer(int8), allocatable :: types(:), cell_ghost_type(:), node_ghost_type(:)
    real(r8) :: scalar_mold

    call this%file%create(trim(env%output_dir)//'/out.vtkhdf', &
        env%comm%mpi_val, stat, errmsg, mode=UG_FIXED_MESH)
    if (stat /= 0) return

    this%mesh => mesh
    allocate(x(3, mesh%nnode))
    x = 0.0_r8
    x(:2,:) = mesh%x(:, :mesh%nnode)
    xcnode = mesh%cstart
    cnode = mesh%cnode
    allocate(types(mesh%ncell))
    do j = 1, mesh%ncell
      nnode = xcnode(j+1) - xcnode(j)
      select case (nnode)
      case (3)
        types(j) = VTK_TRIANGLE
      case (4)
        types(j) = VTK_QUAD
      case default
        stat = 1
        errmsg = 'unsupported 2D cell topology'
        call this%file%close()
        return
      end select
    end do
    call this%file%write_mesh(x, cnode, xcnode, types)
    global_cell_ids = [(mesh%cell_imap%global_index(j), j=1,mesh%ncell)]
    call this%file%write_cell_data('GlobalCellIds', &
        global_cell_ids, attribute='GlobalIds')
    call this%file%write_cell_data('ExternalCellIds', mesh%xcell, &
        attribute='PedigreeIds')
    allocate(cell_ghost_type(mesh%ncell), source=0_int8)
    cell_ghost_type(mesh%ncell_onP+1:) = 1_int8
    call this%file%write_cell_data('vtkGhostType', cell_ghost_type)
    global_node_ids = [(mesh%node_imap%global_index(j), j=1,mesh%nnode)]
    call this%file%write_point_data('GlobalNodeIds', &
        global_node_ids, attribute='GlobalIds')
    call this%file%write_point_data('ExternalNodeIds', mesh%xnode, &
        attribute='PedigreeIds')
    allocate(node_ghost_type(mesh%nnode), source=0_int8)
    node_ghost_type(mesh%nnode_onP+1:) = 1_int8
    call this%file%write_point_data('vtkGhostType', node_ghost_type)

    this%enthalpy = this%file%register_temporal_cell_data('H', scalar_mold)
    this%temperature = this%file%register_temporal_cell_data('T', scalar_mold)
    if (matl_model%nmatl > 1) then
      allocate(this%volume_fraction(matl_model%nmatl))
      do j = 1, matl_model%nmatl
        name = 'vf_' // normalize_material_name(matl_model%matl_name(j))
        this%volume_fraction(j) = &
            this%file%register_temporal_cell_data(name, 0.0_r8)
      end do
    end if
    call register_temporal_fields(this, temporal_output, stat, errmsg)
    if (stat /= 0) then
      call this%file%close()
      return
    end if
    this%is_open = .true.

  end subroutine open

  pure function normalize_material_name(name) result(normalized)
    character(*), intent(in) :: name
    character(:), allocatable :: normalized
    integer :: i

    allocate(character(len(name)) :: normalized)
    normalized = name
    do i = 1, len(name)
      select case (normalized(i:i))
      case ('+', '-', '*', '^', '%')
        normalized(i:i) = '_'
      end select
    end do
  end function normalize_material_name

  subroutine write_solution(this, time, enthalpy, temperature, volume_fraction, temporal_output)

    class(ht_2d_vtkhdf_writer), intent(inout) :: this
    real(r8), intent(in) :: time
    real(r8), intent(in) :: enthalpy(:), temperature(:)
    real(r8), intent(in) :: volume_fraction(:,:)
    type(parameter_list), intent(inout) :: temporal_output
    real(r8), allocatable :: H(:), T(:), v(:)
    integer :: j

    call this%file%start_time_step(time)
    allocate(H(this%mesh%ncell), T(this%mesh%ncell))
    H(:this%mesh%ncell_onP) = enthalpy
    T(:this%mesh%ncell_onP) = temperature
    call this%mesh%cell_imap%gather_offp(H)
    call this%mesh%cell_imap%gather_offp(T)
    call this%file%write_cell_data(this%enthalpy, H)
    call this%file%write_cell_data(this%temperature, T)
    if (allocated(this%volume_fraction)) then
      allocate(v(this%mesh%ncell))
      do j = 1, size(this%volume_fraction)
        v(:this%mesh%ncell_onP) = volume_fraction(j,:)
        call this%mesh%cell_imap%gather_offp(v)
        call this%file%write_cell_data(this%volume_fraction(j), v)
      end do
    end if
    call write_temporal_fields(this, temporal_output)
    call this%file%finalize_time_step()
    call this%file%flush()

  end subroutine write_solution


  !! Register a fixed set of scalar temporal field data from TEMPORAL_OUTPUT.
  !! VTKHDF requires this registration before the first time step is started.
  subroutine register_temporal_fields(this, temporal_output, stat, errmsg)

    class(ht_2d_vtkhdf_writer), intent(inout) :: this
    type(parameter_list), target, intent(in) :: temporal_output
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: iter
    class(*), pointer :: value
    integer :: n

    stat = 0
    n = temporal_output%count()
    if (n == 0) return
    allocate(this%temporal_fields(n))
    iter = parameter_list_iterator(temporal_output)
    n = 0
    do while (.not.iter%at_end())
      n = n + 1
      if (.not.iter%is_scalar()) then
        stat = 1
        errmsg = 'temporal output field "' // iter%name() // '" is not scalar'
        return
      end if
      this%temporal_fields(n)%name = iter%name()
      value => iter%scalar()
      select type (value)
      type is (integer(int32))
        this%temporal_fields(n)%value_type = int32_field
        this%temporal_fields(n)%handle = &
            this%file%register_temporal_field_data(this%temporal_fields(n)%name, 0_int32)
      type is (integer(int64))
        this%temporal_fields(n)%value_type = int64_field
        this%temporal_fields(n)%handle = &
            this%file%register_temporal_field_data(this%temporal_fields(n)%name, 0_int64)
      type is (real(r8))
        this%temporal_fields(n)%value_type = real64_field
        this%temporal_fields(n)%handle = &
            this%file%register_temporal_field_data(this%temporal_fields(n)%name, 0.0_r8)
      class default
        stat = 1
        errmsg = 'temporal output field "' // this%temporal_fields(n)%name // &
            '" must be integer(int32), integer(int64), or real(real64)'
        return
      end select
      call iter%next()
    end do
  end subroutine


  !! Write the current values of the registered scalar temporal field data.
  subroutine write_temporal_fields(this, temporal_output)

    class(ht_2d_vtkhdf_writer), intent(inout) :: this
    type(parameter_list), intent(inout) :: temporal_output

    integer :: j
    integer(int32) :: i32
    integer(int64) :: i64
    real(r8) :: r64

    if (.not.allocated(this%temporal_fields)) then
      if (temporal_output%count() /= 0) error stop 'unregistered temporal output field'
      return
    end if
    if (temporal_output%count() /= size(this%temporal_fields)) then
      error stop 'temporal output field set changed after VTKHDF registration'
    end if
    do j = 1, size(this%temporal_fields)
      select case (this%temporal_fields(j)%value_type)
      case (int32_field)
        call temporal_output%get(this%temporal_fields(j)%name, i32)
        call this%file%write_field_data(this%temporal_fields(j)%handle, i32)
      case (int64_field)
        call temporal_output%get(this%temporal_fields(j)%name, i64)
        call this%file%write_field_data(this%temporal_fields(j)%handle, i64)
      case (real64_field)
        call temporal_output%get(this%temporal_fields(j)%name, r64)
        call this%file%write_field_data(this%temporal_fields(j)%handle, r64)
      end select
    end do
  end subroutine

  subroutine close(this)
    class(ht_2d_vtkhdf_writer), intent(inout) :: this
    if (this%is_open) call this%file%close()
    this%is_open = .false.
    if (allocated(this%volume_fraction)) deallocate(this%volume_fraction)
    if (allocated(this%temporal_fields)) deallocate(this%temporal_fields)
    nullify(this%mesh)
  end subroutine close

end module ht_2d_vtkhdf_writer_type
