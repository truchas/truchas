!!
!! NS_HT_2D_VTKHDF_WRITER_TYPE
!!
!! This module defines NS_HT_2D_VTKHDF_WRITER, which writes the fixed mesh,
!! coupled two-dimensional Navier--Stokes/thermal solution, and solver-
!! published temporal scalar field data to VTKHDF.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module ns_ht_2d_vtkhdf_writer_type

  use,intrinsic :: iso_fortran_env, only: int8, int32, int64, r8 => real64
  use,intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  use parameter_list_type
  use simulation_environment_type
  use unstr_2d_mesh_type
  use material_model_type
  use vtkhdf_ug_file_type, only: vtkhdf_ug_file, vtkhdf_cell_data_handle, &
      vtkhdf_field_data_handle, UG_FIXED_MESH
  implicit none
  private

  integer, parameter :: int32_field = 1, int64_field = 2, real64_field = 3

  type :: temporal_field
    character(:), allocatable :: name
    integer :: value_type
    type(vtkhdf_field_data_handle) :: handle
  end type

  type, public :: ns_ht_2d_vtkhdf_writer
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(vtkhdf_ug_file) :: file
    type(vtkhdf_cell_data_handle) :: pressure, velocity, enthalpy, temperature
    type(vtkhdf_cell_data_handle), allocatable :: vfrac(:)
    type(temporal_field), allocatable :: temporal_fields(:)
    logical :: is_open = .false.
  contains
    procedure :: open
    procedure :: write_solution
    procedure :: close
  end type

contains

  subroutine open(this, env, mesh, matl_model, temporal_output, stat, errmsg)
    use vtkhdf_vtk_cell_types, only: VTK_TRIANGLE, VTK_QUAD

    class(ns_ht_2d_vtkhdf_writer), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(material_model), intent(in) :: matl_model
    type(parameter_list), intent(in) :: temporal_output
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: x(:,:)
    integer, allocatable :: cnode(:), xcnode(:), global_cell_ids(:), global_node_ids(:)
    integer :: c, m, nnode, nfield, nphase, first, last, pid
    character(:), allocatable :: name
    integer(int8), allocatable :: cell_type(:), cell_ghost_type(:), node_ghost_type(:)
    real(r8) :: scalar_mold, vector_mold(3)

    call this%file%create(trim(env%output_dir)//'/out.vtkhdf', env%comm%mpi_val, stat, errmsg, mode=UG_FIXED_MESH)
    if (stat /= 0) return
    this%mesh => mesh
    allocate(x(3,mesh%nnode), cell_type(mesh%ncell))
    x = 0.0_r8
    x(:2,:) = mesh%x(:,:mesh%nnode)
    xcnode = mesh%cstart
    cnode = mesh%cnode
    do c = 1, mesh%ncell
      nnode = xcnode(c+1) - xcnode(c)
      select case (nnode)
      case (3); cell_type(c) = VTK_TRIANGLE
      case (4); cell_type(c) = VTK_QUAD
      case default
        stat = 1
        errmsg = 'unsupported 2D cell topology'
        call this%file%close()
        return
      end select
    end do
    call this%file%write_mesh(x, cnode, xcnode, cell_type)
    global_cell_ids = [(mesh%cell_imap%global_index(c), c=1,mesh%ncell)]
    call this%file%write_cell_data('GlobalCellIds', global_cell_ids, attribute='GlobalIds')
    call this%file%write_cell_data('ExternalCellIds', mesh%xcell, attribute='PedigreeIds')
    allocate(cell_ghost_type(mesh%ncell), source=0_int8)
    cell_ghost_type(mesh%ncell_onP+1:) = 1_int8
    call this%file%write_cell_data('vtkGhostType', cell_ghost_type)
    global_node_ids = [(mesh%node_imap%global_index(c), c=1,mesh%nnode)]
    call this%file%write_point_data('GlobalNodeIds', global_node_ids, attribute='GlobalIds')
    call this%file%write_point_data('ExternalNodeIds', mesh%xnode, attribute='PedigreeIds')
    allocate(node_ghost_type(mesh%nnode), source=0_int8)
    node_ghost_type(mesh%nnode_onP+1:) = 1_int8
    call this%file%write_point_data('vtkGhostType', node_ghost_type)
    this%pressure = this%file%register_temporal_cell_data('pressure', scalar_mold)
    this%velocity = this%file%register_temporal_cell_data('velocity', vector_mold)
    this%enthalpy = this%file%register_temporal_cell_data('H', scalar_mold)
    this%temperature = this%file%register_temporal_cell_data('T', scalar_mold)
    nfield = 0
    do m = 1, matl_model%nmatl
      if (matl_model%have_void .and. m == matl_model%nmatl) then
        nphase = 1
      else
        nphase = matl_model%num_matl_phase(m)
      end if
      if (matl_model%nmatl > 1 .or. nphase > 1) nfield = nfield + nphase
    end do
    if (nfield > 0) then
      allocate(this%vfrac(nfield))
      nfield = 0
      do m = 1, matl_model%nmatl
        if (matl_model%have_void .and. m == matl_model%nmatl) then
          nphase = 1
        else
          nphase = matl_model%num_matl_phase(m)
        end if
        if (matl_model%nmatl == 1 .and. nphase == 1) cycle
        if (nphase == 1) then
          nfield = nfield + 1
          name = 'vf_' // normalize_material_name(matl_model%matl_name(m))
          this%vfrac(nfield) = this%file%register_temporal_cell_data(name, scalar_mold)
        else
          call matl_model%get_matl_phase_index_range(m, first, last)
          do pid = first, last
            nfield = nfield + 1
            name = 'vf_' // normalize_material_name(matl_model%matl_name(m)) // '_' // &
                normalize_material_name(matl_model%phase_name(pid))
            this%vfrac(nfield) = this%file%register_temporal_cell_data(name, scalar_mold)
          end do
        end if
      end do
    end if
    call register_temporal_fields(this, temporal_output, stat, errmsg)
    if (stat /= 0) then
      call this%file%close()
      return
    end if
    this%is_open = .true.
  end subroutine


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
  end function

  subroutine write_solution(this, time, pressure, velocity, enthalpy, temperature, matl_model, vfrac, &
      temporal_output, flow_active)
    class(ns_ht_2d_vtkhdf_writer), intent(inout) :: this
    real(r8), intent(in) :: time, pressure(:), velocity(:,:), enthalpy(:), temperature(:)
    type(material_model), intent(in) :: matl_model
    real(r8), intent(in) :: vfrac(:,:)
    type(parameter_list), intent(inout) :: temporal_output
    logical, intent(in) :: flow_active(:)
    real(r8), allocatable :: p(:), v(:,:), H(:), T(:), vf(:), beta(:,:)
    real(r8) :: qnan
    integer :: c, m, nphase, phase, field

    allocate(p(this%mesh%ncell), v(3,this%mesh%ncell), H(this%mesh%ncell), T(this%mesh%ncell))
    ASSERT(size(flow_active) == this%mesh%ncell)
    p = pressure(:this%mesh%ncell)
    v = 0.0_r8
    v(:2,:) = velocity(:2,:this%mesh%ncell)
    qnan = ieee_value(0.0_r8, ieee_quiet_nan)
    where (.not.flow_active) p = qnan
    where (spread(.not.flow_active, dim=1, ncopies=3)) v = qnan
    H(:this%mesh%ncell_onP) = enthalpy
    T(:this%mesh%ncell_onP) = temperature
    call this%mesh%cell_imap%gather_offp(H)
    call this%mesh%cell_imap%gather_offp(T)
    call this%file%start_time_step(time)
    call this%file%write_cell_data(this%pressure, p)
    call this%file%write_cell_data(this%velocity, v)
    call this%file%write_cell_data(this%enthalpy, H)
    call this%file%write_cell_data(this%temperature, T)
    if (allocated(this%vfrac)) then
      allocate(vf(this%mesh%ncell))
      field = 0
      do m = 1, matl_model%nmatl
        if (matl_model%have_void .and. m == matl_model%nmatl) then
          nphase = 1
        else
          nphase = matl_model%num_matl_phase(m)
        end if
        if (matl_model%nmatl == 1 .and. nphase == 1) cycle
        if (nphase == 1) then
          field = field + 1
          vf(:this%mesh%ncell_onP) = vfrac(m,:)
          call this%mesh%cell_imap%gather_offp(vf)
          call this%file%write_cell_data(this%vfrac(field), vf)
        else
          allocate(beta(nphase, this%mesh%ncell_onP))
          do c = 1, this%mesh%ncell_onP
            if (vfrac(m,c) > 0.0_r8) then
              call matl_model%get_matl_phase_frac(m, temperature(c), beta(:,c))
            else
              beta(:,c) = 0.0_r8
            end if
          end do
          do phase = 1, nphase
            field = field + 1
            vf(:this%mesh%ncell_onP) = vfrac(m,:) * beta(phase,:)
            call this%mesh%cell_imap%gather_offp(vf)
            call this%file%write_cell_data(this%vfrac(field), vf)
          end do
          deallocate(beta)
        end if
      end do
    end if
    call write_temporal_fields(this, temporal_output)
    call this%file%finalize_time_step()
    call this%file%flush()
  end subroutine


  !! Register a fixed set of scalar temporal field data from TEMPORAL_OUTPUT.
  !! VTKHDF requires this registration before the first time step is started.
  subroutine register_temporal_fields(this, temporal_output, stat, errmsg)

    class(ns_ht_2d_vtkhdf_writer), intent(inout) :: this
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

    class(ns_ht_2d_vtkhdf_writer), intent(inout) :: this
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
    class(ns_ht_2d_vtkhdf_writer), intent(inout) :: this
    if (this%is_open) call this%file%close()
    this%is_open = .false.
    if (allocated(this%vfrac)) deallocate(this%vfrac)
    if (allocated(this%temporal_fields)) deallocate(this%temporal_fields)
    nullify(this%mesh)
  end subroutine

end module ns_ht_2d_vtkhdf_writer_type
