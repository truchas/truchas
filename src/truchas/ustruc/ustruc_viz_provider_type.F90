!! Microstructure package adapter for the VTKHDF provider layer.
!!
!! USTRUC output depends on package configuration, so this provider moves
!! ownership of microstructure VTKHDF output into the provider lifecycle.  A
!! selected stream can opt into USTRUC output by naming the meta-field
!! "ustruc".  The stream-local state owns the per-block VTKHDF handles and all
!! VTKHDF-specific registration and write logic.  The driver is used only as a
!! source of microstructure model data.
module ustruc_viz_provider_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use string_utilities, only: i_to_c
  use ustruc_driver, only: ustruc_enabled, ustruc_num_models, ustruc_ncell_onP, &
      ustruc_get_component_flags, ustruc_has_field, ustruc_get_scalar_field, &
      ustruc_get_vector_field
  use viz_field_registry_type, only: viz_field_registry
  use viz_provider_class, only: viz_provider
  use viz_provider_state_class, only: viz_provider_state
  use vtkhdf_mb_file_type, only: vtkhdf_mb_file, vtkhdf_block_handle, &
      vtkhdf_cell_data_handle
  use truchas_timers, only: start_timer, stop_timer
  implicit none
  private

  type :: field_spec
    character(:), allocatable :: data_name
    character(:), allocatable :: vtk_name
    logical :: is_vector = .false.
  end type

  type :: component_spec
    type(field_spec), allocatable :: fields(:)
  end type

  type :: component_data
    type(vtkhdf_cell_data_handle), allocatable :: fields(:)
  end type

  type :: block_data
    type(component_data), allocatable :: comp(:)
  end type

  type, extends(viz_provider), public :: ustruc_viz_provider
  contains
    procedure :: is_active
    procedure :: register_fields
    procedure :: create_state
  end type

  type, extends(viz_provider_state), public :: ustruc_viz_provider_state
    private
    type(component_spec), allocatable :: spec(:)
    type(block_data), allocatable :: blocks(:)
  contains
    procedure :: register_mesh_block_temporal_data
    procedure :: write_mesh_block_timestep
  end type

contains

  logical function is_active(this)

    class(ustruc_viz_provider), intent(in) :: this

    is_active = ustruc_enabled()

  end function is_active

  subroutine register_fields(this, registry, provider_id)

    class(ustruc_viz_provider), intent(in) :: this
    type(viz_field_registry), intent(inout) :: registry
    integer, intent(in) :: provider_id

    call registry%register_field('ustruc', provider_id)

  end subroutine register_fields

  subroutine create_state(this, nblock, state, field_names)

    class(ustruc_viz_provider), intent(in) :: this
    integer, intent(in) :: nblock
    class(viz_provider_state), allocatable, intent(out) :: state
    character(*), intent(in), optional :: field_names(:)
    type(ustruc_viz_provider_state), allocatable :: new_state
    integer :: iblock, n

    allocate(new_state)
    call build_output_spec(new_state%spec)
    allocate(new_state%blocks(nblock))
    do iblock = 1, size(new_state%blocks)
      allocate(new_state%blocks(iblock)%comp(size(new_state%spec)))
      do n = 1, size(new_state%spec)
        if (allocated(new_state%spec(n)%fields)) then
          allocate(new_state%blocks(iblock)%comp(n)%fields(size(new_state%spec(n)%fields)))
        end if
      end do
    end do
    call move_alloc(new_state, state)

  end subroutine create_state

  subroutine register_mesh_block_temporal_data(this, file, iblock, block, block_cells, &
      block_nodes)

    class(ustruc_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)

    integer :: k, n
    do n = 1, size(this%spec)
      if (.not.allocated(this%spec(n)%fields)) cycle
      do k = 1, size(this%spec(n)%fields)
        if (this%spec(n)%fields(k)%is_vector) then
          this%blocks(iblock)%comp(n)%fields(k) = file%register_temporal_cell_data(block, &
              this%spec(n)%fields(k)%vtk_name, [real(r8)::0, 0, 0])
        else
          this%blocks(iblock)%comp(n)%fields(k) = file%register_temporal_cell_data(block, &
              this%spec(n)%fields(k)%vtk_name, 0.0_r8)
        end if
      end do
    end do
  end subroutine register_mesh_block_temporal_data

  subroutine write_mesh_block_timestep(this, file, iblock, block, block_cells, block_nodes)

    class(ustruc_viz_provider_state), intent(inout) :: this
    type(vtkhdf_mb_file), intent(inout) :: file
    integer, intent(in) :: iblock
    type(vtkhdf_block_handle), intent(in) :: block
    integer, intent(in) :: block_cells(:), block_nodes(:)

    integer :: k, n
    real(r8), allocatable :: scalar_out(:), vector_out(:,:)

    if (.not.allocated(this%blocks(iblock)%comp)) return

    call start_timer('Microstructure')
    allocate(scalar_out(ustruc_ncell_onP()), vector_out(3,ustruc_ncell_onP()))

    do n = 1, size(this%spec)
      if (.not.allocated(this%spec(n)%fields)) cycle
      do k = 1, size(this%spec(n)%fields)
        associate (field => this%spec(n)%fields(k), handle => this%blocks(iblock)%comp(n)%fields(k))
          if (field%is_vector) then
            call ustruc_get_vector_field(n, field%data_name, vector_out)
            call file%write_cell_data(block, handle, vector_out(:,block_cells))
          else
            call ustruc_get_scalar_field(n, field%data_name, scalar_out)
            call file%write_cell_data(block, handle, scalar_out(block_cells))
          end if
        end associate
      end do
    end do

    call stop_timer('Microstructure')

  end subroutine write_mesh_block_timestep

  function vtkhdf_component_field_name(label, component, suffix) result(name)

    character(*), intent(in) :: label, component, suffix
    character(:), allocatable :: name

    name = label // '_' // component // '_' // suffix

  end function vtkhdf_component_field_name

  subroutine build_output_spec(spec)

    type(component_spec), allocatable, intent(out) :: spec(:)

    integer :: n
    character(:), allocatable :: label
    logical :: has_gl, has_ldrd

    allocate(spec(ustruc_num_models()))
    do n = 1, size(spec)
      label = 'ustruc' // i_to_c(n)
      call ustruc_get_component_flags(n, has_gl, has_ldrd)

      if (has_gl) then
        call maybe_append_field(spec(n), n, data_name='gl-G', &
            vtk_name=vtkhdf_component_field_name(label, 'gl', 'g'), is_vector=.true.)
        call maybe_append_field(spec(n), n, data_name='gl-L', &
            vtk_name=vtkhdf_component_field_name(label, 'gl', 'l'), is_vector=.false.)
        call maybe_append_field(spec(n), n, data_name='gl-t_sol', &
            vtk_name=vtkhdf_component_field_name(label, 'gl', 't_sol'), is_vector=.false.)
      end if

      if (has_ldrd) then
        call maybe_append_field(spec(n), n, data_name='ldrd-type', &
            vtk_name=vtkhdf_component_field_name(label, 'ldrd', 'type'), is_vector=.false.)
        call maybe_append_field(spec(n), n, data_name='ldrd-lambda1', &
            vtk_name=vtkhdf_component_field_name(label, 'ldrd', 'lambda1'), is_vector=.false.)
        call maybe_append_field(spec(n), n, data_name='ldrd-lambda2', &
            vtk_name=vtkhdf_component_field_name(label, 'ldrd', 'lambda2'), is_vector=.false.)
        call maybe_append_field(spec(n), n, data_name='ldrd-G', &
            vtk_name=vtkhdf_component_field_name(label, 'ldrd', 'g'), is_vector=.false.)
        call maybe_append_field(spec(n), n, data_name='ldrd-V', &
            vtk_name=vtkhdf_component_field_name(label, 'ldrd', 'v'), is_vector=.false.)
        call maybe_append_field(spec(n), n, data_name='ldrd-t_sol', &
            vtk_name=vtkhdf_component_field_name(label, 'ldrd', 't_sol'), is_vector=.false.)
      end if
    end do

  contains

    subroutine maybe_append_field(comp, model_id, data_name, vtk_name, is_vector)
      type(component_spec), intent(inout) :: comp
      integer, intent(in) :: model_id
      character(*), intent(in) :: data_name, vtk_name
      logical, intent(in) :: is_vector
      if (.not.ustruc_has_field(model_id, data_name)) return
      call append_field(comp, data_name, vtk_name, is_vector)
    end subroutine maybe_append_field

    subroutine append_field(comp, data_name, vtk_name, is_vector)
      type(component_spec), intent(inout) :: comp
      character(*), intent(in) :: data_name, vtk_name
      logical, intent(in) :: is_vector
      type(field_spec), allocatable :: tmp(:)
      integer :: n
      if (allocated(comp%fields)) then
        n = size(comp%fields)
        allocate(tmp(n+1))
        tmp(:n) = comp%fields
        call move_alloc(tmp, comp%fields)
      else
        allocate(comp%fields(1))
      end if
      n = size(comp%fields)
      comp%fields(n)%data_name = data_name
      comp%fields(n)%vtk_name = vtk_name
      comp%fields(n)%is_vector = is_vector
    end subroutine append_field

  end subroutine build_output_spec

end module ustruc_viz_provider_type
