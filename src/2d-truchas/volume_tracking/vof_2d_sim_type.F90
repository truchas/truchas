!!
!! VOF_2D_SIM_TYPE
!!
!! This module defines VOF_2D_SIM, a JSON-configured developer simulation for
!! two-dimensional advection of volume fractions with a prescribed velocity.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module vof_2d_sim_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use unstr_2d_mesh_type
  use vector_func_class
  use vector_func_factories, only: alloc_vector_func
  use region_func_type
  use vol_frac_init_procs, only: compute_volume_fractions
  use vof_2d_solver_type
  use vof_2d_vtkhdf_writer_type
  use simulation_environment_type
  implicit none
  private

  type, public :: vof_2d_sim
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    class(vector_func), allocatable :: velocity
    type(vof_2d_solver) :: solver
    type(vof_2d_vtkhdf_writer) :: output
    real(r8), allocatable :: vfrac(:,:), output_times(:)
    real(r8) :: time_step
  contains
    final :: delete
    procedure :: init
    procedure :: run
  end type

contains

  subroutine delete(this)
    type(vof_2d_sim), intent(inout) :: this
    call this%output%close()
    if (associated(this%mesh)) deallocate(this%mesh)
  end subroutine


  subroutine init(this, env, params, stat, errmsg)

    use unstr_2d_mesh_factory
    use geom_axisymmetric, only: mesh_axisymmetry_mod

    class(vof_2d_sim), intent(out) :: this
    type(simulation_environment), intent(in) :: env
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist
    type(region_func) :: regions
    character(:), allocatable :: algorithm, context
    real(r8), allocatable :: owned_vfrac(:,:)
    logical :: axisymmetric
    integer :: refinement_level, material, nmat

    stat = 0
    if (.not.params%is_sublist('mesh')) then
      stat = 1
      errmsg = 'missing "mesh" sublist parameter'
      return
    end if
    plist => params%sublist('mesh')
    this%mesh => new_unstr_2d_mesh(env, plist, stat, errmsg)
    if (stat /= 0) then
      errmsg = 'processing ' // plist%path() // ': ' // errmsg
      return
    end if
    call this%mesh%init_face_centroid

    if (.not.params%is_sublist('vof-model')) then
      stat = 1
      errmsg = 'missing "vof-model" sublist parameter'
      return
    end if
    plist => params%sublist('vof-model')
    context = 'processing ' // plist%path() // ': '
    call plist%get('algorithm', algorithm, stat, errmsg, default='geometric')
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    call plist%get('axisymmetric', axisymmetric, stat, errmsg, default=.false.)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (axisymmetric) call mesh_axisymmetry_mod(this%mesh)
    call alloc_vector_func(plist, 'velocity', this%velocity, stat, errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (this%velocity%dim /= 2) then
      stat = 1
      errmsg = context // '"velocity" must be a two-component vector function'
      return
    end if

    if (.not.params%is_sublist('initial-regions')) then
      stat = 1
      errmsg = 'missing "initial-regions" sublist parameter'
      return
    end if
    plist => params%sublist('initial-regions')
    context = 'processing ' // plist%path() // ': '
    call regions%init(env, this%mesh, plist, stat, errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    nmat = regions%num_region()
    call params%get('initial-region-refinement-level', refinement_level, stat, errmsg, default=6)
    if (stat /= 0) return
    if (refinement_level < 0) then
      stat = 1
      errmsg = '"initial-region-refinement-level" must be >= 0'
      return
    end if
    allocate(owned_vfrac(nmat,this%mesh%ncell_onP), this%vfrac(nmat,this%mesh%ncell))
    call compute_volume_fractions(this%mesh, regions, refinement_level, owned_vfrac, stat)
    if (stat /= 0) then
      errmsg = context // 'unable to compute volume fractions'
      return
    end if
    do material = 1, nmat
      this%vfrac(material,:this%mesh%ncell_onP) = owned_vfrac(material,:)
      call this%mesh%cell_imap%gather_offp(this%vfrac(material,:))
    end do

    call this%solver%init(env, this%mesh, nmat, algorithm, axisymmetric, stat, errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (.not.params%is_sublist('sim-control')) then
      stat = 1
      errmsg = 'missing "sim-control" sublist parameter'
      return
    end if
    plist => params%sublist('sim-control')
    context = 'processing ' // plist%path() // ': '
    call plist%get('time-step', this%time_step, stat, errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    call plist%get('output-times', this%output_times, stat, errmsg)
    if (stat /= 0) then
      errmsg = context // errmsg
      return
    end if
    if (this%time_step <= 0.0_r8 .or. size(this%output_times) == 0 .or. &
        this%output_times(1) /= 0.0_r8 .or. any(this%output_times(2:) <= this%output_times(:size(this%output_times)-1))) then
      stat = 1
      errmsg = context // 'require time-step > 0 and strictly increasing output-times beginning at 0'
      return
    end if
    call this%output%open(env, this%mesh, nmat, stat, errmsg)
    if (stat /= 0) errmsg = 'opening VTKHDF output: ' // errmsg
  end subroutine


  subroutine run(this, env, stat, errmsg)

    class(vof_2d_sim), intent(inout) :: this
    type(simulation_environment), intent(inout) :: env
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n, nmat
    real(r8) :: time, dt

    stat = 0
    nmat = size(this%vfrac,1)
    time = this%output_times(1)
    call write_output(this, time, nmat)
    do n = 2, size(this%output_times)
      do while (time < this%output_times(n))
        dt = min(this%time_step, this%output_times(n)-time)
        call this%solver%step(env, time, dt, this%velocity, this%vfrac)
        time = time + dt
      end do
      call write_output(this, time, nmat)
    end do
    call this%output%close()
  end subroutine


  subroutine write_output(this, time, nmat)
    class(vof_2d_sim), intent(inout) :: this
    real(r8), intent(in) :: time
    integer, intent(in) :: nmat
    call this%output%write_solution(time, this%vfrac)
  end subroutine

end module vof_2d_sim_type
