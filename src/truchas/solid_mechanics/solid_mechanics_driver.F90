!!
!! Zach Jibben <zjibben@lanl.gov>
!! August 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module solid_mechanics_driver

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use unstr_mesh_type
  use solid_mechanics_type
  implicit none
  private

  public :: read_solid_mechanics_namelists
  public :: solid_mechanics_enabled
  public :: solid_mechanics_init
  public :: solid_mechanics_write_checkpoint
  public :: solid_mechanics_read_checkpoint
  public :: solid_mechanics_step
  public :: solid_mechanics_strain_view
  public :: solid_mechanics_stress_view
  public :: solid_mechanics_displacement_view
  public :: solid_mechanics_compute_viz_fields
  public :: solid_mechanics_viscoplasticity_enabled
  public :: solid_mechanics_get_plastic_strain
  public :: solid_mechanics_get_plastic_strain_rate

  !! Bundled up all the driver state
  type :: solid_mechanics_data
    type(unstr_mesh), pointer :: mesh => null() ! reference only -- not owned
    type(solid_mechanics) :: sm

    ! matl_phase indirection
    integer :: nphase
    integer, allocatable :: matl_phase(:)

    real(r8), allocatable :: vof(:,:) ! volume fractions for solids

    ! The SM driver shouldn't need this but temperature is currently stored in
    ! Zone%Temp so we need to keep a version on the new mesh as well
    real(r8), allocatable :: temperature_cc(:)
  end type solid_mechanics_data

  type(solid_mechanics_data), target :: this

contains

  logical function solid_mechanics_enabled()
    use physics_module, only: solid_mechanics
    solid_mechanics_enabled = solid_mechanics
  end function solid_mechanics_enabled

  subroutine read_solid_mechanics_namelists(lun)
    use solid_mechanics_namelist
    use solid_mechanics_bc_namelist
    use viscoplastic_model_namelist
    use parameter_list_type
    integer, intent(in) :: lun
    type(parameter_list), pointer :: plist
    call read_solid_mechanics_namelist(lun)
    plist => params%sublist('model')
    plist => plist%sublist('bc')
    call read_solid_mechanics_bc_namelists(lun, plist)
    plist => params%sublist('viscoplastic-material-models')
    call read_viscoplastic_model_namelists(lun, plist)
  end subroutine read_solid_mechanics_namelists


  subroutine solid_mechanics_init()

    use mesh_manager, only: unstr_mesh_ptr
    use material_model_driver, only: matl_model
    use solid_mechanics_namelist, only: params
    use viscoplastic_material_model_types
    use scalar_func_class
    use scalar_func_containers
    use tm_density
    use parameter_list_type

    integer :: p, m, stat
    real(r8) :: ref_temp
    real(r8), allocatable :: ref_dens(:)
    character(:), allocatable :: errmsg, phase_name
    type(parameter_list), pointer :: vpm_plist => null(), phase_name_plist => null()
    type(viscoplastic_material_model_box), allocatable :: vp(:)
    type(scalar_func_box), allocatable :: lame1(:), lame2(:), density(:)
    class(scalar_func), allocatable :: cte

    this%mesh => unstr_mesh_ptr('MAIN')
    vpm_plist => params%sublist('viscoplastic-material-models')
    call init_phase_table

    allocate(this%vof(this%nphase,this%mesh%ncell), this%temperature_cc(this%mesh%ncell))

    allocate(ref_dens(this%nphase))
    allocate(lame1(this%nphase), lame2(this%nphase), density(this%nphase), vp(this%nphase))
    do p = 1, this%nphase
      m = this%matl_phase(p)
      ref_dens(p) = matl_model%const_phase_prop(m, 'tm-ref-density')
      ref_temp = matl_model%const_phase_prop(m, 'tm-ref-temp')
      call matl_model%get_phase_prop(m, 'tm-linear-cte', cte)
      ASSERT(allocated(cte))
      call alloc_tm_density_func(density(p)%f, ref_dens(p), ref_temp, cte, stat, errmsg)
      if (stat /= 0) call tls_fatal("SOLID MECHANICS ERROR: tm-linear-cte -- " // errmsg)

      call matl_model%get_phase_prop(m, 'tm-lame1', lame1(p)%f)
      call matl_model%get_phase_prop(m, 'tm-lame2', lame2(p)%f)
      ASSERT(allocated(lame1(p)%f))
      ASSERT(allocated(lame2(p)%f))

      phase_name = matl_model%phase_name(m)
      !if (vpm_plist%is_sublist(phase_name)) &
      phase_name_plist => vpm_plist%sublist(phase_name)
      call alloc_viscoplastic_material_model(phase_name_plist, vp(p)%m)
    end do

    call this%sm%init(this%mesh, params, this%nphase, lame1, lame2, density, ref_dens, vp)

    call compute_initial_state()

  end subroutine solid_mechanics_init


  subroutine solid_mechanics_step(t, dt)

    use zone_module, only: zone

    real(r8), intent(in) :: t, dt

    integer :: stat
    character(:), allocatable :: errmsg

    this%temperature_cc(1:this%mesh%ncell_onP) = Zone%Temp
    call this%mesh%cell_imap%gather_offp(this%temperature_cc)
    call get_vof(this%vof)

    call this%sm%step(t, dt, this%vof, this%temperature_cc, stat, errmsg)
    if (stat /= 0) call tls_fatal(errmsg)

  end subroutine solid_mechanics_step


  function solid_mechanics_displacement_view() result(view)
    real(r8), pointer :: view(:,:)
    view => this%sm%displacement_view()
  end function

  function solid_mechanics_strain_view() result(view)
    real(r8), pointer :: view(:,:)
    view => this%sm%strain_view()
  end function

  function solid_mechanics_stress_view() result(view)
    real(r8), pointer :: view(:,:)
    view => this%sm%stress_view()
  end function


  subroutine compute_initial_state()

    use zone_module, only: zone

    integer :: stat
    character(:), allocatable :: errmsg

    this%temperature_cc(1:this%mesh%ncell_onP) = Zone%Temp
    call this%mesh%cell_imap%gather_offp(this%temperature_cc)
    call get_vof(this%vof)
    call this%sm%compute_initial_state(this%vof, this%temperature_cc, stat, errmsg)
    if (stat /= 0) call tls_fatal(errmsg)

  end subroutine compute_initial_state


  subroutine solid_mechanics_compute_viz_fields(displ, thermal_strain, total_strain, &
      elastic_stress, rotation_magnitude, gap_displacement, gap_normal_traction, &
      plastic_strain, plastic_strain_rate)
    real(r8), intent(out), allocatable :: displ(:,:), thermal_strain(:,:), total_strain(:,:), &
        elastic_stress(:,:), rotation_magnitude(:), gap_displacement(:), gap_normal_traction(:), &
        plastic_strain(:,:), plastic_strain_rate(:)
    call this%sm%compute_viz_fields(displ, thermal_strain, total_strain, elastic_stress, &
        rotation_magnitude, gap_displacement, gap_normal_traction, &
        plastic_strain, plastic_strain_rate)
  end subroutine

  !! PHASE_TABLE routines.
  !! These routines interface with the matl_model to provide volume fractions
  !! and a mapping between the matl_model and the phases solid mechanics cares
  !! about (solids).

  subroutine init_phase_table()

    use material_model_driver, only: matl_model

    integer :: m, p

    this%nphase = count(.not.matl_model%is_fluid)
    allocate(this%matl_phase(this%nphase))

    p = 1
    do m = 1, matl_model%nphase_real
      if (matl_model%is_fluid(m)) cycle
      this%matl_phase(p) = m
      p = p + 1
    end do

  end subroutine init_phase_table


  subroutine get_vof(vof)

    use legacy_matl_api, only: gather_vof

    real(r8), intent(out) :: vof(:,:)
    integer :: p

    ASSERT(size(vof,dim=1) == this%nphase)

    do p = 1, this%nphase
      call gather_vof(this%matl_phase(p), this%vof(p,:this%mesh%ncell_onP))
    end do
    call this%mesh%cell_imap%gather_offp(vof)

  end subroutine get_vof


  logical function solid_mechanics_viscoplasticity_enabled()
    solid_mechanics_viscoplasticity_enabled = this%sm%viscoplasticity_enabled()
  end function

  subroutine solid_mechanics_get_plastic_strain(plastic_strain)
    real(r8), intent(out), allocatable :: plastic_strain(:,:)
    call this%sm%get_plastic_strain(plastic_strain)
  end subroutine

  subroutine solid_mechanics_get_plastic_strain_rate(plastic_strain_rate)
    real(r8), intent(out), allocatable :: plastic_strain_rate(:,:)
    call this%sm%get_plastic_strain_rate(plastic_strain_rate)
  end subroutine

  subroutine solid_mechanics_write_checkpoint(seq)
    use truchas_h5_outfile, only: th5_seq_group
    class(th5_seq_group), intent(in) :: seq
    call this%sm%write_checkpoint(seq)
  end subroutine

  subroutine solid_mechanics_read_checkpoint(unit, version)
    integer, intent(in) :: unit, version
    call this%sm%read_checkpoint(unit)
  end subroutine

end module solid_mechanics_driver
