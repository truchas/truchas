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
  use parameter_list_type
  use truchas_logging_services
  use solid_mechanics_type
  implicit none
  private

  public :: solid_mechanics_init
  public :: solid_mechanics_step
  public :: solid_mechanics_strain_view
  public :: solid_mechanics_stress_view
  public :: solid_mechanics_displacement_view

  type(solid_mechanics), target :: this
  type(parameter_list) :: params

contains

  subroutine solid_mechanics_build_parameter_list(lun)

    use physics_module, only: body_force_density
    use sm_namelist

    integer, intent(in) :: lun

    type(parameter_list), pointer :: plist => null()

    call read_sm_namelist(lun, params)

    plist => params%sublist('model')
    call plist%set('body-force-density', body_force_density)

  end subroutine solid_mechanics_build_parameter_list


  subroutine solid_mechanics_init()

    use mesh_manager, only: unstr_mesh_ptr
    use material_model_driver, only: matl_model
    use scalar_func_class
    use scalar_func_containers
    use tm_density

    integer :: m, nmat, stat
    real(r8) :: ref_dens, ref_temp
    character(:), allocatable :: errmsg
    type(scalar_func_box), allocatable :: lame1(:), lame2(:), density(:)
    class(scalar_func), allocatable :: cte

    !! NB: Might need material indirection here.
    allocate(lame1(nmat), lame2(nmat), density(nmat))
    do m = 1, nmat
      ref_dens = matl_model%const_phase_prop(m, 'tm-ref-density')
      ref_temp = matl_model%const_phase_prop(m, 'tm-ref-temp')
      call matl_model%alloc_phase_prop(m, 'tm-linear-cte', cte)
      ASSERT(allocated(cte))
      call alloc_tm_density_func(density(m)%f, ref_dens, ref_temp, cte, stat, errmsg)
      if (stat /= 0) call tls_fatal("SOLID MECHANICS ERROR: tm-linear-cte -- " // errmsg)

      call matl_model%alloc_phase_prop(m, 'tm-lame1', lame1(m)%f)
      call matl_model%alloc_phase_prop(m, 'tm-lame2', lame2(m)%f)
      ASSERT(allocated(lame1(m)%f))
      ASSERT(allocated(lame2(m)%f))
    end do

    call this%init(unstr_mesh_ptr('MAIN'), params, nmat, lame1, lame2, density)
    
  end subroutine solid_mechanics_init


  subroutine solid_mechanics_step(t, dt, vof, temperature_cc)

    real(r8), intent(in) :: t, dt
    real(r8), intent(in) :: vof(:,:), temperature_cc(:)

    integer :: stat
    character(:), allocatable :: errmsg

    call this%step(t, dt, vof, temperature_cc, stat, errmsg)
    if (stat /= 0) call tls_fatal(errmsg)
    
  end subroutine solid_mechanics_step


  function solid_mechanics_displacement_view() result(view)
    real(r8), pointer :: view(:)
    view => this%displacement_view()
  end function

  function solid_mechanics_strain_view() result(view)
    real(r8), pointer :: view(:,:)
    view => this%strain_view()
  end function

  function solid_mechanics_stress_view() result(view)
    real(r8), pointer :: view(:,:)
    view => this%stress_view()
  end function

end module solid_mechanics_driver
