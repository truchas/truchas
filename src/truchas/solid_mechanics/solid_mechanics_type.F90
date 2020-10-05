!!
!! This module provides a type to calculate the thermomechanical
!! response of all relevant solid materials present in a given problem.
!!
!! References:
!!
!! - Truchas Physics & Algorithms handbook.
!!
!! - C. Bailey and M. Cross. A finite volume procedure to solve elastic solid
!! mechanics problems in three dimensions on an unstructured mesh. International
!! Journal for Numerical Methods in Engineering, 38:1757–1776, 1995.
!!
!! - Y.D. Fryer, C. Bailey, M. Cross, and C.H. Lai. A control volume procedure
!! for solving the elastic stress-strain equations on an unstructured mesh.
!! Applied Mathematical Modelling, 15:639–645, 1991.
!!
!! - P.S. Follansbee and U.F. Kocks. A constitutive description of the
!! deformation of copper based on the use of the mechanical threshold stress as
!! an internal state variable. Acta Metallurgica, 36:81–93, 1988.
!!
!! - O. C. Zienkiewicz. The Finite Element Method. McGraw-Hill, New York, NY,
!! 1977.
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

module solid_mechanics_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use sm_model_type
  use sm_ds_precon_type
  use sm_nlsol_model_type
  use nlsol_type
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  type, public :: solid_mechanics
    private
    type(sm_model) :: model
    type(sm_ds_precon) :: precon
    type(sm_nlsol_model) :: solver_model
    type(nlsol) :: solver

    real(r8), allocatable :: displacement(:,:), strain(:,:), stress(:,:)

    integer, public :: thermoelastic_niter = 0 ! linear iteration count
    integer, public :: viscoplastic_niter = 0  ! nonlinear iteration count
  contains
    procedure :: init
    procedure :: step
    procedure :: strain_view
    procedure :: stress_view
    procedure :: displacement_view
  end type solid_mechanics

contains

  subroutine init(this, mesh, params, nmat, lame1f, lame2f, densityf)

    use parameter_list_type
    use unstr_mesh_type
    use scalar_func_containers

    class(solid_mechanics), intent(out), target :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(in) :: nmat
    type(scalar_func_box), allocatable, intent(inout) :: lame1f(:), lame2f(:), densityf(:)

    integer :: stat
    character(:), allocatable :: errmsg
    type(parameter_list), pointer :: plist => null()

    call start_timer("solid mechanics")

    plist => params%sublist("model")
    call this%model%init(mesh, plist, nmat, lame1f, lame2f, densityf)

    plist => params%sublist("preconditioner")
    call this%precon%init(this%model, plist)

    call this%solver_model%init(this%model, this%precon)

    plist => params%sublist("nonlinear-solver")
    call this%solver%init(this%solver_model, plist, stat, errmsg)
    if (stat /= 0) call tls_fatal("SOLID MECHANICS INIT: " // errmsg)

    allocate(this%displacement(3,mesh%nnode_onP), this%stress(6,mesh%nnode_onP), &
        this%strain(6,mesh%nnode_onP))
    this%displacement = 0 ! TODO: initial displacement

    call stop_timer("solid mechanics")

  end subroutine init


  subroutine step(this, t, dt, vof, temperature_cc, stat, errmsg)

    class(solid_mechanics), intent(inout), target :: this
    real(r8), intent(in) :: t, dt, vof(:,:), temperature_cc(:)
    integer, intent(out) :: stat
    character(:), intent(out), allocatable :: errmsg

    !real(r8) :: displ(size(this%displacement)), displ0(size(this%displacement))
    real(r8) :: displ0(size(this%displacement))
    real(r8), pointer :: displ(:)

    ASSERT(dt > 0)

    call start_timer("solid mechanics")

    displ(1:size(this%displacement)) => this%displacement
    displ0 = displ
    call this%model%update_properties(vof, temperature_cc)
    call this%solver%solve(t, dt, displ0, displ, stat)
    if (stat /= 0) then
      errmsg = "NLK-SM did not converge"
      return
    end if
    this%viscoplastic_niter = this%solver%itr

    call stop_timer("solid mechanics")

  end subroutine step


  function displacement_view(this) result(view)
    class(solid_mechanics), intent(in), target :: this
    real(r8), pointer :: view(:,:)
    view => this%displacement
  end function

  function strain_view(this) result(view)
    class(solid_mechanics), intent(in), target :: this
    real(r8), pointer :: view(:,:)
    view => this%strain
  end function

  function stress_view(this) result(view)
    class(solid_mechanics), intent(in), target :: this
    real(r8), pointer :: view(:,:)
    view => this%stress
  end function

end module solid_mechanics_type
