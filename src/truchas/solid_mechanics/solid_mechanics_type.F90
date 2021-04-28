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
    procedure :: compute_initial_state
    procedure :: step
    procedure :: strain_view
    procedure :: stress_view
    procedure :: displacement_view
    procedure :: compute_viz_fields
  end type solid_mechanics

contains

  subroutine init(this, mesh, params, nmat, lame1f, lame2f, densityf, reference_density)

    use parameter_list_type
    use unstr_mesh_type
    use scalar_func_containers

    class(solid_mechanics), intent(out), target :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(in) :: nmat
    type(scalar_func_box), allocatable, intent(inout) :: lame1f(:), lame2f(:), densityf(:)
    real(r8), intent(in) :: reference_density(:)

    integer :: stat
    character(:), allocatable :: errmsg
    type(parameter_list), pointer :: plist => null()

    call start_timer("solid mechanics")

    plist => params%sublist("model")
    call this%model%init(mesh, plist, nmat, lame1f, lame2f, densityf, reference_density)

    plist => params%sublist("preconditioner")
    call this%precon%init(this%model, plist)

    call this%solver_model%init(this%model, this%precon)

    plist => params%sublist("nonlinear-solver")
    if (.not.plist%is_parameter("nlk-tol")) call plist%set("nlk-tol", 1.0_r8) ! default
    call this%solver%init(this%solver_model, plist, stat, errmsg)
    if (stat /= 0) call tls_fatal("SOLID MECHANICS INIT: " // errmsg)

    allocate(this%displacement(3,mesh%nnode), this%stress(6,mesh%nnode_onP), &
        this%strain(6,mesh%nnode_onP))
    this%displacement = 0

    call stop_timer("solid mechanics")

  end subroutine init


  !! Computes the initial state of the solver
  !! Right now, this isn't any different than a step. It computes the reference
  !! density before calling step(), but the current version does this during
  !! update_properties() anyways. In the future, update_properties() will rely
  !! on the density already having a valid value.
  subroutine compute_initial_state(this, vof, temperature_cc, stat, errmsg)
    class(solid_mechanics), intent(inout), target :: this
    real(r8), intent(in) :: vof(:,:), temperature_cc(:)
    integer, intent(out) :: stat
    character(:), intent(out), allocatable :: errmsg
    call start_timer("solid mechanics")
    call TLS_info("SM: Computing initial state...")
    call this%model%compute_reference_density(vof)
    call stop_timer("solid mechanics")
    call this%step(0.0_r8, 1e-6_r8, vof, temperature_cc, stat, errmsg)
    if (stat /= 0) errmsg = "SM initialization failure: " // errmsg
    call TLS_info("") ! line break
  end subroutine compute_initial_state


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
  end function displacement_view

  function strain_view(this) result(view)
    class(solid_mechanics), intent(in), target :: this
    real(r8), pointer :: view(:,:)
    view => this%strain
  end function strain_view

  function stress_view(this) result(view)
    class(solid_mechanics), intent(in), target :: this
    real(r8), pointer :: view(:,:)
    view => this%stress
  end function stress_view


  !! Computes the fields centered on cells for visualization (stresses,
  !! strains, gap displacements, and gap normal tractions). The internal
  !! fields are centered either on nodes or integration points, so
  !! specially-build fields for vizualization must be generated.
  subroutine compute_viz_fields(this, displ, thermal_strain, total_strain, elastic_stress, &
      rotation_magnitude, gap_displacement, gap_normal_traction)

    use sm_bc_utilities, only: compute_gradient_node_to_cell
    use index_partitioning, only: gather_boundary

    class(solid_mechanics), intent(in) :: this
    real(r8), intent(out), allocatable :: displ(:,:), thermal_strain(:,:), total_strain(:,:), &
        elastic_stress(:,:), rotation_magnitude(:), gap_displacement(:), gap_normal_traction(:)

    integer :: j
    real(r8) :: grad_displ(3,3), elastic_strain(6), rotation(3)

    associate (mesh => this%model%mesh)
      allocate(total_strain(6,mesh%ncell_onP), elastic_stress(6,mesh%ncell_onP), &
          rotation_magnitude(mesh%ncell_onP))
      thermal_strain = this%model%thermal_strain(:,:mesh%ncell_onP)
      displ = this%displacement
      call gather_boundary(mesh%node_ip, displ)

      do j = 1, mesh%ncell_onP
        associate (cn => mesh%cnode(mesh%xcnode(j):mesh%xcnode(j+1)-1))
          call compute_gradient_node_to_cell(mesh%x(:,cn), mesh%volume(j), displ(1,cn), grad_displ(:,1))
          call compute_gradient_node_to_cell(mesh%x(:,cn), mesh%volume(j), displ(2,cn), grad_displ(:,2))
          call compute_gradient_node_to_cell(mesh%x(:,cn), mesh%volume(j), displ(3,cn), grad_displ(:,3))
        end associate

        total_strain(1,j) = grad_displ(1,1) ! exx
        total_strain(2,j) = grad_displ(2,2) ! eyy
        total_strain(3,j) = grad_displ(3,3) ! ezz
        total_strain(4,j) = (grad_displ(1,2) + grad_displ(2,1)) / 2 ! exy
        total_strain(5,j) = (grad_displ(1,3) + grad_displ(3,1)) / 2 ! exz
        total_strain(6,j) = (grad_displ(2,3) + grad_displ(3,2)) / 2 ! eyz

        rotation(1) = (grad_displ(1,2) - grad_displ(2,1)) / 2
        rotation(2) = (grad_displ(1,3) - grad_displ(3,1)) / 2
        rotation(3) = (grad_displ(2,3) - grad_displ(3,2)) / 2
        rotation_magnitude(j) = norm2(rotation)

        elastic_strain = total_strain(:,j) - thermal_strain(:,j)
        call this%model%compute_stress(this%model%lame1(j), this%model%lame2(j), &
            elastic_strain, elastic_stress(:,j))
      end do

      displ = displ(:,:mesh%nnode_onP) ! realloc for viz
    end associate

    call this%model%bc%compute_viz_fields(gap_displacement, gap_normal_traction)

  end subroutine compute_viz_fields

end module solid_mechanics_type
