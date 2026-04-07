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
  use sm_precon_class
  use sm_nlsol_model_type
  use nlsol_type
  use truchas_logging_services
  use truchas_timers
  implicit none
  private

  type, public :: solid_mechanics
    private
    type(sm_model) :: model
    class(sm_precon), allocatable :: precon
    type(sm_nlsol_model) :: solver_model
    type(nlsol) :: solver

    real(r8), allocatable :: displacement(:,:), strain(:,:), stress(:,:)

    integer :: max_outer_iterations = 1
    integer, public :: thermoelastic_niter = 0 ! linear iteration count
    integer, public :: viscoplastic_niter = 0  ! nonlinear iteration count
  contains
    procedure :: init
    procedure :: compute_initial_state
    procedure :: write_checkpoint
    procedure :: read_checkpoint
    procedure :: step
    procedure :: strain_view
    procedure :: stress_view
    procedure :: displacement_view
    procedure :: eval_point_fields
    procedure :: eval_cell_fields
    procedure :: get_contact_set_ids
    procedure :: viscoplasticity_enabled
    procedure :: get_plastic_strain
    procedure :: get_plastic_strain_rate
  end type solid_mechanics

  type, public :: solid_mechanics_point_field_request
    logical :: displacement = .false.
    logical :: gap_displacement = .false.
    logical :: gap_normal_traction = .false.
  end type

  type, public :: solid_mechanics_point_field_values
    real(r8), allocatable :: displacement(:,:)
    real(r8), allocatable :: gap_displacement(:)
    real(r8), allocatable :: gap_normal_traction(:)
  end type

  type, public :: solid_mechanics_cell_field_request
    logical :: total_strain = .false.
    logical :: thermal_strain = .false.
    logical :: stress = .false.
    logical :: rotation = .false.
    logical :: plastic_strain = .false.
    logical :: plastic_strain_rate = .false.
  end type

  type, public :: solid_mechanics_cell_field_values
    real(r8), allocatable :: total_strain(:,:)
    real(r8), allocatable :: thermal_strain(:,:)
    real(r8), allocatable :: stress(:,:)
    real(r8), allocatable :: rotation(:)
    real(r8), allocatable :: plastic_strain(:,:)
    real(r8), allocatable :: plastic_strain_rate(:)
  end type

contains

  ! Note: Contact is nonlinear and it's very important for convergence for the
  ! preconditioner to be accurate. We'll use "outer" iterations to drive NKA
  ! multiple times, recomputing the preconditioner as we drive closer to the
  ! solution. Without contact, the problem is linear so there's no need to do
  ! this.
  subroutine init(this, mesh, params, nmat, lame1f, lame2f, densityf, reference_density, vp)

    use parameter_list_type
    use unstr_mesh_type
    use scalar_func_containers
    use viscoplastic_material_model_types
    use sm_ds_precon_type
    use sm_hypre_precon_type

    class(solid_mechanics), intent(out), target :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(in) :: nmat
    type(scalar_func_box), allocatable, intent(inout) :: lame1f(:), lame2f(:), densityf(:)
    real(r8), intent(in) :: reference_density(:)
    type(viscoplastic_material_model_box), allocatable, intent(inout) :: vp(:)

    integer :: stat
    character(:), allocatable :: errmsg, preconditioner_method
    type(parameter_list), pointer :: plist => null()

    call start_timer("solid mechanics")

    plist => params%sublist("model")
    call this%model%init(mesh, plist, nmat, lame1f, lame2f, densityf, reference_density, vp)

    plist => params%sublist("preconditioner")
    call plist%get("method", preconditioner_method, default="boomeramg")
    select case(preconditioner_method)
    case("ds")
      allocate(sm_ds_precon :: this%precon)
      ! The scaling factor may vary by row for DS.
      this%model%use_uniform_scaling_factor = .false.
    case("boomeramg", "ssor")
      allocate(sm_hypre_precon :: this%precon)
    case default
      call TLS_fatal("Invalid selection for solid mechanics preconditioner_type")
    end select
    call this%precon%init(this%model, plist)

    call this%solver_model%init(this%model, this%precon)

    plist => params%sublist("nonlinear-solver")
    if (.not.plist%is_parameter("nlk-tol")) call plist%set("nlk-tol", 1.0_r8) ! default
    if (.not.plist%is_parameter("nlk-max-iter")) call plist%set("nlk-max-iter", 500) ! default
    if (this%model%bc%contact_active) &
        call plist%get("max-outer-iter", this%max_outer_iterations, default=5) ! see note
    call this%solver%init(this%solver_model, plist, stat, errmsg)
    if (stat /= 0) call tls_fatal("SOLID MECHANICS INIT: " // errmsg)

    allocate(this%displacement(3,this%solver_model%size()/3), this%stress(6,mesh%nnode_onP), &
        this%strain(6,mesh%nnode_onP))
    this%displacement = 0

    call stop_timer("solid mechanics")

  end subroutine init


  !! Computes the initial state of the solver
  !! Right now, this isn't any different than a step. It computes the reference
  !! density before calling step(), but the current version does this during
  !! update_properties() anyways. In the future, update_properties() will rely
  !! on the density already having a valid value.
  !!
  !! Note: VOF is needed as a target to hold an internal reference
  !! between update_properties and update_plastic_strain (through step).
  subroutine compute_initial_state(this, vof, temperature_cc, stat, errmsg)
    class(solid_mechanics), intent(inout), target :: this
    real(r8), intent(in), target :: vof(:,:), temperature_cc(:)
    integer, intent(out) :: stat
    character(:), intent(out), allocatable :: errmsg
    call start_timer("solid mechanics")
    call TLS_info("SM: Computing initial state...")
    call this%model%compute_reference_density(vof)
    call stop_timer("solid mechanics")
    call this%step(0.0_r8, 0.0_r8, vof, temperature_cc, stat, errmsg)
    if (stat /= 0) errmsg = "SM initialization failure: " // errmsg
    call this%model%init_viscoplastic_rate()
    call TLS_info("") ! line break
  end subroutine compute_initial_state


  !! Note: VOF is needed as a target to hold an internal reference
  !! between update_properties and update_plastic_strain (through
  !! solve).
  subroutine step(this, t, dt, vof, temperature_cc, stat, errmsg)

    class(solid_mechanics), intent(inout), target :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), target :: vof(:,:), temperature_cc(:)
    integer, intent(out) :: stat
    character(:), intent(out), allocatable :: errmsg

    !real(r8) :: displ(size(this%displacement)), displ0(size(this%displacement))
    integer :: i
    real(r8) :: displ0(size(this%displacement))
    real(r8), pointer :: displ(:)

    call start_timer("solid mechanics")

    displ(1:size(this%displacement)) => this%displacement
    displ0 = displ
    call this%model%update_properties(vof, temperature_cc)

    do i = 1, this%max_outer_iterations
      call this%solver%solve(t+dt, dt, displ0, displ, stat)
      if (stat == 0) exit
      if (i < this%max_outer_iterations) &
          call TLS_info("SM: retrying with updated preconditioner...")
    end do
    if (stat /= 0) then
      errmsg = "NLK-SM did not converge"
      return
    end if

    call this%model%accept_state
    call this%model%mesh%node_imap%gather_offp(this%displacement)
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


  subroutine eval_point_fields(this, nodes, request, values)

    class(solid_mechanics), intent(inout) :: this
    integer, intent(in) :: nodes(:)
    type(solid_mechanics_point_field_request), intent(in) :: request
    type(solid_mechanics_point_field_values), intent(out) :: values

    if (request%displacement) values%displacement = this%displacement(:,nodes)

    if (request%gap_displacement .or. request%gap_normal_traction) then
      if (request%gap_displacement .and. request%gap_normal_traction) then
        call this%model%bc%get_gap_fields(nodes, values%gap_displacement, &
            values%gap_normal_traction)
      elseif (request%gap_displacement) then
        call this%model%bc%get_gap_fields(nodes, gap_displacement=values%gap_displacement)
      else
        call this%model%bc%get_gap_fields(nodes, gap_traction=values%gap_normal_traction)
      end if
    end if

  end subroutine eval_point_fields

  subroutine eval_cell_fields(this, cells, request, values)

    class(solid_mechanics), intent(inout) :: this
    integer, intent(in) :: cells(:)
    type(solid_mechanics_cell_field_request), intent(in) :: request
    type(solid_mechanics_cell_field_values), intent(out) :: values

    real(r8), allocatable :: grad_displ(:,:,:), total_strain(:,:), thermal_strain(:,:), &
        plastic_strain(:,:), stress(:,:)
    logical :: need_grad_displ, need_total_strain, need_thermal_strain, need_plastic_strain, &
        need_stress

    need_grad_displ = request%total_strain .or. request%rotation .or. request%stress .or. &
        request%plastic_strain_rate
    need_total_strain = request%total_strain .or. request%stress .or. request%plastic_strain_rate
    need_thermal_strain = request%thermal_strain .or. request%stress .or. request%plastic_strain_rate
    need_plastic_strain = request%plastic_strain .or. request%stress .or. request%plastic_strain_rate
    need_stress = request%stress .or. request%plastic_strain_rate

    !! NB: Some subsidiary procedures below may not be safe for off-process cells.
    INSIST(all(cells <= this%model%mesh%ncell_onP))

    if (need_grad_displ) call compute_cell_grad_displacement(this, cells, grad_displ)
    if (need_total_strain) call compute_cell_total_strain(this, grad_displ, total_strain)
    if (need_thermal_strain) thermal_strain = this%model%strain_thermal(:,cells)
    if (request%rotation) call compute_cell_rotation(this, grad_displ, values%rotation)
    if (need_plastic_strain) call get_cell_plastic_strain(this, cells, plastic_strain)
    if (need_stress) call compute_cell_stress(this, cells, total_strain, thermal_strain, &
        plastic_strain, stress)

    if (request%total_strain) call move_alloc(total_strain, values%total_strain)
    if (request%thermal_strain) call move_alloc(thermal_strain, values%thermal_strain)
    if (request%plastic_strain) call move_alloc(plastic_strain, values%plastic_strain)

    if (request%plastic_strain_rate) &
        call get_cell_plastic_strain_rate(this, cells, stress, values%plastic_strain_rate)
    if (request%stress) call move_alloc(stress, values%stress)

  end subroutine eval_cell_fields

  subroutine get_contact_set_ids(this, setids)
    class(solid_mechanics), intent(in) :: this
    integer, allocatable, intent(out) :: setids(:)
    call this%model%bc%get_contact_set_ids(setids)
  end subroutine

  subroutine compute_cell_grad_displacement(this, cells, grad_displ)
    use sm_bc_utilities, only: compute_gradient_node_to_cell
    class(solid_mechanics), intent(in) :: this
    integer, intent(in) :: cells(:)
    real(r8), allocatable, intent(out) :: grad_displ(:,:,:)
    integer :: i, j
    allocate(grad_displ(3,3,size(cells)))
    associate (mesh => this%model%mesh)
      do i = 1, size(cells)
        j = cells(i)
        associate (cn => mesh%cnode(mesh%xcnode(j):mesh%xcnode(j+1)-1))
          call compute_gradient_node_to_cell(mesh%x(:,cn), mesh%volume(j), &
              this%displacement(1,cn), grad_displ(:,1,i))
          call compute_gradient_node_to_cell(mesh%x(:,cn), mesh%volume(j), &
              this%displacement(2,cn), grad_displ(:,2,i))
          call compute_gradient_node_to_cell(mesh%x(:,cn), mesh%volume(j), &
              this%displacement(3,cn), grad_displ(:,3,i))
        end associate
      end do
    end associate
  end subroutine

  subroutine compute_cell_total_strain(this, grad_displ, total_strain)
    class(solid_mechanics), intent(in) :: this
    real(r8), intent(in) :: grad_displ(:,:,:)
    real(r8), allocatable, intent(out) :: total_strain(:,:)
    integer :: i
    allocate(total_strain(6,size(grad_displ,dim=3)))
    do i = 1, size(grad_displ,dim=3)
      total_strain(:,i) = this%model%strain_tensor(grad_displ(:,:,i))
    end do
  end subroutine

  subroutine compute_cell_rotation(this, grad_displ, rotation_magnitude)
    class(solid_mechanics), intent(in) :: this
    real(r8), intent(in) :: grad_displ(:,:,:)
    real(r8), allocatable, intent(out) :: rotation_magnitude(:)
    integer :: i
    real(r8) :: rotation(3)
    allocate(rotation_magnitude(size(grad_displ,dim=3)))
    do i = 1, size(grad_displ,dim=3)
      rotation(1) = (grad_displ(1,2,i) - grad_displ(2,1,i)) / 2
      rotation(2) = (grad_displ(1,3,i) - grad_displ(3,1,i)) / 2
      rotation(3) = (grad_displ(2,3,i) - grad_displ(3,2,i)) / 2
      rotation_magnitude(i) = norm2(rotation)
    end do
  end subroutine

  subroutine get_cell_plastic_strain(this, cells, plastic_strain)
    class(solid_mechanics), intent(in) :: this
    integer, intent(in) :: cells(:)
    real(r8), allocatable, intent(out) :: plastic_strain(:,:)
    integer :: i
    allocate(plastic_strain(6,size(cells)))
    do i = 1, size(cells)
      plastic_strain(:,i) = this%model%plastic_strain_cell(cells(i))
    end do
  end subroutine

  subroutine compute_cell_stress(this, cells, total_strain, thermal_strain, plastic_strain, stress)
    use sm_bc_utilities, only: compute_stress
    class(solid_mechanics), intent(in) :: this
    integer, intent(in) :: cells(:)
    real(r8), intent(in) :: total_strain(:,:), thermal_strain(:,:), plastic_strain(:,:)
    real(r8), allocatable, intent(out) :: stress(:,:)
    integer :: i
    real(r8) :: elastic_strain(6)
    allocate(stress(6,size(cells)))
    do i = 1, size(cells)
      elastic_strain = total_strain(:,i) - thermal_strain(:,i) - plastic_strain(:,i)
      call compute_stress(this%model%lame1(cells(i)), this%model%lame2(cells(i)), &
          elastic_strain, stress(:,i))
    end do
  end subroutine

  subroutine get_cell_plastic_strain_rate(this, cells, stress, plastic_strain_rate)
    class(solid_mechanics), intent(inout) :: this
    integer, intent(in) :: cells(:)
    real(r8), intent(in) :: stress(:,:)
    real(r8), allocatable, intent(out) :: plastic_strain_rate(:)
    integer :: i
    allocate(plastic_strain_rate(size(cells)))
    do i = 1, size(cells)
      plastic_strain_rate(i) = this%model%plastic_strain_rate_cell(cells(i), stress(:,i))
    end do
  end subroutine


  logical function viscoplasticity_enabled(this)
    class(solid_mechanics), intent(in) :: this
    viscoplasticity_enabled = this%model%viscoplasticity_enabled()
  end function

  subroutine get_plastic_strain(this, plastic_strain)
    class(solid_mechanics), intent(in) :: this
    real(r8), intent(out), allocatable :: plastic_strain(:,:)
    call this%model%get_plastic_strain(plastic_strain)
  end subroutine

  subroutine get_plastic_strain_rate(this, plastic_strain_rate)
    class(solid_mechanics), intent(in) :: this
    real(r8), intent(out), allocatable :: plastic_strain_rate(:,:)
    call this%model%get_plastic_strain_rate(plastic_strain_rate)
  end subroutine

  subroutine write_checkpoint(this, seq)
    use truchas_h5_outfile, only: th5_seq_group
    class(solid_mechanics), intent(in) :: this
    class(th5_seq_group), intent(in) :: seq
    call this%model%write_checkpoint(seq)
  end subroutine

  subroutine read_checkpoint(this, unit)
    class(solid_mechanics), intent(inout) :: this
    integer, intent(in) :: unit
    call this%model%read_checkpoint(unit)
  end subroutine

end module solid_mechanics_type
