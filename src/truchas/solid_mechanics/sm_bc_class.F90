!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module sm_bc_class

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use sm_bc_node_list_type
  use sm_bc_list_type
  use pcsr_matrix_type
  implicit none
  private

  type, abstract, public :: sm_bc
    integer, allocatable :: index(:)
    logical :: enabled

    ! visualization arrays
    real(r8), allocatable :: displacement(:), traction(:)
  contains
    procedure(init), deferred :: init
    procedure :: add_graph_links
    procedure(apply), deferred :: apply
    procedure(compute_deriv_diag), deferred :: compute_deriv_diag
    procedure(compute_deriv_full), deferred :: compute_deriv_full
  end type sm_bc

  type, public :: sm_bc_box
    class(sm_bc), allocatable :: p
  end type sm_bc_box

  abstract interface
    subroutine init(this, mesh, nodebc, bc, penalty, distance, traction)
      import sm_bc, r8, unstr_mesh, sm_bc_node_list, sm_bc_list
      class(sm_bc), intent(out) :: this
      type(unstr_mesh), intent(in), target :: mesh
      type(sm_bc_node_list), intent(in) :: nodebc
      type(sm_bc_list), intent(in), target :: bc
      real(r8), intent(in) :: penalty, distance, traction
    end subroutine

    subroutine apply(this, time, displ, ftot, stress_factor, r)
      import sm_bc, r8
      class(sm_bc), intent(inout) :: this
      real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
      real(r8), intent(inout) :: r(:,:)
    end subroutine

    subroutine compute_deriv_diag(this, time, displ, ftot, stress_factor, F, diag)
      import sm_bc, r8
      class(sm_bc), intent(inout) :: this
      real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:), F(:,:,:)
      real(r8), intent(inout) :: diag(:,:)
    end subroutine

    subroutine compute_deriv_full(this, time, displ, ftot, stress_factor, Aforce, A)
      import sm_bc, r8, pcsr_matrix
      class(sm_bc), intent(inout) :: this
      real(r8), intent(in) :: time, displ(:,:), ftot(:,:), stress_factor(:)
      type(pcsr_matrix), intent(in) :: Aforce
      type(pcsr_matrix), intent(inout) :: A
    end subroutine
  end interface

contains

  subroutine add_graph_links(this, gforce, g)
    use pcsr_matrix_type
    class(sm_bc), intent(in) :: this
    type(pcsr_graph), intent(in) :: gforce
    type(pcsr_graph), intent(inout) :: g
    ! no-op for displacement BCs
  end subroutine add_graph_links

end module sm_bc_class
