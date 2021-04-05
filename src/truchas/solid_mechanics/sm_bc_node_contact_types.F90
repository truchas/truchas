!!
!! SM_BC_NODE_CONTACT_TYPES
!!
!! TODO
!!
!! Zach Jibben <zjibben@lanl.gov>
!! March 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_bc_node_contact_types

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use unstr_mesh_type
  use sm_bc_list_type
  use sm_bc_node_list_type
  implicit none
  private

  type, public :: sm_bc_c1
    private
    integer, allocatable, public :: index(:,:)
    real(r8), allocatable, public :: value(:,:,:), dvalue(:,:,:)
    logical, public :: enabled

    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    real(r8) :: penalty, distance, normal_traction
    real(r8), allocatable :: area(:), normal(:,:)
  contains
    procedure :: init => c1_init
    procedure :: compute => c1_compute
    procedure :: compute_deriv => c1_compute_deriv
  end type sm_bc_c1

contains

  !! Initialize the single-normal displacement BC nodes.
  !!
  !! The main difficulty here, and in the other below routines for initializing
  !! node-BCs is identifying the appropriate nodes. Here we take a
  !! face_displacement object which holds BCs on faces, along with some helpful
  !! mapping arrays for looking at nodes connected to those faces, and quickly
  !! searching BC faces neighboring each given node. Thus one private procedure
  !! uses faces surrounding a given node to decide whether that node meets the
  !! conditions for this type.
  !!
  !! Here we're looking for nodes touching only faces with a displacement-BC in
  !! a single normal direction. That normal direction may change smoothly over
  !! a curved surface, so there is a small tolerance for deciding what qualifies
  !! as 'distinct' normal vectors between faces. Currently the tolerance is set
  !! to about 25 degrees.
  subroutine c1_init(this, mesh, nodebc, bc, penalty, distance, traction)

    !use integration_geometry_type

    class(sm_bc_c1), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    !type(integration_geometry), intent(in) :: ig
    type(sm_bc_node_list), intent(in) :: nodebc
    type(sm_bc_list), intent(in), target :: bc
    real(r8), intent(in) :: penalty, distance, traction

    character(32) :: msg
    integer :: nnode, li, ni, xbcid

    this%mesh => mesh
    this%penalty = penalty
    this%distance = distance
    this%normal_traction = traction

    nnode = 0
    do li = 1, size(nodebc%link, dim=2)
      if (is_1contact_node(li)) nnode = nnode + 1
    end do
    allocate(this%index(2,nnode), this%value(3,2,nnode), this%dvalue(3,2,nnode), &
        this%normal(3,nnode), this%area(nnode))

    nnode = 0
    do li = 1, size(nodebc%link, dim=2)
      if (.not.is_1contact_node(li)) cycle
      nnode = nnode + 1
      this%index(:,nnode) = nodebc%node(nodebc%link(:,li))
      ni = nodebc%link(1,li)
      xbcid = nodebc%xbcid(ni)
      this%area(nnode) = norm2(nodebc%normal(:,xbcid))
      this%normal(:,nnode) = nodebc%normal(:,xbcid) / norm2(nodebc%normal(:,xbcid))
    end do

    this%enabled = size(this%index) > 0
    if (this%enabled) then
      write(msg,"('SM-1N nodes: ',i6)") size(this%index)
      call TLS_info(trim(msg))
    end if

#ifndef NDEBUG
    block
      integer :: i, n1, n2
      ! TODO-WARN: Need halo node displacements and stresses/residuals?
      !            For now just get everything working in serial.
      do i = 1, size(this%index, dim=2)
        n1 = this%index(1,i)
        n2 = this%index(2,i)
        if (n1 <= mesh%nnode_onP .or. n2 <= mesh%nnode_onP) then
          INSIST(n1 <= mesh%nnode_onP .and. n2 <= mesh%nnode_onP)
        end if
      end do
    end block
#endif

  contains

    !! To qualify as a 1-contact node, the nodes on both sides of the link may
    !! only have a single BC: a contact. If either has another contact or a
    !! displacement BC, it's out.
    !!
    !! Note that there only can be one contact BC for a given link pair. All
    !! contact BCs are identical, so there's no reason to support multiple
    !! contacts for a pair.
    pure logical function is_1contact_node(li)
      integer, intent(in) :: li
      integer :: i, ni, xfi, bc1

      is_1contact_node = .true.
      do i = 1, 2
        ni = nodebc%link(i,li)
        xfi = nodebc%xbcid(ni)
        bc1 = nodebc%bcid(xfi)

        ! If a non-contact bc is here, skip.
        ! If there is any other BC here at all, skip.
        is_1contact_node = is_1contact_node .and. &
            bc1 > bc%dsz &
            .and. nodebc%xbcid(ni+1)-1 == nodebc%xbcid(ni)
      end do

    end function is_1contact_node

  end subroutine c1_init


  subroutine c1_compute(this, displ, ftot, stress_factor)

    class(sm_bc_c1), intent(inout) :: this
    real(r8), intent(in) :: displ(:,:), ftot(:,:), stress_factor(:)

    integer :: i, n1, n2
    real(r8) :: stress1, stress2, x1, x2, s, tn, l, v(2)

    do i = 1, size(this%index, dim=2)
      n1 = this%index(1,i)
      n2 = this%index(2,i)
      stress1 = dot_product(this%normal(:,i), ftot(:,n1))
      stress2 = dot_product(this%normal(:,i), ftot(:,n2))
      x1 = dot_product(this%normal(:,i), displ(:,n1))
      x2 = dot_product(this%normal(:,i), displ(:,n2))

      s = x2 - x1
      tn = - stress1 / this%area(i) ! TODO-WARN: is the sign right?
      l = contact_factor(s, tn, this%distance, this%normal_traction)

      v(1) = stress2 + this%penalty*(x2 - x1) * stress_factor(n1)
      v(2) = stress1 + this%penalty*(x1 - x2) * stress_factor(n2)

      this%value(:,1,i) = this%normal(:,i) * l * v(1)
      this%value(:,2,i) = this%normal(:,i) * l * v(2)
    end do

  end subroutine c1_compute


  subroutine c1_compute_deriv(this, displ, ftot, stress_factor, diag)

    class(sm_bc_c1), intent(inout) :: this
    real(r8), intent(in) :: displ(:,:), ftot(:,:), stress_factor(:), diag(:,:)

    integer :: i, n1, n2
    real(r8) :: stress1, stress2, x1, x2, dldu1(3), dldu2(3), diag1(3), diag2(3)
    real(r8) :: s, tn, l, dl(2), v(2)

    do i = 1, size(this%index, dim=2)
      n1 = this%index(1,i)
      n2 = this%index(2,i)
      stress1 = dot_product(this%normal(:,i), ftot(:,n1))
      stress2 = dot_product(this%normal(:,i), ftot(:,n2))
      x1 = dot_product(this%normal(:,i), displ(:,n1))
      x2 = dot_product(this%normal(:,i), displ(:,n2))
      diag1 = diag(:,n1) * this%normal(:,i) !* stress_factor(n1)
      !diag2 = diag(:,n2) * this%normal(:,i) !* stress_factor(n2)

      s = x2 - x1
      tn = - stress1 / this%area(i)
      l = contact_factor(s, tn, this%distance, this%normal_traction)
      dl = derivative_contact_factor(s, tn, this%distance, this%normal_traction)
      dldu1 = -dl(1)*this%normal(:,i) - dl(2)*diag1 / this%area(i)
      dldu2 =  dl(1)*this%normal(:,i)

      v(1) = stress2 + this%penalty*(x2 - x1) * stress_factor(n1)
      v(2) = stress1 + this%penalty*(x1 - x2) * stress_factor(n2)

      this%dvalue(:,1,i) = this%normal(:,i) * (-l*this%penalty*this%normal(:,i) * stress_factor(n1) + dldu1*v(1))
      this%dvalue(:,2,i) = this%normal(:,i) * (-l*this%penalty*this%normal(:,i) * stress_factor(n2) + dldu2*v(2))
    end do

  end subroutine c1_compute_deriv


  !! PRIVATE ROUTINES

  ! Given the difference in normal displacements s, and the tensile force normal
  ! to the surface tn, compute the contact factor.
  !
  ! If s <= 0, then the nodes are in contact or inside one another.
  ! If tn <= 0, the normal traction is compressive
  real(r8) function contact_factor(s, tn, distance, traction)

    real(r8), intent(in) :: s, tn, distance, traction

    real(r8) :: ls, lt, x

    if (s <= 0) then
      ls = 1
    else if (s >= distance) then
      ls = 0
    else
      x = s / distance - 1
      ls = x**2 * (2*x + 3)
    end if

    if (tn <= 0) then
      lt = 1
    else if (tn >= traction) then
      lt = 0
    else
      x = tn / traction - 1
      lt = x**2 * (2*x + 3)
    end if

    contact_factor = ls * lt
    ASSERT(contact_factor >= 0 .and. contact_factor <= 1)

  end function contact_factor


  pure function derivative_contact_factor(s, tn, distance, traction) result(dl)

    real(r8), intent(in) :: s, tn, distance, traction
    real(r8) :: dl(2)

    real(r8) :: ls, lt, x

    dl = 0

    if (s <= 0) then
      ls = 1
    else if (s >= distance) then
      ls = 0
    else
      x = s / distance - 1
      ls = x**2 * (2*x + 3)
      dl(1) = 6 * x * (x + 1) / distance
    end if

    if (tn <= 0) then
      lt = 1
    else if (tn >= traction) then
      lt = 0
    else
      x = tn / traction - 1
      lt = x**2 * (2*x + 3)
      dl(2) = 6 * x * (x + 1) / traction
    end if

    dl(1) = dl(1) * lt
    dl(2) = dl(2) * ls

  end function derivative_contact_factor

end module sm_bc_node_contact_types
