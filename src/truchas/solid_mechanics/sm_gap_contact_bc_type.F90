!!
!! SM_GAP_CONTACT_BC
!!
!! This module defines a type that implements the solid mechanics gap contact
!! condition on a subset of the interface faces of a mesh type, mapped to
!! appropriate interface nodes. See the htc_intfc_func_type for the HTPC
!! analogue.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! November 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_gap_contact_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use scalar_func_containers
  use intfc_link_group_builder_type
  use integration_geometry_type
  implicit none
  private

  type, public :: sm_gap_contact_bc
    private
    integer, allocatable, public :: index(:,:)
    real(r8), allocatable, public :: value(:,:,:), dvalue(:,:,:)
    real(r8), allocatable, public :: rotation_matrix(:,:,:)
    logical, public :: enabled

    real(r8) :: contact_distance = 1e-7_r8
    real(r8) :: contact_normal_traction = 1e4_r8
    real(r8) :: contact_penalty = 1e3_r8

    real(r8), allocatable :: area(:)

    type(unstr_mesh), pointer :: mesh => null() ! reference only - do not own
    type(integration_geometry), pointer :: ig => null() ! reference only - do not own
    integer :: ngroup
    integer, allocatable :: xgroup(:)
    type(scalar_func_box), allocatable :: displacement(:)
    ! temporaries used during construction
    type(intfc_link_group_builder), allocatable :: builder
    type(scalar_func_list) :: flist
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute
    procedure :: compute_value => compute
    procedure :: compute_deriv
    procedure :: contact_factor
    procedure :: derivative_contact_factor
  end type sm_gap_contact_bc

contains

  subroutine init(this, mesh, ig)
    class(sm_gap_contact_bc), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(integration_geometry), intent(in), target :: ig
    this%mesh => mesh
    this%ig => ig
    allocate(this%builder)
    call this%builder%init(mesh)
  end subroutine init


  subroutine add(this, func, setids, stat, errmsg)
    class(sm_gap_contact_bc), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: func
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    call this%builder%add_link_group(setids, stat, errmsg)
    if (stat /= 0) return
    call this%flist%append(func)
  end subroutine add


  !! This BC is along nodes, specifically linked nodes. We have a list of face
  !! links, from which we must get the node links. We also must group the node
  !! links appropriately, and compute normal rotation matrices.
  subroutine add_complete(this)

    use parallel_communication, only: global_any
    use sm_bc_utilities

    class(sm_gap_contact_bc), intent(inout) :: this

    integer :: n
    integer, allocatable :: lnode_index(:), lface_index(:), fini(:), xfini(:), xgroup_face(:)

    INSIST(allocated(this%builder))
    call this%builder%get_link_id_groups(this%ngroup, xgroup_face, lface_index)
    deallocate(this%builder)
    call scalar_func_list_to_box_array(this%flist, this%displacement)

    call compute_link_index_connectivity(this%ig%lflnode, this%ig%xlflnode, lface_index, &
        fini, xfini, lnode_index)

    this%index = this%ig%lnode(:,lnode_index)

    n = size(this%index,dim=2)
    this%enabled = global_any(n > 0)
    allocate(this%value(3,2,n), this%dvalue(3,2,n), this%rotation_matrix(3,3,n), this%area(n))

    call compute_node_groups(fini, xfini, xgroup_face, this%xgroup)
    call compute_rotation_matrix

  contains

    ! !! TODO: This BC is along nodes, specifically linked nodes. We have a list of
    ! !! face links, from which we must get the node links. Such a structure
    ! !! should probably go into integration geometry.
    ! subroutine compute_index()

    !   integer :: count, lfi, lf, xlfi, ln
    !   logical :: touched(this%ig%nlnode)

    !   count = 0
    !   touched = .false.
    !   do lfi = 1, size(this%lface_index)
    !     lf = this%lface_index(lfi)
    !     do xlfi = this%ig%xlnlface(lf), this%ig%xlnlf(lf+1)
    !       ln = this%ig%lnlface(xlfi)
    !       if (.not.touched(ln)) then
    !         count = count + 1
    !         touched(ln) = .true.
    !       end if
    !     end do
    !   end do

    !   allocate(this%index(2,count))

    !   count = 0
    !   touched = .false.
    !   do lfi = 1, size(this%lface_index)
    !     lf = this%lface_index(lfi)
    !     do xlfi = this%ig%xlnlface(lf), this%ig%xlnlf(lf+1)
    !       ln = this%ig%lnlface(xlfi)
    !       if (.not.touched(ln)) then
    !         count = count + 1
    !         touched(ln) = .true.
    !         this%index(:,count) = this%ig%lnode(:,ln)
    !       end if
    !     end do
    !   end do

    ! end subroutine compute_index


    subroutine compute_rotation_matrix()

      integer :: nip, ni
      integer, allocatable :: face_index(:)
      real(r8), allocatable :: normal_ip(:,:), normal_node(:,:)

      if (size(xfini) <= 1) return
      nip = xfini(size(xfini))-1

      face_index = this%mesh%lface(1,lface_index)
      allocate(normal_ip(3,nip), normal_node(3,size(this%index,dim=2)))

      call compute_ip_normals(face_index, xfini, this%mesh, this%ig, normal_ip)
      call compute_node_normals(fini, xfini, normal_ip, normal_node)

      do ni = 1, size(this%index,dim=2)
        this%rotation_matrix(:,:,ni) = rotation_matrix(normal_node(:,ni))
        this%area(ni) = norm2(normal_node(:,ni))
      end do

    end subroutine compute_rotation_matrix


    !! From the face groups, compute the node groups. The nodes have already
    !! been arranged by compute_link_index_connectivity to be listed in an
    !! appropriate order: first nodes from this face group, then nodes from the
    !! next face group that weren't in the previous groups, and so on. From this
    !! assumption, each node group starts at the node following the maximum node
    !! ID from the previous group.
    subroutine compute_node_groups(fini, xfini, xgroup_face, xgroup)

      integer, intent(in) :: fini(:), xfini(:), xgroup_face(:)
      integer, intent(out), allocatable :: xgroup(:)

      integer :: n, fi, k, last_node

      allocate(xgroup(size(xgroup_face)))
      if (size(xgroup) == 0) return

      xgroup(1) = 1
      last_node = 1
      do n = 1, this%ngroup
        do fi = xgroup_face(n), xgroup_face(n+1)-1
          k = maxval(fini(xfini(fi):xfini(fi+1)-1))
          last_node = max(last_node, k)
        end do
        xgroup(n+1) = last_node+1
        ASSERT(xgroup(n+1) >= xgroup(n))
      end do

    end subroutine compute_node_groups

  end subroutine add_complete


  subroutine compute(this, t, displ, ftot, stress)

    class(sm_gap_contact_bc), intent(inout) :: this
    real(r8), intent(in) :: t, displ(:,:), ftot(:,:), stress(:,:)

    ! integer :: n, lni, j
    ! real(r8) :: args(0:size(this%mesh%x,dim=1))
    integer :: i, n1, n2
    real(r8) :: stress1, stress2, x1, x2, s, tn, l, dl(2), dldu1, dldu2, v(2), normal(3)

    ! do n = 1, this%ngroup
    !   do lni = this%xgroup(n), this%xgroup(n+1)-1
    !     j = this%index(1,lni)
    !     args(1:) = this%mesh%x(:,j)
    !     this%value(lni) = this%displacement(n)%eval(args)
    !   end do
    ! end do

    do i = 1, size(this%index, dim=2)
      n1 = this%index(1,i)
      n2 = this%index(2,i)
      normal = this%rotation_matrix(3,:,i)
      stress1 = dot_product(normal, ftot(:,n1)) !+ this%rhs(:,n1)
      stress2 = dot_product(normal, ftot(:,n2)) !+ this%rhs(:,n2)
      x1 = dot_product(normal, displ(:,n1))
      x2 = dot_product(normal, displ(:,n2))

      s = x2 - x1
      tn = dot_product(normal, ftot(:,n1)) ! TODO-WARN: is this the right idea?
      l = this%contact_factor(s, tn)
      ! dl = this%derivative_contact_factor(s, tn)
      ! dldu1 =  dl(1) + dl(2)*diag1
      ! dldu2 = -dl(1) + dl(2)*diag2

      v(1) = stress2 * this%area(i) + this%contact_penalty*(x2 - x1)
      v(2) = stress1 * this%area(i) + this%contact_penalty*(x1 - x2)

      this%value(:,1,i) = normal * l * v(1)
      this%value(:,2,i) = normal * l * v(2)

      ! this%dvalue(:,1,i) = normal * (-l*this%contact_penalty + dldu1*v(1))
      ! this%dvalue(:,2,i) = normal * (-l*this%contact_penalty + dldu2*v(2))

      !r(1:2,n1) = stress1(1:2) ! If there is a sliding constraint... TODO: is this right?
      !r(3,n1) = stress1(3) + stress2(3)
      !r(3,n1) = stress1(3) + stress2(3) + 1d3*(x1(3) - x2(3) + values(i))
      !r(3,n1) = x1(3) - x2(3) + values(i)
      !r(3,n1) = r(3,n1) + l*(stress2(3) + 1d3*(x2(3) - x1(3) + values(i)))
      !r(1:2,n2) = stress2(1:2) ! If there is a sliding constraint... TODO: is this right?
      !r(3,n2) = x1(3) - x2(3) + values(i)
      !r(3,n2) = stress1(3) + stress2(3) + 1d3*(x2(3) - x1(3) + values(i))
      !r(3,n2) = x1(3) - x2(3) + values(i)
    end do

  end subroutine compute


  subroutine compute_deriv(this, t, displ, ftot, stress, diag)

    class(sm_gap_contact_bc), intent(inout) :: this
    real(r8), intent(in) :: t, displ(:,:), ftot(:,:), stress(:,:), diag(:)

    ! integer :: n, lni, j
    ! real(r8) :: args(0:size(this%mesh%x,dim=1))
    integer :: i, n1, n2
    real(r8) :: stress1, stress2, x1, x2, dldu1, dldu2, diag1, diag2
    real(r8) :: s, tn, l, dl(2), v(2), normal(3)

    ! do n = 1, this%ngroup
    !   do lni = this%xgroup(n), this%xgroup(n+1)-1
    !     j = this%index(1,lni)
    !     args(1:) = this%mesh%x(:,j)
    !     this%value(lni) = this%displacement(n)%eval(args)
    !   end do
    ! end do

    do i = 1, size(this%index, dim=2)
      n1 = this%index(1,i)
      n2 = this%index(2,i)
      normal = this%rotation_matrix(3,:,i)
      stress1 = dot_product(normal, ftot(:,n1)) !+ this%rhs(:,n1)
      stress2 = dot_product(normal, ftot(:,n2)) !+ this%rhs(:,n2)
      x1 = dot_product(normal, displ(:,n1))
      x2 = dot_product(normal, displ(:,n2))
      diag1 = dot_product(normal, diag(3*(n1-1)+1:3*(n1-1)+3))
      diag2 = dot_product(normal, diag(3*(n2-1)+1:3*(n2-1)+3))

      s = x2 - x1
      tn = dot_product(normal, ftot(:,n1)) ! TODO-WARN: is this the right idea?
      l = this%contact_factor(s, tn)
      dl = this%derivative_contact_factor(s, tn)
      dldu1 = -dl(1) + dl(2)*diag1
      dldu2 =  dl(1) !+ dl(2)*diag2

      v(1) = stress2 * this%area(i) + this%contact_penalty*(x2 - x1)
      v(2) = stress1 * this%area(i) + this%contact_penalty*(x1 - x2)

      ! this%value(:,1,i) = normal * l*v(1)
      ! this%value(:,2,i) = normal * l*v(2)

      this%dvalue(:,1,i) = normal * (-l*this%contact_penalty * this%area(i) + dldu1*v(1))
      this%dvalue(:,2,i) = normal * (-l*this%contact_penalty * this%area(i) + dldu2*v(2))

      !r(1:2,n1) = stress1(1:2) ! If there is a sliding constraint... TODO: is this right?
      !r(3,n1) = stress1(3) + stress2(3)
      !r(3,n1) = stress1(3) + stress2(3) + 1d3*(x1(3) - x2(3) + values(i))
      !r(3,n1) = x1(3) - x2(3) + values(i)
      !r(3,n1) = r(3,n1) + l*(stress2(3) + 1d3*(x2(3) - x1(3) + values(i)))
      !r(1:2,n2) = stress2(1:2) ! If there is a sliding constraint... TODO: is this right?
      !r(3,n2) = x1(3) - x2(3) + values(i)
      !r(3,n2) = stress1(3) + stress2(3) + 1d3*(x2(3) - x1(3) + values(i))
      !r(3,n2) = x1(3) - x2(3) + values(i)
    end do

  end subroutine compute_deriv


  ! Given the difference in normal displacements s, and the tensile force normal
  ! to the surface tn, compute the contact factor.
  !
  ! If s <= 0, then the nodes are in contact or inside one another.
  ! If tn <= 0, the normal traction is compressive
  real(r8) function contact_factor(this, s, tn)

    class(sm_gap_contact_bc), intent(in) :: this
    real(r8), intent(in) :: s, tn

    real(r8) :: ls, lt, x

    if (s <= 0) then
      ls = 1
    else if (s >= this%contact_distance) then
      ls = 0
    else
      x = s / this%contact_distance - 1
      ls = 2 * x**3 + 3 * x**2
    end if

    if (tn <= 0) then
      lt = 1
    else if (tn >= this%contact_normal_traction) then
      lt = 0
    else
      x = tn / this%contact_normal_traction - 1
      lt = 2 * x**3 + 3 * x**2
    end if

    contact_factor = ls * lt
    ASSERT(contact_factor >= 0 .and. contact_factor <= 1)
    !contact_factor = 1

  end function contact_factor


  pure function derivative_contact_factor(this, s, tn) result(dl)

    class(sm_gap_contact_bc), intent(in) :: this
    real(r8), intent(in) :: s, tn
    real(r8) :: dl(2)

    real(r8) :: ls, lt, x

    dl = 0

    if (s <= 0) then
      ls = 1
    else if (s >= this%contact_distance) then
      ls = 0
    else
      x = s / this%contact_distance - 1
      ls = 2 * x**3 + 3 * x**2
      dl(1) = 6 * x * (x + 1) / this%contact_distance
    end if

    if (tn <= 0) then
      lt = 1
    else if (tn >= this%contact_normal_traction) then
      lt = 0
    else
      x = tn / this%contact_normal_traction - 1
      lt = 2 * x**3 + 3 * x**2
      dl(2) = 6 * x * (x + 1) / this%contact_normal_traction
    end if

    dl(1) = dl(1) * lt
    dl(2) = dl(2) * ls
    !dl = 0

  end function derivative_contact_factor

end module sm_gap_contact_bc_type
