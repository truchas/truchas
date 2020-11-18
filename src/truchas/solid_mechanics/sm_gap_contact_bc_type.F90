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
    real(r8), allocatable, public :: value(:), rotation_matrix(:,:,:)
    logical, public :: enabled

    integer, allocatable :: lface_index(:)
    ! type(bndry_face_func), allocatable :: bff
    ! real(r8), allocatable :: area_ip(:)
    integer, allocatable :: fini(:), xfini(:)

    type(unstr_mesh), pointer :: mesh => null() ! reference only - do not own
    type(integration_geometry), pointer :: ig => null() ! reference only - do not own
    integer :: ngroup
    integer, allocatable :: xgroup(:), xgroup_face(:)
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
    procedure :: compute_deriv => compute
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


  subroutine add_complete(this)

    use parallel_communication, only: global_any
    use sm_bc_utilities

    class(sm_gap_contact_bc), intent(inout) :: this

    integer, allocatable :: lnode_index(:)

    INSIST(allocated(this%builder))
    call this%builder%get_link_id_groups(this%ngroup, this%xgroup_face, this%lface_index)
    deallocate(this%builder)
    call scalar_func_list_to_box_array(this%flist, this%displacement)

    call compute_link_index_connectivity(this%ig%xlflnode, this%ig%lflnode, this%lface_index, &
        this%xfini, this%fini, lnode_index)
    this%index = this%ig%lnode(:,lnode_index)

    !call compute_index
    allocate(this%value(size(this%index,dim=2)))
    allocate(this%rotation_matrix(3,3,size(this%index,dim=2)))
    call compute_rotation_matrix

    this%enabled = global_any(size(this%index) > 0)

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

      face_index = this%mesh%lface(1,this%lface_index)
      nip = this%xfini(size(this%xfini))-1
      allocate(normal_ip(3,nip), normal_node(3,size(this%index,dim=2)))

      call compute_ip_normals(face_index, this%xfini, this%mesh, this%ig, normal_ip)
      call compute_node_normals(this%xfini, this%fini, normal_ip, normal_node)

      do ni = 1, size(this%index,dim=2)
        this%rotation_matrix(:,:,ni) = rotation_matrix(normal_node(:,ni))
      end do

    end subroutine compute_rotation_matrix

  end subroutine add_complete


  ! TODO-WARN
  subroutine compute_xgroup(xgroup_face, xgroup)

    integer, intent(in) :: xgroup_face(:)
    integer, intent(out), allocatable :: xgroup(:)

    integer :: n

    allocate(xgroup(size(xgroup_face)))

    xgroup(1) = 1
    !xgroup(size(xgroup_face)) = size(index)+1
    do n = 1, size(xgroup_face)
      !xgroup(n+1)
    end do

  end subroutine compute_xgroup


  subroutine compute(this, t)

    class(sm_gap_contact_bc), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: n, lni, j
    real(r8) :: args(0:size(this%mesh%x,dim=1)), x(3)

    do n = 1, this%ngroup
      do lni = this%xgroup(n), this%xgroup(n+1)-1
        j = this%index(1,lni)
        args(1:) = this%mesh%x(:,j)
        this%value(lni) = this%displacement(n)%eval(args)
      end do
    end do

  end subroutine compute

end module sm_gap_contact_bc_type
