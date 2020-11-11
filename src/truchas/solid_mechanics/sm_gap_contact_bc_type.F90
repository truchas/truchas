!!
!! Zach Jibben <zjibben@lanl.gov>
!! October 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_gap_contact_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use intfc_func2_class
  use unstr_mesh_type
  use integration_geometry_type
  implicit none
  private

  type, extends(intfc_func2), public :: sm_gap_contact_bc
    private
    real(r8), public, allocatable :: rotation_matrix(:,:,:)
    logical, public :: enabled

    integer, allocatable :: lface_index(:)
    ! type(bndry_face_func), allocatable :: bff
    ! real(r8), allocatable :: area_ip(:)
    ! integer, allocatable :: fini(:), xfini(:)

    type(unstr_mesh), pointer :: mesh => null() ! reference only - do not own
    type(integration_geometry), pointer :: ig => null() ! reference only - do not own
    integer, allocatable :: xgroup(:)
    type(scalar_func_box), allocatable :: f(:)
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

    class(sm_gap_contact_bc), intent(inout) :: this

    ASSERT(allocated(this%builder))
    call this%builder%get_link_groups(this%ngroup, this%xgroup, this%lface_index)
    deallocate(this%builder)
    call scalar_func_list_to_box_array(this%flist, this%f)

    call compute_index
    allocate(this%value(size(this%index,dim=2)))
    allocate(this%rotation_matrix(3,3,size(this%index,dim=2)))

    this%enabled = global_any(size(this%index) > 0)

  contains

    !! TODO: This BC is along nodes, specifically linked nodes. We have a list of
    !! face links, from which we must get the node links. Such a structure
    !! should probably go into integration geometry.
    subroutine compute_index()

      integer :: count, lfi, lf, xlfi, ln
      logical :: touched(this%ig%nlnode)

      count = 0
      touched = .false.
      do lfi = 1, size(this%lface_index)
        lf = this%lface_index(lfi)
        do xlfi = this%ig%xlnlface(lf), this%ig%xlnlf(lf+1)
          ln = this%ig%lnlface(xlfi)
          if (.not.touched(ln)) then
            count = count + 1
            touched(ln) = .true.
          end if
        end do
      end do

      allocate(this%index(2,count))

      count = 0
      touched = .false.
      do lfi = 1, size(this%lface_index)
        lf = this%lface_index(lfi)
        do xlfi = this%ig%xlnlface(lf), this%ig%xlnlf(lf+1)
          ln = this%ig%lnlface(xlfi)
          if (.not.touched(ln)) then
            count = count + 1
            touched(ln) = .true.
            this%index(:,count) = this%ig%lnode(:,ln)
          end if
        end do
      end do

    end subroutine compute_index

  end subroutine add_complete


  !! TODO-WARN
  subroutine compute(this, t)

    class(sm_gap_contact_bc), intent(inout) :: this
    real(r8), intent(in) :: t

    !! TODO-WARN

  end subroutine compute

end module sm_gap_contact_bc_type
