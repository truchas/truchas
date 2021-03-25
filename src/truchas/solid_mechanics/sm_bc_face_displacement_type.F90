!!
!! SM_BC_FACE_DISPLACEMENT_TYPE
!!
!! This module presents a type for collecting face conditions for solid
!! mechanics. Displacement and tension boundary conditions can be supplied as
!! "face normal" conditions, or as pointing in a direction x/y/z. This type
!! collects all of these into one list of vector-valued conditions. This vector
!! is either calculated from the given x/y/z direction or calculated from the
!! face normal direction. The type will also expose a collection of face sets
!! associated with each face so that a user can decide how to combine these
!! conditions at nodes.
!!
!! A user of this type then doesn't need to know how this value was computed,
!! and instead can focus on how conditions of this sort are used, or interact
!! with other conditions, or interact with each other at nodes which share
!! conditions from multiple neighboring faces.
!!
!! TODO: Generalize this type to work with both displacements and tractions.
!!       This really only affects the parameter lists values are pulled from
!!       during initialization.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! February 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_bc_face_displacement_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use bndry_face_func_type
  implicit none
  private

  type, public :: sm_bc_face_displacement
    private
    integer, allocatable, public :: index(:)
    real(r8), allocatable, public :: value(:,:)
    real(r8), allocatable :: normal_node(:,:)

    real(r8), allocatable :: normal_face(:,:)
    type(bndry_face_func) :: bff(4)
  contains
    procedure :: init
    procedure :: compute
  end type sm_bc_face_displacement

contains

  subroutine init(this, mesh, params, type, stat, errmsg)

    class(sm_bc_face_displacement), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params
    character(*), intent(in) :: type
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, n, j, norm(3)

    call alloc_bc(type, stat, errmsg)
    if (stat /= 0) return

    this%offset(1) = 0
    do i = 1, 4
      this%offset(i+1) = this%offset(i) + size(this%bff(i)%index)
    end do
    n = this%offset(5)-1
    allocate(this%index(n), this%value(3,n), this%normal(3,n))

    do n = 1, 3
      norm = 0
      norm(n) = 1
      do i = 1, size(this%bff(n)%index)
        j = this%offset(n)+i
        this%index(j) = this%bff(n)%index(i)
        this%normal(:,j) = norm
      end do
    end do

    n = 4
    do i = 1, size(this%bff(n)%index)
      f = this%bff(n)%index(i)
      j = this%offset(n)+i
      this%index(j) = this%bff(n)%index(i)
      this%normal(:,j) = this%mesh%normal(:,f) / norm2(this%mesh%normal(:,f))
    end do

    call compute_index_connectivity(mesh, face_displacement%index, fini, xfini, ni_index)
    call compute_inverted_connectivity(size(ni_index), fini, xfini, nifi, xnifi)

  contains

    subroutine alloc_bc(prefix, stat, errmsg)

      character(*), intent(in) :: prefix
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      character(1), parameter :: dirstr(3) = ['x','y','z']
      integer :: d

      stat = 0

      do d = 1, 3
        call this%bff(d)%init(mesh, bndry_only=.false.)
        call iterate_list(params, prefix//'-'//dirstr(d), prefix, this%bff(d), stat, errmsg)
        if (stat /= 0) return
        call this%bff(d)%add_complete
      end do

      d = 4 ! face-normal condition
      call this%bff(d)%init(mesh, bndry_only=.false.)
      call iterate_list(params, prefix//'-n', prefix, this%bff(d), stat, errmsg)
      if (stat /= 0) return
      call this%bff(d)%add_complete

    end subroutine alloc_bc

  end subroutine init


  subroutine compute(this, t)

    class(sm_normal_traction_bc), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: i, j, n

    do n = 1, 4
      call this%bff(n)%compute(t)
      do i = 1, size(this%bff(n)%index)
        j = this%offset(n)+i
        this%value(:,j) = this%bff(n)%value(i) * this%normal(:,j)
      end do
    end do

  end subroutine compute


  !! This auxiliary subroutine iterates over the parameter list and for each
  !! BC sublist that matches the given TYPE, it calls the supplied subroutine
  !! PROC, which is expected to incrementally construct the BC using that data.
  subroutine iterate_list(params, type, data_label, bff, stat, errmsg)

    use scalar_func_class
    use scalar_func_factories, only: alloc_scalar_func
    use string_utilities, only: lower_case

    type(parameter_list), intent(inout) :: params
    character(*), intent(in) :: type, data_label
    type(bndry_face_func), intent(inout) :: bff
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    integer, allocatable :: setids(:)
    character(:), allocatable :: this_type
    class(scalar_func), allocatable :: f

    stat = 0
    piter = parameter_list_iterator(params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      call plist%get('type', this_type, stat=stat, errmsg=errmsg)
      if (stat /= 0) exit
      if (lower_case(this_type) == type) then  ! use this sublist
        call TLS_info('  using SM_BC[' // piter%name() // ']')
        call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
        if (stat /= 0) exit
        call alloc_scalar_func(plist, data_label, f, stat, errmsg)
        if (stat /= 0) exit
        call bff%add(f, setids, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next
    end do
    if (stat /= 0) errmsg = 'SM_BC[' // piter%name() // ']: ' // errmsg

  end subroutine iterate_list

end module sm_bc_face_displacement_type
