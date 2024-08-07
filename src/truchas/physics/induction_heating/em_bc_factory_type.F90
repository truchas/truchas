!!
!! EM_BC_FACTORY_TYPE
!!
!! This module defines a derived type that serves as a factory for creating
!! abstract boundary condition objects used by the electromagnetic models.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!! February 2024
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module em_bc_factory_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use simpl_mesh_type
  use ih_source_factory_type
  use parameter_list_type
  use scalar_func_class
  use vector_func_class
  use string_utilities, only: lower_case
  use truchas_logging_services
  implicit none
  private

  type, public :: em_bc_factory
    private
    type(simpl_mesh), pointer :: mesh => null() ! unowned reference
    type(ih_source_factory), pointer :: src_fac => null() ! unowned reference
    type(parameter_list), pointer :: params => null() ! unowned reference
    logical :: use_legacy_bc = .true.
    integer, allocatable :: pec_setid(:), nxH_setid(:)
  contains
    procedure :: init
    procedure :: alloc_nxE_bc
    procedure :: alloc_nxH_bc
    procedure, private :: iterate_list
  end type

  !! Interface for the boundary condition call back subroutine
  abstract interface
    subroutine bc_cb(plist, setids, stat, errmsg)
      import parameter_list
      implicit none
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
  end interface

contains

  subroutine init(this, mesh, src_fac, params)

    class(em_bc_factory), intent(out) :: this
    type(simpl_mesh), intent(inout), target :: mesh
    type(ih_source_factory), intent(in), target :: src_fac
    type(parameter_list), intent(inout), target :: params

    this%mesh => mesh
    this%src_fac => src_fac
    this%params => params%sublist('bc') !TODO: eliminate the need to dig deeper

    call params%get('use-legacy-bc', this%use_legacy_bc, default=.true.)

    if (this%use_legacy_bc) then  ! legacy data in the top level parameter list !FIXME
      block
        use ih_legacy_bc, only: create_ih_face_sets
        call create_ih_face_sets(this%mesh, params, this%pec_setid, this%nxH_setid)
      end block
    end if

  end subroutine

  !FIXME: this is really just a PEC type, not a general nxE type
  subroutine alloc_nxE_bc(this, bc, stat, errmsg)

    use bndry_func1_class
    use pec_bndry_func_type

    class(em_bc_factory), intent(in) :: this
    class(bndry_func1), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(pec_bndry_func), allocatable :: nxE_bc
    type(parameter_list) :: unused_plist

    call TLS_info('  generating "nxE" electromagnetic boundary condition')

    if (this%use_legacy_bc) then
      call proc(unused_plist, this%pec_setid, stat, errmsg)
    else
      call this%iterate_list('PEC', proc, stat, errmsg)
    end if
    if (stat /= 0) return
    if (allocated(nxE_bc)) then
      call nxE_bc%add_complete(stat)
      if (stat /= 0) then
        stat = 0  ! all edges assigned the same value -- no need to caution
        !call TLS_info('NOTE: duplicate edges in PEC boundary condition')
      end if
    else
      call TLS_info('    none specified')
    end if
    call move_alloc(nxE_bc, bc)

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to the nxE BC specification and incrementally builds the BC object
    !! accordingly. NB: The NXE_BC and MESH objects are accessed from the parent
    !! subroutine through host association.

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      if (.not.allocated(nxE_bc)) then
        allocate(nxE_bc)
        call nxE_bc%init(this%mesh)
      end if
      call nxE_bc%add(setids, stat, errmsg)
      if (stat /= 0) return ! these are fatal errors
    end subroutine

  end subroutine alloc_nxE_bc

  subroutine alloc_nxH_bc(this, bc, stat, errmsg, scale_factor, omit_edge_list)

    use bndry_func1_class
    use nxH_bndry_func_type
    class(scalar_func), allocatable :: f
    class(vector_func), allocatable :: g

    class(em_bc_factory), intent(in) :: this
    class(bndry_func1), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    real(r8), intent(in), optional :: scale_factor
    integer, intent(in), optional :: omit_edge_list(:)

    type(nxH_bndry_func), allocatable :: nxH_bc

    call TLS_info('  generating "nxH" electromagnetic boundary condition')

    if (this%use_legacy_bc) then
      call this%src_fac%alloc_H_waveform_func(f)
      call this%src_fac%alloc_H_profile_func(g, scale_factor)
      allocate(nxH_bc)
      call nxH_bc%init(this%mesh)
      call nxH_bc%add(f, g, this%nxH_setid, stat, errmsg)
      if (stat /= 0) return !TODO: does info need to be added to errmsg?
      call nxH_bc%add_complete(omit_edge_list)
      call move_alloc(nxH_bc, bc)
    else
      call this%iterate_list('ih-hfield', proc, stat, errmsg)
      if (stat /= 0) return
      if (allocated(nxH_bc)) then
        call nxH_bc%add_complete(omit_edge_list)
      else
        call TLS_info('    none specified')
      end if
      call move_alloc(nxh_bc, bc)
    end if

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to the ih-hfield BC specification and incrementally builds the BC object
    !! accordingly. NB: The NXH_BC and MESH objects are accessed from the parent
    !! subroutine through host association.

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      if (.not.allocated(nxH_bc)) then
        allocate(nxH_bc)
        call nxH_bc%init(this%mesh)
      end if
      call this%src_fac%alloc_H_waveform_func(f)
      call this%src_fac%alloc_H_profile_func(g, scale_factor)
      call nxH_bc%add(f, g, setids, stat, errmsg)
      if (stat /= 0) return
    end subroutine

  end subroutine alloc_nxH_bc

  !! This auxiliary subroutine iterates over the parameter list and for each BC
  !! sublist that matches the given TYPE, it calls the supplied subroutine
  !! PROC, which is expected to incrementally construct the BC using that data.

  subroutine iterate_list(this, type, proc, stat, errmsg)

    class(em_bc_factory), intent(in) :: this
    character(*), intent(in) :: type
    procedure(bc_cb) :: proc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    integer, allocatable :: setids(:)
    character(:), allocatable :: this_type

    stat = 0
    piter = parameter_list_iterator(this%params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      call plist%get('type', this_type, stat, errmsg)
      if (stat /= 0) exit
      if (lower_case(this_type) == lower_case(type)) then  ! use this sublist
        call TLS_info('    using ELECTROMAGNETICS_BC[' // piter%name() // ']')
        call plist%get('face-set-ids', setids, stat, errmsg)
        if (stat /= 0) exit
        call proc(plist, setids, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next
    end do
    if (stat /= 0) errmsg = 'ELECTROMAGNETICS_BC[' // piter%name() // ']: ' // errmsg

  end subroutine iterate_list

end module em_bc_factory_type
