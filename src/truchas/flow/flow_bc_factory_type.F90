!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module flow_bc_factory_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use parameter_list_type
  use scalar_func_class
  use vector_func_class
  use truchas_logging_services
  implicit none
  private

  type, public :: flow_bc_factory
    private
    type(unstr_mesh),     pointer :: mesh   => null() ! reference only - do not own
    type(parameter_list), pointer :: params => null() ! reference only - do not own
  contains
    procedure :: init
    procedure :: alloc_dir_vel_bc
    procedure :: alloc_zero_vn_bc
    procedure :: alloc_dir_prs_bc
    procedure :: alloc_neu_prs_bc
    procedure, private :: iterate_list
  end type flow_bc_factory

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

  subroutine init(this, mesh, params)
    class(flow_bc_factory), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(in), target :: params
    this%mesh => mesh
    this%params => params
  end subroutine init

  !! Allocate a BNDRY_VFUNC object that defines all faces where a prescribed
  !! velocity is to be imposed together with the corresponding velocity data.

  subroutine alloc_dir_vel_bc(this, bc, stat, errmsg)

    use bndry_vfunc_class
    use bndry_face_vfunc_type

    class(flow_bc_factory), intent(in) :: this
    class(bndry_vfunc), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_face_vfunc), allocatable :: bff

    allocate(bff)
    call bff%init(this%mesh, bndry_only=.false.)
    call TLS_info('Generating velocity boundary condition for "velocity" type')
    call this%iterate_list('velocity', proc1, stat, errmsg)
    if (stat /= 0) return
    call TLS_info('Generating velocity boundary condition for "no-slip" type')
    call this%iterate_list('no-slip', proc2, stat, errmsg)
    if (stat /= 0) return
    call bff%add_complete
    call move_alloc(bff, bc)

  contains

    !! These call-back subroutines process parameter list data that is specific
    !! to velocity and no-slip BC specifications and incrementally builds the
    !! BC object accordingly. NB: The BFF and MESH objects are accessed from
    !! the parent subroutine through host association.

    subroutine proc1(plist, setids, stat, errmsg)
      use vector_func_class
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(vector_func), allocatable :: f
      call get_vector_func(plist, 'velocity', f, stat, errmsg)
      if (stat /= 0) return
      call bff%add(f, setids, stat, errmsg)
    end subroutine proc1

    subroutine proc2(plist, setids, stat, errmsg)
      use vector_func_factories, only: alloc_const_vector_func
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(vector_func), allocatable :: f
      call alloc_const_vector_func(f, [0.0_r8, 0.0_r8, 0.0_r8])
      call bff%add(f, setids, stat, errmsg)
    end subroutine proc2

  end subroutine alloc_dir_vel_bc

  !! Allocate a BNDRY_FUNC object that identifies all faces where zero normal
  !! velocity is to be imposed (corresponding data is ignored?)

  subroutine alloc_zero_vn_bc(this, bc, stat, errmsg)

    use bndry_func_class
    use bndry_face_func_type

    class(flow_bc_factory), intent(in) :: this
    class(bndry_func), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_face_func), allocatable :: bff

    allocate(bff)
    call bff%init(this%mesh, bndry_only=.false.)

    call TLS_info('Generating velocity boundary condition for "free-slip" type')
    call this%iterate_list('free-slip', proc, stat, errmsg)
    if (stat /= 0) return

    call TLS_info('Generating velocity boundary condition for "marangoni" type')
    call this%iterate_list('marangoni', proc, stat, errmsg)
    if (stat /= 0) return

    call bff%add_complete
    call move_alloc(bff, bc)

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to a pressure BC specifications and incrementally builds the BC object
    !! accordingly. NB: The BFF and MESH objects are accessed from the parent
    !! subroutine through host association.

    subroutine proc(plist, setids, stat, errmsg)
      use scalar_func_factories, only: alloc_const_scalar_func
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      call alloc_const_scalar_func(f, 0.0_r8)
      call bff%add(f, setids, stat, errmsg)
    end subroutine proc

  end subroutine alloc_zero_vn_bc

  !! Allocate a BNDRY_FUNC object that defines all faces where a prescribed
  !! pressure is to be imposed, together with the corresponding pressure data.

  subroutine alloc_dir_prs_bc(this, bc, stat, errmsg)

    use bndry_func_class
    use bndry_face_func_type

    class(flow_bc_factory), intent(in) :: this
    class(bndry_func), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_face_func), allocatable :: bff

    allocate(bff)
    call bff%init(this%mesh, bndry_only=.false.)
    call TLS_info('Generating pressure boundary condition for "pressure" type')
    call this%iterate_list('pressure', proc, stat, errmsg)
    if (stat /= 0) return
    call bff%add_complete
    call move_alloc(bff, bc)

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to a pressure BC specifications and incrementally builds the BC object
    !! accordingly. NB: The BFF and MESH objects are accessed from the parent
    !! subroutine through host association.

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      call get_scalar_func(plist, 'pressure', f, stat, errmsg)
      if (stat /= 0) return
      call bff%add(f, setids, stat, errmsg)
    end subroutine proc

  end subroutine alloc_dir_prs_bc

  !! Allocate a BNDRY_FUNC object that defines all faces where a homogeneous
  !! Neumann pressure condition is to be imposed (data values ignored?)

  subroutine alloc_neu_prs_bc(this, bc, stat, errmsg)

    use bndry_func_class
    use bndry_face_func_type

    class(flow_bc_factory), intent(in) :: this
    class(bndry_func), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_face_func), allocatable :: bff

    allocate(bff)
    call bff%init(this%mesh, bndry_only=.false.)
    call TLS_info('Generating pressure boundary condition for "velocity" type')
    call this%iterate_list('velocity',  proc, stat, errmsg)
    if (stat /= 0) return
    call TLS_info('Generating pressure boundary condition for "no-slip" type')
    call this%iterate_list('no-slip',   proc, stat, errmsg)
    if (stat /= 0) return
    call TLS_info('Generating pressure boundary condition for "free-slip" type')
    call this%iterate_list('free-slip', proc, stat, errmsg)
    if (stat /= 0) return
    call TLS_info('Generating pressure boundary condition for "marangoni" type')
    call this%iterate_list('marangoni', proc, stat, errmsg)
    if (stat /= 0) return
    call bff%add_complete
    call move_alloc(bff, bc)

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to a pressure BC specifications and incrementally builds the BC object
    !! accordingly. NB: The BFF and MESH objects are accessed from the parent
    !! subroutine through host association.

    subroutine proc(plist, setids, stat, errmsg)
      use scalar_func_factories, only: alloc_const_scalar_func
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      call alloc_const_scalar_func(f, 0.0_r8)
      call bff%add(f, setids, stat, errmsg)
    end subroutine proc

  end subroutine alloc_neu_prs_bc

  !! This auxiliary subroutine iterates over the parameter list and for each
  !! BC sublist that matches the given TYPE, it calls the supplied subroutine
  !! PROC, which is expected to incrementally construct the BC using that data.

  subroutine iterate_list(this, type, proc, stat, errmsg)

    use string_utilities, only: lower_case

    class(flow_bc_factory), intent(in) :: this
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
      call plist%get('type', this_type, stat=stat, errmsg=errmsg)
      if (stat /= 0) exit
      if (lower_case(this_type) == type) then  ! use this sublist
        call TLS_info('  using FLOW_BC[' // piter%name() // ']')
        call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
        if (stat /= 0) exit
        call proc(plist, setids, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next
    end do
    if (stat /= 0) errmsg = 'FLOW_BC[' // piter%name() // ']: ' // errmsg

  end subroutine iterate_list

  !! This auxiliary subroutine gets the scalar function specified by the value
  !! of the parameter PARAM in the parameter list PLIST. The parameter value is
  !! either a real acalar or a character string that is the name of a function
  !! in the function table.

  subroutine get_scalar_func(plist, param, f, stat, errmsg)

    use scalar_func_factories
    use scalar_func_table, only: lookup_func

    type(parameter_list), intent(inout) :: plist
    character(*), intent(in) :: param
    class(scalar_func), allocatable, intent(out) :: f
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: const
    character(:), allocatable :: fname

    call plist%get(param, fname, stat=stat)
    if (stat == 0) then ! name of a function
      call lookup_func(fname, f)
      if (.not.allocated(f)) then
        stat = 1
        errmsg = 'unknown function name: ' // fname
        return
      end if
    else  ! it must be a constant value
      call plist%get(param, const, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call alloc_const_scalar_func(f, const)
    end if

  end subroutine get_scalar_func

  !! This auxiliary subroutine gets the vector function specified by the value
  !! of the parameter PARAM in the parameter list PLIST. The parameter value is
  !! either a real array or a character string that is the name of a vector
  !! function in the vector function table.

  subroutine get_vector_func(plist, param, f, stat, errmsg)

    use vector_func_factories
    use vector_func_table, only: lookup_func

    type(parameter_list), intent(inout) :: plist
    character(*), intent(in) :: param
    class(vector_func), allocatable, intent(out) :: f
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable :: const(:)
    character(:), allocatable :: fname

    call plist%get(param, fname, stat=stat)
    if (stat == 0) then ! name of a function
      call lookup_func(fname, f)
      if (.not.allocated(f)) then
        stat = 1
        errmsg = 'unknown function name: ' // fname
        return
      end if
    else  ! it must be a constant value
      call plist%get(param, const, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call alloc_const_vector_func(f, const)
    end if

  end subroutine get_vector_func

end module flow_bc_factory_type
