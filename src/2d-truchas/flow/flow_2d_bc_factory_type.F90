!!
!! FLOW_2D_BC_FACTORY_TYPE
!!
!! This module defines FLOW_2D_BC_FACTORY, a concrete factory for the
!! old-style sparse boundary function objects used by two-dimensional flow.
!! It is patterned after THERMAL_BC_FACTORY1, but has no abstract factory
!! base class.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module flow_2d_bc_factory_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use parameter_list_type
  use scalar_func_class
  use vector_func_class
  use scalar_func_factories, only: alloc_scalar_func, alloc_const_scalar_func
  use vector_func_factories, only: alloc_vector_func, alloc_const_vector_func
  use string_utilities, only: lower_case
  use simulation_environment_type
  implicit none
  private

  type, public :: flow_2d_bc_factory
    private
    type(unstr_2d_mesh), pointer :: mesh => null()  ! unowned reference
    type(parameter_list), pointer :: params => null()  ! unowned reference
  contains
    procedure :: init
    procedure :: alloc_dir_vel_bc
    procedure :: alloc_zero_vn_bc
    procedure :: alloc_dir_prs_bc
    procedure :: alloc_neu_prs_bc
    procedure, private :: iterate_list
  end type

  abstract interface
    subroutine bc_cb(plist, setids, stat, errmsg)
      import parameter_list
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
  end interface

contains

  subroutine init(this, mesh, params)
    class(flow_2d_bc_factory), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(parameter_list), target, intent(in) :: params

    this%mesh => mesh
    this%params => params
  end subroutine


  subroutine alloc_dir_vel_bc(this, bc, stat, errmsg, env)
    use bndry_vfunc_class
    use bndry_face_vfunc_type

    class(flow_2d_bc_factory), intent(in) :: this
    class(bndry_vfunc), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(simulation_environment), intent(in), optional :: env

    type(bndry_face_vfunc), allocatable :: bff

    if (present(env)) call env%simlog%info('  generating "velocity" flow boundary condition')
    allocate(bff)
    call bff%init(this%mesh)
    call this%iterate_list('velocity', velocity, stat, errmsg, env)
    if (stat /= 0) return
    call this%iterate_list('no-slip', no_slip, stat, errmsg, env)
    if (stat /= 0) return
    call bff%add_complete()
    call move_alloc(bff, bc)

  contains

    subroutine velocity(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(vector_func), allocatable :: func

      call alloc_vector_func(plist, 'velocity', func, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(bff)) then
        allocate(bff)
        call bff%init(this%mesh)
      end if
      call bff%add(func, setids, stat, errmsg)
    end subroutine

    subroutine no_slip(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(vector_func), allocatable :: func

      call alloc_const_vector_func(func, [0.0_r8, 0.0_r8])
      if (.not.allocated(bff)) then
        allocate(bff)
        call bff%init(this%mesh)
      end if
      call bff%add(func, setids, stat, errmsg)
    end subroutine

  end subroutine


  subroutine alloc_zero_vn_bc(this, bc, stat, errmsg, env)
    use bndry_func1_class
    use bndry_face_func_type

    class(flow_2d_bc_factory), intent(in) :: this
    class(bndry_func1), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(simulation_environment), intent(in), optional :: env

    type(bndry_face_func), allocatable :: bff

    if (present(env)) call env%simlog%info('  generating "free-slip" flow boundary condition')
    allocate(bff)
    call bff%init(this%mesh)
    call this%iterate_list('free-slip', zero, stat, errmsg, env)
    if (stat /= 0) return
    call bff%add_complete()
    call move_alloc(bff, bc)

  contains

    subroutine zero(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: func

      call alloc_const_scalar_func(func, 0.0_r8)
      if (.not.allocated(bff)) then
        allocate(bff)
        call bff%init(this%mesh)
      end if
      call bff%add(func, setids, stat, errmsg)
    end subroutine

  end subroutine


  subroutine alloc_dir_prs_bc(this, bc, stat, errmsg, env)
    use bndry_func1_class
    use bndry_face_func_type

    class(flow_2d_bc_factory), intent(in) :: this
    class(bndry_func1), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(simulation_environment), intent(in), optional :: env

    type(bndry_face_func), allocatable :: bff

    if (present(env)) call env%simlog%info('  generating "pressure" flow boundary condition')
    allocate(bff)
    call bff%init(this%mesh)
    call this%iterate_list('pressure', pressure, stat, errmsg, env)
    if (stat /= 0) return
    call bff%add_complete()
    call move_alloc(bff, bc)

  contains

    subroutine pressure(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: func

      call alloc_scalar_func(plist, 'pressure', func, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(bff)) then
        allocate(bff)
        call bff%init(this%mesh)
      end if
      call bff%add(func, setids, stat, errmsg)
    end subroutine

  end subroutine


  subroutine alloc_neu_prs_bc(this, bc, stat, errmsg, env)
    use bndry_func1_class
    use bndry_face_func_type

    class(flow_2d_bc_factory), intent(in) :: this
    class(bndry_func1), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(simulation_environment), intent(in), optional :: env

    type(bndry_face_func), allocatable :: bff

    allocate(bff)
    call bff%init(this%mesh)
    call this%iterate_list('velocity', zero, stat, errmsg, env)
    if (stat /= 0) return
    call this%iterate_list('no-slip', zero, stat, errmsg, env)
    if (stat /= 0) return
    call this%iterate_list('free-slip', zero, stat, errmsg, env)
    if (stat /= 0) return
    call bff%add_complete()
    call move_alloc(bff, bc)

  contains

    subroutine zero(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: func

      call alloc_const_scalar_func(func, 0.0_r8)
      if (.not.allocated(bff)) then
        allocate(bff)
        call bff%init(this%mesh)
      end if
      call bff%add(func, setids, stat, errmsg)
    end subroutine

  end subroutine


  subroutine iterate_list(this, type, proc, stat, errmsg, env)
    class(flow_2d_bc_factory), intent(in) :: this
    character(*), intent(in) :: type
    procedure(bc_cb) :: proc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    type(simulation_environment), intent(in), optional :: env

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    integer, allocatable :: setids(:)
    character(:), allocatable :: bc_type

    stat = 0
    piter = parameter_list_iterator(this%params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      call plist%get('type', bc_type, stat, errmsg)
      if (stat /= 0) exit
      if (lower_case(bc_type) == type) then
        if (present(env)) call env%simlog%info('    using FLOW_2D_BC[' // piter%name() // ']')
        call plist%get('face-set-ids', setids, stat, errmsg)
        if (stat /= 0) exit
        call proc(plist, setids, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next()
    end do
    if (stat /= 0) errmsg = 'FLOW_2D_BC[' // piter%name() // ']: ' // errmsg
  end subroutine

end module flow_2d_bc_factory_type
