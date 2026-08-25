!!
!! THERMAL_BC_FACTORY_TYPE
!!
!! This module defines THERMAL_BC_FACTORY, which creates the boundary
!! condition functions used by the two-dimensional heat-transfer model from a
!! parameter list.  It currently supports temperature, flux, heat-transfer
!! coefficient, and radiation boundary-condition definitions.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, November 2018; updated August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module thermal_bc_factory_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use parameter_list_type
  use scalar_func_class
  use scalar_func_factories, only: alloc_scalar_func, alloc_const_scalar_func
  use string_utilities, only: lower_case
  use simulation_environment_type
  implicit none
  private

  type, public :: thermal_bc_factory
    private
    type(unstr_2d_mesh), pointer :: mesh => null()  ! unowned reference
    type(parameter_list), pointer :: params => null()  ! unowned reference
    real(r8) :: sigma, abszero
  contains
    procedure :: init
    procedure :: alloc_dir_bc
    procedure :: alloc_flux_bc
    procedure :: alloc_htc_bc
    procedure :: alloc_rad_bc
    procedure :: alloc_inflow_bc
    procedure :: alloc_outflow_bc
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

  subroutine init(this, mesh, sigma, abszero, params)

    class(thermal_bc_factory), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    real(r8), intent(in) :: sigma, abszero
    type(parameter_list), target, intent(in) :: params

    this%mesh => mesh
    this%sigma = sigma
    this%abszero = abszero
    this%params => params

  end subroutine


  subroutine alloc_dir_bc(this, bc, env, stat, errmsg)

    use bndry_func1_class
    use bndry_face_func_type

    class(thermal_bc_factory), intent(inout) :: this
    class(bndry_func1), allocatable, intent(out) :: bc
    type(simulation_environment), intent(in) :: env
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_face_func), allocatable :: bff

    call env%simlog%info('  generating "temperature" thermal boundary condition')
    call this%iterate_list(env, 'temperature', proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(bff)) call env%simlog%info('    none specified')

    if (allocated(bff)) then
      call bff%add_complete
      call move_alloc(bff, bc)
    end if

  contains

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f

      call alloc_scalar_func(plist, 'temp', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(bff)) then
        allocate(bff)
        call bff%init(this%mesh)
      end if
      call bff%add(f, setids, stat, errmsg)
    end subroutine

  end subroutine


  subroutine alloc_flux_bc(this, bc, env, stat, errmsg)

    use bndry_func1_class
    use bndry_face_func_type

    class(thermal_bc_factory), intent(inout) :: this
    class(bndry_func1), allocatable, intent(out) :: bc
    type(simulation_environment), intent(in) :: env
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_face_func), allocatable :: bff

    call env%simlog%info('  generating "flux" thermal boundary condition')
    call this%iterate_list(env, 'flux', proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(bff)) call env%simlog%info('    none specified')

    if (allocated(bff)) then
      call bff%add_complete
      call move_alloc(bff, bc)
    end if

  contains

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f

      call alloc_scalar_func(plist, 'flux', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(bff)) then
        allocate(bff)
        call bff%init(this%mesh)
      end if
      call bff%add(f, setids, stat, errmsg)
    end subroutine

  end subroutine


  subroutine alloc_htc_bc(this, bc, env, stat, errmsg)

    use bndry_func2_class
    use htc_bndry_func_type

    class(thermal_bc_factory), intent(inout) :: this
    class(bndry_func2), allocatable, intent(out) :: bc
    type(simulation_environment), intent(in) :: env
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(htc_bndry_func), allocatable :: htc

    call env%simlog%info('  generating "htc" thermal boundary condition')
    call this%iterate_list(env, 'htc', proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(htc)) call env%simlog%info('    none specified')

    if (allocated(htc)) then
      call htc%add_complete
      call move_alloc(htc, bc)
    end if

  contains

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f1, f2

      call alloc_scalar_func(plist, 'htc', f1, stat, errmsg)
      if (stat /= 0) return
      call alloc_scalar_func(plist, 'ambient-temp', f2, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(htc)) then
        allocate(htc)
        call htc%init(this%mesh)
      end if
      call htc%add(f1, f2, setids, stat, errmsg)
    end subroutine

  end subroutine


  subroutine alloc_rad_bc(this, bc, env, stat, errmsg)

    use bndry_func2_class
    use rad_bndry_func_type

    class(thermal_bc_factory), intent(inout) :: this
    class(bndry_func2), allocatable, intent(out) :: bc
    type(simulation_environment), intent(in) :: env
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(rad_bndry_func), allocatable :: rad

    call env%simlog%info('  generating "radiation" thermal boundary condition')
    call this%iterate_list(env, 'radiation', proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(rad)) call env%simlog%info('    none specified')

    if (allocated(rad)) then
      call rad%add_complete
      call move_alloc(rad, bc)
    end if

  contains

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f1, f2

      call alloc_scalar_func(plist, 'emissivity', f1, stat, errmsg)
      if (stat /= 0) return
      call alloc_scalar_func(plist, 'ambient-temp', f2, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(rad)) then
        allocate(rad)
        call rad%init(this%mesh, this%sigma, this%abszero)
      end if
      call rad%add(f1, f2, setids, stat, errmsg)
    end subroutine

  end subroutine


  subroutine alloc_inflow_bc(this, bc, env, stat, errmsg)

    use bndry_func1_class
    use bndry_face_func_type

    class(thermal_bc_factory), intent(inout) :: this
    class(bndry_func1), allocatable, intent(out) :: bc
    type(simulation_environment), intent(in) :: env
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_face_func), allocatable :: bff

    call env%simlog%info('  generating "inflow" thermal boundary condition')
    call this%iterate_list(env, 'inflow', proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(bff)) call env%simlog%info('    none specified')

    if (allocated(bff)) then
      call bff%add_complete
      call move_alloc(bff, bc)
    end if

  contains

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f

      call alloc_scalar_func(plist, 'temp', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(bff)) then
        allocate(bff)
        call bff%init(this%mesh)
      end if
      call bff%add(f, setids, stat, errmsg)
    end subroutine

  end subroutine


  subroutine alloc_outflow_bc(this, bc, env, stat, errmsg)

    use bndry_func1_class
    use bndry_face_func_type

    class(thermal_bc_factory), intent(inout) :: this
    class(bndry_func1), allocatable, intent(out) :: bc
    type(simulation_environment), intent(in) :: env
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_face_func), allocatable :: bff

    call env%simlog%info('  generating "outflow" thermal boundary condition')
    call this%iterate_list(env, 'outflow', proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(bff)) call env%simlog%info('    none specified')

    if (allocated(bff)) then
      call bff%add_complete
      call move_alloc(bff, bc)
    end if

  contains

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      call alloc_const_scalar_func(f, 0.0_r8)
      if (.not.allocated(bff)) then
        allocate(bff)
        call bff%init(this%mesh)
      end if
      call bff%add(f, setids, stat, errmsg)
    end subroutine

  end subroutine


  subroutine iterate_list(this, env, type, proc, stat, errmsg)

    class(thermal_bc_factory), intent(in) :: this
    type(simulation_environment), intent(in) :: env
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
      if (lower_case(this_type) == lower_case(type)) then
        call env%simlog%info('    using thermal BC [' // piter%name() // ']')
        call plist%get('face-set-ids', setids, stat, errmsg)
        if (stat /= 0) exit
        call proc(plist, setids, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next
    end do
    if (stat /= 0) errmsg = 'thermal BC [' // piter%name() // ']: ' // errmsg

  end subroutine

end module thermal_bc_factory_type
