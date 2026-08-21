!!
!! HT_2D_SOURCE_FACTORY_TYPE
!!
!! This module defines HT_2D_SOURCE_FACTORY, which creates volumetric thermal
!! source functions for the two-dimensional heat-transfer model from a
!! parameter list.
!!
!! David Neill-Asanza <davidhneill@gmail.com>, August 2020
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2020; updated August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

module ht_2d_source_factory_type

  use unstr_2d_mesh_type
  use parameter_list_type
  use scalar_func_class
  use scalar_func_factories, only: alloc_scalar_func
  use simulation_environment_type
  implicit none
  private

  type, public :: ht_2d_source_factory
    private
    type(unstr_2d_mesh), pointer :: mesh => null()  ! unowned reference
    type(parameter_list), pointer :: params => null()  ! unowned reference
  contains
    procedure :: init
    procedure :: alloc_source_funcs
    procedure, private :: alloc_source_func1
    procedure, private :: alloc_source_func2
    procedure, private :: iterate_list
  end type

  abstract interface
    subroutine src_cb(plist, stat, errmsg)
      import parameter_list
      type(parameter_list), intent(inout) :: plist
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
  end interface

  integer, parameter :: TYPE_NONE = 0
  integer, parameter :: TYPE_SCF1 = 1
  integer, parameter :: TYPE_SCF2 = 2

contains

  subroutine init(this, mesh, params)
    class(ht_2d_source_factory), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(parameter_list), target, intent(in) :: params
    this%mesh => mesh
    this%params => params
  end subroutine


  subroutine alloc_source_funcs(this, env, src, stat, errmsg)

    use scalar_mesh_multifunc_type
    use scalar_cell_func1_type
    use scalar_cell_func2_type

    class(ht_2d_source_factory), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    type(scalar_mesh_multifunc), allocatable, intent(out) :: src
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(scalar_cell_func1), allocatable :: scf1
    type(scalar_cell_func2), allocatable :: scf2

    call this%alloc_source_func1(env, scf1, stat, errmsg)
    if (stat /= 0) return
    call this%alloc_source_func2(env, scf2, stat, errmsg)
    if (stat /= 0) return

    if (allocated(scf1) .or. allocated(scf2)) then
      allocate(src)
      allocate(src%value(this%mesh%ncell))
    end if
    if (allocated(scf1)) then
      call scf1%assemble
      call move_alloc(scf1, src%f1)
    end if
    if (allocated(scf2)) then
      call scf2%assemble
      call move_alloc(scf2, src%f2)
    end if

  end subroutine


  subroutine alloc_source_func1(this, env, src, stat, errmsg)

    use scalar_cell_func1_type

    class(ht_2d_source_factory), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    type(scalar_cell_func1), allocatable, intent(inout) :: src
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    call env%simlog%info('  generating "scalar_cell_func1" thermal source')
    call this%iterate_list(env, TYPE_SCF1, proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(src)) call env%simlog%info('    none specified')

  contains

    subroutine proc(plist, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      character(:), allocatable :: file

      call plist%get('data-file', file, stat, errmsg)
      if (stat /= 0) return
      call alloc_scalar_func(plist, 'prefactor', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(src)) then
        allocate(src)
        call src%init(this%mesh)
      end if
      call src%add(file, f, stat, errmsg)
    end subroutine

  end subroutine


  subroutine alloc_source_func2(this, env, src, stat, errmsg)

    use scalar_cell_func2_type

    class(ht_2d_source_factory), intent(inout) :: this
    type(simulation_environment), intent(in) :: env
    type(scalar_cell_func2), allocatable, intent(inout) :: src
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    call env%simlog%info('  generating "scalar_cell_func2" thermal source')
    call this%iterate_list(env, TYPE_SCF2, proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(src)) call env%simlog%info('    none specified')

  contains

    subroutine proc(plist, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      integer, allocatable :: setids(:)

      call plist%get('cell-set-ids', setids, stat, errmsg)
      if (stat /= 0) return
      call alloc_scalar_func(plist, 'source', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(src)) then
        allocate(src)
        call src%init(this%mesh)
      end if
      call src%add(f, setids, stat, errmsg)
    end subroutine

  end subroutine


  integer function sublist_type(plist)
    type(parameter_list), intent(in) :: plist
    if (plist%is_parameter('data-file')) then
      sublist_type = TYPE_SCF1
    else if (plist%is_parameter('source')) then
      sublist_type = TYPE_SCF2
    else
      sublist_type = TYPE_NONE
    end if
  end function


  subroutine iterate_list(this, env, type, proc, stat, errmsg)

    class(ht_2d_source_factory), intent(in) :: this
    type(simulation_environment), intent(in) :: env
    integer, intent(in) :: type
    procedure(src_cb) :: proc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist

    stat = 0
    piter = parameter_list_iterator(this%params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      if (type == sublist_type(plist)) then
        call env%simlog%info('    using thermal source [' // piter%name() // ']')
        call proc(plist, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next
    end do
    if (stat /= 0) errmsg = 'thermal source [' // piter%name() // ']: ' // errmsg

  end subroutine

end module ht_2d_source_factory_type
